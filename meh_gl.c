#include "meh_gl.h"
#include "memory.h"


/*--------------------------------------------------------------------
  Shaders

  -------------------------------------------------------------------*/
// GL shader
typedef struct vfs_s {
    struct vfs_s *prev, *next;
    char	*nam;
    time_t	time;
	
    int	id;
    uint	p; // program index
    uint	vs,fs;

    //--------- vert attribs
    //	char	**attrodr;
    //	int	nattr;
    int	narray; // number of vertex attrib arrays (can be different with nattr -- one attrib like matrix == multiple arrays)


    char	*dim;// dimension of each vector (attribute)
    short	*acount; // attrib count for each name (may be the name of array, matrix => have multiple indices)
    short	*stride; // for interleaved vertex array
    short	*ofs;   //  for interleaved VA
    short	*stream;
    int	nbytes;

    void	*amem;


    //--------- uniforms
    char	*multitex; // multitex[uniform_index] == multitex_i
    int	ntex;
    int	nu; // num of uniforms
    char	*ucount;
    PFNGLUNIFORM1FVPROC	*vectorfunc;
    PFNGLUNIFORMMATRIX2FVPROC	*matrixfunc;
    void	*umem;  // mem used to store ucount[], vectorfunc[] matrixfunc[] and multitex[]

    //--------- render targets
    int nrts;
    int depthrt;
} vfs_t;




/*-------------------------------------
  Parse vertex/fragment shader program file
  ------------------------------------*/

/*--------------------------------------
  Load vertex/fragment shaders

  The shader named "nam" will be loaded twice:
  Once the string "#define GL_FRAGMENT_SHADER" is inserted at the top of the shader; remove the string after loading.
  Then the string "#define GL_VERTEX_SHADER" is inserted ...

  -------------------------------------*/

static void vfs_printerr(char *nam, int shader)
{
    int infolen;
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infolen);

    if(!infolen) return;
    ++infolen; // null terminated

	
    int shaderlen;
    int trash;
    char *lines[10*1024];
    int n=0;

    glGetShaderiv(shader, GL_SHADER_SOURCE_LENGTH, &shaderlen);
    char *shadersource=alloca(shaderlen);
    glGetShaderSource(shader, shaderlen, &trash, shadersource);

    char *c=shadersource;
    while(1){
	lines[n++]=c;
	while(*c && *c!='\n') ++c;
	if(*c=='\n'){
	    *c++=0;
	}
	else break;

	assert(n<10*1024); // its quite possible that the user reads in some quite large, non-shader files.
    }

    char *info = alloca(infolen);
    glGetShaderInfoLog(shader, infolen, &trash, info);
    //		cs_printf("Shader compiling error:\n%s\n", info);
    info[infolen-1]=0;


    c=info;
    char errstr[1024];
    int sz=0;

    while(*c && *c!='(') ++c; // skip the info header
	
    while(*c){
	char numstr[8];
		
	while(!isdigit(*c)) ++c;

	char *ch=numstr;
	while(isdigit(*c)){
	    *ch++=*c++;
	}
	*ch=0;

	int line=atoi(numstr);

	while(!isalnum(*c)) ++c;
	char *bgn=c;
	while(*c && *c!='\n') ++c;
	if(*c=='\n') *c++=0;

	cs_printf("%s:%i: %s\n---%s\n", nam, line, lines[line-1], bgn);
    }
}


#define VFS_INCLUDE_MAX_RECURSION 256

static int vfs_createshader(char *nam, char **snippets, int count, int shader, int type)
{
    if(!shader){
	shader = glCreateShader(type);
    }

    glShaderSource(shader, count, snippets, NULL);
    glCompileShader(shader);
	
    int ok;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
    if(!ok){
	vfs_printerr(nam, shader);
	glDeleteShader(shader);
	return 0;
    }

    return shader;
}

static int vfs_createprog(int p, int vs, int fs)
{
    int p0=p;
    if(!p){
	p = glCreateProgram();
    }

    int error=glGetError();
    if (error) {
	printf("%s\n", gluErrorString(error));
	assert(!error+111*0);
    }

    glAttachShader(p, vs);

    error=glGetError();
    if (error) {
	printf("p0=%i, p=%i, vs=%i%s\n", p0,p,vs, gluErrorString(error));
	assert(!error+222*0);
    }

    glAttachShader(p, fs);
    glLinkProgram(p);
    int ok;
    glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if(!ok){
	int len;
	glGetProgramiv(p, GL_INFO_LOG_LENGTH, &len);
	if(len){
	    char *info=alloca(len);
	    int trash;
	    glGetProgramInfoLog(p, len, &trash, info);
	    printf("GLSL linking error:\n%s\n", info);
	}
	glDeleteProgram(p);
	return 0;
    }

    return p;
}


/*-----------------------------------------------------------
  Each glProgram is just ONE file -- all vertex and fragment shaders are combined into a single file.

  Q: How do we distinguish between vertex and fragment shaders in one file? 
  I mean, how do we tell whether  e.g.  float shit(vec2 abc){ // stuff } is vertex or fragment shader?

  A: The file is parsed TWICE:
  first time with #define GL_VERTEX_SHADER inserted at the beginning of the file;
  second time with #define GL_FRAGMENT_SHADER ...

  So:

  #ifdef VERTEX_SHADER
  //
  // your vertex shader code here
  //
  #endif // VERTEX_SHADER
 

  #ifdef FRAGMENT_SHADER
  //
  // your fragment shader code here
  //
  #endif // FRAGMENT_SHADER

  And a awesome side-effect is that you can place varying variables outside these #ifdef #endif --
  so it will be compiled into both vertex and fragment shaders -- which is the case!

  You can put any code you want to share between vertex and fragment shader outside the #ifdef #endif pairs.

  ----------------------------------------------------------*/

// file's last modified time
static time_t fmodtime(const char *nam)
{
    struct tm* clock;				// create a time structure
    struct stat attrib;			// create a file attribute structure
    stat(nam, &attrib);		// get the attributes of afile.txt
    return attrib.st_mtime;
}


/*
  Shader should be able to be queried by name, so that in-game shader editting will be possible.

*/


/*
  If the shader is already loaded, then reload

*/




static void vfs_init()
{
    id1k_gen(&G.vfs_idmngr);
}

#define VFS_ID(vptr) ((int)((vptr)-G.vfs))
#define VFS_PTR(id) (G.vfs+(id))


int vfs_load(const char *nam)
{
    /////////////////////// Read in & process shader source codes ///////////////////////////////////


    // quit if the shader has already been loaded, and hasn't been modified since then
    vfs_t *oldp = cvar_getp__(nam, VFS_CKSUM);
    vfs_t *p = alloca(sizeof(vfs_t));

    // People may misuse cvar_seti() unintentionally.

    if(oldp){// the shader under the same name already exists.
	time_t t = fmodtime(nam);
	if(t == oldp->time) return oldp->id;
	p->p = oldp->p;
	p->vs = oldp->vs;
	p->fs = oldp->fs;
    }
    else{
	p->p = p->vs = p->fs = 0;
	p->time = fmodtime(nam); // DEBUGGED: didn't set p->time
    }

    char *source=loadtext_alloca(nam);

    if(!source){ // error -- no such shader file
	printf("Error loading shader %s\n", nam);
	return -1;
    }


    // parse #include stuff
    typedef struct snippet_s{
	struct snippet_s *next;
	char *code;
    } snippet_t;

    snippet_t *head = alloca(sizeof(snippet_t));
    head->next = NULL;
    head->code = source;

    int ninc=0;
    int ns=1;
    for(snippet_t *s=head; s; s=s->next){
	char *c=s->code;

	while(c=strstr(c, "#include")){ // skip some "false alarm", such as  "ohyeah#include_hoax"
	    if(ns > VFS_INCLUDE_MAX_RECURSION){ // error -- too many include, maybe indirect self-inclusion?
		return;
	    }


	    if(c==s->code || isspace(c[-1]) && isspace(c[8])){
		*c=0;
		c+=9;
		break;
	    }
			

	    c+=8;
	}

	// if source doesn't contain #include, done.
	if(!c) continue;


	// parse the remaining part of the line.
	// if the #include isn't valid, then just ignore it

	// skip white
	while(isspace(*c)){
	    if(*c=='\n'){
		if(c[-1] != '\\'){ // error -- #include should be one-line
		    return -1;
		}
	    }
	}


	if(*c != '\"'){		// error -- #include invalid stuff (not included in "")
	    return -1;
	}

	char *begin=++c;
	while(*c!='\"' || c[-1]=='\\'){
	    if(!*c){ // error -- unexpected EOF
	    }
		
	    if(*c == '\n'){
		if(c[-1] != '\\'){ // error -- unexpected newline
		}
	    }
		
	    ++c;
	}
	*c=0; // " -> 0


	// begin of snippet3
	char *code3 = c+1;

	// char *begin stores the name of the file to be included
	if(!strcmp(begin, nam)){ // error -- self-inclusion
	    return -1;
	}


	char *code2 = loadtext_alloca(begin);
	if(!code2){ // error -- cannot find the file
	    return -1;
	}

	snippet_t *s2 = alloca(sizeof(snippet_t)*2);
	snippet_t *s3 = s2+1;
	
	s2->code = code2;
	s3->code = code3;

	s2->next = s3;
	s3->next = s->next;
	s->next = s2;

	ns+=2;
    }


    char **snippets = alloca(ns*sizeof(char*)+1);
    snippet_t *s=head;
    for(int i=1; i<=ns; i++, s=s->next){
	snippets[i] = s->code;
    }
	
    // #include found
    // Source code between and after #include become two code snippets (before=snippet1, after=snippet3)
    // The included file becomes the third code snippet (include=snippet2)

    // recursively call vfs_include_r() on the snippet2 and snippet3


    // tail recursion


    /////////////////////// Analysis the shaders -- how many attributes, uniforms ... ////////////////////////////////


    snippets[0] = "#define VERTEX_SHADER\n";
    p->vs = vfs_createshader(nam, snippets, ns+1, p->vs, GL_VERTEX_SHADER);
    snippets[0] = "#define FRAGMENT_SHADER\n";
    p->fs = vfs_createshader(nam, snippets, ns+1, p->fs, GL_FRAGMENT_SHADER);
    p->p = vfs_createprog(p->p, p->vs, p->fs);
    if(!p->p){
	return NULL;
    }

    //---------------------- How many render targets? -----------------------
    {
	int sz;
	glGetShaderiv(p->fs, GL_SHADER_SOURCE_LENGTH, &sz);
	char mem[sz];
	glGetShaderSource(p->fs, sz, NULL, mem);
	mem[sz-1]=0;

	char *code=strstr(mem, "gl_Frag");
	if (strstr(code, "gl_FragDepth")) { // depth render target
	    p->nrts=1;
	    p->depthrt=1;
	}
	else {
	    p->nrts=0;
	    p->depthrt=0;
	}

	if (strstr(code, "gl_FragColor")) {
	    ++p->nrts;
	}
	else {
	    char *c=code;
	    int n=-1;
	    while (1) {
		c=strstr(c, "gl_FragData");
		if (!c) break;
		c+=11;
		while (isspace(*c)) ++c;
		assert(*c);
		int i=atoi(c);
		if (i>n) n=i;
	    }

	    //			assert(n!=-1);
	    p->nrts+=n+1;
	}
    }

    char buf[256];
    int len, type;

    //----------------------- attributes ---------------------------
    int na;
    glGetProgramiv(p->p, GL_ACTIVE_ATTRIBUTES, &na);
	
    p->amem = malloc(na + na*sizeof(short)*3);
    p->dim = p->amem;
    p->acount = p->dim + na;
    p->ofs = p->acount+na;
    p->stride = p->ofs+na;

    p->narray = 0;

    int realna=0;
    int size;
    int totalfloat=0;
    int acount=0;

    for(int i=0; i<na; i+=acount) {
	glGetActiveAttrib(p->p, i, 256, &len, &size, &type, buf);

	int j=glGetAttribLocation(p->p, buf);  // shit! i2 may be different from i!
	if(j==-1) continue; // TODO:  error handling

	// skip other columns of matrix and array elements
	switch(type){
	case GL_FLOAT:
	    p->dim[j] = 1;
	    p->acount[j] = size;
	    break;

	case GL_FLOAT_VEC2:
	case GL_FLOAT_VEC3:
	case GL_FLOAT_VEC4:
	    p->dim[j] = type-GL_FLOAT_VEC2+2;
	    p->acount[j] = size;
	    break;

	case GL_FLOAT_MAT2:
	case GL_FLOAT_MAT3:
	case GL_FLOAT_MAT4:
	    p->dim[j] = type-GL_FLOAT_MAT2+2;
	    p->acount[j] = size * p->dim[j];
	    break;

	case GL_FLOAT_MAT2x3:
	case GL_FLOAT_MAT2x4:
	    p->dim[j] = type-GL_FLOAT_MAT2x3+3;
	    p->acount[j] = size * 2;
	    break;

	case GL_FLOAT_MAT3x2:
	case GL_FLOAT_MAT3x4:
	    p->dim[j] = (type-GL_FLOAT_MAT3x2+1)<<1;
	    p->acount[j] = size * 3;
	    break;

	case GL_FLOAT_MAT4x2:
	case GL_FLOAT_MAT4x3:
	    p->dim[j] = type-GL_FLOAT_MAT4x2+2;
	    p->acount[j] = size * 4;
	    break;
	}

	totalfloat += p->dim[j]*p->acount[j];
	p->narray += p->acount[j];

	acount = p->acount[j];
    }

    p->nbytes = totalfloat * sizeof(float);

    //-------------- uniforms ----------------------------
    int nu;
    glGetProgramiv(p->p, GL_ACTIVE_UNIFORMS, &nu);


    // find all samplers

    p->umem = malloc(nu + nu + nu*sizeof(PFNGLUNIFORM1FVPROC) + nu*sizeof(PFNGLUNIFORMMATRIX2FVPROC));
    p->multitex = p->umem;
    p->ucount = p->umem + nu;
    p->vectorfunc = p->ucount + nu;
    p->matrixfunc = p->vectorfunc + nu;

    memset(p->multitex, -1, nu);
    int ntex=0;

    for(int i=0; i<nu; i++){
	glGetActiveUniform(p->p, i, 256, &len, p->ucount+i, &type, buf);

	int j=glGetUniformLocation(p->p, buf);
	if(j==-1) continue;

	if(type>=GL_SAMPLER_1D && type<=GL_SAMPLER_2D_SHADOW){
	    p->multitex[i] = ntex++;
	    p->vectorfunc[i] = glUniform1fv;
	    p->matrixfunc[i] = NULL;
	}
	else{
	    switch(type){
	    case GL_FLOAT:
		p->vectorfunc[i] = glUniform1fv;
		p->matrixfunc[i] = NULL;
		break;

	    case GL_FLOAT_VEC2:
		p->vectorfunc[i] = glUniform2fv;
		p->matrixfunc[i] = NULL;
		break;

	    case GL_FLOAT_VEC3:
		p->vectorfunc[i] = glUniform3fv;
		p->matrixfunc[i] = NULL;
		break;

	    case GL_FLOAT_VEC4:
		p->vectorfunc[i] = glUniform4fv;
		p->matrixfunc[i] = NULL;
		break;

	    case GL_FLOAT_MAT2:
		p->vectorfunc[i] = NULL;
		p->matrixfunc[i] = glUniformMatrix2fv;
		break;

	    case GL_FLOAT_MAT3:
		p->vectorfunc[i] = NULL;
		p->matrixfunc[i] = glUniformMatrix3fv;
		break;

	    case GL_FLOAT_MAT4:
		p->vectorfunc[i] = NULL;
		p->matrixfunc[i] = glUniformMatrix4fv;
		break;

	    case GL_FLOAT_MAT2x3:
		p->vectorfunc[i] = NULL;
		p->matrixfunc[i] = glUniformMatrix2x3fv;
		break;

	    case GL_FLOAT_MAT2x4:
		p->vectorfunc[i] = NULL;
		p->matrixfunc[i] = glUniformMatrix2x4fv;
		break;

	    case GL_FLOAT_MAT3x2:
		p->vectorfunc[i] = NULL;
		p->matrixfunc[i] = glUniformMatrix3x2fv;
		break;

	    case GL_FLOAT_MAT3x4:
		p->vectorfunc[i] = NULL;
		p->matrixfunc[i] = glUniformMatrix3x4fv;
		break;

	    case GL_FLOAT_MAT4x2:
		p->vectorfunc[i] = NULL;
		p->matrixfunc[i] = glUniformMatrix4x2fv;
		break;

	    case GL_FLOAT_MAT4x3:
		p->vectorfunc[i] = NULL;
		p->matrixfunc[i] = glUniformMatrix4x3fv;
		break;
	    }
	}
    }
    p->ntex=ntex;

    // finally, store it in svar
    if(oldp){
	return oldp->id;
    }
    else{
	int id=id1k_alloc(G.vfs_idmngr);
	if(id==-1){ // TODO:
	}

	if(!G.vfs[id]) G.vfs[id] = malloc(sizeof(vfs_t));

	cvar_setp__(nam, G.vfs[id], VFS_CKSUM);
	p->id = id;
	*G.vfs[id] = *p;
	return id;
    }
}


static int vfs_isprog(int id)
{
    return id1k_allocated(G.vfs_idmngr, id);
}

#define BUFFER_OFFSET(i)	((char *)NULL + (i))

/*------------------------------------------
  Use OpenGL shader program
  -----------------------------------------*/
void vfs_use(int id)
{
    if(id<0){
	if(G.vfs_cur){		
	    glUseProgram(0);
	    for(int i=0; i<G.vfs_cur->narray; i++)
		glDisableVertexAttribArray(i);
			
	    G.vfs_cur = NULL;
	}
	return;
    }

    if(!vfs_isprog(id)){ // error
	cs_printf("Shader error\n");
	return;
    }


    vfs_t *p = G.vfs[id];

    if (p==G.vfs_cur) return;

    glUseProgram(p->p);

    // disable unused/enable useful vertex attribs
    vfs_t *old = G.vfs_cur;
    if(!old){
	for(int i=0; i<p->narray; i++)
	    glEnableVertexAttribArray(i);
    }
    else{
	if(old->narray > p->narray) // disable unused attrib arrays
	    for(int i=old->narray-1; i>= p->narray; i--)
		glDisableVertexAttribArray(i);
	else if (old->narray < p->narray)
	    for(int i=old->narray; i<p->narray; i++)
		glEnableVertexAttribArray(i);
    }

    G.vfs_cur = p;	
}



void vfs_reload(const char *nam)
{
    vfs_t *p = cvar_getp__(nam, VFS_CKSUM);
    if(!p) return;
	
    vfs_load(nam);
}



/*--------------------------------------
  GLSL uniforms
  -------------------------------------*/


/*---------------------------------------
  Handy uniform setup routines. Alternatively, you can use unifv() instead.
  --------------------------------------*/
void uni1f(const char *nam, float v0)
{
    if (!G.vfs_cur->p) return;
    int loc = glGetUniformLocation(G.vfs_cur->p, nam);
    glUniform1f(loc, v0);
}

void uni2f(const char *nam, float v0, float v1)
{
    if (!G.vfs_cur->p) return;
    int loc = glGetUniformLocation(G.vfs_cur->p, nam);
    glUniform2f(loc, v0, v1);	
}

void uni3f(const char *nam, float v0, float v1, float v2)
{
    if (!G.vfs_cur->p) return;
    int loc = glGetUniformLocation(G.vfs_cur->p, nam);
    glUniform3f(loc, v0, v1, v2);	
}

void uni4f(const char *nam, float v0, float v1, float v2, float v3)
{
    if (!G.vfs_cur->p) return;
    int loc = glGetUniformLocation(G.vfs_cur->p, nam);
    glUniform4f(loc, v0, v1, v2, v3);	
}


/*------------------------------------------------
  Can process all kind of uniforms -- vector, array of vector, matrix, array of matrix, without specifing data types and sizes.
  All info are parsed at shader loading time.
  -----------------------------------------------*/
void unifv(const char *nam, const float *dat)
{
    if (!dat) return;

    vfs_t *p = G.vfs_cur;
    if (!p) return;

    int i = glGetUniformLocation(p->p, nam);
    if(p->vectorfunc[i]) 
	p->vectorfunc[i](i, p->ucount[i], dat);
    else
	p->matrixfunc[i](i, p->ucount[i], 0, dat);
}


/*-----------------------------------------------
  You don't have to specify which channel the texture is bound to -- the info is parsed at shader loading time.
  ----------------------------------------------*/
void samp2d(const char *nam, uint tex)
{
    vfs_t *p = G.vfs_cur;
    if (!p) return;

    int i = glGetUniformLocation(p->p, nam);

    if(i<0){
	//		cs_printf("Cannot find sampler2D \"%s\"\n", nam);
	return;
    }

    glUniform1i(i, p->multitex[i]);
    glActiveTexture(GL_TEXTURE0 + p->multitex[i]);
    glBindTexture(GL_TEXTURE_2D, tex);
    glActiveTexture(GL_TEXTURE0);
}


/*
  An object can be:

  On disk; -- in-core = 0
  In RAM;  -- in-core = 1; cached = 0
  In VRAM; -- cached = 1


  Object can be modified only in RAM. If it's modified, its instance on disk and in VRAM should be updated.

  Hash table should be used to manage objects.

  Object should be stored in binary format on disk. Some tool should be used to translate between binary and text.

  
*/



/*--------------------------------------------------
  Several vertex buffers are allocated at startup, and then controlled by a binary buddy system.

  Each object has a timestamp. Everytime the object is updated (deformed/translated/rotated etc.), the timestamp is updated too.
  If an object isn't used, its VB is released (put into something like "trash bin"). When the object is used again and the trash bin
  isn't clearred, then the VB can be "recoverred".


  Multi-stream attributes can be stored in the same or different VBs. Just let the engine determine automatically.

  -------------------------------------------------*/


/*-----------------------------------------------------------------------
  Vertex buffer

  The minimum block size is 2k. Total vertex mem is 64m == 32k blocks

  Each vertex buffer object is 4m => at most 16 VBOs.
  Each VBO 2k blocks.

  => address / 2k = VBO index
  (address % 2k) * 2k = VBO offset

  ------------------------------------------------------------

  Element buffer

  EB is managed by another buddy system, too.

  You have to store the number of elements in each EB

  4m element buffer is enough.

  The minimum element buffer can be 128 bytes => 4m == 32k *128

  -----------------------------------------------------------------------*/

#define BLKS_PER_VBO (2*1024)
#define VB_BLK_SIZE (2*1024)
#define EB_BLK_SIZE 128
#define MAX_NUM_BLKS (32*1024)
#define VBO_SIZE (4*1024*1024)
#define EBO_SIZE (4*1024*1024)
static void vb_init()
{
    // vertex buffer
    bd_gen(&G.vb_bdsys, 1, MAX_NUM_BLKS); // the memory "allocated" by the bdsys is useless; only the "address" is useful.
    glGenBuffers(1, G.vb_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, G.vb_vbo[0]);
    glBufferData(GL_ARRAY_BUFFER, VBO_SIZE, NULL, GL_STATIC_DRAW);
    G.vb_nvbos = 1;

    // element buffer
    bd_gen(&G.eb_bdsys, 1, MAX_NUM_BLKS);
    glGenBuffers(1, &G.eb_ebo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, G.eb_ebo);
    glBufferData(GL_ARRAY_BUFFER, EBO_SIZE, NULL, GL_STATIC_DRAW);

}

static void vb_shutdown()
{
    glDeleteBuffers(G.vb_nvbos, G.vb_vbo);
    glDeleteBuffers(1, &G.eb_ebo);

    bd_del(&G.vb_bdsys);
    bd_del(&G.eb_bdsys);
}

/*-------------------------------------------
  Return: some kind of "pointer", can only be used with vb_*() functions.

  -------------------------------------------*/
int vb_alloc(int n, int sz, void *dat)
{
    int realsize;
    int size = n*sz/VB_BLK_SIZE;
    if(n*sz%VB_BLK_SIZE) ++size;

    int addr = bd_alloc(G.vb_bdsys, size, &realsize);
    int vb = addr / BLKS_PER_VBO;

    if(vb>=G.vb_nvbos){
	glGenBuffers(1, G.vb_vbo+vb);
	glBindBuffer(GL_ARRAY_BUFFER, G.vb_vbo[vb]);
	glBufferData(GL_ARRAY_BUFFER, VBO_SIZE, GL_STATIC_DRAW, dat);
	G.vb_nvbos = vb;
    }
    else {
	glBindBuffer(GL_ARRAY_BUFFER, G.vb_vbo[vb]);
	int ofs = addr % BLKS_PER_VBO * VB_BLK_SIZE;
	glBufferSubData(GL_ARRAY_BUFFER, ofs, n*sz, dat);
    }

    //	G.vb_numverts[addr] = n;

    G.vb_cap[addr]=realsize*VB_BLK_SIZE;
    G.vb_nv[addr]=n;

    return addr;
}



/*-------------------------------
  To allocate a element buffer, 


  Element buffers don't have something like glElementArrayPointer, but they source
  their indices from that buffer object, using their <indices>
  parameters as offsets into the buffer object

 
  ------------------------------*/
int eb_alloc(int nelems, int elemsz, void *dat)
{
    int realsize;
    int addr = bd_alloc(G.eb_bdsys, nelems * elemsz, &realsize);
    G.eb_numelements[addr] = nelems;
}


void vb_free(int vb)
{
    bd_free(G.vb_bdsys, vb);
}

void eb_free(int eb)
{
    bd_free(G.eb_bdsys, eb);
}

void vb_update(int addr, int n, int sz, void *dat)
{
    int vb=addr / BLKS_PER_VBO;
    int ofs= addr % BLKS_PER_VBO * VB_BLK_SIZE;
    glBindBuffer(GL_ARRAY_BUFFER, G.vb_vbo[vb]);
    glBufferSubData(GL_ARRAY_BUFFER, ofs, n*sz, dat);
    G.vb_nv[addr]=n;
}

void eb_update(int addr, int nbytes, void *dat)
{
    int ofs = addr * EB_BLK_SIZE;
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, G.eb_ebo);
    glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, ofs, nbytes, dat);
}

int vb_getcap(int addr)
{
    return G.vb_cap[addr];
}


/*----------------------------------------
  An object may have multiple streams. 
  
  Multiple VBOs, each VBO may have multiple interleaved attributes.

  rend("vertex {T B N} texcoord", vert, TBN, texcoord);


  Group interleaved attribs with {}. 
  

  Grammar:

  <vb_desc> -> <attrib> {<vb_desc>} |
  <attrib_group> {<vb_desc>} |
  $


  <attrib> -> identifier
  <attrib_group> -> "{" <attrib_list> "}"
  <attrib_list> -> <attrib_interleaved> {<attrib_list>}
  <attrib_interleaved> -> identifier

  ---------------------------------------*/

static void ps_vb_desc(const char *desc)
{
    char *c=desc;
    char buf[256];

    vfs_t *p = G.vfs_cur;
    if(!p) return;	
	
    memset(p->ofs, 0, p->narray*sizeof(p->ofs[0]));
    memset(p->stride, 0, p->narray*sizeof(p->stride[0]));

    int nstreams=0;
    int *attr = G.vb_streams;

    while(*c){
	if(isspace(*c)){
	    ++c;
	    continue;
	}

	if(*c == '{'){
	    ++c;
	    int ns=0;
	    short *strides[256];
	    int ofs=0;
			
	    G.vb_stream_nattr[nstreams] = 0;
	    while(*c){
		if(isspace(*c)){
		    ++c;
		    continue;
		}

		if(*c=='}') break;

		char *b=buf;
		while(*c=='_' || isalnum(*c)) *b++=*c++;
		*b = '\0';

		int a = glGetAttribLocation(p->p, buf);
		if(a==-1){ // TODO: error processing
		}

		p->ofs[a] = ofs;
		ofs += p->dim[a]*p->acount[a]*sizeof(float);
		strides[ns++] = p->stride + a;

		G.vb_stream_nattr[nstreams]++;
		*attr++ = a;
	    }

	    ++nstreams;

	    if(ns>1){
		for(int i=0; i<ns; i++) 
		    *strides[i] = ofs;
	    }
			

	    ++c; // '}'
	    continue;
	}

	while(*c){
	    if(isspace(*c++)) continue;

	    char *b=buf;
	    while(*c=='_' || isalnum(*c)) *b++=*c++;
	    *b = '\0';

	    int a = glGetAttribLocation(p->p, b);
	    p->ofs[a] = p->stride[0] = 0;
	    G.vb_stream_nattr[nstreams++]++;
	    *attr++ = a;
	}
		
    }

    G.vb_nstreams = nstreams;
	
}



void vb_setnv(int addr, int nv)
{
    //	G.vb_numverts[addr]=nv;
    G.vb_nv[addr]=nv;
}

/*
  "fmt" is vertex attribute description

*/
void rend(const char *fmt, int eb,  ...)
{
    va_list va;

    if((G.vb_curvfs != G.vfs_cur || (G.vb_curfmt != fmt && strcmp(G.vb_curfmt, fmt)) ))
	{
	    G.vb_curvfs = G.vfs_cur;
	    G.vb_curfmt = fmt;

	    ps_vb_desc(fmt);
	}

    vfs_t *p = G.vfs_cur;
    int *attr=G.vb_streams;
    int addr;
    va_start(va, eb);
    for(int s=0; s<G.vb_nstreams; s++){
	addr = va_arg(va, int);
	int vb = addr / BLKS_PER_VBO;
	int ofs = addr % BLKS_PER_VBO * VB_BLK_SIZE;

	int vbo=G.vb_vbo[vb];
	if(vbo != G.vbo_cur){
	    G.vbo_cur=vbo;
	    glBindBuffer(GL_ARRAY_BUFFER, G.vb_vbo[vb]);
	}

	for(int i=0; i<G.vb_stream_nattr[s]; i++)
	    {
		int a = *attr++;
		for(int j=0; j<p->acount[a]; j++){
		    glVertexAttribPointer(a+j, p->dim[a], GL_FLOAT, 0, p->stride[a], ofs+p->ofs[a]);
		}
	    }
    }
    va_end(va);

    int nv = G.vb_nv[addr]; //G.vb_numverts[addr];
    if(eb!=-1){
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, G.eb_ebo);
	int type = nv<=0xff ? GL_UNSIGNED_BYTE : nv<=0xffff ? GL_UNSIGNED_SHORT : GL_UNSIGNED_INT; // TODO
	int ofs = eb * EB_BLK_SIZE;
	glDrawRangeElements(GL_TRIANGLES, 0, nv-1, G.eb_numelements[eb], type, ofs);
    }
    else{ // you can render without element buffer.
	glDrawArrays(GL_TRIANGLES, 0, nv);
    }

    ++G.rnd_counter;
}


void rendq(const char *fmt, int eb,  ...)
{
    va_list va;

    if((G.vb_curvfs != G.vfs_cur || (G.vb_curfmt != fmt && strcmp(G.vb_curfmt, fmt)) ))
	{
	    G.vb_curvfs = G.vfs_cur;
	    G.vb_curfmt = fmt;

	    ps_vb_desc(fmt);
	}

    vfs_t *p = G.vfs_cur;
    int *attr=G.vb_streams;
    int addr;
    va_start(va, eb);
    for(int s=0; s<G.vb_nstreams; s++){
	addr = va_arg(va, int);
	int vb = addr / BLKS_PER_VBO;
	int ofs = addr % BLKS_PER_VBO * VB_BLK_SIZE;
	int vbo=G.vb_vbo[vb];
	if(vbo != G.vbo_cur){
	    G.vbo_cur=vbo;
	    glBindBuffer(GL_ARRAY_BUFFER, vbo);
	}

	for(int i=0; i<G.vb_stream_nattr[s]; i++)
	    {
		int a = *attr++;
		for(int j=0; j<p->acount[a]; j++){
		    glVertexAttribPointer(a+j, p->dim[a], GL_FLOAT, 0, p->stride[a], ofs+p->ofs[a]); // DEBUGGED: p->ofs[a] => ofs+p->ofs[a]
		}
	    }
    }
    va_end(va);

    int nv = G.vb_nv[addr]; //G.vb_numverts[addr];
    if(eb!=-1){
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, G.eb_ebo);
	int type = nv<=0xff ? GL_UNSIGNED_BYTE : nv<=0xffff ? GL_UNSIGNED_SHORT : GL_UNSIGNED_INT; // TODO
	int ofs = eb * EB_BLK_SIZE;
	glDrawRangeElements(GL_QUADS, 0, nv-1, G.eb_numelements[eb], type, ofs);
    }
    else{ // you can render without element buffer.
	glDrawArrays(GL_QUADS, 0, nv);
    }

    ++G.rnd_counter;
}


/*
  Immediate mode

  No multi-streams: just interleaved (four vertices)

  => The vertex attrib description is much simpler:
  just write down the names of the attributes, in the same order as in the vertex struct

  rendrect("co texco", co[4], texco[4])

  co[]={x1,y1,x2,y2}

*/
void rendrect(const char *desc, ...)
{
    va_list va;
    vfs_t *p=G.vfs_cur;

    if (!p) return;

    // parse vertex description
    char *c=desc;
    char buf[128];
    int attr[128];
    int na=0;
    int stride=0;
    float *dat[128];

    va_start(va,desc);
    while(1){
	// skip space
	while(isspace(*c)) ++c;

	// *c could be 0
	if(!*c) break;

	// now *c != 0 and isn't space neither
	char *b=buf;
	do{
	    *b++=*c++;
	}while(isalnum(*c) || *c=='_');
	*b=0;

	int a=glGetAttribLocation(p->p, buf);
	if(a==-1){ // error
	}

	attr[na]=a;
	dat[na++]=va_arg(va,float*);
	stride+=p->dim[a]*p->acount[a]*sizeof(float);
    }
    va_end(va);

    glBindBuffer(GL_ARRAY_BUFFER,0);
    for(int i=0; i<na; i++){
	int a=attr[i];
	int dim=p->dim[a];
	int sz=4*dim*p->acount[a]*sizeof(float);
	float *vert=alloca(sz);
	memset(vert,0,sz);

	float *v=vert;
	float *d=dat[i];
	for(int j=0; j<p->acount[a]; j++){
	    float *p=v;
	    v[0]=d[0];
	    v[1]=d[1];
	    v+=dim;

	    v[0]=d[2];
	    v[1]=d[1];
	    v+=dim;

	    v[0]=d[2];
	    v[1]=d[3];
	    v+=dim;

	    v[0]=d[0];
	    v[1]=d[3];
	    v+=dim;

	    glVertexAttribPointer(a+j, dim, GL_FLOAT, 0, 0, p);

	    d+=4;
	}
    }
	
    // plain vertex array (not VBO)
    glDrawArrays(GL_QUADS, 0, 4);
}



uint rnd_getcounter()
{
    return G.rnd_counter;
}

/*------------------------------------------
  Parse render target description

  Grammar:

  <rtdesc> -> <type>{<rtdesc>}   // no white space between types

  <type> -> t | d | D | S | x | X | y | Y | z | Z

  T|t	TEXTURE_2D
  D	DEPTH_COMPONENT16, texture
  d	DEPTH_COMPONENT24, renderbuffer
  S|s	STENCIL_INDEX, renderbuffer
  x|y|z	TEXTURE_CUBE_MAP_NEGATIVE_X|Y|Z
  X|Y|Z	TEXTURE_CUBE_MAP_POSITIVE_X|Y|Z

  -----------------------------------------*/


static void CHECK_FRAMEBUFFER_STATUS(){
    GLenum status;
    status=(GLenum)glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    switch(status) {
    case GL_FRAMEBUFFER_COMPLETE_EXT:
	//		cs_printf("Fine!\n");
	return;
    case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:
	cs_printf("Framebuffer incomplete,incomplete attachment\n");
	return;
    case GL_FRAMEBUFFER_UNSUPPORTED_EXT:
	cs_printf("Unsupported framebuffer format\n");
	return;
    case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
	cs_printf("Framebuffer incomplete,missing attachment\n");
	return;
    case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
    {
	int tex=-1;
	glGetFramebufferAttachmentParameterivEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, 
						 GL_FRAMEBUFFER_ATTACHMENT_OBJECT_NAME_EXT, &tex);
	cs_printf("Framebuffer incomplete,attached images must have same dimensions\n");
	cs_printf("depth attachment=%i",tex);
    }
    return;
    case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
	cs_printf("Framebuffer incomplete,attached images must have same format\n");
	return;
    case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
	cs_printf("Framebuffer incomplete,missing draw buffer\n");
	return;
    case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
	cs_printf("Framebuffer incomplete,missing read buffer\n");
	return;
    }
    return;
}




/*------------------------------------------
  Use render target

  e.g.
  rnd_tg("ttds", nam1, nam2, nam3, nam4);

  rnd_tg(0); // Render to frame buffer

  rnd_tg("td", color, depth);
  -----------------------------------------*/
int frm=0;

void rnd_tg(const char *format, ...)
{
    va_list	ap;
	
    if(!format) {	// "Real" frame buffer
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
	return;
    }

    if(!G.rndtg_fbo){
	glGenFramebuffersEXT(1, &G.rndtg_fbo);
    }


#if 0
    // TESTING
    glBindTexture(GL_TEXTURE_2D, 0);
    // ~TESTING
#endif

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, G.rndtg_fbo);

    int n=strlen(format);
    int clrattach=GL_COLOR_ATTACHMENT0_EXT;
    int depthattach=0;
    va_start(ap, format);

    for(int i=0; i<n; i++){
	int rndtg = va_arg(ap, int);

	switch(format[i]){
	case 't':
	case 'T':
	    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
				      clrattach++, GL_TEXTURE_2D, rndtg, 0);
	    assert(!glGetError()+0);
	    break;
	case 'd':
	case 'D':
	    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
				      GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, rndtg, 0); // DEBUGGED: GL_DEPTH_COMPONENT=>GL_TEXTURE_2D

	    assert(!glGetError()+00);
	    depthattach=1;
	    break;
	case 's':
	case 'S':
	    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
				      GL_STENCIL_ATTACHMENT_EXT,GL_TEXTURE_2D, rndtg, 0);
	    break;
	default:
	    assert(0);
	    break;
	}
    }
    va_end(ap);



    // MRT
    int nrt=clrattach-GL_COLOR_ATTACHMENT0_EXT;

#if 1
    if (G.rndtg_nca>nrt) {
	for (int i=nrt; i<G.rndtg_nca; ++i ){
	    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
				      GL_COLOR_ATTACHMENT0_EXT+i, GL_TEXTURE_2D, 0,0);
	}
    }
    G.rndtg_nca=nrt;

    if (!depthattach){
	if (G.rndtg_nda) {
	    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
				      GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, 0, 0); // DEBUGGED: GL_DEPTH_COMPONENT=>GL_TEXTURE_2D
	    G.rndtg_nda=0;
	}
    }
    else {
	G.rndtg_nda=1;
    }

#endif

    if(nrt==1){
	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
	glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
    }
    else if(nrt>1){
	uint *mrt=alloca(sizeof(uint)*nrt);
	for(int i=0; i<nrt; i++){
	    mrt[i]=GL_COLOR_ATTACHMENT0_EXT+i;
	}
	glDrawBuffers(nrt, mrt);
    }
    else{	// No color attachment, only render to depth buffer
	if(format[0]=='d' || format[0]=='D'){
	    glDrawBuffer(GL_NONE);
	    glReadBuffer(GL_NONE);
	}
    }

    CHECK_FRAMEBUFFER_STATUS();
}


static void rnd_tgn(const char *format, uint *rndtg)
{
    if(!format) {	// "Real" frame buffer
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
	return;
    }

    if(!G.rndtg_fbo){
	glGenFramebuffersEXT(1, &G.rndtg_fbo);
    }

#if 0
    // TESTING
    glBindTexture(GL_TEXTURE_2D, 0);
    // ~TESTING
#endif

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, G.rndtg_fbo);

    int n=strlen(format);
    int clrattach=GL_COLOR_ATTACHMENT0_EXT;
    int depthattach=0;

    for(int i=0; i<n; i++){
	switch(format[i]){
	case 't':
	case 'T':
	    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
				      clrattach++, GL_TEXTURE_2D, rndtg[i], 0);
	    assert(!glGetError()+1*0);
	    break;
	case 'd':
	case 'D':
	    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
				      GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, rndtg[i], 0); // DEBUGGED: GL_DEPTH_COMPONENT=>GL_TEXTURE_2D
	    depthattach=1;
	    assert(!glGetError()+11*0);
	    break;
	case 's':
	case 'S':
	    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
				      GL_STENCIL_ATTACHMENT_EXT,GL_TEXTURE_2D, rndtg[i], 0);
	    assert(!glGetError()+111*0);
	    break;
	default:
	    assert(0);
	    break;
	}
    }

    // MRT
    int nrt=clrattach-GL_COLOR_ATTACHMENT0_EXT;

#if 1
    if (G.rndtg_nca>nrt) {
	for (int i=nrt; i<G.rndtg_nca; ++i ){
	    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
				      GL_COLOR_ATTACHMENT0_EXT+i, GL_TEXTURE_2D, 0,0);
	}
    }
    G.rndtg_nca=nrt;

    if (!depthattach){
	if (G.rndtg_nda) {
	    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
				      GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, 0, 0); // DEBUGGED: GL_DEPTH_COMPONENT=>GL_TEXTURE_2D
	    G.rndtg_nda=0;
	}
    }
    else {
	G.rndtg_nda=1;
    }
#endif

    if(nrt==1){
	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
	glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
    }
    else if(nrt>1){
	uint *mrt=alloca(sizeof(uint)*nrt);
	for(int i=0; i<nrt; i++){
	    mrt[i]=GL_COLOR_ATTACHMENT0_EXT+i;
	}
	glDrawBuffers(nrt, mrt);
    }
    else{	// No color attachment, only render to depth buffer
	if(format[0]=='d' || format[0]=='D'){
	    glDrawBuffer(GL_NONE);
	    glReadBuffer(GL_NONE);
	}
    }

    CHECK_FRAMEBUFFER_STATUS();
}




/*------------------------------------
  Common texture filters:
  mipminmag[0][1][1]

  Usage:
  uint tex=tex3b(512,512,000,data);
  -----------------------------------*/
static void fltrtyp_nomip(int f, int *minf, int *magf)
{
    static int init=0;
    static char m[212]={[000]=0, [001]=1, [010]=2, [011]=3, 
			[100]=0, [101]=1, [110]=2, [111]=3,
			[200]=0, [201]=1, [210]=2, [211]=3}; // mipminmag[]

    static int minfltrtb[6]={GL_NEAREST,
			     GL_LINEAR,
			     GL_NEAREST_MIPMAP_NEAREST,
			     GL_LINEAR_MIPMAP_NEAREST,
			     GL_NEAREST_MIPMAP_LINEAR,
			     GL_LINEAR_MIPMAP_LINEAR	};
    static int magfltrtb[2]={GL_NEAREST, GL_LINEAR};

    int typ=m[f % 212]; // So that m[] won't overflow. And all unknown filter types are treated 
    *minf=minfltrtb[typ/2];
    *magf=magfltrtb[typ%2];
}

static void fltrtyp(int f, int *minf, int *magf)
{
    static int init=0;
    static char m[212]={[000]=0, [001]=1, [010]=2, [011]=3, 
			[100]=4, [101]=5, [110]=6, [111]=7,
			[200]=8, [201]=9, [210]=10, [211]=11}; // mipminmag[]

    static int minfltrtb[6]={GL_NEAREST,
			     GL_LINEAR,
			     GL_NEAREST_MIPMAP_NEAREST,
			     GL_LINEAR_MIPMAP_NEAREST,
			     GL_NEAREST_MIPMAP_LINEAR,
			     GL_LINEAR_MIPMAP_LINEAR	};
    static int magfltrtb[2]={GL_NEAREST, GL_LINEAR};
    int typ=m[f % 212];
    *minf=minfltrtb[typ/2];
    *magf=magfltrtb[typ%2];
}

uint tex1b(int w, int h, int filter, void *data)
{
    uint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);

#if 0
    glTexImage2D(GL_TEXTURE_2D, 0/*level*/, GL_LUMINANCE8_ALPHA8, w, h, 0/*border*/,
		 GL_LUMINANCE_ALPHA,GL_UNSIGNED_BYTE, data);
#endif

    glTexImage2D(GL_TEXTURE_2D, 0/*level*/, GL_ALPHA8, w, h, 0/*border*/,
		 GL_ALPHA, GL_UNSIGNED_BYTE, data);


    int minf,magf;
    fltrtyp(filter, &minf, &magf);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minf);	
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magf);
    return tex;
}

uint tex3b(int w, int h, int filter, void *data)
{
    uint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0/*level*/, GL_RGB8, w, h, 0/*border*/,
		 GL_RGB,	GL_UNSIGNED_BYTE, data);

    int minf,magf;
    fltrtyp(filter, &minf, &magf);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minf);	
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magf);
    return tex;
}


uint tex4b(int w, int h, int filter, void *data)
{
    uint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGB,
		 GL_UNSIGNED_BYTE, data);

    int minf,magf;
    fltrtyp(filter, &minf, &magf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magf);
    return tex;
}


uint tex1h(int w, int h, int filter, void *data)
{
    uint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA16F_ARB, w, h, 0, GL_ALPHA,
		 GL_FLOAT, data);
	
    int minf, magf;
    fltrtyp_nomip(filter, &minf, &magf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magf);
    return tex;
}


uint tex3h(int w, int h, int filter, void *data)
{
    uint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F_ARB, w,h,0, GL_RGB, 
		 GL_FLOAT, data);

    int minf, magf;
    fltrtyp_nomip(filter, &minf, &magf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magf);
    return tex;	
}

uint tex4h(int w, int h, int filter, void *data)
{
    uint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F_ARB, w,h,0, GL_RGBA, 
		 GL_FLOAT, data);

    int minf, magf;
    fltrtyp_nomip(filter, &minf, &magf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magf);
    return tex;	
}

uint tex1f(int w, int h, int filter, void *data)
{
    uint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA32F_ARB, w, h, 0, GL_ALPHA,
		 GL_FLOAT, data);
	
    int minf, magf;
    fltrtyp_nomip(filter, &minf, &magf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magf);
    return tex;
}

uint tex3f(int w, int h, int filter, void *data)
{
    uint tex;
    glGenTextures(1, &tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F_ARB, w,h,0, GL_RGB, 
		 GL_FLOAT, data);

    int minf, magf;
    fltrtyp_nomip(filter, &minf, &magf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magf);
    return tex;	
}


uint tex3bdxt1(int w, int h, int filter, void *data)
{
    uint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);

    if (data) {
	int sz=8*(w/4)*(h/4);
	glCompressedTexImage2D(GL_TEXTURE_2D, 0,
			       GL_COMPRESSED_RGB_S3TC_DXT1_EXT,
			       w, h, 0, sz, data);
    }
    else {  // DEBUGGED: should use glTexImage2D() to reserve space for texture
	glTexImage2D(GL_TEXTURE_2D, 0, 
		     GL_COMPRESSED_RGB_S3TC_DXT1_EXT, 
		     w,h,0,GL_RGB,GL_UNSIGNED_BYTE, NULL);
    }

    int minf, magf;
    fltrtyp_nomip(filter, &minf, &magf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magf);
    return tex;
}


uint tex4f(int w, int h, int filter, void *data)
{
    uint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, w,h,0, GL_RGBA, 
		 GL_FLOAT, data);

    int minf, magf;
    fltrtyp_nomip(filter, &minf, &magf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magf);
    return tex;	
}

uint dtex(int w, int h, int filter, void *data)
{
    uint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, w,h,0, 
		 GL_DEPTH_COMPONENT, GL_UNSIGNED_SHORT, data);

    int minf, magf;
    fltrtyp_nomip(filter, &minf, &magf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magf);
    return tex;	
}

uint dtex2b(int w, int h, int filter, void *data)
{
	
    uint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT16, w,h,0, 
		 GL_DEPTH_COMPONENT, GL_UNSIGNED_SHORT, data);

    int minf, magf;
    fltrtyp_nomip(filter, &minf, &magf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magf);
    return tex;	
}


#if 0
uint dtex3b(int w, int h, int filter, void *data)
{
    uint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, w,h,0,
		 GL_DEPTH_COMPONENT, GL_UNSIGNED_SHORT, data);

    int minf, magf;
    fltrtyp_nomip(filter, &minf, &magf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magf);
    return tex;	
}
#endif

uint dtex3b(int w, int h, int filter, void *data)
{
    uint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, w,h,0,
		 GL_DEPTH_COMPONENT, GL_UNSIGNED_SHORT, data);
    //		     GL_LUMINANCE, GL_FLOAT, data);

    int minf, magf;
    fltrtyp_nomip(filter, &minf, &magf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magf);

#if 0
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
    glTexParameteri(GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE, GL_LUMINANCE);
#endif

    return tex;	
}


uint dstex(int w, int h, int filter, void *data)
{
    uint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH24_STENCIL8_EXT, w,h,0, 
		 GL_DEPTH_STENCIL_EXT, GL_UNSIGNED_INT_24_8_EXT, data);

    int minf, magf;
    fltrtyp_nomip(filter, &minf, &magf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minf);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magf);
    return tex;	
}




const uint (*texfunc[RT_MAX])(int,int,int,void*)={[RT_4B]=tex4b, [RT_4H]=tex4h, [RT_DEPTH]=dtex};

void rt_init()
{
    int nlevels=lg2(min(G.screenw,G.screenh))+1;

    for (int level=0; level<nlevels; ++level) {
	int w=G.screenw>>level;
	int h=G.screenh>>level;
	const int nrts[RT_MAX]={[RT_4B]=4, [RT_4H]=2, [RT_DEPTH]=2};

	for (int type=0; type<RT_MAX; ++type) {
	    int bit=1;
	    for (int i=0; i<nrts[type]; ++i) {
		// texture
		uint t=texfunc[type](w,h,011,0);
		arr_push(G.rt_tex[type][level], t);
		//				arr_pushn(G.rt_tex[type][level], &t, 1);

		// bitmap
		G.rt_bm[type][level] |= bit<<i;

		// info
		int n=arr_len(G.rt_info); // rt_info[tex_name]
		if (t>=n) {
		    rtinfo_t info[t-n+1];
		    for (int k=0; k<t-n+1; ++k) {
			info[k].type=-1;
		    }
					
		    arr_pushn(G.rt_info, info, t-n+1);
		}
		rtinfo_t info={.type=type, .level=level, .index=i};

		if (type==0 && level==0) {
		    printf("texture=%u\n", t);
		}

		assert(t<1000); /////////////////////////// dddddddddddddddddddddddd
		G.rt_info[t]=info;
	    }
	}
    }
}

void rt_shutdown()
{
}

/*
  gets garbage collected every frame
*/
uint rt_alloca(int type, int level)
{
}


/*
  
 */
void rt_lock(uint tex)
{
}

void rt_unlock(uint tex)
{
}


uint rt_alloc(int type, int level)
{
    type %= RT_MAX;

    if (G.rt_bm[type][level]) {
	int n=arr_len(G.rt_tex[type][level]);
	for (int i=0; i<n; i++) {
	    if (G.rt_bm[type][level] & (1<<i) ) {
		G.rt_bm[type][level] &= (~(1<<i));
		uint t=G.rt_tex[type][level][i];
		rtinfo_t info={.type=type, .level=level, .index=i};

		if (t>1000) {
		    printf("type=%i, level=%i\n", type, level);
		    assert(t<1000); ////////////////////// dddddddddddddddddddd
		}
				
		G.rt_info[t]=info;
		return t;
	    }
	}			
    }
    else {
	int sz=arr_len(G.rt_tex[type][level]);
	if (sz>10) return 0; // TODO:

	uint t=texfunc[type](G.screenw>>level, G.screenh>>level, 011, 0);
	arr_push(G.rt_tex[type][level], t);
	//		arr_pushn(G.rt_tex[type][level], &t, 1);
	if (t>=sz) {
	    rtinfo_t info[t-sz+1];
	    for (int k=0; k<t-sz+1; ++k) {
		info[k].type=-1;
	    }
			
	    arr_pushn(G.rt_info, info, t-sz+1);
	}

	rtinfo_t info={.type=type, .level=level, .index=arr_len(G.rt_tex[type][level])-1};		
	G.rt_info[t]=info;
	return t;
    }
}

static int rt_valid(uint tex)
{
    if (tex>=arr_len(G.rt_info)) return 0;			// access violation

    rtinfo_t *info=G.rt_info+tex;			// not a (known) RT texture
    if (info->type==-1) return 0;

    return 1;
}

int rt_fmt(uint tex)
{
    if (rt_valid(tex)) {
	return G.rt_info[tex].type;
    }

    return -1;
}

int rt_getinfo(uint tex, rtinfo_t *info)
{
    if (rt_valid(tex)) {
	*info=G.rt_info[tex];
	return 0;
    }

    return 1;
}

void rt_free(uint tex)
{
    if (rt_valid(tex)) {
	rtinfo_t *info=G.rt_info+tex;
	G.rt_bm[info->type][info->level] |= (1<<info->index);
    }
}

