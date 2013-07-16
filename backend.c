
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>
#include <alloca.h>
#include <dlfcn.h>

#include <sys/stat.h>
#include <unistd.h>
#include <sys/time.h>

#include <math.h>

#define GL_GLEXT_PROTOTYPES

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>

#include <SDL.h>

#include "backend.h"

typedef unsigned int	uint;
typedef unsigned short	ushort;
typedef unsigned char	uchar;
typedef uchar		byte;
typedef ushort		word;
typedef uint		dword;





typedef struct cvar_s {
    uint cksum;
    struct cvar_s *next;
    struct {int i; float f;} value;
    char isnum;
    short datbgn;
    char mem[256]; // name start from mem[0];  data start from mem[datbgn]
} cvar_t;


#define CVAR_BUCKET_SIZE 997


typedef struct fntch_s{
    short x,y,w,h;
    short advance;
    float u, v; // texcoord
    float texw, texh;
} fntch_t;

typedef struct fnt_s{
    fntch_t ch[128];
    char nam[96];
    short size;
    uint tex;
    char occupied;
    short height, ascender;
    float avgw; // average char width
} fnt_t;

typedef struct {
    float x,y,u,v; // 4x4=16   ---\____ total 20 bytes
    //	char r,g,b,a; // 4         ---/
    float rgba[4];
} fnt_vbfmt_t; 

typedef struct mblist_s {
    struct mblist_s *prev, *next;
    char mb[256];
} mblist_t;

typedef struct guibtn_s{
    struct guibtn_s *prev, *next;
    char *nam;
    char *txt;
    char *cmd;
    int x,y,w,h;
    float fg[4];
    float bg[4];
    short txtnol;
    float txtscale;
} guibtn_t;

typedef struct vu_s{
    struct vu_s *prev, *next;
    char *kcmd[SDLK_LAST];
    char *cmd;
    char *nam;
    guibtn_t *btnlist;

    int g; // grid id, 1 to 9
    int x,y,w,h;
    //	int center[2];
    float fg[4], bg[4];

    int btnid;
} rndvu_t;

typedef struct {
    //	float x,y,z,u,v; // 4x5=20   ---\____ total 24 bytes
    //	char r,g,b,a; // 4         ---/

    float x,y,z,u,v;
    float rgba[4];
} btn_vbfmt_t;

typedef struct rtinfo_s{
    char type, level;
    int index;
    int o; // is output?
} rtinfo_t;


typedef struct kw_s{
    char	*str;
    int	tok;
} kw_t;

typedef struct cveclut_s{
    char *nam;
    float *v;
} cveclut_t;

typedef struct lex_s{
    char	*buf;
    char	*p;
    kw_t	*kw;
    int nkw;
    int	tok;
    union {
	int i;
	float f;
	char str[1024];	/* the content of sem.str and sem.id will be overwritten after next gettok ()*/
	char id[1024];
    } sem;
} lex_t;



enum {SHADER_DS=1, SHADER_PACK=2};
#define FILTER_CKSUM 0x8da23798
#define FILTER_INPUT 0
#define FILTER_OUTPUT 1


typedef struct {
    int type;
    short indeg, outdeg;
    short inout[64];
    int valid;

    struct {
	int type;
	int level;
	uint tex;
    } buffer;

    struct {
	struct {
	    char nam[64];
	    uint u; // uniform
	    int b; // buffer idx -- -1 if input texture is absent
	} *in;

	struct {
	    int b;
	} *out;

	uint p;
	uint p_dscolor;
	uint p_dsdepth;
	short ni, no;
	char colormask[4];
	char *nam;
	short type;
	int level; // texture level of its RTs
    } shader;


    struct {
	float x,y,w,h;
	int rank;
    } disp; // how to render the node in ppedit mode.

} filternode_t;


typedef struct {
    char *nam;
    int id;
    float *tex;
} filterio_t;


typedef struct {
    char *nam;
    filterio_t *in;
    filterio_t *out;
    filternode_t *nodes;
    short *validshaders;
    int ns;
    int hasfreenode;
    time_t time;
    struct { // for filter editor
	float x,y;
	float zoom;
	int maxrank;
	float mousex, mousey;
    } disp;
} filter_t;


typedef struct {
    char *nam;
    char *text;
    char *cmd;
    char *submenu;
    struct menu_s *link;
} menuitem_t;

typedef struct menu_s{
    char *nam;
    menuitem_t *items;
    struct {
	int x,y,w,h;
	int rhs; // pop submenus in right hand side?
    } disp;

    char *mem;
} menu_t;

typedef struct {
    menu_t *menus;
    int *disps;
    float textscale;
    int x2;
    int itemh;
    int curmenu;
    int curitem;
} menuenv_t;

typedef struct{
    v3_t  pos,oldpos;
    v3_t  axis;
    float angle;
    float u, v;			/* coordinates in a spherical coord system. horizental:u, vertical:v */
    v3_t  rotaxis[3];
    float fov;
    float znear;
    float proj[16];
    float rotmat[16];
    struct {
	float u, v, a;
    } rotmin, rotmax;
    struct {
	int u, v, a;		/* if clamprot.xx == 0, then user provided rotmin and rotmax are ignored and set to 0 and 360. */
    } clamprot;
    int updated;
    float w, h;
    float acc[3];
    float angacc[3];
    float vel[3];
    float angvel[3];
    float friction,angfriction,maxvel,maxangvel;
    short sid;
} cam1_t;

typedef struct {
    float x,y;
    float zoom;
} cam2_t;

typedef struct{
    v3_t pos;
    v3_t velocity;
    v3_t target;
    v3_t spherical;
    float d0;
    float Ks, Kd;
    float Kslarge, Kdlarge;
    float Kssmall, Kdsmall;
    float angvelocity;
    float angKd;
    v3_t view;
    v3_t strafe;
    v3_t spherical_prediction;
    double modelview[16];
    double proj[16];
    int viewport[4];
} cam3_t;


typedef struct {
    short sid, eid;
    int flag;
    float t1,t2;
} pvs2d_wall_t;

typedef struct {
    short sid;
    short nv;
    float t[5],y;
    short eid[5];
} pvs2d_flor_t;

typedef struct {
    pvs2d_wall_t *walls;
    pvs2d_flor_t *flors;
    short *nwalls;
    struct {short sid,fid;} *flor_sectlut;
    int wall_ns,nf;
} pvs2d_t;

typedef struct {
    short id;
    short type;
} rayinfo_t;

typedef struct {
    cam2_t cam;
    char *cmd;
    char *kcmd[SDLK_LAST];
    int level;
    float mousex,mousey;
    int cur_sect;
    cam1_t cam1;
    int dim;
    pvs2d_t pvs;
    rayinfo_t rayinfo;
} editor_t;

typedef struct {
    int x,y;
} nook_t;

#define MAXNWALLS 0x7fff
#define MAXNSECTS 0X3fff
#define MAXNLIGHTS 0x3fff

#define IS_SOLID(flag) ((flag)>=0x7fff0000)

typedef struct {
    short adjw, adjs, nextco; // nextwall != nextco if the sector has hole
    int x,y,floorz,ceilingz;
    // A wall can be split into multiple parts if the object is self-intersecting.
    // All these sub-walls have the same adjacent sector as the original wall.

    int type; // wall types: solid, shared, invalid, 
    uint tex;
//    commonwallvert_t *wallvert;
//    sharedwall_t *sharedwall;

    short *elist[2]; // each wall has two tess triangle lists--the triangles that's above and blow it
    short ne[2];
    int reverse;

    // texture info

} wall_t;

typedef struct {
    short firstw, nw;
    int bnd[4];
    float *tess;
    int *e; // edge data: 
    // 1) sold wall--0x7fff0000 | wall_id (>=0x7fff0000)
    // 2) internal--id to adj tess edge (>=0 && <0x7fff)
    // 3) shared--id to the shared wall +1 (short), and larger index (>0) to the edge list entries on the shared wall (short).
    //    to avoid being confused with 2)
    
    int nt;
    int floorz,ceilingz;
    short floorgrad, ceilinggrad;
    short hinge;

    // plane equations
    float floorpln[4], ceilingpln[4];

    // accelerate ray casting

    uint ts; // time stamp

    short visited;
    short hittype; // 0 -- not hit by the ray; 1--floor; 2--ceiling
    float hitt;

    // texture info
    
    short *lightids;
} sect_t;


static uint sect_ts=1;
#define sect_tick() (++sect_ts)
#define sect_copytick(s) ((s)->ts = sect_ts)
#define sect_matchtick(s) ((s)->ts == sect_ts)

typedef struct {
    short sid;
    short e1,e2;
    float t1,t2;
} rayseg_t;

// Each tess triangle stores begin/end points
typedef struct {
    rayseg_t *segments;
    float dst[2];
    int ns;
} rayintx_t;

typedef struct {
    float	r,g,b,a;
} color_t;

typedef struct {
    v3_t	pos;
    pvs2d_t	pvs;
    color_t	color;
    float	radius;
} light_t;

typedef struct {
} model_t;

struct G_s{
    idmngr1k_t *id1ks;
    pool_t *pools;
    bdsys_t *bds;

    SDL_Surface *surface;

    vfs_t *vfs_list;
    vfs_t *vfs_cur;
    int vfs_idmngr;
    vfs_t *vfs[1024];

    char *vb_curfmt;
    vfs_t *vb_curvfs; 	// It's possible that different VFS can share the same VB format.
    int vb_nstreams;
    int vb_stream_nattr[256];
    int vb_streams[256];
    int vbo_cur;

    int cvar_pool;
    cvar_t *cvar_buckets[CVAR_BUCKET_SIZE];
    //	cvar_t *var_buckets[CVAR_BUCKET_SIZE];
	
    int vb_bdsys;
    int eb_bdsys;
    uint vb_vbo[16];
    int vb_nvbos;
	
    int vb_nv[32*1024];
    int vb_cap[32*1024];


    uint eb_ebo;
    uint eb_numelements[32*1024];

    int rndtg_fbo;
    int rndtg_nca; // num of color attachments
    int rndtg_nda; // num of depth attachment (0 or 1)

    /*
      uint *rt_4b[12];
      uint *rt_4h[12];
      uint *rt_depth[12];
      uint rt_4bbm[12]; // bitmap to record the render target usages
      uint rt_4hbm[12];
      uint rt_depthbm[12];
    */
    uint *rt_tex[RT_MAX][12];
    uint rt_bm[RT_MAX][12];
    rtinfo_t *rt_info;

    fnt_t fnt_slots[8];
    fnt_t *fnt_cur;
    fnt_vbfmt_t *fnt_vb;
    uint fnt_tex;
    uint fnt_shader;
    //	uchar fnt_r, fnt_g, fnt_b, fnt_a;
    float fnt_rgba[4];

    int fnt_align;
    float fnt_scale;
    uint fnt_vbo;
    int fnt_vbocap;

    char *cmd_buf;
    char *cmd_cur;

    short cs_enabled;
    char **cs_lines;
    int cs_w, cs_h;
    int cs_outh;
    int cs_out;
    int cs_in;
    int cs_nol;
	 
    mblist_t *cs_mbmem, *cs_mbhead, *cs_mbcur;
    int cs_nmb;
    int cs_scroll;
    int cs_mbpos;
    int cs_csrblink:6; // [-31,31]

    float cs_fg[4];
    float cs_bg[4];
    float cs_mbfg[4];
    float cs_mbbg[4];


    int cs_winh; // height of the console window

    rndvu_t *vu_cur;
    rndvu_t *vu_edit; // the current view selected by command, used by vu_bgn() vu_end() kind of functions
    rndvu_t *vulist;
    guibtn_t *btn_edit;
    btn_vbfmt_t *btn_vb;
    int btn_shader;
    int vu_id;
	

    short guied_LMBdown;
    short guied9;
    short guied_op;
    guibtn_t *guied_btn; // the button currently under the mouse cursor
    rndvu_t *guied_vu; // not the same as G.vu_edit.
    rndvu_t *guied_vuall[9];
    int guied_vucenter[2];
    short guied_changed;
    uint guied_tntex; // thumb nail texture
    short guied9_dragpnt[2]; // relative position inside the selected view
    float guied_statusbarfg[4], guied_statusbarbg[4];
    char *guied_kcmd[SDLK_LAST];	

    int screenw, screenh;
    int mt;
    int mt_counter;
    float mt_tween, mt_dt;
    rndvu_t *mt_vuvu[2];
    uint mt_vuvu_shader;

    uint tex_screensz, tex_screensz1;
    uint dtex_screensz;
    uint tex_128x128;
    uint dtex_128x128;
	
    char keydown[SDLK_LAST]; // DEBUGGED: keydown[256] => keydown[SDLK_LAST] -- SDLK_* may be greater than 255!
    int mousex, mousey;

    int mode, prevmode;
    time_t timebase, timecur;
    int pause;
    int updthumb;
    short vu_tnsz;

    uint rnd_counter;

    lex_t *lex_cur;

    char *filter_curnam;
    int filter_curnodeid;
    int filter_curchannel;
    filter_t *filter_cur;
    int filter_ctxmenu_on;
    menuenv_t filter_menuenv;

    float cur_z;

    void *dll;
    cveclut_t *cvec_lut;

    editor_t ed;

    wall_t walls[MAXNWALLS];
    sect_t sects[MAXNSECTS];
    short nwalls, nsects;
    wall_t *ed_wallstk;

    int nktype;
    nook_t *nkstk;

    rayintx_t lray, rray;
    pvs2d_t tracepvs;

    light_t lights[MAXNLIGHTS];
    short nlights;
} G;


void mymemcpy(void *dest, void *src, int size)
{
    char *d=dest;
    char *s=src;
    int n=0;
    while(n<=size){
	*d++ = *s++;
	n++;
    }
}

#define VFS_CKSUM 0x7eb608cf // cksum shader


/*------------------------------------
  Just copy the value of *sv to the svar named "nam" (if not exist add it)
  ------------------------------------*/
static void cvar_init()
{
    pool_gen(&G.cvar_pool, 128, sizeof(cvar_t));
}

static void cvar_shutdown()
{
    pool_del(&G.cvar_pool);
}


static uint cvar_hash(char *s)
{
    uint h=5381;
    char c;

    while(c = *s++){
	h = (h<<5)+h+c;
    }

    return h;
}


#define CVAR_ITER	int bkt; for(cvar_t *var=G.cvar_buckets[bkt=cvar_hash(nam) % CVAR_BUCKET_SIZE]; var; var=var->next)

void cvar_set(const char *nam, const char *str)
{
    CVAR_ITER{
	if(!strcmp(var->mem, nam)){
	    var->isnum = 0;
	    strcpy(var->mem + var->datbgn, str);
	    return;
	}
    }

    cvar_t *var = pool_alloc(G.cvar_pool);
    var->next = G.cvar_buckets[bkt];
    G.cvar_buckets[bkt] = var;
    var->isnum = 0;
	
    char *s=nam;
    char *d=var->mem;
    while(*s){
	*d++ = *s++;
    }

    *d++=0;

    s=str;
    var->datbgn = d - var->mem;
    while(*s){
	*d++ = *s++;
    }
	
}

char *cvar_get(const char *nam, int *ok)
{
    CVAR_ITER{
	if(!strcmp(var->mem, nam)){
	    if(ok) *ok = 1;
	    return var->mem + var->datbgn;
	}
    }

    if(ok) *ok = 0;
    return NULL;
}



void cvar_seti(const char *nam, int i)
{
    CVAR_ITER{
	if(!strcmp(var->mem, nam)){
	    var->isnum = 1;
	    var->value.i = i;
	    var->value.f = i;
	    sprintf(var->mem + var->datbgn, "%i", i);
	    return;
	}
    }

    cvar_t *var = pool_alloc(G.cvar_pool);
    var->value.i = i;
    var->value.f = i;
    var->isnum = 1;

    char *src=nam, *dest=var->mem;
    while(*src){
	*dest++ = *src++;
    }

    *dest++ = 0;
    var->datbgn = dest-var->mem;

    sprintf(var->mem + var->datbgn, "%i", i);

    var->next = G.cvar_buckets[bkt];
    G.cvar_buckets[bkt] = var;
}


int cvar_geti(const char *nam, int *ok)
{
    CVAR_ITER{
	if(!strcmp(var->mem, nam)){
	    if(var->isnum){
		if(ok) *ok = 1;
		return var->value.i;
	    }
			
	    break;
	}
    }

    if(ok) *ok = 0;
    return 0;
}


static void cvar_setp__(const char *nam, void *ptr, int cksum)
{
    CVAR_ITER{
	if(!strcmp(var->mem, nam)){
	    var->isnum = 0;
	    var->value.i = ptr;
	    var->cksum = cksum;
	    var->mem[var->datbgn] = 0;
	}
    }

    cvar_t *var = pool_alloc(G.cvar_pool);
    var->isnum = 0;
    var->value.i = ptr;
    var->cksum = cksum;

    char *src=nam, *dest=var->mem;
    while(*src){
	*dest++ = *src++;
    }

    *dest++ = 0;
    var->datbgn = dest-var->mem;
    var->mem[var->datbgn] = 0;

    var->next = G.cvar_buckets[bkt];
    G.cvar_buckets[bkt] = var;
}


void *cvar_getp__(const char *nam, int cksum)
{
    int b=cvar_hash(nam) % CVAR_BUCKET_SIZE;
    for (cvar_t *v=G.cvar_buckets[b]; v; v=v->next) {
	if(!strcmp(v->mem, nam)){
	    if(v->cksum == cksum){
		return v->value.i;
	    }
	    break;
	}
    }
    return NULL;	
}



void cvar_setf(const char *nam, float f)
{
    CVAR_ITER{
	if(!strcmp(var->mem, nam)){
	    var->isnum = 1;
	    var->value.i = f;
	    var->value.f = f;
	    sprintf(var->mem + var->datbgn, "%f", f);
	    return;
	}
    }

    cvar_t *var = pool_alloc(G.cvar_pool);
    var->value.i = f;
    var->value.f = f;
    var->isnum = 1;

    char *src=nam, *dest=var->mem;
    while(*src){
	*dest++ = *src++;
    }

    *dest++ = 0;
    var->datbgn = dest-var->mem;

    sprintf(var->mem + var->datbgn, "%f", f);

    var->next = G.cvar_buckets[bkt];
    G.cvar_buckets[bkt] = var;
	
}

float cvar_getf(const char *nam, int *ok)
{
    CVAR_ITER{
	if(!strcmp(var->mem, nam)){
	    if(var->isnum){
		if(ok) *ok = 1;
		return var->value.f;
	    }
	    break;
	}
    }
	
    if(ok) *ok = 0;
    return 0;
}


void cvar_del(char *nam)
{
    int bkt = cvar_hash(nam) % CVAR_BUCKET_SIZE;
    cvar_t *prev = NULL;
    for(cvar_t *var=G.cvar_buckets[bkt]; var; var=var->next){
	if(!strcmp(var->mem, nam)){
	    if(prev) prev->next = var->next;
	    else G.cvar_buckets[bkt] = var->next;
	    pool_free(G.cvar_pool, var);
	    return;
	}
    }	
}


/*
  cvec -- console floating piont vectors

*/

typedef struct cvec_padding_s{
    char *nam;
    float v[16]; // 4x16+8=72 + 24 = 96
    int id;
    float *padding;
} cvec_padding_t;

static void cvec_init()
{
    char *c=loadtext_alloca("cvec.cfg");
    if (!c) return;

    lex_t lex;

    lex_bgn(&lex, c, 0);
    lex_gettok();
    while (lex.tok!=TOK_END) {
	if (lex.tok!=TOK_ID) { // ignore all non-identifier tokens
	    lex_gettok();
	    continue;
	}

	char nam[128]="cvec_init_";
	strcat(nam, lex.sem.id);
	void (*init)(cvec_padding_t**,int*)=dlsym(G.dll, nam);
	if (!dlerror()) {
	    cvec_padding_t *foo;
	    int n;
	    init(&foo,&n);
			
	    cveclut_t lut[n];
	    for (int i=0; i<n; ++i) {
		lut[i].nam=foo[i].nam;
		lut[i].v=foo[i].v;
	    }
	    arr_pushn(G.cvec_lut, lut, n);
	}

	lex_gettok();
    }

    lex_end();
}

static void cvec_shutdown()
{
    arr_kill(G.cvec_lut);
}

float *cvec_get(char *nam)
{
    int sz=arr_len(G.cvec_lut);
    for (int i=0; i<sz; ++i) {
	if (!strcmp(G.cvec_lut[i].nam, nam)) return G.cvec_lut[i].v;
    }

    return NULL;
}



#include <emmintrin.h>  // sse2
#include <string.h>

/*
  gcc -msse2

*/


#define C565_5_MASK 0xF8 // 0xFF minus last three bits
#define C565_6_MASK 0xFC // 0xFF minus last two bits

#define INSET_SHIFT 4 // inset the bounding box with ( range >> shift )

#if !defined(MAX_INT)
#define       MAX_INT         2147483647      /* max value for an int 32 */
#define       MIN_INT         (-2147483647-1) /* min value for an int 32 */
#endif


#if defined(__GNUC__)
#define   ALIGN16(_x)   _x __attribute((aligned(16)))
#else
#define   ALIGN16( x ) __declspec(align(16)) x
#endif


static void ExtractBlock_Intrinsics( const byte *inPtr, int width, byte *colorBlock );
static void GetMinMaxColors_Intrinsics( const byte *colorBlock, byte *minColor, byte *maxColor );
static void EmitColorIndices_Intrinsics( const byte *colorBlock, const byte *minColor, const byte *maxColor, byte **outData);
static word ColorTo565( const byte *color );
static void EmitWord( word s, byte**);


static word ColorTo565( const byte *color )
{
    return ( ( color[ 0 ] >> 3 ) << 11 ) |
	( ( color[ 1 ] >> 2 ) <<  5 ) |
	(   color[ 2 ] >> 3 );
}


static void EmitWord( word s, byte **outData)
{
    (*outData)[0] = ( s >> 0 ) & 255;
    (*outData)[1] = ( s >> 8 ) & 255;
    *outData += 2;
}


static void ExtractBlock_Intrinsics( const byte *inPtr, int width, byte *colorBlock ) 
{
    __m128i t0, t1, t2, t3;
    register int w = width << 2;  // width*4

    t0 = _mm_load_si128 ( (__m128i*) inPtr );
    _mm_store_si128 ( (__m128i*) &colorBlock[0], t0 );   // copy first row, 16bytes

    t1 = _mm_load_si128 ( (__m128i*) (inPtr + w) );
    _mm_store_si128 ( (__m128i*) &colorBlock[16], t1 );   // copy second row

    t2 = _mm_load_si128 ( (__m128i*) (inPtr + 2*w) );
    _mm_store_si128 ( (__m128i*) &colorBlock[32], t2 );   // copy third row

    inPtr = inPtr + w;     // add width, intead of *3

    t3 = _mm_load_si128 ( (__m128i*) (inPtr + 2*w) );
    _mm_store_si128 ( (__m128i*) &colorBlock[48], t3 );   // copy last row
}


#define R_SHUFFLE_D( x, y, z, w ) (( (w) & 3 ) << 6 | ( (z) & 3 ) << 4 | ( (y) & 3 ) << 2 | ( (x) & 3 ))

ALIGN16( static byte SIMD_SSE2_byte_0[16] ) = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

static void GetMinMaxColors_Intrinsics( const byte *colorBlock, byte *minColor, byte *maxColor )
{
    __m128i t0, t1, t3, t4, t6, t7;

    // get bounding box
    // ----------------
    
    // load the first row
    t0 = _mm_load_si128 ( (__m128i*) colorBlock );
    t1 = _mm_load_si128 ( (__m128i*) colorBlock );

    __m128i t16 = _mm_load_si128 ( (__m128i*) (colorBlock+16) );
    // Minimum of Packed Unsigned Byte Integers
    t0 = _mm_min_epu8 ( t0, t16);
    // Maximum of Packed Unsigned Byte Integers
    t1 = _mm_max_epu8 ( t1, t16);
    
    __m128i t32 = _mm_load_si128 ( (__m128i*) (colorBlock+32) );
    t0 = _mm_min_epu8 ( t0, t32);
    t1 = _mm_max_epu8 ( t1, t32);
    
    __m128i t48 = _mm_load_si128 ( (__m128i*) (colorBlock+48) );
    t0 = _mm_min_epu8 ( t0, t48);
    t1 = _mm_max_epu8 ( t1, t48);
    
    // Shuffle Packed Doublewords
    t3 = _mm_shuffle_epi32( t0, R_SHUFFLE_D( 2, 3, 2, 3 ) );
    t4 = _mm_shuffle_epi32( t1, R_SHUFFLE_D( 2, 3, 2, 3 ) );
    
    t0 = _mm_min_epu8 ( t0, t3);
    t1 = _mm_max_epu8 ( t1, t4);
    
    // Shuffle Packed Low Words
    t6 = _mm_shufflelo_epi16( t0, R_SHUFFLE_D( 2, 3, 2, 3 ) );
    t7 = _mm_shufflelo_epi16( t1, R_SHUFFLE_D( 2, 3, 2, 3 ) );
    
    t0 = _mm_min_epu8 ( t0, t6);
    t1 = _mm_max_epu8 ( t1, t7);
    
    // inset the bounding box
    // ----------------------
    
    // Unpack Low Data
    //__m128i t66 = _mm_set1_epi8( 0 );
    __m128i t66 = _mm_load_si128 ( (__m128i*) SIMD_SSE2_byte_0 );
    t0 = _mm_unpacklo_epi8(t0, t66);
    t1 = _mm_unpacklo_epi8(t1, t66);
    
    // copy (movdqa)
    //__m128i t2 = _mm_load_si128 ( &t1 );
    __m128i t2 = t1;
    
    // Subtract Packed Integers
    t2 = _mm_sub_epi16(t2, t0);
    
    // Shift Packed Data Right Logical 
    t2 = _mm_srli_epi16(t2, INSET_SHIFT);
    
    // Add Packed Integers
    t0 = _mm_add_epi16(t0, t2);
    
    t1 = _mm_sub_epi16(t1, t2);
    
    // Pack with Unsigned Saturation
    t0 = _mm_packus_epi16(t0, t0);
    t1 = _mm_packus_epi16(t1, t1);
    
    // store bounding box extents
    // --------------------------
    _mm_store_si128 ( (__m128i*) minColor, t0 );
    _mm_store_si128 ( (__m128i*) maxColor, t1 );
}


ALIGN16( static word SIMD_SSE2_word_0[8] ) = { 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000 };
ALIGN16( static word SIMD_SSE2_word_1[8] ) = { 0x0001, 0x0001, 0x0001, 0x0001, 0x0001, 0x0001, 0x0001, 0x0001 };
ALIGN16( static word SIMD_SSE2_word_2[8] ) = { 0x0002, 0x0002, 0x0002, 0x0002, 0x0002, 0x0002, 0x0002, 0x0002 };
ALIGN16( static word SIMD_SSE2_word_div_by_3[8] ) = { (1<<16)/3+1, (1<<16)/3+1, (1<<16)/3+1, (1<<16)/3+1, (1<<16)/3+1, (1<<16)/3+1, (1<<16)/3+1, (1<<16)/3+1 };
ALIGN16( static byte SIMD_SSE2_byte_colorMask[16] ) = { C565_5_MASK, C565_6_MASK, C565_5_MASK, 0x00, 0x00, 0x00, 0x00, 0x00, C565_5_MASK, C565_6_MASK, C565_5_MASK, 0x00, 0x00, 0x00, 0x00, 0x00 };

static void EmitColorIndices_Intrinsics( const byte *colorBlock, const byte *minColor, const byte *maxColor, byte **outData )
{
    ALIGN16( byte color0[16] );
    ALIGN16( byte color1[16] );
    ALIGN16( byte color2[16] );
    ALIGN16( byte color3[16] );
    ALIGN16( byte result[16] );
	
    // mov esi, maxColor
    // mov edi, minColor

    __m128i t0, t1, t2, t3, t4, t5, t6, t7;

    t7 = _mm_setzero_si128();
    //t7 = _mm_xor_si128(t7, t7);
    _mm_store_si128 ( (__m128i*) &result, t7 );


    //t0 = _mm_load_si128 ( (__m128i*)  maxColor );
    t0 = _mm_cvtsi32_si128( *(int*)maxColor);

    // Bitwise AND
    __m128i tt = _mm_load_si128 ( (__m128i*) SIMD_SSE2_byte_colorMask );
    t0 = _mm_and_si128(t0, tt);

    t0 = _mm_unpacklo_epi8(t0, t7);

    t4 = _mm_shufflelo_epi16( t0, R_SHUFFLE_D( 0, 3, 2, 3 ));
    t5 = _mm_shufflelo_epi16( t0, R_SHUFFLE_D( 3, 1, 3, 3 ));

    t4 = _mm_srli_epi16(t4, 5);
    t5 = _mm_srli_epi16(t5, 6);

    // Bitwise Logical OR
    t0 = _mm_or_si128(t0, t4);
    t0 = _mm_or_si128(t0, t5);   // t0 contains color0 in 565




    //t1 = _mm_load_si128 ( (__m128i*)  minColor );
    t1 = _mm_cvtsi32_si128( *(int*)minColor);

    t1 = _mm_and_si128(t1, tt);

    t1 = _mm_unpacklo_epi8(t1, t7);

    t4 = _mm_shufflelo_epi16( t1, R_SHUFFLE_D( 0, 3, 2, 3 ));
    t5 = _mm_shufflelo_epi16( t1, R_SHUFFLE_D( 3, 1, 3, 3 ));

    t4 = _mm_srli_epi16(t4, 5);
    t5 = _mm_srli_epi16(t5, 6);

    t1 = _mm_or_si128(t1, t4);
    t1 = _mm_or_si128(t1, t5);  // t1 contains color1 in 565



    t2 = t0;

    t2 = _mm_packus_epi16(t2, t7);

    t2 = _mm_shuffle_epi32( t2, R_SHUFFLE_D( 0, 1, 0, 1 ));

    _mm_store_si128 ( (__m128i*) &color0, t2 );

    t6 = t0;
    t6 = _mm_add_epi16(t6, t0);
    t6 = _mm_add_epi16(t6, t1);

    // Multiply Packed Signed Integers and Store High Result
    __m128i tw3 = _mm_load_si128 ( (__m128i*) SIMD_SSE2_word_div_by_3 );
    t6 = _mm_mulhi_epi16(t6, tw3);
    t6 = _mm_packus_epi16(t6, t7);

    t6 = _mm_shuffle_epi32( t6, R_SHUFFLE_D( 0, 1, 0, 1 ));

    _mm_store_si128 ( (__m128i*) &color2, t6 );

    t3 = t1;
    t3 = _mm_packus_epi16(t3, t7);
    t3 = _mm_shuffle_epi32( t3, R_SHUFFLE_D( 0, 1, 0, 1 ));

    _mm_store_si128 ( (__m128i*) &color1, t3 );

    t1 = _mm_add_epi16(t1, t1);
    t0 = _mm_add_epi16(t0, t1);

    t0 = _mm_mulhi_epi16(t0, tw3);
    t0 = _mm_packus_epi16(t0, t7);

    t0 = _mm_shuffle_epi32( t0, R_SHUFFLE_D( 0, 1, 0, 1 ));
    _mm_store_si128 ( (__m128i*) &color3, t0 );

    __m128i w0 = _mm_load_si128 ( (__m128i*) SIMD_SSE2_word_0);
    __m128i w1 = _mm_load_si128 ( (__m128i*) SIMD_SSE2_word_1);
    __m128i w2 = _mm_load_si128 ( (__m128i*) SIMD_SSE2_word_2);

    // mov eax, 32
    // mov esi, colorBlock
    int x = 32;
    //const byte *c = colorBlock;
    while (x >= 0)
	{
	    t3 = _mm_loadl_epi64( (__m128i*) (colorBlock+x+0));
	    t3 = _mm_shuffle_epi32( t3, R_SHUFFLE_D( 0, 2, 1, 3 ));
	    
	    t5 = _mm_loadl_epi64( (__m128i*) (colorBlock+x+8));
	    t5 = _mm_shuffle_epi32( t5, R_SHUFFLE_D( 0, 2, 1, 3 ));

	    t0 = t3;
	    t6 = t5;
	    // Compute Sum of Absolute Difference
	    __m128i c0 = _mm_load_si128 ( (__m128i*)  color0 );
	    t0 = _mm_sad_epu8(t0, c0);
	    t6 = _mm_sad_epu8(t6, c0);
	    // Pack with Signed Saturation 
	    t0 = _mm_packs_epi32 (t0, t6);

	    t1 = t3;
	    t6 = t5;
	    __m128i c1 = _mm_load_si128 ( (__m128i*)  color1 );
	    t1 = _mm_sad_epu8(t1, c1);
	    t6 = _mm_sad_epu8(t6, c1);
	    t1 = _mm_packs_epi32 (t1, t6);

	    t2 = t3;
	    t6 = t5;
	    __m128i c2 = _mm_load_si128 ( (__m128i*)  color2 );
	    t2 = _mm_sad_epu8(t2, c2);
	    t6 = _mm_sad_epu8(t6, c2);
	    t2 = _mm_packs_epi32 (t2, t6);

	    __m128i c3 = _mm_load_si128 ( (__m128i*)  color3 );
	    t3 = _mm_sad_epu8(t3, c3);
	    t5 = _mm_sad_epu8(t5, c3);
	    t3 = _mm_packs_epi32 (t3, t5);


	    t4 = _mm_loadl_epi64( (__m128i*) (colorBlock+x+16));
	    t4 = _mm_shuffle_epi32( t4, R_SHUFFLE_D( 0, 2, 1, 3 ));
	    
	    t5 = _mm_loadl_epi64( (__m128i*) (colorBlock+x+24));
	    t5 = _mm_shuffle_epi32( t5, R_SHUFFLE_D( 0, 2, 1, 3 ));

	    t6 = t4;
	    t7 = t5;
	    t6 = _mm_sad_epu8(t6, c0);
	    t7 = _mm_sad_epu8(t7, c0);
	    t6 = _mm_packs_epi32 (t6, t7);
	    t0 = _mm_packs_epi32 (t0, t6);  // d0

	    t6 = t4;
	    t7 = t5;
	    t6 = _mm_sad_epu8(t6, c1);
	    t7 = _mm_sad_epu8(t7, c1);
	    t6 = _mm_packs_epi32 (t6, t7);
	    t1 = _mm_packs_epi32 (t1, t6);  // d1

	    t6 = t4;
	    t7 = t5;
	    t6 = _mm_sad_epu8(t6, c2);
	    t7 = _mm_sad_epu8(t7, c2);
	    t6 = _mm_packs_epi32 (t6, t7);
	    t2 = _mm_packs_epi32 (t2, t6);  // d2

	    t4 = _mm_sad_epu8(t4, c3);
	    t5 = _mm_sad_epu8(t5, c3);
	    t4 = _mm_packs_epi32 (t4, t5);
	    t3 = _mm_packs_epi32 (t3, t4);  // d3

	    t7 = _mm_load_si128 ( (__m128i*) result );

	    t7 = _mm_slli_epi32( t7, 16);

	    t4 = t0;
	    t5 = t1;
	    // Compare Packed Signed Integers for Greater Than
	    t0 = _mm_cmpgt_epi16(t0, t3); // b0
	    t1 = _mm_cmpgt_epi16(t1, t2); // b1
	    t4 = _mm_cmpgt_epi16(t4, t2); // b2
	    t5 = _mm_cmpgt_epi16(t5, t3); // b3
	    t2 = _mm_cmpgt_epi16(t2, t3); // b4
	      
	    t4 = _mm_and_si128(t4, t1); // x0
	    t5 = _mm_and_si128(t5, t0); // x1
	    t2 = _mm_and_si128(t2, t0); // x2

	    t4 = _mm_or_si128(t4, t5);
	    t2 = _mm_and_si128(t2, w1);
	    t4 = _mm_and_si128(t4, w2);
	    t2 = _mm_or_si128(t2, t4);

	    t5 = _mm_shuffle_epi32( t2, R_SHUFFLE_D( 2, 3, 0, 1 ));

	    // Unpack Low Data
	    t2 = _mm_unpacklo_epi16 ( t2, w0);
	    t5 = _mm_unpacklo_epi16 ( t5, w0);

	    //t5 = _mm_slli_si128 ( t5, 8);
	    t5 = _mm_slli_epi32( t5, 8);

	    t7 = _mm_or_si128(t7, t5);
	    t7 = _mm_or_si128(t7, t2);

	    _mm_store_si128 ( (__m128i*) &result, t7 );

	    x -=32;
	}

    t4 = _mm_shuffle_epi32( t7, R_SHUFFLE_D( 1, 2, 3, 0 ));
    t5 = _mm_shuffle_epi32( t7, R_SHUFFLE_D( 2, 3, 0, 1 ));
    t6 = _mm_shuffle_epi32( t7, R_SHUFFLE_D( 3, 0, 1, 2 ));

    t4 = _mm_slli_epi32 ( t4, 2);
    t5 = _mm_slli_epi32 ( t5, 4);
    t6 = _mm_slli_epi32 ( t6, 6);

    t7 = _mm_or_si128(t7, t4);
    t7 = _mm_or_si128(t7, t5);
    t7 = _mm_or_si128(t7, t6);

    //_mm_store_si128 ( (__m128i*) outData, t7 );

    int r = _mm_cvtsi128_si32 (t7);
    memcpy(*outData, &r, 4);   // Anything better ?

    *outData += 4;
}

void img_todxt1( const byte *inBuf, byte *outBuf,
		 int width, int height, int *outputBytes )
{
    ALIGN16( byte *outData );
    ALIGN16( byte block[64] );
    ALIGN16( byte minColor[4] );
    ALIGN16( byte maxColor[4] );

    outData = outBuf;
    for ( int j = 0; j < height; j += 4, inBuf += width * 4*4 ) {
	for ( int i = 0; i < width; i += 4 ) {
	    ExtractBlock_Intrinsics( inBuf + i * 4, width, block );
	    GetMinMaxColors_Intrinsics( block, minColor, maxColor );
	    EmitWord( ColorTo565( maxColor ), &outData );
	    EmitWord( ColorTo565( minColor ), &outData );
	    EmitColorIndices_Intrinsics( block, minColor, maxColor, &outData );
	}
    }
    *outputBytes = (int) ( outData - outBuf );
}



static SDL_Surface *vid_init()
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE,  24);	
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);

    int w,h;
    w=G.screenw=800;
    h=G.screenh=600;
    G.surface=SDL_SetVideoMode(w, h, 0, SDL_OPENGL);	
}



#include <ft2build.h>
#include FT_FREETYPE_H



/*--------------------------------------------
  Load a font. A texture atlas is generated.

  You can have at most 8 textures; A character is at most 32x32 pixels

  => Each font has 128 x 32 x 32 pixels
  => 8 fonts have 128 x 32 x 32 x 8 = 1024x1024 = 1M


  Wanna use additional fonts? Higher resolutions? Then you have to DIY.

  There should be one default font -- console/minibuffer/when no font can be found.

  => only 7 custom fonts.

  You can read 0-7, but can write only 1-7 -- font[0] is reserved.


  -------------------------------------------*/
static void fnt_force_load(const char *name, fnt_t *fnt)
{
    FT_Library	library;
    FT_Face		face;

    if(FT_Init_FreeType(&library)) return NULL;
    if(FT_New_Face(library,name,0,&face)) return NULL;

    fnt->occupied = 1;

    // first pass: collect fnt info

    // The font can be at any size. We have to scale the text so that it can fit 32x32 better.
    FT_Set_Pixel_Sizes(face, 18,18); // TODO: how to set font size?

    short max_w=0;
    short max_h=0;
    float totalw=0;
    for(uchar c=0; c<128; c++){
	int err=FT_Load_Char(face, c, FT_LOAD_RENDER);
	int w = face->glyph->metrics.width>>6;
	if(w > max_w) max_w=w;
	int h = face->glyph->metrics.height>>6;
	if(h > max_h) max_h = h;

	fnt->height = face->size->metrics.height>>6; // TODO: need more investigation
	fnt->ascender = face->size->metrics.ascender>>6;
	fnt->ch[c].w = w;
	fnt->ch[c].h = h;
	fnt->ch[c].x = face->glyph->metrics.horiBearingX>>6;
	fnt->ch[c].y = (face->glyph->metrics.horiBearingY>>6) - h;
	fnt->ch[c].advance = face->glyph->metrics.horiAdvance>>6;

	totalw+=w;
    }

    fnt->avgw=totalw/128;

    //	int size = nextp2(max_w>max_h ? max_w : max_h);

    //	glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, G.fnt_tex);  // TODO

    int sz = 32;

    // second pass: store all charactors in a 16x8 grid
    // if we store all characters in one texture, then texture filtering becomes a problem.

    int sz2=sz*sz;
    uchar *emptybmp = alloca(sz2);
    memset(emptybmp, 0,sz2);
    int x = 0;
    int y = (fnt-G.fnt_slots)*4*sz;
    int texw = 32*sz;
    for(uchar c=0; c<128; c++){
	FT_Load_Char(face, c, FT_LOAD_RENDER);
	FT_Bitmap bitmap = face->glyph->bitmap;

	// DEBUGGED: the bitmap buffer is upsidedown => (fnt->ch[c].h-1-j) to correct it
	for(int j=0; j<32; j++){
	    for(int i=0; i<32; i++){
		if(i<fnt->ch[c].w && j<fnt->ch[c].h) 
		    emptybmp[j*32+i] = bitmap.buffer[(fnt->ch[c].h-1-j)*fnt->ch[c].w+i];
		else 
		    emptybmp[j*32+i]=0;
	    }
	}

	fnt->ch[c].u = (float)x/(float)texw;
	fnt->ch[c].v = (float)y/(float)texw;
	fnt->ch[c].texw= (float)fnt->ch[c].w/(float)texw;
	fnt->ch[c].texh= (float)fnt->ch[c].h/(float)texw;

	// upload character image
	glTexSubImage2D(GL_TEXTURE_2D, 0, x, y, 32,32,
			GL_ALPHA, GL_UNSIGNED_BYTE, emptybmp); // DEBUGGED: GL_LUMINANCE_ALPHA => GL_ALPHA

	if(x==texw-sz){
	    y+=sz;
	}
	x = (x+sz)%texw;
    }

    FT_Done_Face(face);
    FT_Done_FreeType(library);
}


#define FNT_VB_NVERTS 2048
static void fnt_init()
{
    G.fnt_tex = tex1b(1024, 1024, 011, NULL);
    fnt_force_load("courbd.ttf", G.fnt_slots);
    //	fnt_force_load("Test.ttf", G.fnt_slots);

    G.fnt_vbo = vb_alloc(FNT_VB_NVERTS, sizeof(fnt_vbfmt_t), &G.fnt_vbocap);
    G.fnt_shader = vfs_load("shaders/font.c");
    G.fnt_cur = G.fnt_slots;
    G.fnt_scale = 1;
}

static void fnt_shutdown()
{
    glDeleteTextures(1, &G.fnt_tex);
}

/*
  
 */
int fnt_load(const char *nam)
{
    // find an empty slot
    int i;
    for(i=1; i<8; i++){
	if(!G.fnt_slots[i].occupied){
	    fnt_force_load(nam, G.fnt_slots+i);
	    break;
	}
    }
	
    return G.fnt_slots[i%8].occupied * (i%8);
}



void fnt_loadslot(const char *nam, int i)
{
    i %= 8;
    if(!i) return; // cannot overwrite font[0]

    fnt_force_load(nam, G.fnt_slots+i);
}



/*
  set the current font to "nam"
*/
void fnt_use(const char *nam)
{
    for(int i=0; i<8; i++){
	fnt_t *fnt = G.fnt_slots + i;
	if(fnt->occupied){
	    if(!strcmp(fnt->nam, nam)){
		G.fnt_cur = fnt;
		return;
	    }
	}
    }

    G.fnt_cur = G.fnt_slots;
}


void fnt_useslot(int i)
{
    i %= 8;
    fnt_t *fnt = G.fnt_slots + i;
    G.fnt_cur = fnt->occupied ? fnt : G.fnt_slots;
}


void fnt_color3f(float r, float g, float b)
{
    float *rgba = G.fnt_rgba;
    rgba[0]=r;
    rgba[1]=g;
    rgba[2]=b;
    rgba[3]=1;
}

void fnt_color3fv(float c[3])
{
    float *rgba=G.fnt_rgba;
    rgba[0]=c[0];
    rgba[1]=c[1];
    rgba[2]=c[2];
    rgba[3]=1;
}

void fnt_color4f(float r, float g, float b, float a)
{
    float *rgba = G.fnt_rgba;
    rgba[0]=r;
    rgba[1]=g;
    rgba[2]=b;
    rgba[3]=a;
}

void fnt_color4fv(float c[4])
{
    float *rgba=G.fnt_rgba;
    rgba[0]=c[0];
    rgba[1]=c[1];
    rgba[2]=c[2];
    rgba[3]=c[3];
}

void fnt_getcolor(float c[4])
{
    float *rgba=G.fnt_rgba;
    c[0]=rgba[0];
    c[1]=rgba[1];
    c[2]=rgba[2];
    c[3]=rgba[3];
}

void fnt_scale(float s)
{
    G.fnt_scale = s;
}


fnt_t *fnt_getcur()
{
    return G.fnt_cur;
}


void rnd_text(float bnd[4], char *fmt, ...)
{
    char str[128];

    va_list va;
    va_start(va, fmt);
    vsnprintf(str, 128, fmt, va);
    va_end(va);
	
    fnt_t *fnt=G.fnt_cur;
    float factor = (bnd[3]-bnd[1])/G.fnt_cur->height;;
    char *c=str;
    float x=bnd[0];
    float y=bnd[3]-fnt->ascender*factor;
    float z=G.cur_z;
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, G.fnt_tex);
    glBegin(GL_QUADS);
    while (*c) {
	fntch_t *ch=&fnt->ch[*c];
	float w=factor*ch->w;
	float h=factor*ch->h;
	float texw=ch->texw;
	float texh=ch->texh;
	float dx=factor*ch->x;
	float dy=factor*ch->y;

	if (x+dx+w>bnd[2]) break;
	glTexCoord2f(ch->u,ch->v); glVertex3f(x+dx, y+dy, z);
	glTexCoord2f(ch->u+texw,ch->v); glVertex3f(x+dx+w, y+dy, z);
	glTexCoord2f(ch->u+texw,ch->v+texh); glVertex3f(x+dx+w, y+dy+h, z);
	glTexCoord2f(ch->u,ch->v+texh); glVertex3f(x+dx, y+dy+h, z);
	x+=factor*ch->advance;

	++c;
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);
}


void fnt_printf(int posx, int posy, int alignment, const char *fmt, ...)
{
    char string[2048];

    va_list va;

    va_start(va, fmt);
    vsnprintf(string, 2048, fmt, va);
    va_end(va);
	
    float scalefactor = G.fnt_scale;

    char *line=string;
    int i=0;
    int totalw=0;
    char *bgn = string;
    fnt_vbfmt_t f[4];

    for(int i=0; i<4; i++){
	memcpy(f[i].rgba, G.fnt_rgba, 16);
    }

    fnt_t *fnt=G.fnt_cur;
    float y=posy;

    // render each line
    while(*line){
	if(alignment){
	    char *c=line;
	    while(*c && *c!='\n'){
		totalw += fnt->ch[*c].advance;
		++c;
	    }
	}


	float ofs=0;
	if(alignment == 1){
	    ofs = -totalw * 0.5f * scalefactor;

	    // character VB: each character = 4 x {float2 position, float2 texcoord, byte4 rgba}
	}
	else if(alignment == 2){
	    ofs = -totalw * scalefactor;
	}
		

	float x0=posx + ofs;

	/*
	  3 2
	  0 1
	*/

	char *c;
	int kk=0;
	float x=x0;
	for(c=line; *c && *c!='\n'; ++c){
	    fntch_t *ch = &fnt->ch[*c];
	    float w=scalefactor * ch->w;
	    float h=scalefactor * ch->h;
	    float texw = ch->texw;
	    float texh = ch->texh;

	    f[0].x = f[3].x = x+ch->x;
	    f[1].x = f[2].x = f[0].x+w;
	    f[0].y = f[1].y = y+ch->y;
	    f[2].y = f[3].y = f[0].y+h;

	    f[0].u = f[3].u = ch->u;
	    f[1].u = f[2].u = ch->u+texw;
	    f[0].v = f[1].v = ch->v;
	    f[2].v = f[3].v = ch->v+texh;

	    x+=ch->advance*scalefactor;
	    arr_pushn(G.fnt_vb, f, 4);
	}

	if(*c=='\n') ++c;
	line = c;
	y-=fnt->height*scalefactor;
    }
}


static void fnt_flush()
{
    rnd_scrncoord();
    vfs_use(G.fnt_shader);
    samp2d("fnt_tex", G.fnt_tex);

    glEnable(GL_BLEND);
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

    int n = arr_len(G.fnt_vb);
    int vbo=G.fnt_vbo;
    int sz=sizeof(fnt_vbfmt_t);
    int cap = vb_getcap(vbo)/sz;

    fnt_vbfmt_t *vb=G.fnt_vb;
    while(n>0){
	if(n/cap){
	    vb_update(vbo,cap,sz,vb);
	    vb+=cap;
	}
	else{
	    vb_update(vbo,n,sz,vb);
	}

	n-=cap;
	rendq("{co texco color}", -1, vbo);
    }

    arr_popall(G.fnt_vb);
    glDisable(GL_BLEND);
}





/*
  console
*/
#define MINIBUF_SIZE 255
#define MBHIST_RECALL_LIMIT	64
#define CS_MAX_W	157
#define CS_NOL	512

/*--------------------------------------
  w - number of chars each line
  h - number of lines of the output buffer
  -------------------------------------*/
static void cs_init()
{
    int w=70, h=CS_NOL;
    G.cs_w = w;
    G.cs_h = h;
    G.cs_lines = malloc(sizeof(char*)*h);
    char *mem = malloc((CS_MAX_W+1)*h);
    for(int i=0; i<h; i++, mem+=(CS_MAX_W+1)){
	mem[0] = 0;
	G.cs_lines[i]=mem;
    }

    // mini-buffer cmd history
    G.cs_mbmem = G.cs_mbhead = malloc(sizeof(mblist_t) * MBHIST_RECALL_LIMIT);
    mblist_t *end = G.cs_mbhead + MBHIST_RECALL_LIMIT;
    for(mblist_t *mb=G.cs_mbhead; mb<end; mb++){
	mb->prev = mb-1;
	mb->next = mb+1;	
	mb->mb[0] = 0;
    }
    end[-1].next = G.cs_mbhead;
    G.cs_mbhead->prev = end-1;
    G.cs_mbcur = G.cs_mbhead;
    G.cs_outh = 13; // TODO: shouldn't be hardcoded

    for(int i=0; i<4; i++){
	G.cs_fg[i]=1.0f;
	G.cs_bg[i]=.5f;
	G.cs_mbfg[i]=1.0;
	G.cs_mbbg[i]=.3f;
    }
}

static void cs_shutdown()
{
    free(G.cs_lines[0]);
    free(G.cs_lines);
    free(G.cs_mbmem);
}

void cs_printf(const char *fmt, ...)
{
    char buf[512];
    va_list	ap;

    va_start(ap, fmt);
    vsnprintf(buf, 512, fmt, ap);
    vprintf(fmt, ap); // good old printf
    va_end(ap);

    int pos=0;
    int in=G.cs_in;
    int w=G.cs_w;
    int h=G.cs_h;
    int n=G.cs_nol;

    char **lines = G.cs_lines;
    char *line = lines[in];
    char *c = buf;
    while(*c){
	if(pos==w){
	    line[w] = '\0';
	    in = (in+1)%h;
	    line = lines[in];
	    pos = 0;
	    n = n+1;
	    if(n>h) n=h;
	}			

	if(*c=='\n'){
	    line[pos] = '\0';
	    in = (in+1)%h;
	    line = lines[in];
	    pos = 0;
	    n=n+1;
	    if(n>h) n=h;
	    ++c;
	    continue;
	}
	line[pos++] = *c++;
    }

    G.cs_nol=n;
    G.cs_in=in;

    int out=G.cs_out;
    int outh=G.cs_outh;
    if(h < outh){
	if(n<h){
	    out=0;
	}
	else{
	    out=in;
	}
    }
    else{
	if(n<outh)
	    out=0;
	else
	    out=(in-outh+h)%h;
    }

    G.cs_out=out;
}

static void cs_scrollup()
{
    int n=G.cs_nol;
    int h=G.cs_h;
    int in=G.cs_in;
    int out=G.cs_out;

    if (n<h) {
	if (!out) return;
    }
    else{
	if (out==in) return;
    }

    G.cs_out = (out-1+h)%h;
}

static void cs_scrolldown()
{
    int n=G.cs_nol;
    int h=G.cs_h;
    int in=G.cs_in;
    int out=G.cs_out;
    int outh=G.cs_outh;
    if(outh>h) outh=h;

    if (n<h) {
	if (out==in) return;
    }
    else {
	if ((out+outh)%h == in) return;
    }

    G.cs_out = (out+1)%h;
}


static void cs_end()
{
    int n=G.cs_nol;
    int h=G.cs_h;
    int in=G.cs_in;
    int outh=G.cs_outh;
    if(outh>h) outh=h;

    G.cs_out = n<h ? 0 : (in-outh+h)%h;
}


void rnd_scrncoord()
{
    glViewport(0,0,G.screenw,G.screenh);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0,G.screenw, 0, G.screenh, -1, 1);
}


void rnd_scrncoordsz(int w, int h)
{
    glViewport(0,0,w,h);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0,w, 0, h, -1, 1);
}


void rnd_quad(float rgba[4], float qbnd[4])
{
    rnd_scrncoord();

    glColor4fv(rgba);
    glBegin(GL_QUADS);
    glVertex2f(qbnd[0], qbnd[1]);
    glVertex2f(qbnd[2], qbnd[1]);
    glVertex2f(qbnd[2], qbnd[3]);
    glVertex2f(qbnd[0], qbnd[3]);
    glEnd();

    ++G.rnd_counter;
}


void rnd_quadsz(int w, int h, float rgba[4], float qbnd[4])
{
    rnd_scrncoordsz(w,h);

    glColor4fv(rgba);
    glBegin(GL_QUADS);
    glVertex2f(qbnd[0], qbnd[1]);
    glVertex2f(qbnd[2], qbnd[1]);
    glVertex2f(qbnd[2], qbnd[3]);
    glVertex2f(qbnd[0], qbnd[3]);
    glEnd();

    ++G.rnd_counter;
}



/*------------------------------------
  bnd[] = {xmin,ymin,xmax,ymax}

  single texture;
  quad is in screen space

  -----------------------------------*/
void rnd_texquad(uint tex, float tcbnd[4], float qbnd[4])
{
    rnd_scrncoord();
    glColor4f(1,1,1,1);
    glEnable(GL_TEXTURE_2D); //TODO
    glBindTexture(GL_TEXTURE_2D, tex);
    glBegin(GL_QUADS);
    glTexCoord2f(tcbnd[0], tcbnd[1]);	glVertex2f(qbnd[0], qbnd[1]);
    glTexCoord2f(tcbnd[2], tcbnd[1]);	glVertex2f(qbnd[2], qbnd[1]);
    glTexCoord2f(tcbnd[2], tcbnd[3]);	glVertex2f(qbnd[2], qbnd[3]);
    glTexCoord2f(tcbnd[0], tcbnd[3]);	glVertex2f(qbnd[0], qbnd[3]);
    glEnd();
    glDisable(GL_TEXTURE_2D); //TODO

    ++G.rnd_counter;	
}


void rnd_texquadsz(uint tex, int w, int h, float tcbnd[4], float qbnd[4])
{
    rnd_scrncoordsz(w,h);
    glEnable(GL_TEXTURE_2D); //TODO
    glBindTexture(GL_TEXTURE_2D, tex);
    glBegin(GL_QUADS);
    glTexCoord2f(tcbnd[0], tcbnd[1]);	glVertex2f(qbnd[0], qbnd[1]);
    glTexCoord2f(tcbnd[2], tcbnd[1]);	glVertex2f(qbnd[2], qbnd[1]);
    glTexCoord2f(tcbnd[2], tcbnd[3]);	glVertex2f(qbnd[2], qbnd[3]);
    glTexCoord2f(tcbnd[0], tcbnd[3]);	glVertex2f(qbnd[0], qbnd[3]);
    glEnd();
    glDisable(GL_TEXTURE_2D); //TODO

    ++G.rnd_counter;	
}

void rnd_texrect(uint tex, float tcbnd[4], float qbnd[4])
{
    glColor4f(1,1,1,1);
    glEnable(GL_TEXTURE_2D); //TODO
    glBindTexture(GL_TEXTURE_2D, tex);
    glBegin(GL_QUADS);
    glTexCoord2f(tcbnd[0], tcbnd[1]);	glVertex2f(qbnd[0], qbnd[1]);
    glTexCoord2f(tcbnd[2], tcbnd[1]);	glVertex2f(qbnd[2], qbnd[1]);
    glTexCoord2f(tcbnd[2], tcbnd[3]);	glVertex2f(qbnd[2], qbnd[3]);
    glTexCoord2f(tcbnd[0], tcbnd[3]);	glVertex2f(qbnd[0], qbnd[3]);
    glEnd();
    glDisable(GL_TEXTURE_2D); //TODO
}



/*
  command buffer

*/
#define CMD_CKSUM 0x17cb315b

void cmd_addtxt(char *t)
{
    int n=strlen(t)+1;
    arr_pushn(G.cmd_buf, t, n);
}


int cmd_geti(int *i)
{
    char *end;
    *i = strtol(G.cmd_cur, &end, 0);
	
    char *d=end-G.cmd_cur;
    G.cmd_cur = end;
    return d;
}

int cmd_getf(float *f)
{
    char *end;
    *f = strtof(G.cmd_cur, &end);

    char *d=end-G.cmd_cur;
    G.cmd_cur = end;
    return d;
}

/*----------------------------------
  Strings don't need to quoted, although you can do that.
  ---------------------------------*/
int cmd_gets(char *s)
{
    char *s0 = s;
    char *c = G.cmd_cur;
	
    while (isspace(*c)) ++c;

    if(*c=='\"') ++c;
    while(*c){
	if(*c=='\"' || *c==';' || *c=='\n'){
	    G.cmd_cur = c + (*c=='\"');
	    *s='\0';
	    return s-s0;
	}
	*s++ = *c++;
    }

    *s='\0';
    G.cmd_cur = c;
    return s-s0;
}

void cmd_execnow(char *cmd)
{
    if(!cmd) return;

    char buf[strlen(cmd)];

    while (*cmd) {
	char *c=cmd;
	while(*c && *c =='\b' || *c =='\r' || *c ==' '  || *c =='\t')
	    ++c;

	if(!*c) return;

	char *b=buf;
	char sign=0;
	if(*c=='+' || *c=='-') *b++=*c++;

	if(!isalnum(*c)){
	    cs_printf("Invalid command: %s\n", cmd);
	    return;
	}
	do{
	    *b++=*c++; 
	} while(isalnum(*c) || *c=='_');
	*b = '\0';

	G.cmd_cur = c;  // func() may use cmd. So update it
	void (*func)() = cvar_getp__(buf, CMD_CKSUM);

	if (func) func();
	else if (buf[0]=='+' || buf[0]=='-') {
	    func=cvar_getp__(buf+1, CMD_CKSUM);
	    if (!func) {goto error;}
	    if (buf[0]=='+') func();
	}
	else {goto error;}

	// ignore all unused cmd arguments (skip to '\n' or ';')
	c = G.cmd_cur;
	while(*c && *c!='\n' && *c!=';')
	    ++c; // skip '\n', ';' and '\0'
		
	if(*c) ++c;

	cmd = c;
    }

    return;

  error:
    cs_printf("Invalid command: %s\n", cmd);

}

void cmd_exec()
{
    // parse the command buffer	
    cmd_execnow(G.cmd_buf);		

    if(G.cmd_buf){
	arr_popall(G.cmd_buf);
	G.cmd_buf[0] = 0;
	G.cmd_cur = G.cmd_buf;
    }
}


void cmd_add(const char *nam, void (*func)())
{
    if(!nam || !func) return;

    cvar_setp__(nam, func, CMD_CKSUM);
}


/*-------------------------------------
  Output the minibuf to output buffer.
  ------------------------------------*/

static void mb_enter()
{
    char *mb = G.cs_mbcur->mb;
    if(!*mb) return;

    cs_printf("]$ %s\n", mb);
    cmd_execnow(G.cs_mbcur->mb);

    G.cs_mbhead = G.cs_mbhead->next;

    strcpy(G.cs_mbhead->mb, mb);
    G.cs_mbcur = G.cs_mbhead;
    G.cs_mbcur->mb[0] = 0;

    ++G.cs_nmb;
    if(G.cs_nmb > MBHIST_RECALL_LIMIT) G.cs_nmb = MBHIST_RECALL_LIMIT;
    G.cs_scroll = 0;
    G.cs_mbpos = 0;
}

static void mb_ins(char ins)
{
    if (G.cs_mbpos==MINIBUF_SIZE) return;
    char *c, *c0; 
    c = c0 = G.cs_mbcur->mb+G.cs_mbpos;
    while (*c) ++c;

    char *end = c;
    while(c!=c0){
	c[0] = c[-1];
	--c;
    }
	
    if(end - G.cs_mbcur->mb == MINIBUF_SIZE) *end=0;
    else end[1] = 0;
    *c0 = ins;
	
    if(end==c0){
	c[1]=0; // insert at the end of the minibuffer
    }

    ++G.cs_mbpos;
}

static void mb_backspace()
{
    if(!G.cs_mbpos) return;
    char *c=G.cs_mbcur->mb + G.cs_mbpos--;
    c[-1]=*c;
    if(!*c) return;
    while(*c){
	*c++ = c[1];
    }
}

static void mb_del()
{
    for(char *c=G.cs_mbcur->mb + G.cs_mbpos; *c; c++) *c=c[1];
}

static void mb_left()
{
    if(G.cs_mbpos>0) --G.cs_mbpos;
}

static void mb_right()
{
    if(G.cs_mbcur->mb[G.cs_mbpos]) ++G.cs_mbpos;
}

static void mb_prev()
{
    if(G.cs_scroll >= G.cs_nmb) return;
	
    ++G.cs_scroll;

    G.cs_mbcur = G.cs_mbcur->prev;
    G.cs_mbpos = strlen(G.cs_mbcur->mb);
}

static void mb_next()
{
    if(G.cs_scroll <= 0) return;

    --G.cs_scroll;

    G.cs_mbcur = G.cs_mbcur->next;
    G.cs_mbpos = strlen(G.cs_mbcur->mb);
}



/*
  viewport:

  x: [-1,1] => [0, w]
  y: [-1,1] => [0, h]


*/


//#define lrbt(l,r,b,t) glViewport((l), (b), (r)-(l), (t)-(b))


enum {MODE_GUIED=1, MODE_VU, MODE_FILTER, MODE_EDITOR};


/*
  1st person camera control
*/
static void updaterotmat(cam1_t *cam)
{
    float *m=cam->rotmat;

    m[0] = -cam->rotaxis[0].x;
    m[4] = -cam->rotaxis[0].y;
    m[8] = -cam->rotaxis[0].z;

    m[1] = cam->rotaxis[1].x;
    m[5] = cam->rotaxis[1].y;
    m[9] = cam->rotaxis[1].z;

    m[2] = cam->rotaxis[2].x;
    m[6] = cam->rotaxis[2].y;
    m[10] = cam->rotaxis[2].z;
}



enum {
    LOC_WALL=1,
    LOC_NOOK=2,
    LOC_IN=4,
    LOC_OUT=8,
};

enum {
    TYPE_PENDING,
    TYPE_IN,
    TYPE_OUT,
    TYPE_SPLIT,
};


enum {
    HIT_NOTHING,
    HIT_SOLID,
    HIT_FLOOR,
    HIT_CEILING,
    HIT_UPPER,
    HIT_LOWER
};

static char *hitstrings[] = {
    "HIT_NOTHING",
    "HIT_SOLID",
    "HIT_FLOOR",
    "HIT_CEILING",
    "HIT_UPPER",
    "HIT_LOWER",
};

int pvs2d(pvs2d_t *pvs, short sid, v3_t pos, v3_t dir, float fov);
static int chkloc(float x, float y, short *dat, int datsz);
float *interp2(float *c, const float *a, const float *b, float t);

int tracemotion(short s1, short e1, short *s2, short *e2, float co[3], float dst[3], float r);
static void cam1_tracemotion()
{
    
}


int ray2d(rayinfo_t *info, short s1, float *co, float *dst, float *hitt, rayintx_t *intx);


/*
  Move camera in its own coord system

  Does:
  1. Find which sector it's in:
  1) If not inited, find which sector it's currently in
  2) If the cam is previously inside none of the sectors, check which sector it's in
  3) If it's inited but the sector it's currently in has changed, find which sector it's currently in


  2. Trace the camera's motion to resolve collisions

  3. Update which sector it's in according to the result of the tracemotion() function
*/
void cam1_move(cam1_t *cam, float dx, float dy, float dz)
{
    if (dx==0 && dy==0 && dz==0) {
	++cam->updated;
	return;
    }

    // 1.1 & 1.2
    if (cam->sid == -1) {
	short sid;
	int loc=chkloc(cam->pos.x, cam->pos.y, &sid, 2);
	if (loc==LOC_IN) {
	    cam->sid=sid;
	}
	else {
	    return;
	}
    }

    // forget about 1.3 now

    v3_t newpos = add(add(scale(cam->rotaxis[0], dx),cam->pos),
		   add(scale(cam->rotaxis[1], dy),scale(cam->rotaxis[2], dz)));

    float p1[3]={cam->pos.x, cam->pos.y};
    float p2[3]={newpos.x, newpos.y};

    // 2. Tracemotion
    // Try sth. simple

    rayinfo_t ri;
    float t;
    if (ray2d(&ri, cam->sid, p1, p2, &t, NULL) == HIT_SOLID) {
	printf("hit!\n");
	return;
    }


/*
    if (cam->sid!=-1) {    // clip motion if the cam is currently in some sector
	float t;
	rayintx_t intx;
//	cam->sid=ray2d(cam->sid,p1,p2,&t,&intx); //TODO

	// update PVS
//	pvs2d();
    }
    else {

#if 0
	// test which sector it's in
	short sid;
	int loc=chkloc(cam->pos.x, cam->pos.y, &sid, 1); // TODO: user should be able to switch between overlapping sectors
	if (loc==4) { // loc==LOC_IN.   TODO: change this to a separate function
	    cam->sid=sid;
	    pvs2d();
	}

#endif

    }
*/

    if (ri.type == HIT_NOTHING) {
	cam->sid = ri.id;
    }

    cam->pos = newpos;
}



static void clampcamrot(cam1_t *cam, float du, float dv, float da){
    cam->u += du;
    cam->v += dv;
    cam->angle += da;

    if(cam->clamprot.u){
	if(cam->u > cam->rotmax.u){
	    cam->u = cam->rotmax.u;
	}
	else if(cam->u < cam->rotmin.u){
	    cam->u = cam->rotmin.u;
	}
    }
    else {
	if(cam->u > 360.0){
	    cam->u -= 360.0;
	}
	else if(cam->u < 0.0){
	    cam->u += 360.0;
	}
    }


    if(cam->clamprot.v){
	if(cam->v > cam->rotmax.v){
	    cam->v = cam->rotmax.v;
	}
	else if(cam->v < cam->rotmin.v){
	    cam->v = cam->rotmin.v;
	}
    }
    else {
	if(cam->v > 360.0){
	    cam->v -= 360.0;
	}
	else if(cam->v < 0.0){
	    cam->v += 360.0;
	}
    }


    if(cam->clamprot.a){
	if(cam->angle > cam->rotmax.a){
	    cam->angle = cam->rotmax.a;
	}
	else if(cam->angle < cam->rotmin.a){
	    cam->angle = cam->rotmin.a;
	}
    }
    else {
	if(cam->angle > 360.0){
	    cam->angle -= 360.0;
	}
	else if(cam->angle < 0.0){
	    cam->angle += 360.0;
	}
    }
}

void cam1_rot(cam1_t *cam, float du, float dv, float da){
    if (du==0 && dv==0 && da==0) {
	++cam->updated;
	return;
    }

    v3_t view, right, up;

    clampcamrot(cam, du, dv, da);

    /* rotate around z axis */
    view.x = COSD(cam->u);
    view.y = SIND(cam->u);
    view.z = 0.0;

    right = normalize(cross(vec(0.0f, 0.0f, 1.0f), view)); /* right handed coord system */

    /* rotate around x axis */
    view = scale(view, COSD(cam->v));
    view.z = SIND(cam->v);

    up = normalize(cross(view, right));

    /* rotate around z axis */
    right = add(scale(right, COSD(cam->angle)), scale(up, SIND(cam->angle)));
    up = normalize(cross(view, right));

    cam->axis = view;		/* useless */

    cam->rotaxis[0] = right;
    cam->rotaxis[1] = view;
    cam->rotaxis[2] = up;
}



// TODO 3.20 -- clean this mess up!
#if 0
static void cam1_updatepvs(cam1_t *cam)
{
    // find 
    short sid[100];
    int loc=chkloc(cam->pos.x, cam->pos.y, sid, 100);
    
    if (loc==LOC_IN) {
	pvs2d(&G.ed.pvs, sid[0], cam->pos, cam->rotaxis[1], M_PI_2);
	float t;
	float len=20;

	float dst[3]={cam->pos.x + len*cam->rotaxis[1].x,
		      cam->pos.y + len*cam->rotaxis[1].y,
		      cam->pos.z + len*cam->rotaxis[1].z};
//	int hit=ray2d(&G.ed.rayinfo, sid[0], &cam->pos, dst, &t, NULL);
/*
	if (hit!=HIT_NOTHING) {
	    printf("%s %i\n", hitstrings[hit], G.ed.rayinfo.id);
	}
*/
    }
}
#endif

static void cam1_updatepvs(cam1_t *cam)
{
    if (cam->sid != -1) {
	pvs2d(&G.ed.pvs, cam->sid, cam->pos, cam->rotaxis[1], M_PI_2);	
    }
}

static void cam1_update(cam1_t *cam)
{
    for (int i=0;i<3;i++) {
	float a=cam->acc[i];
	float v=cam->vel[i];
	// update velocity
	if (a==0) { // break
	    if (v!=0) {
		float f= cam->friction * (v>0 ? -1 : 1);
		cam->vel[i]+=f;
		if (cam->vel[i]*v<0) cam->vel[i]=0;
	    }
	}
	else {
	    v+=a;
	    v=v>cam->maxvel ? cam->maxvel : v;
	    v=v<-cam->maxvel ? -cam->maxvel : v;
	    cam->vel[i]=v;
	}

	// update angular velocity
	float aa=cam->angacc[i];
	float av=cam->angvel[i];
	if (aa==0) { // break
	    if (av!=0) {
		float f= cam->angfriction * (av>0 ? -1 : 1);
		cam->angvel[i]+=f;
		if (cam->angvel[i]*av<0) cam->angvel[i]=0;
	    }
	}
	else {
	    av+=aa;
	    av=av>cam->maxangvel ? cam->maxangvel : av;
	    av=av<-cam->maxangvel ? -cam->maxangvel : av;
	    cam->angvel[i]=av;
	}
    }

//    cam->updated=0;
    cam1_rot(cam, cam->angvel[0], cam->angvel[1], cam->angvel[2]);
    cam1_move(cam, cam->vel[0], cam->vel[1], cam->vel[2]);

    if (cam->pos.x!=cam->oldpos.x || cam->pos.y!=cam->oldpos.y || cam->angvel[0]!=0) {
	cam1_updatepvs(cam);
    }
    
    cam->oldpos=cam->pos;
}


void cam1_use(cam1_t *cam)
{
    cam1_update(cam);
    static float projmat[] = {
	1, 0, 0, 0,
	0, 1, 0, 0,
	0, 0, -0.98, -1,
	0, 0, -1.98, 0
    };

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
/*      gluPerspective(90, 4.0/3.0, 1, 500);    */
    glMultMatrixf(cam->proj); 
/*     glMultMatrixf(projmat);  */

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    updaterotmat(cam);

    glRotatef(-90,1,0,0);
    glMultMatrixf(cam->rotmat);
    glTranslatef(-cam->pos.x, -cam->pos.y, -cam->pos.z);

    float swapyz[16]={1,0,0,0,
		      0,0,1,0,
		      0,1,0,0,
		      0,0,0,1};

//    glMultMatrixf(swapyz);

}


static void computeprojmat(cam1_t *cam){
    float top;

    top = cam->h / cam->w;

    memset(cam->proj, 0, 16*sizeof(float));
    cam->proj[0] = cam->znear;
    cam->proj[5] = cam->znear / top;
    cam->proj[10] = -0.99;
    cam->proj[11] = -1;
    cam->proj[14] = -1.99 * cam->znear;

}

/*
  size = size of the monitor in inch
  aspect = aspect ratio of the monitor (i.e. w / h)
  dist = distance from eye to monitor
*/
void cam1_frommonitor(cam1_t *cam, float size, float aspect, float dist){
    float halfw;
    cam->clamprot.u = 0;
    cam->clamprot.v = 1;
    cam->clamprot.a = 0;
    cam->rotmin.v = -80;
    cam->rotmax.v = 80;

    cam->w = 640;
    cam->h = 480;

    halfw = aspect * sqrtf(size*size / (1.0 + aspect*aspect)) / 2.0;

    cam->znear = dist / halfw;
    cam->fov = RAD_TO_ANG(2.0 * atanf(halfw / dist));

    computeprojmat(cam);

    cam->u = 0;
    cam->v = 0;
    cam->angle = 0;
    cam->pos = vec(0, 0, 0);

    cam->rotmat[3] = 
	cam->rotmat[7] = 
	cam->rotmat[11] = 
	cam->rotmat[12] =
	cam->rotmat[13] =
	cam->rotmat[14] = 0.0;
    cam->rotmat[15] = 1.0;

    cam1_rot(cam, 90, 0, 0);

    return cam;
}


void cam1_init(cam1_t *cam,float w,float h,float fov)
{
    memset(cam,0,sizeof(*cam));
    cam->clamprot.u = 0;
    cam->clamprot.v = 1;
    cam->clamprot.a = 0;
    cam->rotmin.v = -80;
    cam->rotmax.v = 80;

    cam->w = w;
    cam->h = h;
    cam->znear = 1.0 / tanf(ANG_TO_RAD(fov / 2));
    cam->fov = fov;

    computeprojmat(cam);

    cam->u = 0;
    cam->v = 0;
    cam->angle = 0;
    cam->oldpos = cam->pos = vec(0, 0, 0);

    cam->rotmat[3] = 
	cam->rotmat[7] =
	cam->rotmat[11] =
	cam->rotmat[12] =
	cam->rotmat[13] =
	cam->rotmat[14] = 0.0f;
    cam->rotmat[15] = 1.0f;
    cam1_rot(cam, 90, 0, 0);
    cam->friction=cam->angfriction=0.1;
    cam->maxvel=2;
    cam->maxangvel=2;

    // the camera can start anywhere -- outside all sectors, or inside one or more sectors
    // you have to test to find the answer.

    cam->sid=-1;
    return;

}


void cam1_setres(cam1_t *cam, float w,float h){
    float s;

    cam->w = w;
    cam->h = h;

    computeprojmat(cam);
}

void cam1_setpos(cam1_t *cam,float x,float y,float z)
{
    cam->pos=vec(x,y,z);
}


/*
  camera controller:

  1. Camera should be aware of the environment
  it should be updated when:
  1) cam init -- which sector it's in?
  2) cam movement -- prev & current sector
  3) cam setpos -- same as init
  4) environment changes

  

 */


static v3_t pixel2camvec(cam1_t *cam, float xrel, float yrel){
    v3_t v;
    float centerx, centery;
    v3_t camvec;

    centerx = cam->w*0.5f;
    centery = cam->h*0.5f;

    v.x = (xrel - centerx) / centerx;
    v.y = (yrel - centery) / centerx;
    v.z = -cam->znear;

    camvec = add(scale(cam->rotaxis[0], v.x), scale(cam->rotaxis[1], v.y));
    camvec = add(camvec, scale(cam->rotaxis[2], v.z));

    return normalize(camvec); 
}




/*
  3rd person camera control
*/
#define SWIZZLE 1
/*
  spherical: (rho,theta,phi)
  to cardesian
*/
v3_t s2c(v3_t s)
{
    v3_t c;
    c.x = -s.x*sinf(s.y)*cosf(s.z);
    c.y = -s.x*sinf(s.z);
    c.z = s.x*cosf(s.y)*cosf(s.z);
    return c;
}

static v3_t c2s(v3_t c)
{
    v3_t s;
    s.x = length(c);
    s.y = atan2f(-c.x,c.z);
    //	s.y = -atan2f(c.x,c.z);
    s.z = asinf(-c.y/s.x);

    return s;
}


#define CAM3_CKSUM 0x582db30d

void cam3_update(char *nam, float du, float dt)
{
    cam3_t *cam=cvar_getp__(nam, CAM3_CKSUM);
    if (!cam) return;

#if 1

#if 1

    //	printvec("cam->strafe",&cam->strafe,3);

    float ang0 = atan2f(-(cam->pos.x-cam->target.x),(cam->pos.z-cam->target.z));
    cam->target = add(cam->target, cam->strafe);
    float ang1 = atan2f(-(cam->pos.x-cam->target.x),(cam->pos.z-cam->target.z));
    float strafeang = ang1-ang0;
    //	printf("strafeang=%f\n", strafeang);
    //	printf("ang0=%f ang1=%f\n", ang0, ang1);
#endif

    float K = 2.5;
    float anga = du*K - cam->angvelocity*cam->angKd;

    v3_t s = c2s(sub(cam->pos, cam->target));

    if(du){
	/* how many steps will cam->angvelocity converge to zero? */
	float angv = cam->angvelocity;
	float a=anga;
	float y=s.y;
	while(fabsf(angv)>0.0001){
	    angv += a*dt;
	    y += angv*dt;
	    a = -angv*cam->angKd;
	}
	cam->spherical_prediction = cam->spherical;
	cam->spherical_prediction.y = y;
    }

    //	printvec("prediction",&cam->spherical_prediction,3);


    cam->angvelocity += anga*dt;
    //	printf("s.y=%f\n", s.y);

#if 1
    s.y += cam->angvelocity*dt;
#else
    s.y += anga*dt;
#endif

    v3_t temp = s2c(s);

    //	printf("cam->spherical.y=%f\n", cam->spherical.y);
    //	printvec("s",&s,3);
    //	printvec("s2c",&temp,3);

    cam->pos = add(cam->target,s2c(s));


    float ang2 = atan2f(-(cam->pos.x-cam->target.x),cam->pos.z-cam->target.z);
    cam->spherical.y += ang2 - strafeang -ang0;
	


#if 0
    cam->spherical.y = atan2f(-(cam->pos.x-cam->target.x),
			      (cam->pos.z-cam->target.z));
    printf("cam->spherical.y=%f\n", cam->spherical.y);
    cam->spherical.y -= strafeang;
    printf("cam->angvelocity=%f\n", cam->angvelocity);
    printf("cam->spherical.y=%f\n", cam->spherical.y);
#endif

    v3_t restpos = add(cam->target, s2c(cam->spherical));
    v3_t d = sub(cam->pos, restpos);
    //	printf("d=%f %f %f\n", d.x,d.y,d.z);

    //printf("ks=%f, kd=%f\n", cam->Ks, cam->Kd);

    v3_t a = add(scale(d, -cam->Ks), scale(cam->velocity, -cam->Kd));

    //	printvec("cam->pos", &cam->pos,3);
    //	printvec("d",&d,3);


    cam->velocity = add(cam->velocity, scale(a, dt));
    cam->pos = add(cam->pos, scale(cam->velocity, dt));

    //	printvec("cam->view",&cam->view,3);


    if(cam->angvelocity){
	//		printf("cam->angv=%f\n",cam->angvelocity);
	cam->view = sub(cam->target, cam->pos);
    }
    else{
	//		printf("shit\n");
    }


    memset(&cam->strafe,0,sizeof(cam->strafe));

	
    //	printf("a=%f %f %f\n", a.x,a.y,a.z);
    //	printf("cam->pos=%f %f %f\n",cam->pos.x, cam->pos.y, cam->pos.z);
    //	printf("cam->target=%f %f %f\n", cam->target.x, cam->target.y, cam->target.z);
#endif


#if 0

    cam->oldcam = cam->newcam;
    cam->oldtarget = cam->newtarget;
#endif
}


void cam3_movetg(char *nam, float x, float y, float z)
{
    cam3_t *cam=cvar_getp__(nam, CAM3_CKSUM);
    if (!cam) return;

    //	cam->view = sub(cam->pos, cam->target);
    //	printf("view=%f %f %f\n", view.x, view.y, view.z);
    //	printvec("cam->pos",&cam->pos,3);
    //	printvec("cam->target",&cam->target, 3);


    v3_t forward = s2c(cam->spherical_prediction);
#if 0
    if(z){
	forward = sub(cam->pos, cam->target);
    }
    else{
	forward = cam->view;
    }
#endif

    forward.y = 0;
    forward = normalize(forward);
    //	printf("forward=%f %f %f\n", forward);

    v3_t right = cross(forward, vec(0,-1,0));
    v3_t offset = scale(forward, z);

    //	printvec("right",&right,3);

    offset = add(offset, scale(right, x));
    //	printf("offset=%f %f %f\n", offset.x, offset.y,offset.z);

    offset.y += y;
    //	cam->target = add(cam->target, offset);
    cam->strafe = add(cam->strafe, offset);

    //	printvec("strafe",&cam->strafe,3);

}



void cam3_use(char *nam)
{
    cam3_t *cam=cvar_getp__(nam, CAM3_CKSUM);
    if (!cam){
	cs_printf("Cannot find camera %s\n", nam);
	rnd_scrncoord();
	return;
    }

    static float projmat[] = {
	1, 0, 0, 0,
	0, 1, 0, 0,
	0, 0, -0.98, -1,
	0, 0, -1.98, 0
    };

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity( );
    gluPerspective( 80, 4.0/3.0, 1, 500 );    
    //  glMultMatrixf( cam->proj ); 
    /*     glMultMatrixf( projmat );  */

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity( );

    //  printvec("sph",&cam->spherical,3);
    v3_t target = add(cam->pos,scale(s2c(cam->spherical),-1));

    gluLookAt(cam->pos.x, cam->pos.y, cam->pos.z,
	      //		  cam->target.x, cam->target.y, cam->target.z,
	      //	 	  cam->pos.x, cam->target.y, cam->pos.z+cam->spherical.x+1,
	    
	      target.x,target.y,target.z,
	      0,1,0);

    //	glGetDoublev(GL_PROJECTION_MATRIX, cam->proj);
    //	glGetDoublev(GL_MODELVIEW_MATRIX, cam->modelview);
    //	glGetIntegerv(GL_VIEWPORT, cam->viewport);
}


void cam3_add(char *nam, v3_t target, v3_t pos, float Kd, float angKd)
{
    cam3_t *cam=cvar_getp__(nam, CAM3_CKSUM);
    if (!cam) {
	cam = malloc(sizeof(cam3_t));
	cvar_setp__(nam, cam, CAM3_CKSUM);
    }
    cam->target = target;
    cam->pos = pos;
    cam->spherical = cam->spherical_prediction = c2s(sub(pos,target));
    cam->d0 = cam->spherical.x;

    //	printf("sph0=%f %f %f\n", cam->spherical.x, cam->spherical.y, cam->spherical.z);

    cam->Kd = Kd;
    cam->Ks = Kd*Kd*.25;
    cam->Kdsmall = cam->Kd;
    cam->Kssmall = cam->Ks;
    cam->Kdlarge = 4*Kd;
    cam->Kslarge = cam->Kdlarge*Kd;
	

    //	printf("ks=%f, kd=%f\n", cam->Ks, Kd);
    cam->velocity = vec(0,0,0);
    cam->angvelocity = 0;
    cam->angKd = angKd;
    cam->view = sub(cam->target, cam->pos);
    cam->strafe = vec(0,0,0);
#if 0
    cam->angle = 0;
    cam->targetaxis[0] = normalize(cross(sub(pos, target), vec(0,1,0)));
    cam->targetaxis[1] = vec(0,1,0);
    cam->targetaxis[2] = cross(cam->targetaxis[0], cam->targetaxis[1]);
#endif
    return cam;
}


static void cam2_use(float x, float y, float zoom) 
{
    glViewport(0,0,G.screenw,G.screenh);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(G.screenw*.5, G.screenh*.5, 0);
    glScalef(zoom,zoom,1);
    glTranslatef(-x,-y,0);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //	glOrtho(0,G.screenw, 0, G.screenh, -1, 1);	
    gluOrtho2D(0,G.screenw, 0,G.screenh);
}



///////////////////////////////////////////////////////////////////////////////////////////

enum {FILTER_BUFFER, FILTER_SHADER, FILTER_IN, FILTER_OUT};

#define isnum(ch) ((ch) >= '0' && (ch) <= '9')
#define getcha() (c=*(++lex->p))

/*
  skip white space and comments (c-style)
*/
void lex_skipwhite()
{
    lex_t *lex=G.lex_cur;
    if (!lex) return;

    char *c=lex->p;
    while (1) {
	while(isspace(*c)) ++c;
		
	if (c[0]=='/' && c[1]=='*') c+=2;
	else {
	    break;;
	}
		
	while (*c && (c[0]!='*' || c[1]!='/')) ++c;

	if (!*c) {
	    break;;
	}
		
	c+=2;
    }

    lex->p=c;
}

int lex_gettok() {
    int numdots;
    int i;
    char c0;

    lex_t *lex=G.lex_cur;
    if (!lex) return TOK_END;
    lex_skipwhite(lex);

    char c=*lex->p;
    if (c==0) {
	return lex->tok=TOK_END;
    }

    char buf[32];
    char *b = buf;
    numdots = 0;
    if (isnum(c)||c=='.') { /* numbers */
	do {
	    *b++ = c;
	    if (c== '.') {
		numdots++;
	    }
	    getcha();
	} while (isnum(c)||c=='.');
	*b = 0;
	if (numdots == 0) {
	    lex->tok = TOK_INT_CONST;
	    char *end;
	    lex->sem.i=strtol(buf,&end,0);
	    if (end==buf) { // error
		lex_error();
	    }
	}
	else if (numdots == 1) {
	    lex->tok = TOK_FLOAT_CONST;
	    char *end;
	    lex->sem.f=strtof(buf,&end);
	    if (end=buf) { // error
		lex_error();
	    }
	}
	else {
	    lex_error();
	}
    }
    else if (isalpha(c)||c=='_') {	/* identifier or keyword */
	do {
	    *b++ = c;
	    getcha();
	} while (isalnum(c)||c=='_');
	*b = 0;

	for (i=0; i<lex->nkw; ++i) {
	    if (!strcmp(buf, lex->kw[i].str)) {
		lex->tok=lex->kw[i].tok;
		break;
	    }
	}
	if (i==lex->nkw) {
	    lex->tok= TOK_ID;
	    strcpy(lex->sem.id, buf);
	}
    }
    else if (c == '\"') {
	do {
	    getcha();
	    if (c == '\\') {
		getcha();
		switch (c) {
		case 'f':
		    c='\f';
		    break;
		case 'n':
		    c='\n';
		    break;
		case 'r':
		    c='\r';
		    break;
		case 't':
		    c='\t';
		    break;
		case 'v':
		    c='\v';
		    break;
		case '\\':
		    c='\\';
		    break;
		case '\"':
		    c='\"';
		    break;
		default:
		    lex_error();
		    break;
		}
	    }
	    *b++ = c;
	} while (c != '\"');
	getcha();
	b[-1] = 0;
	lex->tok= TOK_STRING_CONST;
	strcpy(lex->sem.str,buf);
    }
    else {
	switch (c) {
	case '(':
	case ')':
	case '{':
	case '}':
	case '+':
	case '-':
	case '*':
	case '/':
	case '%':
	case ';':
	case ',':
	case '[':
	case ']':
	    lex->tok= c;
	    getcha();
	    break;
	case '>':
	case '<':
	case '=':
	case '!':
	    c0 = c;
	    getcha();
	    if (c == '=') {
		switch (c0) {
		case '>':
		    lex->tok= TOK_GE;
		    break;
		case '<':
		    lex->tok= TOK_LE;
		    break;
		case '=':
		    lex->tok= TOK_EQ;
		    break;
		case '!':
		    lex->tok= TOK_NEQ;
		    break;
		}
		getcha();
	    }
	    else {
		lex->tok= c0;
	    }
	    break;
	default:
	    lex_error();
	    break;
	}
      
    }
    return lex->tok;
}


void lex_expect(int tok)
{
    lex_t *lex=G.lex_cur;
    if (lex->tok!=tok) lex_error();
    lex_gettok();
}

static jmp_buf lex_jmpbuf;


void lex_bgn(lex_t *lex, char *code, kw_t *kw, ...)
{
    va_list va;

    if (!lex || !code){
	lex->buf=lex->p=lex->kw=NULL;
	return;
    }

    lex->buf=lex->p=code;
    lex->kw=kw;
	
    va_start(va,kw);
    if (kw) lex->nkw=va_arg(va,int);
    else lex->nkw=0;

    lex->tok=0;

    G.lex_cur=lex;

    lex_end();	
}

void lex_end()
{
    int a=setjmp(lex_jmpbuf);
    if (!a) return;

    // lex_err() will goto here
    G.lex_cur=NULL;

}

void lex_error()
{
    longjmp(lex_jmpbuf,1);
}

int lex_getint()
{
    lex_t *lex=G.lex_cur;
    if (lex->tok==TOK_INT_CONST) {
	int i=lex->sem.i;
	lex_gettok();
	return i;
    }
    else if (lex->tok==TOK_FLOAT_CONST) {
	float f=lex->sem.f;
	lex_gettok();
	return f;
    }
    else{
	lex_error();
    }
}

void lex_getstr(char *str, int n)
{
    lex_t *lex=G.lex_cur;
    if (lex->tok==TOK_STRING_CONST) {
	strncpy(str, lex->sem.str, n);
	lex_gettok();
	return;
    }
	
    lex_error();
}



/*
  The only API that's exposed to the user is filter_use().

  However there're several commands.

  Plus you can change filter files and reload at runtime.

  If a filter which isn't exist neither in RAM nor on disk, it'll have zero input/output.
  You can add input/output by command:

  filter_addinput "color"
  filter_addinput "depth"
  filter_addoutput "color".

  Or simply write a filter file which only have in{}/out{} nodes

  input vars are stored in filter itself

  
  If in editing mode, 
*/
void rnd_setz(float z)
{
    G.cur_z=z;
}

float rnd_getz()
{
    return G.cur_z;
}

static void rnd_bezier(float x1, float y1, float x2, float y2)
{
    float z=G.cur_z;
    float d=fabs((x1-x2)*.5);

    float p[4][3]={
	{x1,y1,z},
	{x1+d,y1,z},
	{x2-d,y2,z},
	{x2,y2,z}
    };

    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 4, p);
    glEnable(GL_MAP1_VERTEX_3);

    glBegin(GL_LINE_STRIP);
    for (int i = 0; i <= 30; i++) 
	glEvalCoord1f((GLfloat) i/30.0);
    glEnd();

    glDisable(GL_MAP1_VERTEX_3);

    glDisable(GL_LINE_SMOOTH);
}


static void texclamp()
{
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
}


enum {CNR_LB=1, 
      CNR_RB=2, 
      CNR_RT=4, 
      CNR_LT=8,
      CNR_ALL=15};

static void rnd_roundrect(float bnd[4], float r, int s, int corner)
{
    float z=G.cur_z;
    float cx=(bnd[0]+bnd[2])*.5;
    float cy=(bnd[1]+bnd[3])*.5;
    glBegin(GL_TRIANGLE_FAN);
    // center
    glVertex3f(cx,cy,z);

    // LB
    if (corner & CNR_LB) {
	for (int i=0; i<=s; i++ ){
	    float angle=M_PI + i*M_PI_2/s;
	    float x=r*cosf(angle) + bnd[0]+r;
	    float y=r*sinf(angle) + bnd[1]+r;
	    glVertex3f(x,y,z);
	}
    }
    else {
	glVertex3f(bnd[0],bnd[1]+r,z);
	glVertex3f(bnd[0],bnd[1],z);
	glVertex3f(bnd[0]+r,bnd[1],z);
    }

    // RB
    if (corner & CNR_RB) {
	for (int i=0; i<=s; i++ ){
	    float angle=M_PI*1.5f + i*M_PI_2/s;
	    float x=r*cosf(angle) + bnd[2]-r;
	    float y=r*sinf(angle) + bnd[1]+r;
	    glVertex3f(x,y,z);
	}
    }
    else {
	glVertex3f(bnd[2]-r,bnd[1],z);
	glVertex3f(bnd[2],bnd[1],z);
	glVertex3f(bnd[2],bnd[1]+r,z);
    }

    // RT
    if (corner & CNR_RT) {
	for (int i=0; i<=s; i++ ){
	    float angle=i*M_PI_2/s;
	    float x=r*cosf(angle) + bnd[2]-r;
	    float y=r*sinf(angle) + bnd[3]-r;
	    glVertex3f(x,y,z);
	}
    }
    else {
	glVertex3f(bnd[2],bnd[3]-r,z);
	glVertex3f(bnd[2],bnd[3],z);
	glVertex3f(bnd[2]-r,bnd[3],z);
    }

    // LT
    if (corner & CNR_LT) {
	for (int i=0; i<=s; i++ ){
	    float angle=M_PI_2 + i*M_PI_2/s;
	    float x=r*cosf(angle) + bnd[0]+r;
	    float y=r*sinf(angle) + bnd[3]-r;
	    glVertex3f(x,y,z);
	}
    }
    else {
	glVertex3f(bnd[0]+r,bnd[3],z);
	glVertex3f(bnd[0],bnd[3],z);
	glVertex3f(bnd[0],bnd[3]-r,z);
    }

    glVertex3f(bnd[0], bnd[1]+r,z);

    glEnd();
}


static void rnd_roundtexrect(uint tex, float bnd[4], float r, int s, int corner)
{
    float z=G.cur_z;
    float cx=(bnd[0]+bnd[2])*.5;
    float cy=(bnd[1]+bnd[3])*.5;

    glColor4f(1,1,1,1);
    glEnable(GL_TEXTURE_2D); //TODO
    glBindTexture(GL_TEXTURE_2D, tex);

    glBegin(GL_TRIANGLE_FAN);
    // center
    glTexCoord2f(.5,.5);
    glVertex3f(cx,cy,z);

    float trx= r/(bnd[2]-bnd[0]);
    float try= r/(bnd[3]-bnd[1]);

    // LB
    if (corner & CNR_LB) {
	for (int i=0; i<=s; i++ ){
	    float angle=M_PI + i*M_PI_2/s;
	    float x=r*cosf(angle) + bnd[0]+r;
	    float y=r*sinf(angle) + bnd[1]+r;
	    float tx=trx*cosf(angle) + trx;
	    float ty=try*sinf(angle) + try;
	    glTexCoord2f(tx,ty);
	    glVertex3f(x,y,z);
	}
    }
    else {
	glTexCoord2f(0,try);
	glVertex3f(bnd[0],bnd[1]+r,z);
	glTexCoord2f(0,0);
	glVertex3f(bnd[0],bnd[1],z);
	glTexCoord2f(trx,0);
	glVertex3f(bnd[0]+r,bnd[1],z);
    }

    // RB
    if (corner & CNR_RB) {
	for (int i=0; i<=s; i++ ){
	    float angle=M_PI*1.5f + i*M_PI_2/s;
	    float x=r*cosf(angle) + bnd[2]-r;
	    float y=r*sinf(angle) + bnd[1]+r;
	    float tx=trx*cosf(angle) + 1-trx;
	    float ty=try*sinf(angle) + try;
	    glTexCoord2f(tx,ty);
	    glVertex3f(x,y,z);
	}
    }
    else {
	glTexCoord2f(1-trx,0);
	glVertex3f(bnd[2]-r,bnd[1],z);
	glTexCoord2f(1,0);		
	glVertex3f(bnd[2],bnd[1],z);
	glTexCoord2f(1,try);
	glVertex3f(bnd[2],bnd[1]+r,z);
    }

    // RT
    if (corner & CNR_RT) {
	for (int i=0; i<=s; i++ ){
	    float angle=i*M_PI_2/s;
	    float x=r*cosf(angle) + bnd[2]-r;
	    float y=r*sinf(angle) + bnd[3]-r;
	    float tx=trx*cosf(angle) + 1-trx;
	    float ty=try*sinf(angle) + 1-try;
	    glTexCoord2f(tx,ty);
	    glVertex3f(x,y,z);
	}
    }
    else {
	glTexCoord2f(1,1-try);
	glVertex3f(bnd[2],bnd[3]-r,z);
	glTexCoord2f(1,1);
	glVertex3f(bnd[2],bnd[3],z);
	glTexCoord2f(1-trx,1);
	glVertex3f(bnd[2]-r,bnd[3],z);
    }

    // LT
    if (corner & CNR_LT) {
	for (int i=0; i<=s; i++ ){
	    float angle=M_PI_2 + i*M_PI_2/s;
	    float x=r*cosf(angle) + bnd[0]+r;
	    float y=r*sinf(angle) + bnd[3]-r;
	    float tx=trx*cosf(angle) + trx;
	    float ty=try*sinf(angle) + 1-try;
	    glTexCoord2f(tx,ty);			
	    glVertex3f(x,y,z);
	}
    }
    else {
	glTexCoord2f(trx,1);
	glVertex3f(bnd[0]+r,bnd[3],z);
	glTexCoord2f(0,1);
	glVertex3f(bnd[0],bnd[3],z);
	glTexCoord2f(0,1-try);
	glVertex3f(bnd[0],bnd[3]-r,z);
    }

    glTexCoord2f(0,try);
    glVertex3f(bnd[0], bnd[1]+r,z);

    glEnd();

    glDisable(GL_TEXTURE_2D);
}


/*
  Solid circle
*/
static void rnd_circle1(float cx, float cy, float r, int s)
{
    float z=G.cur_z;
    glBegin(GL_TRIANGLE_FAN);
    glVertex3f(cx,cy,z);
    for (int i=0; i<=s; ++i) {
	float angle=2*M_PI*i/s;
	float x=cx + r*cosf(angle);
	float y=cy + r*sinf(angle);
	glVertex3f(x,y,z);
    }
    glEnd();
}


/*
  Hollow circle
*/
static void rnd_circle0(float cx, float cy, float r, int s)
{
    float z=G.cur_z;
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    glBegin(GL_LINE_LOOP);
    for (int i=0; i<s; ++i) {
	float angle=2*M_PI*i/s;
	float x=cx + r*cosf(angle);
	float y=cy + r*sinf(angle);
	glVertex3f(x,y,z);
    }
    glEnd();

    glDisable(GL_LINE_SMOOTH);
}

/*
  Solid square
*/
static void rnd_square1(float cx, float cy, float r)
{
    float z=G.cur_z;
    glBegin(GL_QUADS);
    glVertex3f(cx-r,cy-r,z);
    glVertex3f(cx+r,cy-r,z);
    glVertex3f(cx+r,cy+r,z);
    glVertex3f(cx-r,cy+r,z);
    glEnd();
}

/*
  Hollow square
*/
static void rnd_square0(float cx, float cy, float r)
{
    float z=G.cur_z;
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    glBegin(GL_LINE_LOOP);
    glVertex3f(cx-r,cy-r,z);
    glVertex3f(cx+r,cy-r,z);
    glVertex3f(cx+r,cy+r,z);
    glVertex3f(cx-r,cy+r,z);
    glEnd();

    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_BLEND);
}



static void rnd_rectshadow(float bnd0[4], float zoom)
{
    float z=G.cur_z;
    if (zoom<0.4) zoom=0.4;
    float w=20.0/zoom;
    float d=7;
    float a=0.2;

    float bnd[4]={bnd0[0],bnd0[1],bnd0[2],bnd0[3]};

    bnd[0]+=w;
    bnd[1]+=w;
    bnd[2]-=w;
    bnd[3]-=w;

    bnd[1]-=0.5*w;
    bnd[3]-=0.5*w;

    w*=2.5;

    static uint fallofftex=0;
    if (!fallofftex) {
	// sine off
	uchar img[256];
	for (int i=0; i<256; ++i) {
	    img[i]=sinf(M_PI_2*i/255)*255;
	    printf("img[i]=%i\n", img[i]);
	}

	fallofftex=tex1b(256,1,011,img);
	texclamp();
    }

    vfs_use(-1);

    glColor4f(0,0,0,1);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, fallofftex);

    // sides
    glBegin(GL_QUADS);

    // left
    glTexCoord1f(0.0); glVertex3f(bnd[0]-w, bnd[3],z);
    glTexCoord1f(0.0); glVertex3f(bnd[0]-w, bnd[1],z);
    glTexCoord1f(a); glVertex3f(bnd[0], bnd[1],z);
    glTexCoord1f(a); glVertex3f(bnd[0], bnd[3],z);

    // right
    glTexCoord1f(0.0); glVertex3f(bnd[2]+w, bnd[3],z);
    glTexCoord1f(0.0); glVertex3f(bnd[2]+w, bnd[1],z);
    glTexCoord1f(a); glVertex3f(bnd[2], bnd[1],z);
    glTexCoord1f(a); glVertex3f(bnd[2], bnd[3],z);
	
    // bottom
    glTexCoord1f(0.0); glVertex3f(bnd[0], bnd[1]-w,z);
    glTexCoord1f(0.0); glVertex3f(bnd[2], bnd[1]-w,z);
    glTexCoord1f(a); glVertex3f(bnd[2], bnd[1],z);
    glTexCoord1f(a); glVertex3f(bnd[0], bnd[1],z);

    // top
    glTexCoord1f(0.0); glVertex3f(bnd[0], bnd[3]+w,z);
    glTexCoord1f(0.0); glVertex3f(bnd[2], bnd[3]+w,z);
    glTexCoord1f(a); glVertex3f(bnd[2], bnd[3],z);
    glTexCoord1f(a); glVertex3f(bnd[0], bnd[3],z);
	
    glEnd();

    // corners

    // 0,0
    glBegin(GL_TRIANGLE_FAN);
    glTexCoord1f(a); glVertex3f(bnd[0], bnd[1],z);
    glTexCoord1f(0.0);
    for (int i=0; i<=10; i++) {
	float angle=M_PI + i*M_PI_2/10;
	float x=w*cosf(angle) + bnd[0];
	float y=w*sinf(angle) + bnd[1];
	glVertex3f(x,y,z);
    }
    glEnd();

    // 1,0
    glBegin(GL_TRIANGLE_FAN);
    glTexCoord1f(a); glVertex3f(bnd[2], bnd[1],z);
    glTexCoord1f(0.0);
    for (int i=0; i<=10; i++) {
	float angle=M_PI*1.5 + i*M_PI_2/10;
	float x=w*cosf(angle) + bnd[2];
	float y=w*sinf(angle) + bnd[1];
	glVertex3f(x,y,z);
    }
    glEnd();

    // 1,1
    glBegin(GL_TRIANGLE_FAN);
    glTexCoord1f(a); glVertex3f(bnd[2], bnd[3],z);
    glTexCoord1f(0.0);
    for (int i=0; i<=10; i++) {
	float angle=i*M_PI_2/10;
	float x=w*cosf(angle) + bnd[2];
	float y=w*sinf(angle) + bnd[3];
	glVertex3f(x,y,z);
    }
    glEnd();

    // 0,1
    glBegin(GL_TRIANGLE_FAN);
    glTexCoord1f(a); glVertex3f(bnd[0], bnd[3],z);
    glTexCoord1f(0.0);
    for (int i=0; i<=10; i++) {
	float angle=M_PI_2 + i*M_PI_2/10;
	float x=w*cosf(angle) + bnd[0];
	float y=w*sinf(angle) + bnd[3];
	glVertex3f(x,y,z);
    }
    glEnd();

    glDisable(GL_TEXTURE_2D);

}

static void filter_drawbadnode(filter_t *ft, int id)
{
    filternode_t *node=ft->nodes + id;
    float bnd[4]={node->disp.x, node->disp.y, node->disp.x+node->disp.w, node->disp.y+node->disp.h};
    if (node->type==FILTER_SHADER) {
	if (node->shader.type==SHADER_DS) {
	    glColor3f(.5,.5,.5);
	    rnd_roundrect(bnd, 15, 10, CNR_ALL);
	    glColor3f(0,0,0);
	    float texbnd[4]={bnd[0]+50, bnd[1]+80, bnd[2], bnd[1]+120};
	    rnd_text(texbnd, "1 / 2");

	    float y=(bnd[1]+bnd[3])*.5;
	    glColor3f(1,1,0);
	    rnd_circle1(bnd[0], y, 15, 20);
	    glColor3f(0.3,0.3,0.3);
	    rnd_circle0(bnd[0], y, 15, 20);

	    /*
	      glColor3f(0.6,0.6,0.6);
	      rnd_square1(bnd[2], y, 13.3);
	      glColor3f(0.3,0.3,0.3);
	      rnd_square0(bnd[2], y, 13.3);
	    */
	}
	else if (node->shader.type==SHADER_PACK) {
	    glColor3f(0.5,0.5,0.5);
	    rnd_roundrect(bnd, 15, 10, CNR_ALL);
	    glColor3f(0,0,0);
	    float texbnd[4]={bnd[0]+20, bnd[1]+80, bnd[2], bnd[1]+120};
	    rnd_text(texbnd, "%s", "+");
	}
	else {
	    printf("hahaha\n");
	}
    }
    else if (node->type==FILTER_BUFFER) {
    }
}

/*
  Draw shadows in depth order, disable depth write
*/
static void filter_drawinputshadow(filternode_t *node, float zoom)
{
    float tcbnd[4]={0,0,1,1};
    float qbnd[4];
    qbnd[0]=node->disp.x;
    qbnd[1]=node->disp.y;
    qbnd[2]=qbnd[0]+ node->disp.w;//(G.screenw>>node->buffer.level);
    qbnd[3]=qbnd[1]+ node->disp.h;//(G.screenh>>node->buffer.level);

    float qshadow[4]={qbnd[0], qbnd[1], qbnd[2], qbnd[3]};
    rnd_rectshadow(qshadow, zoom);
}


/*
  Draw opaque images in arbitrary order
*/
static void filter_drawinput(filternode_t *node, char *nam, float zoom)
{
    float tcbnd[4]={0,0,1,1};
    float qbnd[4];
    qbnd[0]=node->disp.x;
    qbnd[1]=node->disp.y;
    qbnd[2]=qbnd[0]+ (G.screenw>>node->buffer.level);
    qbnd[3]=qbnd[1]+ (G.screenh>>node->buffer.level);
	

    //	fnt_printf(qbnd[0], qbnd[3], 0, "%s\n", nam);

    //	fnt_puts();
	
    float bw=30;
	
    // shadow

#if 0
    {
	float qshadow[4]={qbnd[0], qbnd[1], qbnd[2], qbnd[3]+40};
	rnd_rectshadow(qshadow, zoom);
    }
#endif

    // title bar
    {
	glColor3f(.5,.51,.55);
	float qtitle[4]={qbnd[0], qbnd[3], qbnd[2], qbnd[3]+40};
	rnd_roundrect(qtitle, 15, 10, CNR_LT|CNR_RT);
    }

    // title text
    {
	glColor4f(1,1,1,1);
	float textbnd[4]={qbnd[0]+10, qbnd[3], qbnd[2], qbnd[3]+35}; // TODO: don't hard code
	rnd_text(textbnd, nam);
    }

    // texture
    {
	glColor3f(0,0,0);
	rnd_roundtexrect(node->buffer.tex, qbnd, 15, 10, CNR_LB|CNR_RB);
    }


    // output point
    {
	glColor3f(1,1,0);
	rnd_circle1(node->disp.x+node->disp.w, node->disp.y+node->disp.h*.5, 15, 20);
	glColor3f(0.3,0.3,0.3);
	rnd_circle0(node->disp.x+node->disp.w, node->disp.y+node->disp.h*.5, 15, 20);
    }

}

static void filter_drawoutput()
{
}


static void filter_drawbuffer(filter_t *ft, int id)
{
    filternode_t *b=ft->nodes+id;
    float qbnd[4];
    qbnd[0]=b->disp.x;
    qbnd[1]=b->disp.y;
    qbnd[2]=qbnd[0]+ (G.screenw>>b->buffer.level);
    qbnd[3]=qbnd[1]+ (G.screenh>>b->buffer.level);
	
    float z=(float)b->disp.rank/ft->disp.maxrank;
    rnd_setz(z);
    rnd_roundtexrect(b->buffer.tex, qbnd, 15, 10, CNR_ALL);

    glColor3f(1,1,0);
    rnd_circle1(b->disp.x+b->disp.w, b->disp.y+b->disp.h*.5, 15, 20);
    glColor3f(0.3,0.3,0.3);
    rnd_circle0(b->disp.x+b->disp.w, b->disp.y+b->disp.h*.5, 15, 20);
}

static void filter_drawshader(filter_t *ft, int id)
{
}

static void filter_drawds()
{
}

static void filter_drawpack()
{
}


/*
  g -- the desired distance (in pixels) between two lines on screen
*/
void rnd_grid(cam2_t *cam, int g)
{
    float s= g*cam->zoom; // if zoom=100, g=10, then s=10/100=0.1, which means two lines with 10 unit apart with be only
    // 0.1 pixels apart. Obviously we can't use 0.1 as the distance between lines.
    // we have to enlarge the line distance to the (10,100) range
    float sz=g;		// real (enlarged) distance between two lines
    int g2=g*g;
    // find suitable grid size
    if (s<g) {
	while (s<g) {
	    s*=g;
	    sz*=g;
	}
		
    }
    else if (s>g2) {
	while (s>g2) {
	    s/=g;
	    sz*=g;
	}
    }

    // clamp grid lines
    float xmin,xmax,ymin,ymax;
    float hw=G.screenw/cam->zoom*.5; // half width
    float hh=G.screenh/cam->zoom*.5;
    xmin=cam->x-hw;
    xmax=cam->x+hw;
    ymin=cam->y-hh;
    ymax=cam->y+hh;

    xmin = ((int)(xmin/sz)-1)*sz;
    xmax = ((int)(xmax/sz)+1)*sz;
    ymin = ((int)(ymin/sz)-1)*sz;
    ymax = ((int)(ymax/sz)+1)*sz;

    float c=.5-s/g2*.2;
    glColor3f(c,c,c);
    glBegin(GL_LINES);
    for (float x=xmin; x<=xmax; x+=sz) {
	glVertex2f(x,ymin);
	glVertex2f(x,ymax);
    }

    for (float y=ymin; y<=ymax; y+=sz) {
	glVertex2f(xmin,y);
	glVertex2f(xmax,y);
    }
    glEnd();

    glColor3f(.3,.3,.3);
    glBegin(GL_LINES);
    for (float x=xmin; x<=xmax; x+=sz) {
	if ((int)x % ((int)sz*g) == 0) {
	    glVertex2f(x,ymin);
	    glVertex2f(x,ymax);
	}
    }

    for (float y=ymin; y<=ymax; y+=sz) {
	if ((int)y % ((int)sz*g) == 0) {
	    glVertex2f(xmin,y);
	    glVertex2f(xmax,y);
	}
    }
    glEnd();

}

#define WORLD_MAX 16*1024
void rnd_simplegrid(cam2_t *cam, int level)
{
    // clip in world space
    float xmin,xmax,ymin,ymax;
    float hw=G.screenw/cam->zoom*.5;
    float hh=G.screenh/cam->zoom*.5;
    xmin=cam->x - hw;
    if (xmin > WORLD_MAX) return;
    if (xmin < -WORLD_MAX) xmin=-WORLD_MAX;

    ymin=cam->y - hh;
    if (ymin > WORLD_MAX) return;
    if (ymin < -WORLD_MAX) ymin=-WORLD_MAX;

    xmax=cam->x + hw;
    if (xmax < -WORLD_MAX) return;
    if (xmax > WORLD_MAX) xmax=WORLD_MAX;
    
    ymax=cam->y + hh;
    if (ymax < -WORLD_MAX) return;
    if (ymax > WORLD_MAX) ymax=WORLD_MAX;

    float d=(1<<level);
    xmin=((int)(xmin/d)-1)*d;
    xmax=((int)(xmax/d)+1)*d;
    ymin=((int)(ymin/d)-1)*d;
    ymax=((int)(ymax/d)+1)*d;
    
    glColor3f(.5,.5,.5);
    if ((xmax-xmin)/d>G.screenw) { // lines are too dense => just draw a quad
	glBegin(GL_QUADS);
	glVertex2f(xmin,ymin);
	glVertex2f(xmax,ymin);
	glVertex2f(xmax,ymax);
	glVertex2f(xmin,ymax);
	glEnd();
	printf("yeah %f\n", (xmax-xmin)/d);
    }
    else {
	glBegin(GL_LINES);
	for (float x=xmin; x<=xmax; x+=d) {
	    glVertex2f(x,ymin);
	    glVertex2f(x,ymax);
	}
	
	for (float y=ymin; y<=ymax; y+=d) {
	    glVertex2f(xmin,y);
	    glVertex2f(xmax,y);
	}
	glEnd();
    }

    glBegin(GL_LINES);
    glColor3f(1,0,0);
    glVertex2f(-WORLD_MAX,0);
    glVertex2f(WORLD_MAX,0);
    glColor3f(0,1,0);
    glVertex2f(0,-WORLD_MAX);
    glVertex2f(0,WORLD_MAX);
    glEnd();
}


static int cmprank(void *a, void *b)
{
    struct {
	int rank;
	int id;
    } *r1, *r2;

    r1=a;
    r2=b;

    return r1->rank-r2->rank;
}

void filter_use(const char *nam)
{
    if (!nam || !*nam) return;

    filter_t *ft=cvar_getp__(nam, FILTER_CKSUM);
    if (!ft) {
	ft=filter_load(nam, NULL);
	if (!ft) {
	    // create a new filter
	    ft=malloc(sizeof(filter_t));
	    ft->nam=strdup(nam);
	    ft->nodes=0;
	    ft->validshaders=0;
	    ft->ns=0;
	    ft->in=ft->out=0;
	    ft->disp.zoom=.5; // TODO
	    cvar_setp__(nam, ft, FILTER_CKSUM);
	    return;
	}

	cvar_setp__(nam, ft, FILTER_CKSUM);
    }

    int edit=0;
    G.mode=MODE_FILTER;
    G.filter_curnam=nam;
    G.filter_cur=ft; // TODO: just testing
    extern float zoom;
    extern float filterx, filtery;

    if (G.mode==MODE_FILTER && !strcmp(nam, G.filter_curnam)) {
	rnd_tg(0);
	//		rnd_clearall();
	glClearColor(.5,.5,.5,1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,1);

	cam2_use(ft->disp.x, ft->disp.y, ft->disp.zoom);
	//		cam2_use(filterx, filtery, zoom);
	edit=1;

	//		filter_drawbg();

	cam2_t cam={ft->disp.x,ft->disp.y,ft->disp.zoom};
	rnd_grid(&cam,10);
//		rnd_grid(&cam, 8);
//		rnd_simplegrid(&cam,8);
    }

    //	rnd_scrncoord();

    // loop over valid shaders
    int sz=arr_len(ft->nodes);
    if (!sz) return; // TODO: draw some info?

    int outdeg[sz]; // used to recycle RTs -- RT with zero outdeg will be freed if not locked
    int nvalids=0;
    // init outdeg[]

    for (int i=0; i<sz; i++) {
	filternode_t *node=ft->nodes+i;
	outdeg[i]=node->outdeg;
	if (node->indeg!=-1){
	    //			if (node->type==FILTER_IN){  // TODO
	    ++nvalids;
	    //			}
	}
    }

    int drawn[sz]; // array to store whether a node is drawn
    memset(drawn,0,sizeof(drawn));

    // each input buffer has a name and an ID
    // init input nodes
    for (int i=0; i<arr_len(ft->in); ++i) {
	float *f=cvec_get(ft->in[i].nam);
	if (f) {
	    int tex=f[0];
	    rtinfo_t info;
	    int err=rt_getinfo(tex, &info);
	    if (!err) {
		filternode_t *node=ft->nodes + ft->in[i].id;
		assert(node->indeg!=-1);
		node->buffer.level=info.level;
		node->buffer.tex=tex;
		node->buffer.type=info.type;

		// if in edit mode, draw the input texture
		if (edit) {
		    //					printf("%s\n", ft->in[i].nam);
		    glEnable(GL_DEPTH_TEST);
		    glEnable(GL_BLEND);
		    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		    glDepthFunc(GL_LEQUAL);
		    float z=(float)node->disp.rank/ft->disp.maxrank;
		    rnd_setz(z);
		    filter_drawinput(node, ft->in[i].nam, ft->disp.zoom);
		    glDisable(GL_DEPTH_TEST);
		    glDisable(GL_BLEND);

		    drawn[ft->in[i].id]=1;
		}
	    }
	}
    }
	

    // traverse valid shaders
    int ns=arr_len(ft->validshaders);
    for (int i=0; i<ns; ++i) {
	int sid=ft->validshaders[i];
	filternode_t *s= ft->nodes + sid;
	uint intex;
		
	if (s->type==FILTER_OUT) {
	    //			assert(s->shader.ni==1);
	    if (s->shader.ni!=1) continue;
	    filternode_t *b=ft->nodes + s->inout[0];
	    int j;
	    for (j=0; j<arr_len(ft->out); ++j) { // find the outbuffer's texture name
		if (ft->out[j].id==sid) {
		    float *tex=cvec_get(ft->out[j].nam);
		    if (tex) {
			tex[0]=b->buffer.tex;
			b->buffer.tex=0;
			break;
		    }
		}
	    }

	    drawn[sid]=1;
	    continue;
	}


	if (s->shader.type==SHADER_PACK) {
	    // don't draw pack nodes until the end
	    continue; // ??
	}

	/*
	  if (!s->outdeg) { // output the texture ??
	  assert(s->indeg==1);
	  filternode_t *b=ft->nodes + s->inout[0];
	  continue;
	  }
	*/

	if (s->shader.type==SHADER_DS) {
	    assert(s->shader.in);
	    int bid=s->shader.in->b;
	    assert(bid!=-1);
	    intex=ft->nodes[bid].buffer.tex;
	    int type=rt_fmt(intex);
	    uint p = type==RT_DEPTH ? s->shader.p_dsdepth : s->shader.p_dscolor;
	    vfs_use(p);
	}
	else {
	    vfs_use(s->shader.p);
	}

	// shader's I/O (input=multitexture, output=MRT)
	int n=s->shader.no;
	char fmt[n+1];
	uint rt[n];
	int outlevel=0;
	int depthtest=0;
	if (s->shader.type==SHADER_DS) {			
	    // input
	    int bid=s->shader.in[0].b;
	    assert(bid!=-1);
	    filternode_t *in=ft->nodes + bid;
	    --outdeg[bid];

	    // output

	    // find a unpacked buffer
	    filternode_t *b=ft->nodes + s->shader.out[0].b;

	    while (b->outdeg==1) {
		filternode_t *nexts=ft->nodes + b->inout[63];
		if (nexts->shader.type==SHADER_PACK) {
		    b=ft->nodes + nexts->inout[0];
		}
		else {
		    break;
		}
	    }

	    // alloc buffer if necesary
	    if (!b->buffer.tex) {
		// filter is downsample, level++
		int level=in->buffer.level+1;
		b->buffer.level=level;
		b->buffer.tex=rt_alloc(in->buffer.type, level);
		b->buffer.type=in->buffer.type;
	    }
	    fmt[0]=in->buffer.type==RT_DEPTH ? depthtest=1,'d' : 't';
	    fmt[1]=0;
	    rt[0]=b->buffer.tex;

	    char *nam = fmt[0]=='d' ? "depth_tex" : "color_tex";
	    samp2d(nam, intex);

	    outlevel=in->buffer.level+1;
	    b->disp.w=G.screenw>>outlevel;
	    b->disp.h=G.screenh>>outlevel;

	    drawn[s->shader.out[0].b]=1;
	}
	else {
	    // input-- multitexture (setup texture samplers)
	    for (int j=0; j<s->shader.ni; ++j) {
		int bid=s->shader.in[j].b;
		if (bid!=-1) {
		    filternode_t *b=ft->nodes + bid;
		    samp2d(s->shader.in[j].nam, b->buffer.tex);
		    --outdeg[bid];
		}
	    }
			
	    // output-- render targets
	    // size of the render targets are controlled by the shader node
	    for (int j=0; j<n; ++j) {
		// find a unpacked buffer
		filternode_t *b=ft->nodes + s->shader.out[j].b;
		while (b->outdeg==1) {
		    filternode_t *nexts=ft->nodes + b->inout[63-j];
		    if (nexts->shader.type==SHADER_PACK) {
			b=ft->nodes + nexts->inout[0];
		    }
		    else {
			break;
		    }
		}

		// alloc buffer if necesary
		if (!b->buffer.tex) {
		    b->buffer.level=s->shader.level;
		    b->buffer.tex=rt_alloc(b->buffer.type, s->shader.level);
		}
		fmt[j]=b->buffer.type==RT_DEPTH ? depthtest=1,'d' : 't';
		fmt[j+1]=0;
		rt[j]=b->buffer.tex;
	    }
			
	    outlevel=s->shader.level;
	}
		
	glColorMask(s->shader.colormask[0], s->shader.colormask[1],
		    s->shader.colormask[2], s->shader.colormask[3]);


	//		glDisable(GL_DEPTH_TEST);
	rnd_tgn(fmt, rt);
	rnd_scrncoordsz(G.screenw>>outlevel, G.screenh>>outlevel);

	//		glColorMask(1,1,1,1);
		
	// draw a quad
	float co[4]={0,0,G.screenw>>outlevel, G.screenh>>outlevel};
	float texco[4]={0,0,1,1};
	glColor4f(1,1,1,1);

	if (depthtest) {
	    glEnable(GL_DEPTH_TEST);
	    rnd_clearall();
	}
	rendrect("co texco", co, texco);
	if (depthtest) {
	    glDisable(GL_DEPTH_TEST);
	}

	glColorMask(1,1,1,1);

	// if in edit mode, output the RT directly onto screen
	if (edit) {
	    vfs_use(-1);
	    rnd_tg(0);
	    cam2_use(ft->disp.x, ft->disp.y, ft->disp.zoom);

	    glEnable(GL_BLEND);
	    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	    glEnable(GL_DEPTH_TEST);

	    // draw RTs
	    float tcbnd[4]={0,0,1,1};
	    float qbnd[4];
	    for (int i=0; i<n; ++i) {
		filter_drawbuffer(ft, s->inout[63-i]);
	    }

	    // draw shader node

	    float z=(float)s->disp.rank/ft->disp.maxrank;
	    rnd_setz(z);
	    filter_drawshader(ft, sid);
	    //			drawn[sid]=1;
			
	    glDisable(GL_DEPTH_TEST);
	    glDisable(GL_BLEND);
	    // draw edges
	    // TODO:
	}

	// reclaim input buffers
	for (int j=0; j<s->indeg; ++j) {
	    if (!outdeg[s->inout[j]]) {
		filternode_t *b=ft->nodes+s->inout[j];
		rt_free(b->buffer.tex);
		b->buffer.tex=0;
	    }
	}
    }
	
    // reclaim all tmp buffers except output buffers.
    for (int i=0; i<sz; ++i) {
	filternode_t *node=ft->nodes+i;
	if (node->type==FILTER_BUFFER || node->type==FILTER_IN){
	    uint t=node->buffer.tex;
	    if (t>0) {
		rt_free(t);
		node->buffer.tex=0;
	    }
	}
    }
	

    // draw gui stuff
    if (edit) {
	vfs_use(-1);
	rnd_tg(0);
	cam2_use(ft->disp.x, ft->disp.y, ft->disp.zoom);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_DEPTH_TEST);

	// remained nodes
	glDepthFunc(GL_LEQUAL);
	for (int i=0; i<sz; ++i) {
	    if (drawn[i]) continue;

	    filternode_t *node=ft->nodes + i;
	    if (node->indeg==-1) continue;

	    node->disp.w=200;
	    node->disp.h=200;

	    float z=(float)node->disp.rank/ft->disp.maxrank;
	    rnd_setz(z);
	    filter_drawbadnode(ft, i);
	}

	// edges
	rnd_setz(-1);
	for (int i=0; i<sz; ++i) {
	    filternode_t *node=ft->nodes + i;
	    if (node->indeg==-1) continue;
	    if (node->type==FILTER_SHADER || node->type==FILTER_OUT) {
		for (int j=0; j<node->shader.no; ++j) {
		    if (node->shader.out[j].b!=-1) {
			filternode_t *b=ft->nodes + node->shader.out[j].b;
			glColor3f(1,1,1);
			float x1,y1,x2,y2;
			x1=node->disp.x+node->disp.w;
			y1=node->disp.y+node->disp.h*.5;
			x2=b->disp.x;
			y2=b->disp.y+b->disp.h*.5;
			rnd_bezier(x1,y1,x2,y2);
		    }
		}

		for (int j=0; j<node->shader.ni; ++j) {
		    if (node->shader.in[j].b!=-1) {
			int bid=node->shader.in[j].b;
			filternode_t *b=ft->nodes + bid;
			glColor3f(1,1,1);
			float x1,y1,x2,y2;
			float dy=node->disp.h/(node->shader.ni+1);
			x1=b->disp.x+b->disp.w;
			y1=b->disp.y+b->disp.h*.5;
			if (G.filter_curchannel==j && 
			    G.filter_curnodeid==i) {
			    x2=ft->disp.mousex;
			    y2=ft->disp.mousey;
			}
			else {
			    x2=node->disp.x;
			    y2=node->disp.y+ (j+1)*dy;//node->disp.h*.5;
			}

			rnd_bezier(x1,y1,x2,y2);
		    }
		    else {
#if 1
			float x1,y1,x2,y2;
			if (G.filter_curchannel==j &&
			    G.filter_curnodeid==i) {
			    float dy= node->disp.h/(node->shader.ni+1);
			    x1=ft->disp.mousex;
			    y1=ft->disp.mousey;
			    x2=node->disp.x;
			    y2=node->disp.y+(j+1)*dy;
			    rnd_bezier(x1,y1,x2,y2);
			}
#endif
		    }
		}

	    }
	}

	// current edge being edited
	if (G.filter_curchannel>=0) {
	    assert(G.filter_curnodeid>=0);
	    filternode_t *node=ft->nodes + G.filter_curnodeid;
	    if (node->type==FILTER_IN || node->type==FILTER_BUFFER) {
		float x1,y1,x2,y2;
		x1=node->disp.x+node->disp.w;
		y1=node->disp.y+node->disp.h*.5;
		x2=ft->disp.mousex;
		y2=ft->disp.mousey;
		rnd_bezier(x1,y1,x2,y2);
	    }
	}

		
	// draw shadows
	// (black) shadow blending is order-independent: (I=intensity, A=alpha, s=source, d=dest)
	// I_d2 = Is*A_s1 + I_d1*(1-A_s1),    and Is=0
	// => I_d2 = I_d1*(1-A_s1)
	// => I_d3 = I_d2*(1-A_s2) = I_d*(1-A_s1)*(1-A_s2)
	// => I_dn = I_d1*(1-A_s1)*(1-A_s2)*...*(1-A_sn_minus_1) -- multiplication is order-independent	
	glDepthFunc(GL_LESS);
	glDepthMask(0);
	for (int i=0; i<sz; ++i) {
	    filternode_t *node=ft->nodes+i;
	    if (node->indeg==-1) continue;
	    float z=(float)node->disp.rank/ft->disp.maxrank;
	    rnd_setz(z);
	    filter_drawinputshadow(node, ft->disp.zoom);
	}

	glDepthMask(1);

	glClear(GL_DEPTH_BUFFER_BIT);
	// context menu
	if (G.filter_ctxmenu_on) {
	    rnd_scrncoord();
	    rnd_setz(1);
	    glColor3f(1,1,1);
	    filter_drawmenu();
	}


	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
    }

	
}

#undef FILTER_INPUT
#undef FILTER_OUTPUT

static void filter_updaterank(filter_t *ft, int id)
{
    int sz=arr_len(ft->nodes);
    filternode_t *top=ft->nodes+id;

    // move all nodes above the selected node down
    for (int i=0; i<sz; ++i) {
	filternode_t *node=ft->nodes+i;
	if (node->indeg!=-1) {
	    if (node->disp.rank>top->disp.rank) {
		--node->disp.rank;
	    }
	}
    }

    // move the selected node to top
    top->disp.rank=ft->disp.maxrank;
}

static void filter_init()
{
    G.filter_curnodeid=-1;
    G.filter_curchannel=-1;
    menu_init(&G.filter_menuenv);
    menu_add(&G.filter_menuenv, "filter_ctxmenu{hello(Hello World){print hello world} yeah(Filters)[filters]}");
    menu_add(&G.filter_menuenv, "filters{downsampler(DownSampler){print ds} packer(Pack){print pack} blur(Blur)[blurmenu]}");
    menu_add(&G.filter_menuenv, "blurmenu{blurx(BlurX){a} blury(BlurY){a} motionblur(MotionBlur){a}}");

    G.mode=MODE_FILTER;
}

//////////////////////////////////////////////////////////////////
//// GUI CODE




/*-----------------------------------------------------------------
  map the screen coord to filter editor space
  -----------------------------------------------------------------*/
static void filter_updatemousepos(filter_t *ft, short mx, short my)
{
    float x1,y1,x2,y2;
    float w=G.screenw/ft->disp.zoom;
    float h=G.screenh/ft->disp.zoom;
    x1=ft->disp.x - w*.5;
    y1=ft->disp.y - h*.5;
	
    // cursor position filter space
    float x=x1 + w*mx/G.screenw;
    float y=y1 + h*my/G.screenh;
	
    ft->disp.mousex=x;
    ft->disp.mousey=y;
}

float zoom=1;
float filterx, filtery;
static void filter_mmfunc(SDL_MouseMotionEvent motion)
{
    filter_t *ft=G.filter_cur;
    if (!ft) return;

    filter_updatemousepos(ft, motion.x, motion.y);

    int btn=SDL_GetMouseState(NULL, NULL);
    if (btn&SDL_BUTTON(2)) {
	ft->disp.x-= motion.xrel/ft->disp.zoom;
	ft->disp.y-= motion.yrel/ft->disp.zoom;
    }
    else if (btn&SDL_BUTTON(1)) {
	if (G.filter_curnodeid>=0) {
	    filternode_t *node=G.filter_cur->nodes + G.filter_curnodeid;
	    if (G.filter_curchannel==-1) {  // drag node
		node->disp.x += motion.xrel/ft->disp.zoom;
		node->disp.y += motion.yrel/ft->disp.zoom;
	    }
	    else if (node->type==FILTER_SHADER) { // drag shader's input -- deledge + (possibly) addedge
		printf("drag shader\n");
	    }
	    else if (node->type==FILTER_BUFFER) { // drag buffer's output -- addedge
		printf("drag buffer\n");
	    }
	    else if (node->type==FILTER_IN) {
		printf("drag input\n");
	    }
	}
    }
    else {
	if (G.filter_ctxmenu_on) {
	    menu_focus(&G.filter_menuenv, motion.x,motion.y);
	}
    }

}

static int pointincircle(float px, float py, float cx, float cy, float r)
{
    float dx=px-cx;
    float dy=py-cy;
    return dx*dx+dy*dy<=r*r;
}


static int filter_picknode(filter_t *ft)
{
    float x=ft->disp.mousex;
    float y=ft->disp.mousey;
    int maxrank=-1;
    int maxwho=-1;

    int sz=arr_len(ft->nodes);
    for (int i=0; i<sz; ++i) {
	filternode_t *node=ft->nodes + i;
	if (x>=node->disp.x && x<=node->disp.x+node->disp.w &&
	    y>=node->disp.y && y<=node->disp.y+node->disp.h) {
	    if (node->disp.rank>maxrank) {
		maxrank=node->disp.rank;
		maxwho=i;
	    }
	}
    }

    return maxwho;
}


static int filter_pickpoint(filter_t *ft, int *channel)
{
    float x=ft->disp.mousex;
    float y=ft->disp.mousey;
    struct {
	int maxrank;
	int maxwho;
	int channel;
    } point={-1,-1,-1};
	
    int sz=arr_len(ft->nodes);
    for (int i=0; i<sz; ++i) {
	filternode_t *node=ft->nodes+i;
	if (node->type==FILTER_OUT || node->type==FILTER_SHADER) { // shader input
	    filternode_t *s=node;
	    int ni=s->shader.ni;
	    float inx=s->disp.x;
	    float dy=s->disp.h/(ni+1);
	    float iny=s->disp.y+dy;
			
	    for (int j=0; j<ni; ++j, iny+=dy) {							
		if (pointincircle(x,y,inx,iny,15)) {
		    if (s->disp.rank>point.maxrank) {
			point.maxrank=s->disp.rank;
			point.maxwho=i;
			point.channel=j;
		    }
		}
	    }
	}
	else if (node->type==FILTER_IN || node->type==FILTER_BUFFER) { // buffer output
	    filternode_t *b=node;
	    float outx,outy;
	    outx=b->disp.x+b->disp.w;
	    outy=b->disp.y+b->disp.h*.5;
	    if (pointincircle(x,y,outx,outy,15)) {
		if (b->disp.rank>point.maxrank) {
		    point.maxrank=b->disp.rank;
		    point.maxwho=i;
		    point.channel=0;
		}
	    }
	}
    }

    if (point.maxwho!=-1) {
	*channel=point.channel;
    }
	
    return point.maxwho;

}

/* 
   menuname {
   name (caption)  {cmd}
   name (caption)  {cmd}
   name (caption)  [submenuname]
   }

   menuname {
   
   }

*/


static jmp_buf menujmpbuf;
static void menuend()
{
    int a=setjmp(menujmpbuf);
    if (!a) return;
}

static void menubgn()
{
    menuend();
}

static void menuerr()
{
    longjmp(menujmpbuf,1);
}



/*
  return beginning of the next token
*/
static char *getname(char *d, char **source)
{
    char *s=*source;
    while (isspace(*s)) ++s;	
    while (isalnum(*s)||*s=='_'){
	*d++=*s++;
    }
    *d=0;
    *source=s;
    if (!*s) menuerr();

    return d+1;
}

static char *getblock(char *d, char **source, char bgn, char end)
{
    char *d0=d;
    char *s=*source;
    while (isspace(*s)) ++s;
    if (!*s || *s!=bgn) menuerr();

    int cnt=0;
    do {
	cnt += (*s==bgn)-(*s==end);
	*d++=*s++;
    } while (cnt && *s);
    *d=0;
    *source=s;
    if (!*s) menuerr();
	
    *d0=' ';
    *(d-1)=' ';

    return d+1;
}

static int blocklen(char *block, char bgn, char end)
{
    int cnt=0;
    char *b=block;

    // find first occurence of bgn
    while (*b && *b!=bgn) ++b;
    if (!*b) return 0;

    char *start=b;
    do {
	cnt += (*b==bgn)-(*b==end);
	++b;
    } while (cnt && *b);

    return b-start;
}

static int namecmp(char *n1, char *n2)
{
    while (*n1 && isspace(*n1)) ++n1;
    while (*n2 && isspace(*n2)) ++n2;

    while (1) {
	if (*n1==*n2) {
	    if (!*n1) return 0;
	    ++n1;
	    ++n2;
	    continue;
	}

	int end1=!*n1 || isspace(*n1);
	int end2=!*n2 || isspace(*n2);
	if (end1&&end2) return 0;
	else return 1;
    }

}

static menu_t *menuptr(menuenv_t *env, char *nam)
{
    int sz=arr_len(env->menus);
    for (int i=0; i<sz; ++i) {
	menu_t *m=env->menus+i;
	if (!namecmp(m->nam, nam)) return m;
    }
    return NULL;
}

int fnt_strwidth(int f, char *str, int h)
{
    f%=8;
    fnt_t *fnt;
    if (!G.fnt_slots[f].occupied) {
	fnt=G.fnt_slots;
    }
    else {
	fnt=G.fnt_slots+f;
    }

    float factor= (float)h/fnt->height;
    char *c=str;
    float w=0;
    while (*c) {
	w+=fnt->ch[*c++].advance;
    }

    return w*factor;
}

static void updmenusz(menuenv_t *env, int id)
{
    menu_t *menu=env->menus+id;
    int n=arr_len(menu->items);
    int maxw=0;
    for (int i=0; i<n; ++i) {
	int w=fnt_strwidth(0, menu->items[i].text, env->itemh);
	if (w>maxw) maxw=w;
    }

    menu->disp.w=maxw;
    menu->disp.h=n*env->itemh;
    menu->disp.rhs=1;
}

static int menuid(menuenv_t *env, char *nam)
{
    int sz=arr_len(env->menus);
    for (int i=0; i<sz; ++i) {
	menu_t *m=env->menus+i;
	if (!namecmp(m->nam, nam)) return i;
    }

    return -1;
}

// TODO: not finished!
void menu_add(menuenv_t *env, char *code)
{
    char buf[1024];
    char *c=code;
    menubgn();
    while (1) {
	while (isspace(*c)) ++c;
	if (!*c) break;

	// menu name
	getname(buf, &c);

	int namlen=strlen(buf);
	int len=blocklen(c,'{','}') + namlen;
	if (len<=0) menuerr();

	int id=menuid(env,buf);
	char *mem;
	if (id>=0) { // reload all
	    menu_t *m=env->menus+id;
	    m->mem=realloc(m->mem, len);
	    arr_popall(m->items);
	    mem=m->mem;
	    m->nam=mem;
	    strcpy(mem,buf);
	    mem+=namlen+1;
	}
	else {
	    int sz=arr_len(env->menus);
	    arr_pushn(env->menus,NULL,1);
	    menu_t *m=env->menus+sz;
	    m->items=NULL;
	    m->mem=malloc(len+10);
	    mem=m->mem;
	    m->nam=mem;
	    strcpy(mem,buf);
	    mem+=namlen+1;

	    id=sz;
	}

	// item properties
	/// name
	while (*c && *c!='{') ++c;
	if (!*c) menuerr();
	++c; //'{'
		
	while (*c && *c!='}') {
	    menuitem_t item;
	    item.nam=mem;
			
	    mem=getname(mem, &c);
			
	    /// caption
	    item.text=mem;
	    mem=getblock(mem, &c, '(', ')');
			
	    while (isspace(*c)) ++c;
	    if (*c=='[') {
		item.submenu=mem;
		item.cmd=NULL;
		mem=getblock(mem, &c,'[',']');
	    }
	    else if (*c=='{') { // command
		item.cmd=mem;
		item.submenu=NULL;
		mem=getblock(mem, &c, '{','}');
	    }
	    else {
		menuerr();
	    }
			
	    while (isspace(*c)) ++c;

	    arr_push(env->menus[id].items, item);
	}

	if (!*c) break;

	++c; // '}'

	updmenusz(env,id);
    }
    menuend();
}


/* 
   point (x,y) is NOT inside the menu at first
*/
void menu_show(menuenv_t *env, char *menunam, short x, short y)
{
    assert(env && menunam);
    arr_popall(env->disps);
       
    menu_t *m=menuptr(env, menunam);
    if (!m) return;

    // choose a nice position for the menu
    // default pos: the cursor is at the left-top corner
    int x2=x+m->disp.w;
    m->disp.x = x2>env->x2 ? x-(x2-env->x2) : x;
    m->disp.y = y-m->disp.h<0 ? y : y-m->disp.h;
    m->disp.rhs=1;
    int id=m - env->menus;
    arr_push(env->disps, id);
}


void menu_focus(menuenv_t *env, short x, short y)
{
    // traverse each visible menu to see if the cursor is in it
    int sz=arr_len(env->disps);
    int d=-1;
    for (int i=0; i<sz; ++i) {
	int oldd=d;
	d=env->disps[i];
	menu_t *m=env->menus+d;
	if (! (x>m->disp.x && x<m->disp.x+m->disp.w && 
	       y>m->disp.y && y<m->disp.y+m->disp.h) ){
	    d=oldd;
	}
    }

    if (d==-1){
	env->curmenu=-1;
	return;
    }

    // unshow all sub menus of the focused menu
    arr_popn(env->disps, sz-d-1);
	
    menu_t *base=env->menus+d;
    // show submenu
    int i= (base->disp.y+base->disp.h-y)/env->itemh;
    env->curmenu=d;
    env->curitem=i;

    char *sm=base->items[i].submenu;
    if (sm) {
	// reposition the menu to fit the window
	int id=menuid(env, sm);
	if (id==-1) {
	    assert(0);
	    return;
	}

	menu_t *m=env->menus+id;
	m->disp.rhs=base->disp.rhs;
	int x1,y1,x2,y2;
	if (base->disp.rhs) { // pop in RHS
	    x1=base->disp.x+base->disp.w;
	    x2=x1+m->disp.w;
			
	    if (x2>env->x2) {
		m->disp.rhs=base->disp.rhs=0;
		x1=base->disp.x-m->disp.w;
	    }
	}
	else { // pop in LHS
	    x1=base->disp.x-m->disp.w;

	    if (x1<0) {
		m->disp.rhs=base->disp.rhs=1;
		x1=base->disp.x+base->disp.w;
	    }
	}
		
	y1=base->disp.y+base->disp.h-i*env->itemh-m->disp.h;
	x2=x1+m->disp.w;
		
	if (y1<0) {
	    y1=0;
	}

	m->disp.x=x1;
	m->disp.y=y1;

	arr_push(env->disps, id);
    }
}



void menu_hide(menuenv_t *env)
{
    if (env) {
	arr_popall(env->disps);
    }
}


void menu_draw(menuenv_t *env)
{
    int sz=arr_len(env->disps);
    for (int i=0; i<sz; ++i) {
	menu_t *m=env->menus + env->disps[i];
	float bnd[4]={m->disp.x, m->disp.y, m->disp.x+m->disp.w, m->disp.y+m->disp.h};
	glColor3f(1,1,1);

	float z=(float)i/sz;


	rnd_setz(z+0.01);
	rnd_roundrect(bnd,15,10,0);
	rnd_setz(z);
	rnd_rectshadow(bnd,4);
	rnd_setz(z+0.02);
	int n=arr_len(m->items);
	float y=m->disp.y+ m->disp.h-env->itemh;
	for (int j=0; j<n; ++j,y-=env->itemh) {
	    float itembnd[4]={m->disp.x, y, m->disp.x+m->disp.w, y+env->itemh};
	    glColor3f(0,0,0);
	    rnd_text(itembnd,"%s",m->items[j].text);
	}
    }
}


void menu_click()
{
}

void menu_init(menuenv_t *env)
{
    env->menus=NULL;
    env->disps=NULL;
    env->itemh=20;
    env->x2=G.screenw;
}

void filter_drawmenu()
{
    menu_draw(&G.filter_menuenv);
}

static void filter_mbfunc(SDL_MouseButtonEvent button, int down)
{
    filter_t *ft=G.filter_cur;
    if (!ft) return;

    if (G.filter_ctxmenu_on) { // in context menu mode, only LMB & mouse motion is processed
	if (down && button.button==SDL_BUTTON_LEFT) {
	    // find the menu item clicked
	    menuenv_t *env=&G.filter_menuenv;
	    menu_focus(env, button.x, button.y);
	    if (env->curmenu!=-1) {
		menu_t *m=env->menus+env->curmenu;
		char *cmd=m->items[env->curitem].cmd;
		if (cmd) {
		    cmd_execnow(cmd);
		}
	    }
	    else {
		G.filter_ctxmenu_on=0;
		menu_hide(env);
	    }
	}
    }


    if (button.button==SDL_BUTTON_WHEELDOWN) {
	if (down) {
	    ft->disp.zoom /=1.1;
	}
    }
    else if (button.button==SDL_BUTTON_WHEELUP) {
	if (down) {
	    ft->disp.zoom *=1.1;
	}
    }
    if (button.button==SDL_BUTTON_LEFT) {
	filter_updatemousepos(ft, button.x, button.y);
	float x=ft->disp.mousex;
	float y=ft->disp.mousey;
	struct {
	    int maxrank;
	    int maxwho;
	    int channel;
	} point={-1,-1,-1};
	if (down) {
	    // map the screen coord to filter editor space

	    // picking
	    int n=arr_len(ft->nodes);
	    int maxrank=-1;
	    int maxwho=-1;
	    for (int i=0; i<n; ++i) {
		filternode_t *node=ft->nodes + i;
		int rank=node->disp.rank;
		// picking node
		if (node->indeg!=-1) {
		    if (x>node->disp.x && x<node->disp.x+node->disp.w &&
			y>node->disp.y && y<node->disp.y+node->disp.h+40) {
			if (rank>maxrank) {
			    maxrank=rank;
			    maxwho=i;
			}
		    }
		}

		// picking connection point
		if (node->type==FILTER_SHADER || node->type==FILTER_OUT) { // pick shader's input points
		    int ni=node->shader.ni;
		    float inx=node->disp.x;
		    float diny=node->disp.h/(ni+1);
		    float iny=node->disp.y+diny;
		    for (int j=0; j<ni; ++j, iny+=diny) {
			if (pointincircle(x,y,inx,iny,15)) {
			    if (rank>point.maxrank) {
				point.maxrank=rank;
				point.maxwho=i;
				point.channel=j;
			    }
			}

		    }
		}
		else if (node->type==FILTER_BUFFER || node->type==FILTER_IN) { // pick buffer's output points
		    float outx=node->disp.x+node->disp.w;
		    float outy=node->disp.y+node->disp.h*.5;
		    if (pointincircle(x,y,outx,outy,15)) {
			if (rank>point.maxrank) {
			    point.maxrank=rank;
			    point.maxwho=i;
			    point.channel=0;
			}
		    }
		}
	    }

	    // save info about the picked stuff
	    if (maxwho==-1 && point.maxwho==-1) { // nothing is picked
		G.filter_curnodeid=-1;
		G.filter_curchannel=-1;
	    }
	    else if (maxwho==-1 || (point.maxwho != -1 && point.maxrank>=maxrank)) { // point is picked
		G.filter_curchannel=point.channel;
		G.filter_curnodeid=point.maxwho;
		filter_updaterank(ft, point.maxwho);				
	    }
	    else { // node is picked
		G.filter_curchannel=-1;
		G.filter_curnodeid=maxwho;
		filter_updaterank(ft, maxwho);				
	    }
	}
	else {
	    if (G.filter_curchannel!=-1) {
		int id=G.filter_curnodeid;
		filternode_t *node=ft->nodes + id;
		//				int sz=arr_len(ft->nodes);
		int channel=0;
		int sid=filter_pickpoint(ft,&channel);

		// find if the cursor is currently inside some points
		if (sid!=-1) {
		    if (node->type==FILTER_IN || node->type==FILTER_BUFFER && sid!=id) { // from buffer
			filternode_t *s=ft->nodes+sid;;
			if (s->type!=FILTER_IN && s->type!=FILTER_BUFFER) {
			    if (s->shader.in[channel].b==-1) {
				// add new edge
				filter_addedge(ft,id,sid,channel);
			    }
			    else {
				// first remove old edge
				int oldb=ft->nodes[sid].shader.in[channel].b;
				if (oldb!=id) {
				    filter_deledge(ft,oldb,sid);

				    // add new edge
				    int err=filter_addedge(ft,id,sid,channel);
				    if (err) { // error recovery
					filter_addedge(ft,oldb,sid,channel);
				    }
				}
			    }
			}
			else {
			}
		    }
		    else if (node->type==FILTER_OUT || node->type==FILTER_SHADER) {
			if (!(sid==id && channel==G.filter_curchannel)) {
			    int bid=node->shader.in[G.filter_curchannel].b;

			    if (bid==-1) {
#if 1
				filternode_t *b=ft->nodes + sid;
				if (b->type==FILTER_IN || b->type==FILTER_BUFFER) {
				    filter_addedge(ft,sid,id,G.filter_curchannel);
				}
#endif
			    }
			    else {
				assert(bid!=-1);
				filter_deledge(ft,bid,id);
								
				filternode_t *s=ft->nodes+sid;
				if (s->type!=FILTER_IN && s->type!=FILTER_BUFFER) {
				    // first remove old edge
				    int oldb=ft->nodes[sid].shader.in[channel].b;
				    filter_deledge(ft,oldb,sid);
									
				    // add new edge
				    filter_addedge(ft,bid,sid,channel);
				}
			    }
			}
		    }
		}
		else {
		    if (node->type==FILTER_OUT || node->type==FILTER_SHADER) {
			int bid=node->shader.in[G.filter_curchannel].b;
			if (bid!=-1) {
			    filter_deledge(ft,bid,id);
			}
		    }
		}
	    }

	    G.filter_curnodeid=-1;
	    G.filter_curchannel=-1;			
	}
    }
    else if (button.button==SDL_BUTTON_RIGHT) { // pop-up menu
	if (down) {
	    if (G.filter_ctxmenu_on) {
		menu_hide(&G.filter_menuenv);
		G.filter_ctxmenu_on=0;
	    }
	    else {
		menu_show(&G.filter_menuenv, "filter_ctxmenu", button.x, button.y);
		G.filter_ctxmenu_on=1;
	    }
	}
    }
}


static void filter_keyfunc(SDL_keysym keysym, int down)
{
}



typedef struct {
    short width;
    short height;
    uchar type;
    uchar depth;
    uchar *data;
} image_t;

///////////////////// TGA image I/O ///////////////////
uchar *readtga(char *filename, int *w, int *h, int *bpp)
{
    FILE *f;
    uchar idlen;
    uchar trash[255];
    uchar type;
    uchar desc;
    image_t image;
    uchar swapx, swapy;
    int size;
    int i;
    uchar tmp;
    int row, col;
    int x, xmirror;
    int y, ymirror;
    int pos;

    f = fopen(filename, "r" );
    if(!f ){
	return NULL;
    }

    /* Load header */
    fread(&idlen, sizeof(uchar), 1, f);
    fread(trash, sizeof(uchar), 1, f);
    fread(&type, sizeof(uchar), 1, f);
    fread(trash, sizeof(uchar), 5, f);
    fread(trash, sizeof(short), 2, f);
    fread(&image.width, sizeof(short), 1, f);
    fread(&image.height, sizeof(short), 1, f);
    fread(&image.depth, sizeof(uchar), 1, f);
    fread(&desc, sizeof(uchar), 1, f);
    fread(trash, sizeof(uchar), idlen, f);

    /* Only allows uncompressed true color(2) or black-white(3) image. */
    if ((type != 2) && (type != 3)) {
	fclose (f);
	return NULL;
    }

    image.type = image.depth / 8; /* num bytes */

    /* Load image data */
    if (image.type != 1 && image.type != 3 && image.type != 4) {
	fclose (f);
	return NULL;
    }

    swapx = desc & (1 << 4);	/* left to right or right to left */
    swapy = desc & (1 << 5);	/* bottom to top or vise versa */
    size = image.width * image.height * image.type;
    image.data = malloc (size);

    fread (image.data, sizeof (uchar), size, f);
    fclose (f);

    /* tga stores BGR(A). swap B and R */

    if (image.type != 1) {
	for (i=0; i<image.width*image.height; i++) {
	    pos = i * image.type;
	    tmp = image.data[pos];
	    image.data[pos] = image.data[pos+2];
	    image.data[pos+2] = tmp;
	}
    }


    /* OpenGL texture coord centers at bottom-left. */
    if (swapx) {		/* swap columns */
	for (row=0; row<image.height; row++) {
	    for (col=0; col<image.width/2; col++) {
		x = (row*image.width + col) * image.type;
		xmirror = (row*image.width + (image.width - 1 - col)) * image.type;

		memcpy(trash, &(image.data[x]), image.type);
		memcpy(&(image.data[x]), &(image.data[xmirror]), image.type);
		memcpy(&(image.data[xmirror]), trash, image.type);
	    }
	}
    }

    if (swapy) {		/* swap rows */
	for (col=0; col<image.width; col++) {
	    for (row=0; row<image.height/2; row++) {
		y = (row*image.width + col) * image.type;
		ymirror = ((image.height - 1 - row)*image.width + col) * image.type;

		memcpy(trash, &(image.data[y]), image.type);
		memcpy(&(image.data[y]), &(image.data[ymirror]), image.type);
		memcpy(&(image.data[ymirror]), trash, image.type);
	    }
	}
    }

    *w=image.width;
    *h=image.height;
    *bpp=image.depth;
    return image.data;
}

int writetga(char *filename, int w, int h, int bpp)
{
    
}



/////////////////////////////////////////////////////////
///// Editor

enum {
    EDSTATE_IDLE, EDSTATE_NEWLOOP
};

enum {
    EDCMD_ADDNOOK
};

int nearestwall()
{
}

int nearestnook()
{
}



/////////////////////////////////////////////////////////
/*
  Nook possibilities:
  1. in some sectors => inner loop of an existing sector
  2. in none of the sectors => a new sector
  3. on some wall => invalid op
  4. overlap with another nook => new sector/spliting the sector, depending on the next nook
*/

/*
  after a loop being added (new sector or existing sector's new inner loop), you have to update:

  After the following operations:
  1) Adding a loop (outer or inner)
  2) moving nooks

  the sector needs a clean up:
  1) resolve the sector's self-intersections
  2) find shared sectors
  3) find valid/invalid shared walls
  4) group non-intersecting sectors into regions
  5) optimise region boundaries
  6) update lights' visibilities
*/


/*
 * Given M points, returns triangles

 co_in: 
 tri_out: [v1_id v2_id v3_id]
*/

#define dot2(a,b) ((a)[0]*(b)[0]+(a)[1]*(b)[1])
#define cross2(a,b) ((a)[0]*(b)[1]-(a)[1]*(b)[0])
#define proj2(a,b) ((a)[0]*(b)[0]+(a)[1]*(b)[1])
#define length2(a) (sqrtf(dot2((a),(a))))


void normalize2(float b[2], float a[2])
{
    float l=length2(a);
    if (!l) {
	b[0]=b[1]=0;
	return;
    }

    b[0]=a[0]/l;  
    b[1]=a[1]/l;
}

typedef struct {
    float co[2];
    int wid; // vertex index;
} tessvert_t;

static float tess_v[60000];
static int tess_wid[60000];
static int tess_nv;
static char tess_iswall[60000];
static char tess_iswallflag;
static tessvert_t tess_combinev[10000];
static tessvert_t *tess_combinevp=tess_combinev;

static void tess_vertcb(tessvert_t *v)
{
    float *t=tess_v + 2*tess_nv;

    t[0]=v->co[0];
    t[1]=v->co[1];
    tess_iswall[tess_nv]=tess_iswallflag;
    tess_wid[tess_nv]=v->wid;
    ++tess_nv;
    tess_combinevp=tess_combinev;
}

static void tess_edgecb(int iswall)
{
    tess_iswallflag=iswall;    
}

static void tess_begincb()
{
    tess_combinevp=tess_combinev;
    tess_nv=0;
}

static void tess_endcb()
{
}

static int commonedge(int i1, int i2)
{
    if (i1<0xffff) { // vertex
	if (i2<0xffff) {  // vertex
	    if (abs(i1-i2)>1) {
		return i1<i2 ? i2 : i1;
	    }
	    else {
		return i1<i2 ? i1 : i2;
	    }
	}
	else { // intersect
	    int hi=i2>>16;
	    return i1==hi ? hi : i2&0xffff;
	}
    }
    else { // intx
	if (i2<0xffff) { // vertex
	    int hi=i1>>16;
	    return i2==hi ? hi : i1&0xffff;
	}
	else { // intx
	    int hi1=i1>>16;
	    int hi2=i2>>16;
	    int lo1=i1&0xffff;
	    int lo2=i2&0xffff;
	    return (hi1==hi2||hi1==lo2) ? hi1 : lo1;
	}
    }
}

static void tess_combinecb( GLdouble coords[3],
			    tessvert_t *tv[4],
			    GLfloat weight[4], void **data )
{
    tessvert_t *v=tess_combinevp++;
    v->co[0]=coords[0];
    v->co[1]=coords[1];
    int e1=commonedge(tv[0]->wid,tv[1]->wid);
    int e2=commonedge(tv[2]->wid,tv[3]->wid);

    v->wid=e1>e2 ? (e1<<16)|e2 : (e2<<16)|e1;   // to avoid (0<<16)|x == x
    *data = v;
}


/*
  Sort floating point values by radix sort.
  But instead of writing the sorted data directly into the original data, we just store the rank of each element in an array

  rank[0] -- index of minimum element in input[]
  rank[n-1] -- index of maximum element in input[]

  So:
  input[rank[0]] <= input[rank[1]] <= ... <= input[rank[n-1]]

*/
void rank32f(float *input, uint *rank, int n)
{
    // write the (unsorted/incorrect) rank
    uint *src,*dst;
    src=rank;
    dst=alloca(sizeof(uint)*n);

    for(int i=0; i<n; i++) src[i]=i;

    // sort floats just as uints
    int counter[256], offset[256];
    uint *f=input;
    for(int pass=0; pass<4; pass++){
	uint shift=pass<<3;

	memset(counter, 0, sizeof(counter));
	memset(offset, 0, sizeof(offset));

	for(int i=0; i<n; i++){
	    uchar radix=(f[i]>>shift) & 0xff;
	    counter[radix]++;
	}

	if(pass!=3){	
	    offset[0] = 0;
	    for(int i=1; i<256; i++){
		offset[i] = offset[i-1] + counter[i-1];
	    }

	    for(int i=0; i<n; i++){
		uchar radix=(f[src[i]]>>shift) & 0xff;
		dst[offset[radix]++] = src[i];
	    }
			
	}
	else{	// fix floats
	    offset[0] = 0;
	    for(int i=1; i<128; i++){
		offset[i] = offset[i-1] + counter[i-1];				
	    }

	    offset[128] = 0;
	    int nn=counter[128];
	    for(int i=129; i<256; i++){
		offset[i] = offset[i-1] - counter[i-1];
		nn+=counter[i];
	    }

	    dst += nn;
	    for(int i=0; i<n; i++){
		uchar radix=(f[src[i]]>>shift) & 0xff;
		if(radix>127){ // neg
		    int j=--offset[radix];
		    dst[j] = src[i];

//					dst[--offset[radix]] = src[i];
		}
		else{
		    int j=offset[radix]++;
		    dst[j] = src[i];

//					dst[offset[radix]++] = src[i];
		}
	    }
			
	    dst -= nn;
	}

	uint *tmp=src;
	src=dst;
	dst=tmp;
	shift += 8;
    }
}


void rsort32(uint data[], int n)
{
    int counter[256], offset[256];
    uint *dst = alloca(n*sizeof(uint));
    uint *src = data;
    for(int pass=0; pass<4; pass++){
	uint shift = pass<<3;
	
	memset(counter, 0, sizeof(counter));
	memset(offset, 0, sizeof(offset));
	
	for(int i=0; i<n; i++){
	    uchar radix = (src[i]>>shift) & 0xff;
	    counter[radix]++;
	}
	
	offset[0] = 0;
	for(int i=1; i<256; i++){
	    offset[i] = offset[i-1] + counter[i-1];
	}
	
	for(int i=0; i<n; i++){
	    uchar radix = (src[i]>>shift) & 0xff;
	    dst[offset[radix]++] = src[i];
	}
	
	uint *tmp = src;
	src = dst;
	dst = tmp;
	shift += 8;
    }
}


/*
  triangle connectivities:
  1) inner edges -- same sector, one adjacent edge
  2) outer edges -- one different sector, zero or more adjacent triangle edges (zero or one adjacent wall)

  each shared wall should know all the triangle edges it contains

  

 */


/*
  water reflection is easy -- visibility computed for the main view point can be reused for reflection.
  
  EXACT point visibility is a simple polygon without holes (however in 3D this visibility isn't exact -- conservative visibility)
  for each sector, the following info are needed:
    1) visible floor/ceiling polys (should be converted to 3D)
    2) walls
    3) shared walls
  
  first pass:
  simply write down all polys:

  poly:
  N verts: x,y
  N edges: wall/shared/internal; if shared, id of adjacent wall (not edge)

  the data structure should be compact

  many polygon may share a single vertex

  should all tess vertices be stored together? or indexed by (sector_id,vert_id)?
  in each sector, all duplicated vertices should be removed. -- don't!

  In a sector, each tess edge should be classified -- wall, shared, internal.
  For shared edge, the wall it's in should also be known.


  A shared wall should record all tess vertices on it, sorted in increasing distance to the first vertex of the wall.
  
  A shared wall can be tagged as transparent -- sorting them from near to far is done implicitly during the PVS computation.


 */



/*
  point P on edge P1P2 test

  if P is on P1P2, then t is:
  P = P1+(P2-P1)*t
 */
static int pntonedge(float p[2], float p1[2], float p2[2], float *t)
{
    
}


static int edgeonwall()
{
}


/*
  wall--0x7fffffff
  internal--id to the adj tess edge
  shared--negative id to the shared wall (short), and beginning index to the edge list entries on the shared wall (short).

  for each internal tess edge: find shared edge (sort by gradient)
  for each external tess edge: find which wall it's in (sort walls by gradient; sort tess edges by gradient too)

 */

#if 0
typedef struct {
    float *v; // vertex data [x1,y1],[x2,y2],[x3,y3]...
    int *e; // edge data e12,e23,e31,...
} tess_t;
#endif



typedef struct {    
    

} pvs_t;



/*
  each sector is tesselated to a bunch of triangles
 */
typedef struct {
    float *v[2];
    uint id; // (sid<<16) | eid
} tessedge_t;


/*
  e1<e2 if 
 */
static int cmpedges(tessedge_t *e1, tessedge_t *e2)
{
    float dx1=e1->v[1][0]-e1->v[0][0];
    float dy1=e1->v[1][1]-e1->v[0][1];
    float v1[4];
    if (dx1<0 || (dx1==0&&dy1<0)) {
	v1[0]=e1->v[1][0];  v1[1]=e1->v[1][1];
	v1[2]=e1->v[0][0];  v1[3]=e1->v[0][1];
    }
    else {
	v1[0]=e1->v[0][0];  v1[1]=e1->v[0][1];
	v1[2]=e1->v[1][0];  v1[3]=e1->v[1][1];
    }

    float dx2=e2->v[1][0]-e2->v[0][0];
    float dy2=e2->v[1][1]-e2->v[0][1];
    float v2[4];
    if (dx2<0 || (dx2==0&&dy2<0)) {
	v2[0]=e2->v[1][0];  v2[1]=e2->v[1][1];
	v2[2]=e2->v[0][0];  v2[3]=e2->v[0][1];
    }
    else {
	v2[0]=e2->v[0][0];  v2[1]=e2->v[0][1];
	v2[2]=e2->v[1][0];  v2[3]=e2->v[1][1];
    }
    
    for (int i=0;i<4;i++) {
	if (v1[i]<v2[i]) return -1;
	if (v1[i]>v2[i]) return 1;
    }

    return 0;
}


typedef struct {
    float t;
    short w;
    int vid; // tess vertex indices -- hi:vertex with larger t; lo: vertex with smaller t
} commonwallvert_t;

static int cmpwallvert(commonwallvert_t *v1, commonwallvert_t *v2)
{
    if (v1->w = v2->w) {
	float dt=v1->t - v2->t;
	return dt==0 ? v1->vid-v2->vid : (dt<0 ? -1 : (dt>0 ? 1 : 0));
    }
    else return v1->w-v2->w;
}


static int abovewall(float *v, float *v1, float *v2)
{
    if (v1[0]==v2[0]) return v[0]<v1[0];
    return (v[1]-v1[1])*(v2[0]-v1[0]) > (v[0]-v1[0])*(v2[1]-v1[1]);
    
}


void printedge(tessedge_t *e)
{
    printf("(%f %f) (%f %f)\n", e->v[0][0], e->v[0][1], e->v[1][0], e->v[1][1]);
}

/*
  CCW:
  >0 => point is to the left of R1R2
  <0 => ................right ......
  =0 => .........on R1R2

 */
static inline float sign(const float v[2], const float r1[2], const float r2[2])
{
    return (r2[0]-r1[0])*(v[1]-r1[1])-(r2[1]-r1[1])*(v[0]-r1[0]);
}



float *add2(float *c, const float *a, const float *b)
{
    c[0]=a[0]+b[0];
    c[1]=a[1]+b[1];
    return c;
}


float *sub2(float *c, const float *a, const float *b)
{
    c[0]=a[0]-b[0];
    c[1]=a[1]-b[1];
    return c;
}

float *interp2(float *c, const float *a, const float *b, float t)
{
    c[0]=a[0]+(b[0]-a[0])*t;
    c[1]=a[1]+(b[1]-a[1])*t;
    return c;
}


/*
  input: one edge, two clip vectors

  all clip vectors start at the PVS query point 
  => two clip vectors can be represented by the entry points to the triangle
 
  x1->x2 is ray
  x3->x4 is edge

 */

static float rayedgeintx(float v1[2], float v2[2], float v3[2], float v4[2])
{
    float denom=(v4[1]-v3[1])*(v2[0]-v1[0]) - (v4[0]-v3[0])*(v2[1]-v1[1]);
    assert(denom);

    float numer2=(v2[0]-v1[0])*(v1[1]-v3[1]) - (v2[1]-v1[1])*(v1[0]-v3[0]);

    float t=numer2/denom;
//    assert(t>=0 && t<=1 || (printf("t=%f\n",t),0));
//    if (t<0) t=0;
//    if (t>1) t=1;
    return t;//numer/denom;
}


static int segsegintx(float v1[2], float v2[2], float v3[2], float v4[2], float *t1, float *t2)
{
    float denom=(v4[1]-v3[1])*(v2[0]-v1[0]) - (v4[0]-v3[0])*(v2[1]-v1[1]);

    if (!denom) return 0;

    float numer1=(v4[0]-v3[0])*(v1[1]-v3[1]) - (v4[1]-v3[1])*(v1[0]-v3[0]);
    float numer2=(v2[0]-v1[0])*(v1[1]-v3[1]) - (v2[1]-v1[1])*(v1[0]-v3[0]);
    *t1=numer1/denom;
    *t2=numer2/denom;

    return *t1>=0 && *t1<=1 && *t2>=0 && *t2<=1;
}


/*
  t1,t2 should be inited when passed in.

  v1v2 -- line
  v3v4 -- seg

  both t1 and t2 are correspond to the segment (== have nothing to do with the line)
 */
static void lineclipseg(float v1[2], float v2[2], float v3[2], float v4[2], float *t1, float *t2)
{
    float line[2], seg[2];
    sub2(line,v2,v1);
    sub2(seg,v4,v3);
    int infront = cross2(line,seg)>0;

    float denom=(v4[1]-v3[1])*(v2[0]-v1[0]) - (v4[0]-v3[0])*(v2[1]-v1[1]);
    if (!denom) {
	if (!infront) *t2=*t1;
	return 0;
    }

//.    float numer1=(v4[0]-v3[0])*(v1[1]-v3[1]) - (v4[1]-v3[1])*(v1[0]-v3[0]);
    float numer2=(v2[0]-v1[0])*(v1[1]-v3[1]) - (v2[1]-v1[1])*(v1[0]-v3[0]);

    float t=numer2/denom;

    if (infront) { // t1=t, t2=1
	if (t>1) {
	    *t2=*t1;
	    return 1;
	}
	else if (t>0) {
	    if (t>*t1) *t1=t;
	    return 1;
	}
	else
	    return 0;
    }
    else { // t1=0, t2=t
	if (t<0) {
	    *t2=*t1;
	    return 1;
	}
	else if (t<1) {
	    if (t<*t2) *t2=t;
	    return 1;
	}
	else
	    return 0;
    }
}


/*
  Return:
  0 -- the two lines are parallel
  1 -- the two lines intersect

 */
static int linelineintx(float v1[2], float v2[2], float v3[2], float v4[2], float *t1, float *t2)
{
    float denom=(v4[1]-v3[1])*(v2[0]-v1[0]) - (v4[0]-v3[0])*(v2[1]-v1[1]);

    if (!denom) return 0;

    float numer1=(v4[0]-v3[0])*(v1[1]-v3[1]) - (v4[1]-v3[1])*(v1[0]-v3[0]);
    float numer2=(v2[0]-v1[0])*(v1[1]-v3[1]) - (v2[1]-v1[1])*(v1[0]-v3[0]);
    *t1=numer1/denom;
    *t2=numer2/denom;

    return 1;
}


/*
  the point should already be on the line.
 */
static inline int pntonseg(const float p[2], const float s1[2], const float s2[2])
{
    float dx=s2[0]-s1[0], dy=s2[1]-s1[1];
    float t=((p[0]-s1[0])*dx + (p[1]-s1[1])*dy)/(dx*dx+dy*dy);
    return t>0 && t<1;
}

/*
  same as pntonseg
 */
static inline float prj2seg(const float p[2], const float s1[2], const float s2[2])
{
    float dx=s2[0]-s1[0], dy=s2[1]-s1[1];
    float t=((p[0]-s1[0])*dx + (p[1]-s1[1])*dy)/(dx*dx+dy*dy);
    return t;
}


/*
  O(nlogn)
 */
static void inittessedgeinfo(short sid, const char iswall[])
{
    sect_t *s=G.sects+sid;
    int nv=s->nt*3;

    int ex[nv]; // external tess edges
    int nex=0;

    // 1. classify tess edge types
    tessedge_t edges[nv];
    int ne=0;

    // clean up old info
    for (int i=0;i<s->nw;i++) {
	wall_t *w=G.walls + s->firstw + i;
	w->ne[0]=w->ne[1]=0;

	float dx=G.walls[w->nextco].x-w->x;
	float dy=G.walls[w->nextco].y-w->y;
	w->reverse = dx>dy || dx==dy&&dy<0;

	if (w->adjw!=-1) { // shared wall
	    wall_t *aw=&G.walls[w->adjw];
	    sect_t *as=&G.sects[w->adjs];
	    for (int j=0;j<aw->ne[0];j++) {
		as->e[aw->elist[0][j]]= ((w->adjw - as->firstw +1)<<16) | j;
	    }
	    for (int j=0;j<aw->ne[1];j++) {
		as->e[aw->elist[1][j]]= ((w->adjw - as->firstw +1)<<16) | j;
	    }	    
	}
    }

    s->e=realloc(s->e, nv*sizeof(int));
    for (int i=0;i<nv;i++) {
	if (iswall[i]) {
	    int i2=(i%3==2) ? i-2 : i+1;
	    int wid=commonedge(tess_wid[i],tess_wid[i2]); // real_wall_id = wid + s->firstw

	    wall_t *w=&G.walls[s->firstw + wid];
	    wall_t *w2=&G.walls[w->nextco];
	    // add to wall's elist

	    // two shared walls should have compatible definition of 'side'
	    float v1[2]={w->x,w->y}, v2[2]={w2->x,w2->y};
	    int i3=(i2%3==2) ? i2-2: i2+1;
	    float sn=sign(&s->tess[2*i3], v1, v2);
	    if (w->reverse) sn*=-1;
	    int side=sn<0; // left=>0, right=>1
	    w->elist[side]=realloc(w->elist[side], (w->ne[side]+1)*sizeof(short)); // TODO: mem management
	    w->elist[side][w->ne[side]++]=i;

	    if (w->adjw==-1) { // solid wall
		s->e[i]=0x7fff0000 | wid;
	    }
	    else { // shared wall
		s->e[i]= ((wid+1)<<16) | (w->ne[side]-1); 
	    }
	}
	else {
	    tessedge_t *e=edges+ne++;
	    float *t=s->tess+i*2;
	    e->v[0]=t;
	    e->v[1]=(i%3==2) ? t-4 : t+2;
	    e->id=i;
	}
    }

    // sort internal edges
    // each internal edge has exactly one adjacent internal edge
    qsort(edges,ne,sizeof(tessedge_t),cmpedges);

    // now in the sorted array, edges[0] is adjacent to edges[1]; edges[2] adj. edges[3]; ...
    for (int i=0;i<ne;i+=2) {
	s->e[edges[i].id]=edges[i+1].id;
	s->e[edges[i+1].id]=edges[i].id;
	if (!(edges[i].v[0][0]==edges[i+1].v[1][0] && 
	      edges[i].v[0][1]==edges[i+1].v[1][1] && 
	      edges[i].v[1][0]==edges[i+1].v[0][0] &&
	      edges[i].v[1][1]==edges[i+1].v[0][1])) {
	    assert(0);
	};
    }
}


/*
  adj_out[] stores each edge's adjacent edge in the same triangle
 */
int tess(int co_in[], int nci,
	 float **tri_out, int *nt)
{
    GLUtesselator *t=gluNewTess();
    gluTessNormal(t,0,0,1);
    gluTessCallback(t, GLU_TESS_BEGIN, tess_begincb);
    gluTessCallback(t, GLU_TESS_VERTEX, tess_vertcb);
    gluTessCallback(t, GLU_TESS_END, tess_endcb);
    gluTessCallback(t, GLU_TESS_EDGE_FLAG, tess_edgecb);
    gluTessCallback(t, GLU_TESS_COMBINE, tess_combinecb);
    gluBeginPolygon(t);
    
    tessvert_t *verts=malloc(2*nci*sizeof(tessvert_t));
    tessvert_t *v=verts;
    int *co=co_in;
    gluTessBeginContour(t);
    for (int i=0; i<nci; i++) {
	if (*co==0x7fffffff) { // a new contour begins
	    gluTessEndContour(t);
	    gluTessBeginContour(t);
	    ++co;
	}
	GLdouble loc[3];
	loc[0]=v->co[0]=co[0];
	loc[1]=v->co[1]=co[1];
	loc[2]=0;
	v->wid=i; // relative wall id -- real_wall_id = v->wid + sector->firstw
	co+=2;
	gluTessVertex(t,loc,v);
	++v;
    }
    gluTessEndContour(t);

    gluEndPolygon(t);
    gluDeleteTess(t);
    
    int sz=2*sizeof(float)*tess_nv;
    float *tris=malloc(sz);
    memcpy(tris,tess_v,sz);
    free(verts);

    *tri_out=tris;
    *nt=tess_nv/3;

    // for each edge of each triangle, find adjacent edge & triangle

}

/*
  check location (x,y)

  Only the first nook of a loop needs a full check -- 1~4
  Intermedia nooks only need to check walls & internal nooks;
  Can't add new wall over shared wall -- the 2nd nook can't be on the
  same shared wall.

  TYPE_IN/TYPE_OUT -- wall
  TYPE_SPLIT -- wall,nook
  1ST_NOOK/TYPE_PENDING -- wall,nook,in/out

  Wall test is always needed. If the nook is found on a wall, return immediately.
  Nook test is a part of wall test.

  
*/

static int chkloc(float x, float y, short *dat, int datsz)
{
    short c[10000],nc=0; // coincide nooks
    short s[10000],ns=0; // sectors that the nook is inside

    for (int i=0; i<G.nsects; i++) {
	int *b=G.sects[i].bnd;
	if (x<b[0] || x>b[2] || y<b[1] || y>b[3]) continue;
		
	int intx=0;
	int w=G.sects[i].firstw;
	int lastw=w+G.sects[i].nw;
	for (;w<lastw;w++) {
	    int y1=G.walls[w].y;
	    int y2=G.walls[G.walls[w].nextco].y;
	    int x1=G.walls[w].x;
	    int x2=G.walls[G.walls[w].nextco].x;
	    if (x==x1&&y==y1) {
		c[nc++]=w;
		break;
	    }
	    if (x==x2&&y==y2) {
		c[nc++]=G.walls[w].nextco;
		break;
	    }
	    if (y1==y2) {
		if (y==y1 && (x>x1&&x<x2 || x>x2&&x<x1)) return LOC_WALL;
		continue;
	    }
	    if (y<y1&&y<y2 || y>y2&&y>y1) continue;

	    if (y==y2 || y==y1) {
		if (y2>y1) { // y2 should <= y1
		    int tmp;
		    tmp=y2; y2=y1; y1=tmp;
		    tmp=x2; x2=x1; x1=tmp;
		}
		if (y==y1 && y>y2 && x<x1) {
		    intx^=1;
		}
		continue;
	    }

	    int k1=(x-x1)*(y2-y1);
	    int k2=(y-y1)*(x2-x1);
	    if (k1==k2) return LOC_WALL;
	    else if (y2>y1&&k2>k1 || y2<y1&&k2<k1) intx^=1;
	}

	if (intx) s[ns++]=i;
    }    

    if (nc) {
	int sz=nc*sizeof(short);
	if (sz>datsz) sz=datsz;
	memcpy(dat,c,sz);
	return LOC_NOOK;
    }
    else if (ns) {
	int sz=ns*sizeof(short);
	if (sz>datsz) sz=datsz;
	memcpy(dat,s,sz);
	return LOC_IN;
    }
    else {
	return LOC_OUT;
    }

}

#define SAMENOOK(w,w2) (G.walls[(w)].x==G.walls[(w2)].x && G.walls[(w)].y==G.walls[(w2)].y)

static void findsharedwalls(int sid)
{
    sect_t *s=G.sects+sid;
    for (int i=0; i<G.nsects; i++) {
	sect_t *s2=G.sects+i;
	if (s==s2) continue;

	if (s->bnd[0]>s2->bnd[2] || s->bnd[2]<s2->bnd[0] ||
	    s->bnd[1]>s2->bnd[3] || s->bnd[3]<s2->bnd[1]) continue;

	int lastw=s->firstw+s->nw;
	int lastw2=s2->firstw+s2->nw;
	for (int w=s->firstw; w<lastw; w++) {
	    for (int w2=s2->firstw; w2<lastw2; w2++) {
		if ((SAMENOOK(w,w2) && SAMENOOK(G.walls[w].nextco, G.walls[w2].nextco)) ||
		    (SAMENOOK(w,G.walls[w2].nextco) && SAMENOOK(G.walls[w].nextco, w2))) {
		    if (G.walls[w].adjw==-1 && G.walls[w2].adjw==-1) {
			// tests
			G.walls[w].adjs=i;
			G.walls[w].adjw=w2;
			G.walls[w2].adjs=sid;
			G.walls[w2].adjw=w;
		    }
		}
	    }
	}
    }
}

static void closesectloop()
{
    int sz=arr_len(G.nkstk);
    float *tris;
    int nt;
    tess(G.nkstk,sz,&tris,&nt);
    
    // walls
    wall_t *w=G.walls + G.nwalls;
    sect_t *s=G.sects + G.nsects;
    s->tess=tris;
    s->nt=nt;
    s->firstw=G.nwalls;
    s->nw=sz;
    short next=s->firstw+1;
    s->bnd[0]=s->bnd[1]=2*WORLD_MAX;
    s->bnd[2]=s->bnd[3]=-2*WORLD_MAX;
    for (int i=0;i<sz;i++) {
	w->x=G.nkstk[i].x;
	w->y=G.nkstk[i].y;
	w->nextco=next++;
	w->adjw=-1;
	if (w->x<s->bnd[0]) s->bnd[0]=w->x;
	if (w->x>s->bnd[2]) s->bnd[2]=w->x;
	if (w->y<s->bnd[1]) s->bnd[1]=w->y;
	if (w->y>s->bnd[3]) s->bnd[3]=w->y;

	// TEMP!
	w->floorz=-20;
	w->ceilingz=20;

	++w;
    }
    w[-1].nextco=s->firstw;
    G.nwalls+=sz;
    arr_popall(G.nkstk);

    findsharedwalls(G.nsects);

    inittessedgeinfo(G.nsects, tess_iswall);

    ++G.nsects;
}


static void addtoloop()
{
    int sz=arr_len(G.nkstk);
    sect_t *s=G.sects + G.ed.cur_sect;
    int co[2*(2*s->nw+sz)];
    int w=s->firstw;
    int lastw=s->firstw+s->nw;
    int *c=co;
    int prev1st;
    int nc=0; // num of contours
    
    for (w;w<lastw;w++) {
	c[0]=G.walls[w].x;
	c[1]=G.walls[w].y;
	c+=2;
	if (w+1!=G.walls[w].nextco) *c++=0x7fffffff;
    }

    // update links
    for (int i=0;i<lastw;i++) {
	if (G.walls[i].adjw>=lastw) G.walls[i].adjw+=sz;
    }
    for (int i=lastw;i<G.nwalls;i++) {
	if (G.walls[i].adjw>=lastw) G.walls[i].adjw+=sz;
	if (G.walls[i].nextco>=lastw) G.walls[i].nextco+=sz;
    }
    for (int i=0;i<G.nsects;i++) {
	if (G.sects[i].firstw>=lastw) G.sects[i].firstw+=sz;
    }

    // shift walls
    for (int i=G.nwalls-1;i>=lastw;i--) {
	G.walls[i+sz]=G.walls[i];
    }

    w=lastw;
    for (int i=0;i<sz;i++) {
	int x,y;
	x=c[0]=G.walls[w].x=G.nkstk[i].x;
	y=c[1]=G.walls[w].y=G.nkstk[i].y;
	G.walls[w].nextco=w+1;
	G.walls[w].adjw=-1;
	if (x<s->bnd[0]) s->bnd[0]=x;
	if (x>s->bnd[2]) s->bnd[2]=x;
	if (y<s->bnd[1]) s->bnd[1]=y;
	if (y>s->bnd[3]) s->bnd[3]=y;
	++w;
	c+=2;
    }
    G.walls[w-1].nextco=lastw;
    G.nwalls+=sz;
    s->nw+=sz;

    float *tris;
    int nt;
    tess(co,s->nw,&tris,&nt);
    if (s->tess) free(s->tess);
    s->tess=tris;
    s->nt=nt;

    arr_popall(G.nkstk);

    findsharedwalls(G.ed.cur_sect);
    inittessedgeinfo(G.ed.cur_sect, tess_iswall);

    printf("num tris=%i\n",nt);
    float *t=s->tess;
    for (int i=0;i<nt;i++) {
	printf("%f %f, %f %f, %f %f\n", t[0],t[1],t[2],t[3],t[4],t[5]);
	t+=6;
    }
    printf("------------------\n");
}


static void splitloop()
{
}

int pushnook(int x,int y)
{
    int loc;
    short dat[1000];
    int sz=arr_len(G.nkstk);
    if (sz==0) { // first point of a loop, check check!
	loc=chkloc(x,y,dat,sizeof(dat));

	switch (loc) {
	case LOC_WALL:
	    return 1;

	case LOC_NOOK:
	    G.nktype=TYPE_PENDING;
	    break;

	case LOC_IN:
	    G.nktype=TYPE_IN;
	    G.ed.cur_sect=dat[0];
	    break;

	case LOC_OUT:
	    G.nktype=TYPE_OUT;
	    break;
	}
    }
    else {
	int type=G.nktype;
	nook_t *nk=G.nkstk;
	// finish?
	if (nk[0].x==x&&nk[0].y==y) { // loop closed
	    if (sz<3) { // too few nooks
		return 1;
	    }
	    else {
		switch (G.nktype) {
		case TYPE_SPLIT:
		    splitloop();
		    return 1;

		case TYPE_IN:
		    addtoloop();
		    return 0;

		case TYPE_OUT:
		    closesectloop();
		    return 0;
		}
	    }
	}

	for (int i=1;i<sz;i++) { // coinside nook is not allowed when the loop is init
	    if (nk[i].x==x&&nk[i].y==y) {
		return 1;
	    }
	}

	// find conflict with other sectors
	loc=chkloc(x,y,dat,sizeof(dat));
	if (loc==LOC_WALL) {
	    return 1;
	}

	if (type==TYPE_PENDING) {
	    G.nktype = loc==LOC_IN ? TYPE_IN : TYPE_OUT;
	}
	
    }

    nook_t n;
    n.x=x;
    n.y=y;
    arr_push(G.nkstk, n);

    return 0;
}

void popnook()
{
    arr_pop(G.nkstk);
}


void c_editor_pushnk()
{
    if (G.mode!=MODE_EDITOR) return;
    pushnook(G.ed.mousex,G.ed.mousey);
}

void c_editor_popnk()
{
    if (G.mode!=MODE_EDITOR) return;
    popnook();
}

////////////////////////////////////////////////////////


void xlatsect(int i)
{
	
}

void insnook(int x,int y)
{
}

void xlatnook(int dx, int dy)
{
    /* 
       Lots of stuff need to be updated:
       1. Adding/removing self-intersections
       2. Sectors became non-overlap => merge sectors into region
       3. Sectors became overlap => split region
       4. light visibility changes
       5. Objects intersecting with walls/floors will be removed or hided.
       6. Lights outside all sectors will be removed or hided
       7. Shared edges may change shapes
    */
    
}

void delnook()
{
}


/////////////////////////////////////////////////////////


/*
  Each light has its own PVS.

  When it moves or the sector that contains it moves, the PVS should be updated.


  The rendering loop looks like this:
  
  Pass 1:
  Depth only pass for all sectors in PVS

  Pass 2:
  for each sector in PVS:
    render & write stencil buffer
    for each light whose LPVS touches this sector:
      render while stencil test is enabled

  
  There won't be too many lights so iterating over all lights every frame is not necessarily a bad idea:
  it makes things much simpler.


  TODO: 5.13
 */
static void bgnpvs(pvs2d_t *p, float r);
static void endpvs(pvs2d_t *p);

ushort addlight(ushort sid, v3_t pos, color_t color, float radius)
{
    if (G.nlights >= MAXNLIGHTS) return;

    // init
    light_t *light = &G.lights[G.nlights];
    light->pos = pos;
    light->radius = radius;
    light->color = color;

    // LPVS
    bgnpvs(&light->pvs, radius);
    pvs2d(&light->pvs, sid, pos, vec(0,1,0), M_PI*2.0f);
    endpvs(&light->pvs);

    // book-keeping stuff (build a connection between sector & light)
    ushort id = G.nlights;
    pvs2d_wall_t *w = light->pvs.walls;
    for (int i=0; i<light->pvs.wall_ns; i++) {
	arr_push(G.sects[w->sid].lightids, id);
	w += light->pvs.nwalls[i];
    }
    
    return G.nlights++;
}


void setlightcolor(int id, color_t color)
{
    if (id >= G.nlights) return;

    G.lights[id].color = color;
}

void xlatlight()
{
}

void drawgood(uint walltex, uint flortex)
{
    
}

/*
  Deleting a light == move the last light in the light array to the pos of the deleted light
 */
void dellight()
{
}

void drawlight(light_t *light)
{
    
}

/////////////////////////////////////////////////////////


void addobj()
{
}

void xlatobj()
{
}

void delobj()
{
}


///////////////////////////////////////////////////////////
void addbip()
{
}

void xlatbip()
{
}

void delbip()
{
}



#define SOLID_WALL 1

static int cmpwalls(pvs2d_wall_t *w1, pvs2d_wall_t *w2)
{
    if (w1->sid!=w2->sid) return w1->sid-w2->sid;
    if (w1->flag!=w2->flag) return w1->flag-w2->flag; // shared wall < solid wall
    if (w1->eid!=w2->eid) return w1->eid-w2->eid;
    return w1->t1==w2->t1 ? 0 : (w1->t1<w2->t1 ? -1 : 1);
}


float pntsegdist2(float *p, float *s1, float *s2)
{
    float dx=s2[0]-s1[0];
    float dy=s2[1]-s1[1];

    float dx1=p[0]-s1[0];
    float dy1=p[1]-s1[1];
    float prj1=dx1*dx + dy1*dy;
    float d1=dx1*dx1+dy1*dy1;
    if (prj1 <= 0) {
	return d1; 
    }

    float dx2=p[0]-s2[0];
    float dy2=p[1]-s2[1];
    float prj2=dx2*dx + dy2*dy;
    if (prj2 >= 0) {
	return dx2*dx2+dy2*dy2;
    }

    float d=dx*dx+dy*dy;
    return d1-prj1*prj1/d;
}

static int pvs2d_r(pvs2d_t *pvs, short sid, short eid, float t1, float t2, int overloadtype);

static float pvs2d_center[2], pvs2d_dir[2];
static float pvs2d_r2;
static int pvs2d_directional;
enum {PVS2D_TYPE_SOLID=1, PVS2D_TYPE_INTERNAL};

static int frusculledge(float ev1[2], float ev2[2], short sid, short id1, float t[2], pvs2d_t *pvs)
{
    short id2=id1%3==2 ? id1-2 : id1+1;

    sect_t *s=&G.sects[sid];
    float *v1=&s->tess[id1<<1];
    float *v2=&s->tess[id2<<1];

    float t3=0, t4=1;

    // t3 
    if (pvs2d_directional) {
	float dst[2]={ev2[0]+pvs2d_dir[0], ev2[1]+pvs2d_dir[1]};
	lineclipseg(dst,ev2,v1,v2,&t3,&t4);
    }
    else {
	lineclipseg(pvs2d_center,ev2,v1,v2,&t3,&t4);
    }
//    if (t3>=1) return 0;
//    if (t3<0) t3=0;
    
    if (t3==t4) return 0;

    // t4
    if (pvs2d_directional) {
	float dst[2]={ev1[0]+pvs2d_dir[0], ev1[1]+pvs2d_dir[1]};
	lineclipseg(ev1,dst,v1,v2,&t3,&t4);
    }
    else {
	lineclipseg(ev1,pvs2d_center,v1,v2,&t3,&t4);
    }
//    if (t4<=0) return 0;
//    if (t4>1) t4=1;

    if (t3==t4) return 0;

//    float p1[2]={v1[0]+(v2[0]-v1[0])*t3, v1[1]+(v2[1]-v1[1])*t3};
//    float p2[2]={v1[0]+(v2[0]-v1[0])*t4, v1[1]+(v2[1]-v1[1])*t4};

    float p1[2], p2[2];
    interp2(p1,v1,v2,t3);
    interp2(p2,v1,v2,t4);

    if (pntsegdist2(pvs2d_center,p1,p2)>pvs2d_r2) {
	t[0]=t3;
	t[1]=t4;
	return 1; // Without recursion
    }

    int eflag=G.sects[sid].e[id1];
    
    if (eflag<0x7fff) { // Internal
	pvs2d_r(pvs,sid,eflag,1.0f-t4,1.0f-t3,0);
    }
    else { // Solid or shared
	pvs2d_r(pvs,sid,id1,t3,t4,0);
    }

    t[0]=t3;
    t[1]=t4;

    return 1;
}


static int genedges(sect_t *s, wall_t *w, float *ev1, float *ev2, int side, int asinternal, pvs2d_t *pvs)
{
    side^=asinternal; // As internal => different side; As solid => same side
    wall_t *aw=&G.walls[w->adjw];
    sect_t *as=&G.sects[w->adjs];
    for (int i=0;i<aw->ne[side];i++) {
	int i1=aw->elist[side][i];
	int i2=i1%3==2 ? i1-2 : i1+1;

	float *av1=&as->tess[i1*2];
	float *av2=&as->tess[i2*2];

	// project the two input point to adjacent tess edge
	float at1,at2;
	at1=prj2seg(ev2,av1,av2);
	at2=prj2seg(ev1,av1,av2);
	    
	if (at1>=1 || at2<=0) continue;

	if (at1<0) at1=0;
	if (at2>1) at2=1;

	if (asinternal) pvs2d_r(pvs,w->adjs,i1,at1,at2, PVS2D_TYPE_INTERNAL);

	// TODO: shared wall should be two sided
	pvs2d_wall_t wall;
	wall.sid=w->adjs;
	wall.eid=i1;
	wall.flag=!asinternal;
	wall.t1=at1;
	wall.t2=at2;
	arr_push(pvs->walls,wall);
    }

}

static int analyzesharededge(short sid, short eid, float *ev1, float *ev2, float *v3, pvs2d_t *pvs)
{
    // find adjacent tess edges that intx t1t2
    sect_t *s=&G.sects[sid];
    wall_t *w=&G.walls[(s->e[eid]>>16)-1 + s->firstw];
    assert(w->adjw!=-1);
//    wall_t *aw=&G.walls[w->adjw];
//    sect_t *as=&G.sects[w->adjs];

    float wv1[2]={w->x,w->y};
    float wv2[2]={G.walls[w->nextco].x, G.walls[w->nextco].y};

    float sn=sign(v3,wv1,wv2);
    int side=w->reverse ? sn>0 : sn<0;

    genedges(s,w,ev1,ev2,side,0,pvs); // As solid walls    
    genedges(s,w,ev1,ev2,side,1,pvs); // As internal walls
}

static void add2florpvs(int hit1, int hit2, float t1, float t2, float t3[2], float t4[2], short sid, short i1, short i2, short i3, pvs2d_t *pvs)
{
    assert(hit1||hit2);

    pvs2d_flor_t f;
    f.sid=sid;
    f.t[0]=t1;  f.eid[0]=i1;
    f.t[1]=t2;  f.eid[1]=i1;
    int nv;
    if (hit1&&hit2) {
	if (t2==1) {
	    f.t[2]=0;
	    f.eid[2]=i3;
	    nv=3;
	}
	else {
	    f.t[2]=t3[0]; f.eid[2]=i2;
	    f.t[3]=0;  f.eid[3]=i3;
	    nv=4;
	}
	
	if (t1!=0) {
	    f.t[nv]=t4[1];
	    f.eid[nv++]=i3;
	}
	f.nv=nv;
    }
    else if (hit1) { // p1
	if (t2==1) {
	    f.t[2]=t3[1];
	    f.eid[2]=i2;
	    f.nv=3;
	}
	else {
	    f.t[2]=t3[0];
	    f.t[3]=t3[1];
	    f.eid[2]=f.eid[3]=i2;
	    f.nv=4;
	}
    }
    else { // hit2 -- p2
	if (t1==0) {
	    f.t[2]=t4[0];
	    f.eid[2]=i3;
	    f.nv=3;
	}
	else {
	    f.t[2]=t4[0];
	    f.t[3]=t4[1];
	    f.eid[2]=f.eid[3]=i3;
	    f.nv=4;
	}
    }

    // add to floor PVS
    arr_push(pvs->flors,f);
}

static int pvs2d_r(pvs2d_t *pvs, short sid, short eid, float t1, float t2, int overloadtype)
{
    if (t1==t2) return;

    if (t1>=t2) {
	printf("t1=%f, t2=%f\n", t1,t2);
	assert(t1<t2);
	return;
    }

    sect_t *s=G.sects+sid;
    if ((s->e[eid]>=0x7fff0000 && overloadtype==0) || overloadtype==PVS2D_TYPE_SOLID) { // solid wall
	// add wall to wall PVS (wall_id,t1,t2)
	pvs2d_wall_t wal={sid,eid,SOLID_WALL,t1,t2};
	arr_push(pvs->walls,wal);
	return;
    }

    // the other vertices of the triangle
    int id2=eid%3==2 ? eid-2 : eid+1;
    int id3=id2%3==2 ? id2-2 : id2+1;
    
    // project the 3rd vertex onto the edge
    float *v1=s->tess + 2*eid;
    float *v2=s->tess + 2*id2;
    float *v3=s->tess + 2*id3;

    float dx=v2[0]-v1[0], dy=v2[1]-v1[1];
    float ev1[2]={v1[0]+dx*t1, v1[1]+dy*t1};
    float ev2[2]={v1[0]+dx*t2, v1[1]+dy*t2};
    
    if ((s->e[eid]<0x7fff&&overloadtype==0) || overloadtype==PVS2D_TYPE_INTERNAL) { //internal edge/shared edge that there's no need to scan
	float t3[2],t4[2];
	int hit1=frusculledge(ev1,ev2,sid,id2,t3,pvs);
	int hit2=frusculledge(ev1,ev2,sid,id3,t4,pvs);

	add2florpvs(hit1,hit2,t1,t2,t3,t4,sid,eid,id2,id3,pvs);
    }
    else { // shared edge
	analyzesharededge(sid,eid,ev1,ev2,v3,pvs);
    }
}

/*
  returns:
  0 -- not in the triangle
  1 -- inside
  2 -- on one of the edges
  3 -- at one of the corners

 */
static int pntintri(float p[2], float tri[6], int tid, int *data)
{
    int nz=0;
    int zid[3];
    for (int i=0;i<3;i++) {
	float *v1=tri+2*i;
	float *v2=tri+2*((i+1)%3);
	float sn=sign(p,v1,v2);
	if (sn<0) return 0;
	if (sn==0) zid[nz++]=i;
    }
    

    *data=tid*3;

    // On an edge
    if (!nz) return 1;
    if (nz==1) {
	*data+=zid[0];
	return 2;
    }

    // At a corner
    assert(nz==2);
    if (zid[0]==0) {
	if (zid[1]!=2) (*data)+=1;
    }
    else {
	(*data)+=2;
    }
    return 3;
}


/*
  quick sort for unsigned int
 */
static inline void swap32u(uint *a, uint *b)
{
    uint t=*a; *a=*b; *b=t;
}

static inline void swap16(short *a, short *b)
{
    short t=*a; *a=*b; *b=t;
}


static void qsort32u_r(uint a[], int b, int e)
{
    if (e>b+1) {
	int p=a[b], l=b+1, r=e;
	while (l<r) {
	    if (a[l]>p) swap32u(&a[l], &a[--r]); // >pivot element => append it to the end of the array
	    else ++l;
	}

	// now l==r. Any element before r <=pivot
	swap32u(&a[--l], &a[b]);
	qsort32u_r(a,b,l);
	qsort32u_r(a,r,e);
    }
}

static void qsort16_r(short a[], int b, int e)
{
    if (e>b+1) {
	int p=a[b], l=b+1, r=e;
	while (l<r) {
	    if (a[l]>p) swap32u(&a[l], &a[--r]); // >pivot element => append it to the end of the array
	    else ++l;
	}

	// now l==r. Any element before r <=pivot
	swap16(&a[--l], &a[b]);
	qsort16_r(a,b,l);
	qsort16_r(a,r,e);
    }
}


void qsort32u(uint a[], int n)
{
    qsort32u_r(a,0,n);
}

void qsort16(short a[], int n)
{
    qsort16_r(a,0,n);
}


/*
  Ray casting -- Find the hit wall/floor

  return values:
  

 */


/*
  Called by incfloorslope(), incceilingslope(), floorhinge(), ceilinghinge()

  
 */
static void updfloorceilingeqations(short sid)
{
    
}


void incfloorslope(short sid, int deg)
{
}

void incceilingslope(short sid, int deg)
{
}

void floorhinge(short sid, short wid)
{
}

void ceilinghinge(short sid, short wid)
{
}

/*
  Plane-segment intersection
 */
int plnsegintx(float p[4], float s1[3], float s2[3], float *ret)
{
    float denom=p[0]*(s1[0]-s2[0]) + p[1]*(s1[1]-s2[1]) + p[2]*(s1[2]-s2[2]);
    if (denom==0) { // parallel the plane
	return 0;
    }

    float numer=p[0]*s1[0] + p[1]*s1[1] + p[2]*s1[2] + p[3]*s1[3];
    float t=numer/denom;
    *ret=t;
    return t>=0 && t<=1;
}


/*
  Hit solid shared wall?

  Input:
  p[2] -- point on the tess edge
  ref[2] -- reference point
  s1 -- sector id which ref[] is in
  e1 -- tess edge id which p[] is on
  
  Output:
  s2 -- s1's adjacent sector id
  e2 -- e1's adjacent tess edge in s2's adjacent sector
  wid -- e1's wall's adjacent wall's id

  Return:
  0 -- not hit
  1 -- hit

 */

// TODO 4-3
static int hitsolidsw(float p[2], float ref[2], short s1, short e1, short *s2, short *e2, short *wid)
{
    sect_t *s=&G.sects[s1];
    wall_t *w=&G.walls[(s->e[e1]>>16)-1 + s->firstw];
    wall_t *aw=&G.walls[w->adjw]; // adjacent wall
    sect_t *as=&G.sects[w->adjs]; // adjacent sector
    float wv1[2]={w->x,w->y};
    float wv2[2]={G.walls[w->nextco].x, G.walls[w->nextco].y};
    float sn=sign(ref,wv1,wv2);
    int side=w->reverse ? sn>0 : sn<0;

    for (int i=0;i<aw->ne[!side];i++) {
	int i1=aw->elist[!side][i];
	int i2=i1%3==2 ? i1-2 : i1+1;
	float *av1 = &as->tess[i1*2];
	float *av2 = &as->tess[i2*2];

	if (pntonseg(p,av1,av2)) {
	    *s2=w->adjs;
	    *e2=i1;
	    *wid=aw - &G.walls[as->firstw];
	    return 0;
	}
    }

    return 1;

}


static int pntaboveseg(float p[3], float s1[3], float s2[3])
{
    float pln[4]={s2[0]-s1[0], s2[1]-s1[1], 0.0f};
    pln[3]= -p[0]*pln[0] - p[1]*pln[1];
    
    float t;
    plnsegintx(pln,s1,s2,&t);
    float z=s1[2] + (s2[2]-s1[2])*t;

    return z>p[2];
}

/*
  Ray-SharedWall test

  Input:
  p[3] -- point coordinates
  sid -- sector_id the point is currently in
  eid -- tess_edge_id the point is currently on

  Output:
  adjsid -- 


  Return:
  0 -- No intersection
  1 -- Hit solid wall
  2 -- Hit lower (floor) shared wall
  3 -- Hit upper (ceiling) shared wall

 */
/*
  Test whether a point is on the upper shared wall
 */
static int hituppersw(float p[3], short wid)
{
    wall_t *w1=&G.walls[wid];
    wall_t *w2=&G.walls[w1->nextco];

    float dx=w2->x - w1->x;
    float dy=w2->y - w1->y;
    float t= fabsf(dx)>fabsf(dy) ? (p[0] - w1->x)/dx : (p[1] - w1->y)/dy; // project point to the shared wall (a 3d line)
    float z= w1->ceilingz + (w2->ceilingz - w1->ceilingz)*t;

    return p[2]>z;
}

/* 
   Test whether a point is on the lower shared wall
 */
static int hitlowersw(float p[3], short wid)
{
    wall_t *w1=&G.walls[wid];
    wall_t *w2=&G.walls[w1->nextco];

    float dx=w2->x - w1->x;
    float dy=w2->y - w1->y;
    float t= fabsf(dx)>fabsf(dy) ? (p[0] - w1->x)/dx : (p[1] - w1->y)/dy; // project point to the shared wall (a 3d line)
    float z= w1->floorz + (w2->floorz - w1->floorz)*t;

    return p[2]<z;
}

static int rayhitsw(float co[3], float dst[3], short s1, short e1, short *adjs, short *adje, float t)
{
    float p[3]={co[0]+(dst[0]-co[0])*t, co[1]+(dst[1]-co[1])*t, co[2]+(dst[2]-co[2])*t};
    int s2,e2,wid;
    if (hitsolidsw(p,co,s1,e1,&s2,&e2,&wid)) {
	return HIT_SOLID;
    }

    if (hitlowersw(p,wid)) {
	return HIT_LOWER;
    }
    else if (hituppersw(p,wid)) {
	return HIT_UPPER;
    }
    else { // Hit nothing
	*adjs=s2;
	*adje=e2;
	return HIT_NOTHING;
    }
}


/*
  Return: 
  0 -- No intersections
  other -- Offset to the input edge's index

 */
static int rayhitedge(float co[3], float dst[3], short sid, short eid, short *eid2, float tmin, float tmax, float *t, float *tt)
{
    // Shorthands
    sect_t *s=&G.sects[sid];
    int id[3];
    id[0]=eid%3==2 ? eid-2 : eid+1;
    id[1]=id[0]%3==2 ? id[0]-2 : id[0]+1;
    id[2]=eid;
    
    // Intersect the ray with the other two tess edges
    for (int i=0;i<2;i++) {
	float t1,t2;
	float *v1=&s->tess[id[i]*2];
	float *v2=&s->tess[id[i+1]*2];
	
	if (segsegintx(co,dst,v1,v2,&t1,&t2)) { // intesect with tess edge
	    if (t1>=tmin && t1<=tmax) {
		*t=t1;
		*tt=t2;
		*eid2=id[i];
		return i+1;
	    }
	}
    }

    return 0;

}

static int rayhitfloorceiling(float co[3], float dst[3], short sid, float t1, float t2, float *t)
{
    sect_t *s=&G.sects[sid];
    if (sect_matchtick(s)) {
	if (s->hittype==0 || (s->hitt<t1 || s->hitt>t2)) { // Ray will hit nothing
	    return HIT_NOTHING;
	}
	else {
	    *t=s->hitt;
	    return s->hittype;
	}
    }

    float tf,tc;
    int hitf,hitc;
    hitf=plnsegintx(s->floorpln,co,dst,&tf);
    hitc=plnsegintx(s->ceilingpln,co,dst,&tc);
    int firsthit=HIT_NOTHING;
    float tmin=100;
    if (hitf) {
	tmin=tf;
	firsthit=HIT_FLOOR;
    }

    if (hitc) {
	if (tc<tmin) {
	    tmin=tc;
	    firsthit=HIT_CEILING;
	}
    }

    *t=tmin;
    sect_copytick(s);
    s->hittype=firsthit;
    s->hitt=tmin;

    return (firsthit!=HIT_NOTHING && tmin>=t1 && tmin<=t2) ? firsthit : HIT_NOTHING;
}


/*
  1. find which triangle the starting point is in
  2. compute intersection between ray & floor/ceiling
  3. if hit the floor/ceiling inside the triangle, return
     otherwise find a edge

  4. for the chosen edge, find one of the other two edges to continue

  return value:
  -1 -- nothing
  0 -- floor
  1 -- wall
  2 -- ceiling
  3 -- shared wall
  

  hitt -- pos + dir*t is the intersection point
  data -- sector id if floor/ceiling is hit;
          wall id if wall/shared wall is hit;

 */


/*
  Return:


 */
static int raybgntri(float co[3], float dst[3], short sid, short id1, float *t, rayintx_t *intx)
{
    // degenerate case -- v3 is on the ray
    // simply choose tess edge v2v3 to trace

    // 1.1 ray intersects with three edges
    sect_t *s=&G.sects[sid];
    float t1,t2;
    int hite=-1; // hit tess edges
    short id2=id1%3==2 ? id1-2 : id1+1;
    short id3=id2%3==2 ? id2-2 : id2+1;
    float *v[4]={&s->tess[2*id1], &s->tess[2*id2], &s->tess[2*id3], &s->tess[2*id1]};
    for (int i=0;i<3;i++) {
	if (segsegintx(co,dst,v[i],v[i+1],&t1,&t2)) { // start tracing from this tess edge
	    hite=i;
	    break;
	}
    }

    // hit floor/ceiling
    int hitfc=rayhitfloorceiling(co,dst,sid,0,t1,&t2);
    
    if (hite==HIT_NOTHING && hitfc==HIT_NOTHING) { // hit nothing
	return -1;
    }

    if (hitfc>=0) { // hit floor/ceiling
	return hitfc;
    }

    // hit tess edge
    return 1;
}


/*
  TODO: not finished

 */
static int raybgnedge(float co[3], float dst[3], short sid, short id1, float *t, rayintx_t *intx)
{
    sect_t *s=&G.sects[sid];
    int id2=id1%3==2?id1-2:id1+1;
    float sn=sign(dst,&s->tess[id1*2],&s->tess[id2*2]);

    int eflag=s->e[id1];

    if (eflag>=0x7fff0000) { // solid edge
	if (sn>=0) { // facing it => visibility is blocked
	    *t=0; // TODO: t
//	    *data=s->firstw + (eflag&0xffff); // wall id
	    return 1;
	}
	else {	// else goto edge tracing
//	    *data=id1;
	    return 0;
	}
    }
    else if (eflag>=0 && eflag<0x7fff) { // internal edge
	// choose the correct side
//	*data = sn>=0 ? eflag : id1;
    }
    else { // shared
	if (sn>=0) { // facing shared wall
	    // find the interval that contains the hit point
	    float *v3=&s->tess[(id2%3==2? id2-2 : id2+1)*2];
	    short s2,e2,wid;
	    if (hitsolidsw(co,v3,sid,id1,&s2,&e2,&wid)) {
		// *data=wid;
		return 1;
	    }
	    else {
		if (hituppersw(co,wid)) {
		    // *data=wid;
		    return 3;
		}
		else if (hitlowersw(co,wid)) {
		    // *data=wid;
		    return 2;
		}

		//data[0]=s2;
		//data[1]=e2;
		return 0;
	    }
	}
	else { // easy case: 
	    //data[0]=id1;
	    return 0;
	}
    }

}

/*
  see illustrate
 */
static int pntinangle(sect_t *s, float p[3], short i1)
{
    int i2=i1%3==2 ? i1-2 : i1+1;
    int i3=i2%3==2 ? i2-2 : i2+1;

    float *v1=&s->tess[i1*2];
    float *v2=&s->tess[i2*2];
    float *v3=&s->tess[i3*2];

    float sn12=sign(p,v1,v2);
    float sn23=sign(p,v2,v3);
    
    return sn12>=0 && sn23>=0;
}


/*
  Input:
  co[3] -- Ray begin
  dst[3] -- Ray end
  s1 -- Angle's sector
  i1 -- Angle's vertex id

  Output:
  rays -- Sector that contains the ray
  rayi -- Angle that contains the ray

  Return:
  0 -- Not found
  1 -- Found

 */

static int findraytri(float co[3], float dst[3], short s1, short i1, short *rays, short *rayi)
{
    // First search i1i2
    sect_t *s=&G.sects[s1];
    short i12=i1;
    while (1) {
	int in=pntinangle(s,dst,i12);
	
	if (in) { // Gotcha!
	    return 1;
	}

	int eflag=s->e[i12];
	if (eflag>=0 && eflag<0x7fff) { // Internal tess edge
	    i12=eflag%3==2 ? eflag-2 : eflag+1;
	}
	else if (eflag>=0x7fff0000) { // Solid wall
	    break; // Not found, search in the other direction
	}
	else { // Shared wall
	    // Two cases:

	    // 1) co[] is on the adjacent tess edge's edge
	    
	}
    }


    // Then search i2i3
}

/*
  Especially necessary for shared walls

 */
static int raybgnpnt(float co[3], float dst[3], short sid, short id1, float *t,rayintx_t *intx)
{
    // Find first triangle that contains the ray
    sect_t *s=&G.sects[sid];
    int id2=id1%3==2 ? id1-2 : id1+1;
    int id3=id2%3==2 ? id2-2 : id2+1;
	    
    float *v1=&s->tess[id1*2];
    float *v2=&s->tess[id2*2];
    float *v3=&s->tess[id3*2];

    float sn12=sign(dst,v1,v2);
    float sn23=sign(dst,v2,v3);
	    
    int hite;

#if 0
    if (sn12==0) { // degenerate case:
	t=prj2seg(v2,co,dst);
//	eid=id3;
    }
    else if (sn23==0) {
	t=prj2seg(v3,co,dst);
//	eid=id2;
    }
    else {
	if (sn12>0 && sn23>0) {
	    hite=segsegintx(co,dst,v2,v3);
		    
	    // shared edge

//	    eid=id2;
	}
	else {
	    if (sn12<0) {
			
	    }
	}
    }

    if (t>=1) {
	return;
    }
    sid=s1;
#endif

}


int ray2d(rayinfo_t *info, short s1, float *co, float *dst, float *hitt, rayintx_t *intx)
{
    if (intx) {
	arr_popall(intx->segments);
    }

    // ray clipped by walls

    // ray clipped by floor/ceiling

    // degenerate cases

    // when a sector is visited, if .visited==0 than put it in the svisited[] array, and set .visited=1

    // compute t for each encountered sector's floor/ceiling
    
    // 1. find which triangle the start point is in
    sect_t *s = &G.sects[s1];
    float t=0;
    int sid=s1,eid=-1;
    float tt1,tt2;
    int ret=HIT_NOTHING;
    for (int i=0;i<s->nt;i++) {
	int data;
	int type = pntintri(co,&s->tess[i*6],i,&data);
	if (type) {
	    if (type==1) { // triangle
		// cast ray in reverse direction to find the edge to trace
		float invdst[2];
		sub2(invdst,co,dst);
		add2(invdst,co,invdst);
		for (int j=0;j<3;j++) {
		    int k=i*6 + j*2;
//		    float invt=rayedgeintx(co,invdst, &s->tess[k], &s->tess[j==2 ? k-4 : k+2] );

		    float invt, edget;
		    segsegintx(co,invdst,&s->tess[k], &s->tess[j==2 ? k-4 : k+2], &invt, &edget);
		    if (invt>0 && edget>=0 && edget<=1) {
			eid=k>>1;
			break;
		    }
		}
	    }
	    else if (type==2) {
	    }
	    else if (type==3) {
	    }

	    break;
	}
	
    }


    if (eid==-1) {
	ret=HIT_NOTHING;
	goto bye;
    }

    // -------------- Edge tracing -----------------------
    while (t<=1) {
	float t1;
	short eid2;
	int hite=rayhitedge(co,dst,sid,eid,&eid2,t,1,&t1,&tt2); // hite=0,1,2

	if (intx) {
	    // Add intersection
	    rayseg_t seg;
	    seg.sid = sid;
	    seg.e1 = eid;
	    seg.e2 = eid2;
	    seg.t1=tt1;
	    seg.t2=tt2;
	    tt1=tt2;

	    arr_push(intx->segments,seg);
	    if (!hite) { // Done
		intx->ns=arr_len(intx->segments);
		intx->dst[0]=dst[0];
		intx->dst[1]=dst[1];
		ret=HIT_NOTHING;
		goto bye;
	    }
	}
	else {
	    float t2;
	    int hit;
	    if (hite) {
		hit=rayhitfloorceiling(co,dst,sid,t,t1,&t2);
	    }
	    else {
		hit=rayhitfloorceiling(co,dst,sid,t,1,&t2);
		if (hit==HIT_NOTHING) { // Hit nothing --> done
		    ret=HIT_NOTHING;
		    goto bye;
		}
	    }
	    

	    // TODO: 3.24
	    if (hit==HIT_FLOOR || hit==HIT_CEILING) { // Hit floor first --> done
		info->id = sid;
		ret=hit;
		goto bye;
	    }
	}

	// ------------- Hit tess edge first ---------------
	int eflag=s->e[eid2];

	if (eflag>=0x7fff0000) { // solid wall --> done
	    ret=HIT_SOLID;
	    info->id = s->firstw + (eflag&0xffff);
	    goto bye;
	}
	else if (eflag>=0 && eflag<0x7fff) { // internal wall --> continue
	    t=t1;
	    eid=eflag;
	    continue;
	}
	else { // shared wall
	    short newsid,neweid;
	    int hitsw=rayhitsw(co,dst,sid,eid2,&newsid,&neweid,t1);

	    if (hitsw==HIT_SOLID) { // hit solid wall --> done
		info->id = s->firstw + (s->e[eid2]>>16) - 1;
		ret=HIT_SOLID;
		goto bye;
	    }
	    else if (hitsw==HIT_UPPER || hitsw==HIT_LOWER) { // hit lower/upper (floor/ceiling) shared wall --> done
		if (!intx) {
		    info->id = s->firstw + (s->e[eid2]>>16) - 1;		    
		    ret=hitsw;
		    goto bye;
		}
	    }
	    else { // hit nothing
		sid=newsid;
		eid=neweid;
		t=t1;
		s=&G.sects[sid];
		continue;
	    }
	}
    }

  bye:
    sect_tick();
    info->type = ret;
    if (ret==HIT_NOTHING) {
	info->id = sid;
    }
    return ret;
}



static int rayhitcircle(float v1[2], float v2[2], float pos[2], float r, float *t)
{
    // possibly two intersections
    float dx=v2[0]-v1[0];
    float dy=v2[1]-v1[1];

    float d1[2]={v1[0]-pos[0], v1[1]-pos[1]};
    float dd=dot2(d1,d1);
    float rr=r*r;

    if (dd<rr) {    // p1
	*t=0;
	return 1;
    }

    if (dd==rr) {   // p2
	float tangent[2]={-dy,dx};  // p3
	if (cross2(d1,tangent) >= 0.0f) {
	    return 0;
	}
    }


    float a=dx*dx + dy*dy;
    float b=2.0f*(dx*d1[0] + dy*d1[1]);
    float c=pos[0]*pos[0] + pos[1]*pos[1] + v1[0]*v1[0] + v1[1]*v1[1] - 2.0f*(pos[0]*v1[0]+pos[1]*v1[1]) - r*r;

    float k=b*b-4.0*a*c;
    if (k<0) return 0;

    float invdenom=.5f/a;
    if (k==0) { // p4
	return 0;
    }

    k=sqrtf(k);
    float t1=(-b-k)*invdenom;
    float t2=(-b+k)*invdenom;

    if (t1>1.0f || t2<0.0f) { // p5
	return 0;
    }

    *t=t1; // p0 -- valid intersection time
    return 1;
}

/*
  co[], dst[] can be rewritten
 */
static int rayhitcapsule(short sid, short eid, float et1, float et2, float co[3], float dst[3], float r, float *t)
{
    sect_t *s=&G.sects[sid];
    short id2=eid%3==2 ? eid-2 : eid+1;

    float *ev1=&s->tess[eid<<1];
    float *ev2=&s->tess[id2<<1];

    float v1[2]={ev1[0]+(ev2[0]-ev1[0])*et1, ev1[1]+(ev2[1]-ev1[1])*et1};
    float v2[2]={ev1[0]+(ev2[0]-ev1[0])*et2, ev1[1]+(ev2[1]-ev1[1])*et2};

    float v12[2]={v2[0]-v1[0], v2[1]-v1[1]};
    float len=length2(v12);

    float d[2]={dst[0]-co[0], dst[1]-co[1]};

    float wallofs[2];
    float k=r/len;
    if (cross2(d,v12)>0.0f) { // p3.1
	wallofs[0]=v12[1]*k;
	wallofs[1]=-v12[0]*k;
    }
    else { // p3.2
	wallofs[0]=-v12[1]*k;
	wallofs[1]=v12[0]*k;
    }

    float shiftv1[2]={v1[0]+wallofs[0], v1[1]+wallofs[1]};
    float shiftv2[2]={v2[0]+wallofs[0], v2[0]+wallofs[1]};

    float tmin=100;
    float t1,t2;
    if (segsegintx(co,dst,shiftv1,shiftv2,&t1,&t2)) { // Circle vs. edge
	tmin=t1;
    }

    if (rayhitcircle(co,dst,v1,r,&t1)) { // Circle vs. point1
	if (t1<tmin) tmin=t1;
    }

    if (rayhitcircle(co,dst,v2,r,&t2)) { // Circle vs. point2
	if (t2<tmin) tmin=t2;
    }
    
    if (tmin>1.0f) {
	*t=1.0f;
	return 0;
    }

    *t=tmin;
    return 1;
}


/*
  Static collision detection between [player's bounding circle] & [solid walls]

  Can be used to place player in the map
 */
int circlestucked(float *co, float r)
{
}


typedef struct {
    int i;
    float t;
} hitbycircle_t;

int cmphit(hitbycircle_t *a, hitbycircle_t *b)
{
    return a->t==b->t ? 0 : (a->t<b->t ? -1 : 1);
}


/*
  Floor/ceiling can block the motion of the viewer too.

  
 */
static int sweepcircleztest(pvs2d_flor_t *flor, float *co, float *dst, float r, float height, float *t1, float *t2)
{
    
}


/*
  Multiple colliding points
 */
static int sweepcircle(float *co, float *dst, float r, pvs2d_t *pvs, float *t)
{
    float t1;
    float tmin=100;
    int emin=-1, emin2=-1;
    int wid;
    int sz=arr_len(pvs->walls);

    hitbycircle_t hit[1000];
    int h=0;

    // Add all possible colliding edges to a list (which will be sorted later)
    for (int i=0;i<sz;i++) {
	pvs2d_wall_t *w=&pvs->walls[i];
	if (rayhitcapsule(w->sid,w->eid,w->t1,w->t2,co,dst,r,&t1)) {
	    hit[h].i=i;
	    hit[h++].t=t1;
	}
    }

    // Sort the edge list by t
    qsort(hit,h,sizeof(*hit),cmphit);

    // Traverse all edges
    for (int i=0;i<h;i++) {
	pvs2d_wall_t *w=&pvs->walls[hit[i].i];
	sect_t *s=&G.sects[w->sid];

	int eflag = s->e[w->eid];
	if (eflag >= 0x7fff0000) { // Solid edge => done
	}

	// Shared edge
	

    }


    // Find the first stuck tess triangle
    sz=arr_len(pvs->flors);
    for (int i=0;i<sz;i++) {
    }

    


    return wid;
}


static int sweeprays(rayintx_t *ray, float *dir, pvs2d_t *pvs)
{
    for (int i=0;i<ray->ns;i++) {
	rayseg_t *seg = &ray->segments[i];
	sect_t *s=&G.sects[seg->sid];

	if (seg->e1==seg->e2) {
	    pvs2d_r(pvs,seg->sid, seg->e1, seg->t1, seg->t2, 0);
	}
	else {
	    // Compute 3rd vertex's id (p1)
	    short n = seg->e1/3;
	    short e3= 9*n + 3 - seg->e1 - seg->e2;

	    // For each tess edge, cast beam
	    float *v1,*v2;
	    float edge[2];
	    float t1,t2;
	    float xx1[2], xx2[2]; // p4

	    // Edge1
	    v1=&s->tess[seg->e1<<1];
	    v2=&s->tess[seg->e1%3==2 ? (seg->e1-2)<<1 : (seg->e1+1)<<1];
	    sub2(edge,v2,v1);
	    interp2(xx1,v1,v2,seg->t1);
	    if (cross2(dir,edge)>0.0f) { // Cast this edge. See p1
		if (dot2(dir,edge)<0.0f) { // p3
		    t1=0;    
		    t2=seg->t1;
		}
		else {
		    t1=seg->t1;
		    t2=1;
		}
		pvs2d_r(pvs, seg->sid, seg->e1, t1, t2, 0);
	    }

	    // Edge2 (same algorithm as edge1)
	    v1=&s->tess[seg->e2<<1];
	    v2=&s->tess[seg->e2%3==2 ? (seg->e2-2)<<1 : (seg->e2+1)<<1];
	    sub2(edge,v2,v1);
	    interp2(xx2,v1,v2,seg->t2);
	    if (cross2(dir,edge)>0.0f) { // Cast this edge. See p1
		if (dot2(dir,edge)<0.0f) { // p3
		    t1=0;    
		    t2=seg->t2;
		}
		else {
		    t1=seg->t2;
		    t2=1;
		}
		pvs2d_r(pvs, seg->sid, seg->e2, t1, t2, 0);
	    }

	    // Edge3
	    v1=&s->tess[e3<<1];
	    v2=&s->tess[e3%3==2 ? (e3-2)<<1 : (e3+1)<<1];
	    sub2(edge,v2,v1);
	    if (cross2(dir,edge)>0.0f) { // Cast this edge. See p1
		t1=rayedgeintx(xx1,dir,v1,v2);
		t2=rayedgeintx(xx2,dir,v1,v2);

		if (t1>t2) {
		    float tmp=t1; 
		    t1=t2;
		    t2=tmp;
		}

		if (t1<0) t1=0;
		if (t2>1) t2=1;
		pvs2d_r(pvs,seg->sid, e3, t1, t2, 0);
	    }

	}
    }
}


float *scale2(float *b, float *a, float t)
{
    b[0]=a[0]*t;
    b[1]=a[1]*t;
    return b;
}

/*
  See p1
 */
static void computeslidedir(pvs2d_wall_t *w, float *co, float *dst, float *t)
{
    sect_t *s=&G.sects[w->sid];
    float *v1=&s->tess[w->eid*2];
    float *v2=&s->tess[w->eid%3==2 ? (w->eid-2)<<1 : (w->eid+1)<<1];
    float v12[2],e[2];
    float d[2];
    sub2(d,dst,co);
    scale2(d,d,1-*t); // d=(dst-co)*(1-t)

    sub2(v12,v2,v1);
    normalize2(e,v12); // e=normalize(v2-v1) => normalized edge direction

    scale2(d,e,dot2(d,e)); // d = e dot(d,e)  => project d onto e
    interp2(co,co,dst,*t);    // Update co[]:  co = co + (dst-co)*t
    add2(dst,co,d); // Update dst[]
    *t=0;
}

/*
  p SHOULD be global variable!

 */
static void bgnpvs(pvs2d_t *p, float r)
{
    assert(p);
    
    arr_popall(p->walls);
    arr_popall(p->flors);
    arr_popall(p->nwalls);
    arr_popall(p->flor_sectlut);
    p->wall_ns=0;

    pvs2d_r2=r*r;
}

static void endpvs(pvs2d_t *p)
{
    assert(p);

    // process the PVS:
    // sort walls by sectors, edge id & t => can combine neighbouring walls
    int nw=arr_len(p->walls);
    qsort(p->walls,nw,sizeof(pvs2d_wall_t),cmpwalls);

    // count num of walls for each sector & combine connected ones
    pvs2d_wall_t *wal=p->walls;
    int ns = nw>0;
    short nws[1000]={1}; // num of walls per sector
    for (int i=1,j=0;i<nw;i++) {
	if (wal[i].sid!=wal[i-1].sid) {
	    nws[ns++]=1;
	    wal[++j]=wal[i];
	}
	else {
	    if (wal[i].sid==wal[i-1].sid &&
		wal[i].eid==wal[i-1].eid && 
		wal[i].t1==wal[i-1].t2) {
		wal[j].t2=wal[i].t2;
	    }
	    else {
		++nws[ns-1];
		wal[++j]=wal[i];
	    }
	}
    }
    arr_pushn(p->nwalls, nws, ns);
    p->wall_ns=ns;

    int nf=arr_len(p->flors);
    struct {short sid,fid;} flr[nf];
    for (short i=0;i<nf;i++) {
	flr[i].sid=p->flors[i].sid;
	flr[i].fid=i;
    }
    
    qsort32u(flr,nf);
    arr_pushn(p->flor_sectlut, flr, nf);
    p->nf=nf;
    
}


/*
  The function assumes that the current view position is valid.
 */
int tracemotion(short s1, short e1, short *s2, short *e2, float co[3], float dst[3], float r)
{
    int wid;
    for (int step=0; step<2; step++) {
	// Cast two rays perpendicular to the motion vector
	float dx=dst[0]-co[0];
	float dy=dst[1]-co[1];
	float d=sqrtf(dx*dx+dy*dy);
	float k=r/d;
	dx*=k; dy*=k;

	// Left: left_dx=-dy, left_dy=dx
	float left[3]={-dy,dx};
	float right[3]={dy,-dx};

	float lt,rt; // left_time, right_time
	rayinfo_t linfo, rinfo;
	ray2d(&linfo, s1,co,left, &lt, &G.lray);
	ray2d(&rinfo, s1,co,right,&rt, &G.rray);
	
	// Sweep line segment in the direction of the motion vector
	float dir[3]={dx*(d+r), dy*(d+r)};
	bgnpvs(&G.tracepvs,r);
	sweeprays(&G.lray, dir, &G.tracepvs);
	sweeprays(&G.rray, dir, &G.tracepvs);

	// Sweep the circle against the walls found, and find the first hit time
	//    Two or more collision points => get stucked
	float t;
	wid = sweepcircle(co,dst,r,&G.tracepvs,&t);
	if (wid < 0) { // Stucked / Hit nothing
	    endpvs(&G.tracepvs);
	    break;
	}

        // Hit something
	// Change the motion direction (slide) & trace again
	if (step==0) {
	    computeslidedir(&G.tracepvs.walls[wid], co, dst, &t);
	}

	endpvs(&G.tracepvs);
    }

    return wid;
}



/*-----------------------------------
  FOV can be:
  1) 0 -- Directional
  2) (0,PI] && 2PI
  ---------------------------------*/
static int pvs2dbgntri(int sid, int tid, float fov, pvs2d_t *pvs)
{
    // Add the triangle itself to the PVS -- TODO: wrong
    pvs2d_flor_t f;
    f.sid=sid;
//    f.nv=3;

    sect_t *sect=&G.sects[sid];
    float x=pvs2d_dir[0], y=pvs2d_dir[1];

    int nv=0;
    int i0=-1;
    // Three edges
    for (int i=0;i<3;i++) {
	int id=3*tid + i;
	int adj=sect->e[id];

	float t1,t2;
	float d;
	if (fov==0) { // Directional
	}
	else if (fov==M_PI*2.0f) {
	    t1=0;
	    t2=1;

	    f.eid[i]=id;
	    f.t[i]=0;
	}
	else { // <PI
	    float *v1=&sect->tess[id<<1];
	    float *v2=&sect->tess[id%3==2 ? (id-2)<<1 : (id+1)<<1];

	    float c=cosf(fov*0.5f);
	    float s=sinf(fov*0.5f);
	    float bgndst[2]={pvs2d_center[0]+x*c+y*s, pvs2d_center[1]+y*c-x*s}; // Rotate -fov/2
	    float enddst[2]={pvs2d_center[0]+x*c-y*s, pvs2d_center[1]+y*c+x*s}; // Rotate fov/2
	    
	    float d1,d2;
	    int intx1=linelineintx(pvs2d_center,bgndst,v1,v2,&d1,&t1);
	    int intx2=linelineintx(pvs2d_center,enddst,v1,v2,&d2,&t2);

	    if (t1==t2) continue; // ignore zero sized portal

	    int hite = d2>0 && t2>=0 && t2<=1;

	    if (!intx1) d1=-1;
	    if (!intx2) d2=-1;

	    if (d1<0 && d2<0) continue; // p2

	    if (d1<0) t1=0;
	    else if (d2<0) t2=1;

	    if (t1<0) t1=0;
	    if (t2>1) t2=1;

	    if (t2<0 || t1>1) continue;

	    if (t1==t2) continue; // ignore zero sized portal

	    if (!nv || (t1!=0 && t2!=1)) {
		f.eid[nv]=f.eid[nv+1]=id;
		f.t[nv]=t1;
		f.t[nv+1]=t2;
		nv+=2;
	    }
	    else {
		int skip1=0,skip2=0;
		if (t1==0) {
		    int j;
		    for (j=0;j<nv;j++) {
			if (f.t[j]==1 && (f.eid[j] == (id%3==0 ? id+2 : id-1))) {
			    skip1=1;
			    break;
			}
		    }
		}
		
		if (t2==1) {
		    int j;
		    for (j=0;j<nv;j++) {
			if (f.t[j]==0 && (f.eid[j] == (id%3==2 ? id-2 : id+1))) {
			    skip2=1;
			    break;
			}
		    }		    
		}

		if (!skip1) {
		    f.eid[nv]=id;
		    f.t[nv++]=t1;
		}

		if (!skip2) {
		    f.eid[nv]=id;
		    f.t[nv++]=t2;
		}
	    }

	    if (hite) {
		f.eid[nv]=-1; // TODO: 3.19
		f.t[nv++]=pvs2d_center[0];
		f.y=pvs2d_center[1];
	    }

	    d = pntsegdist2(pvs2d_center,v1,v2);
	}

	if (d > pvs2d_r2) continue;

	if (adj<0x7fff) pvs2d_r(pvs,sid,adj,1.0f-t2,1.0f-t1,0); // internal edge
	else pvs2d_r(pvs,sid,id,t1,t2,0); // solid/shared edges

    }

    // Add the point itself to the floor pvs
#if 0
    if (i0==-1) i0=nv;
    f.eid[i0]=-1; // TODO: 3.19
    f.t[i0]=pvs2d_center[0];
    f.y=pvs2d_center[1];
#endif
    f.nv=nv;
    arr_push(pvs->flors,f);
}


static int pvs2dbgnedge(short sid, int eid, float *co, pvs2d_t *pvs)
{
    // two edges
    // ignore solid edges
    sect_t *s=&G.sects[sid];
    if (s->e[eid]>=0x7fff0000) { // wall
	pvs2d_r(pvs,sid,eid,0,1, 0);
	pvs2d_r(pvs,sid,eid,0,1, PVS2D_TYPE_INTERNAL);
    }
    else if (s->e[eid]<0x7fff) { // internal
	pvs2d_r(pvs,sid,eid,0,1, 0);
	pvs2d_r(pvs,sid,s->e[eid]&0xffff,0,1, 0);
    }
	    
    // TODO: 
    else { // shared
	// find the interval the point is in
	wall_t *w=&G.walls[s->firstw + (s->e[eid]>>16) - 1];
	wall_t *w2=&G.walls[w->nextco];
	float wv1[2]={w->x,w->y},wv2[2]={w2->x,w2->y};
	float sn=sign(&s->tess[2*eid],wv1,wv2);
	int side=w->reverse ? sn>0 : sn<0;

	float *v1=&s->tess[2*eid];
	float *v2=&s->tess[2*(eid%3==2 ? eid-2 : eid+1)]; 

	wall_t *aw=&G.walls[w->adjw];
	sect_t *as=&G.sects[w->adjs];

	int blockvis=1;
	for (int j=0;j<aw->ne[side];j++) { // adjacent sector, below
	    int j1=aw->elist[side][j];
	    int j2=j1%3==2 ? j1-2 : j1+1;
	    float *av1=&as->tess[2*j1],  *av2=&as->tess[2*j2];
	    if (pntonseg(co,av1,av2)) {
		blockvis=0;
			
		// treated as two internal walls
		float t1,t2;
		t1=prj2seg(av2,v1,v2);
		t2=prj2seg(av1,v1,v2);
		if (t1<0) t1=0;
		if (t2>1) t2=1;
		pvs2d_r(pvs,sid,eid,t1,t2,PVS2D_TYPE_INTERNAL);
			
		float t3,t4;
		t3=prj2seg(v2,av1,av2);
		t4=prj2seg(v1,av1,av2);
		if (t3<0) t3=0;
		if (t4>1) t4=1;
		pvs2d_r(pvs,w->adjs,j1,t3,t4,PVS2D_TYPE_INTERNAL);
			
		// no need to add to wall PVS here -- you cannot see the walls (you're on the wall plane)
		break;
	    }
	}

	if (blockvis) { // treated as a solid wall
	    // find the interval that blocks it
	    for (int j=0;j<aw->ne[1];j++) { // adjacent sector, above
		int j1=aw->elist[0][j];
		int j2=j1%3==2 ? j1-2 : j1+1;
		float *av1=&as->tess[2*j1],  *av2=&as->tess[2*j2];
		if (pntonseg(co,av1,av2)) {
		    // add to wall PVS
		    float t1,t2;
		    t1=prj2seg(av2,v1,v2);
		    t2=prj2seg(av1,v1,v2);
		    if (t1<0) t1=0;
		    if (t2>1) t2=1;
			    
		    pvs2d_wall_t wall;
		    wall.sid=sid;
		    wall.eid=eid;
		    wall.flag=1;
		    wall.t1=t1;
		    wall.t2=t2;
		    arr_push(pvs->walls,wall);
		    break;
		}
	    }
		    
	    pvs2d_r(pvs,sid,eid,0,1,PVS2D_TYPE_INTERNAL);
	}
    }
}

static int pvs2dbgnpnt(short sid, short id1, pvs2d_t *pvs)
{
    // corner view point cannot see stuff in other sectors
    // find all triangles adjacent to the corner
	    
    // the other two vertices in the same triangle:
    sect_t *s=&G.sects[sid];
    int id2,id3;
    float *v1,*v2,*v3;
    id2=(id1+1)%3==0 ? id1-2 : id1+1;
    id3=(id2+1)%3==0 ? id2-2 : id2+1;
    v1=s->tess + 2*id1;	    v2=s->tess + 2*id2;	    v3=s->tess + 2*id3;
    pvs2d_r(pvs,sid,id2,0,1,0);

    int left,right;
    left= id1%3==0 ? id1+2 : id1-1; // mod==0 => 2; mod==1 => 0; mod==2 => 1
    right=id1%3==2 ? id1-2 : id1+1; // mod==0 => 1; mod==1 => 2; mod==0 => 2
	    
    // all left edges
    while (s->e[left]<0x7fff) {
	id1=s->e[left]&0xffff;
	id2=(id1+1)%3==0 ? id1-2 : id1+1;
	id3=(id2+1)%3==0 ? id2-2 : id2+1;
	v1=s->tess + 2*id1;	    v2=s->tess + 2*id2;	    v3=s->tess + 2*id3;

	// add the triangle to PVS
	// nothing to be added to wall PVS
	// add to floor PVS

	pvs2d_r(pvs,sid,id2,0,1, 0);
	left=id3;
    }

    // all right edges
    while (s->e[right]<0x7fff) {
	id1=s->e[right]&0xffff;
	id2=(id1+1)%3==0 ? id1-2 : id1+1;
	id3=(id2+1)%3==0 ? id2-2 : id2+1;
	v1=s->tess + 2*id1;	    v2=s->tess + 2*id2;	    v3=s->tess + 2*id3;

	// add the triangle to PVS

	pvs2d_r(pvs,sid,id3,0,1, 0);
	right=id2;
    }
}



/*
  Given a point and a sector, compute the __exact__ PVS

  propagate the visibility to other triangles & sectors

  Can be used for:
  1. View frustum culling
  2. Light frustum culling
  3. Find possible colliding set

  
 */
int pvs2d(pvs2d_t *pvs, short sid, v3_t pos, v3_t dir, float fov)
{
    bgnpvs(pvs,1000);

    pvs2d_dir[0]=dir.x;	pvs2d_dir[1]=dir.y;
    pvs2d_directional = fov==0.f;

    // project view frustum to 2D
//    if (fov>=360) fov=360;

    // find which triangle the viewpoint is in
    sect_t *s=G.sects + sid;
    float *t=s->tess;
    float co[2]={pos.x,pos.y};
    pvs2d_center[0]=co[0];
    pvs2d_center[1]=co[1];
    for (int i=0;i<s->nt;i++) {
	int info;
	int type=pntintri(co,t,i,&info);

	if (type) {
	    if (type==1) { // triangle
		pvs2dbgntri(sid,i,fov,pvs);
	    }
	    else if (type==2) { // edge
		pvs2dbgnedge(sid,info,co,pvs);
	    }
	    else if (type==3) { // corner -- FOV can be infinite if visibility can get past shared walls =>
		pvs2dbgnpnt(sid,info,pvs);
	    }
	    break;
	}
	
	// otherwise not in this triangle

	t+=6;
    }

    endpvs(pvs);
}

static void editor_updatemousepos(short mx, short my)
{
    cam2_t *c=&G.ed.cam;
    float w=G.screenw/c->zoom;
    float h=G.screenh/c->zoom;
    float x1=c->x - w*.5;
    float y1=c->y - h*.5;
    
    // cursor position filter space
    float x=x1 + w*mx/G.screenw;
    float y=y1 + h*my/G.screenh;
    
    // snap to grid
    int xmin,xmax,ymin,ymax;
    int lv=G.ed.level;
    xmin=((int)floorf(x))>>lv<<lv;
    ymin=((int)floorf(y))>>lv<<lv;
    xmax=xmin+(1<<lv);
    ymax=ymin+(1<<lv);

    G.ed.mousex= (x-xmin)<(xmax-x) ? xmin : xmax;
    G.ed.mousey= (y-ymin)<(ymax-y) ? ymin : ymax;
}

static void editor_mmfunc(SDL_MouseMotionEvent motion)
{
    editor_updatemousepos(motion.x,motion.y);
    int btn=SDL_GetMouseState(NULL, NULL);
    if (btn&SDL_BUTTON(2)) {
	cam2_t *c=&G.ed.cam;
	c->x -= motion.xrel/c->zoom;
	c->y -= motion.yrel/c->zoom;
    }
}

static void editor_mbfunc(SDL_MouseButtonEvent button, int down)
{
    editor_updatemousepos(button.x,button.y);
    if (button.button==SDL_BUTTON_WHEELDOWN) {
	if (down) {
	    G.ed.cam.zoom /=1.1;
	}
    }
    else if (button.button==SDL_BUTTON_WHEELUP) {
	if (down) {
	    G.ed.cam.zoom *=1.1;
	}
    }
    else if (button.button==SDL_BUTTON_RIGHT) {
	if (down) {
	    cam1_setpos(&G.ed.cam1,G.ed.mousex,G.ed.mousey,0);
	    v3_t pos={G.ed.cam1.pos.x,G.ed.cam1.pos.y,G.ed.cam1.pos.z};
	    short undermouse[1000];
	    int loc=chkloc(G.ed.mousex,G.ed.mousey,undermouse,1000);

	    if (loc==LOC_IN) {
		pvs2d(&G.ed.pvs, undermouse[0], pos, G.ed.cam1.rotaxis[1], M_PI_2);
		printf("nf=%i, ns=%i\n", G.ed.pvs.nf, arr_len(G.ed.pvs.walls));
	    }
	}
    }
}

static void editor_keyfunc(SDL_keysym keysym, int down)
{
    if (G.mode!=MODE_EDITOR) return;

    int k=keysym.sym;
    char *cmd = G.ed.kcmd[k];
    if(cmd){
	char *buf = alloca(strlen(cmd) + 2);
	buf[0] = down ? '+' : '-';
	strcpy(buf+1, cmd);
	cmd_addtxt(buf);
    }
    
}

static void editor_kcmd(char *key, char *cmd)
{
    if (!cmd || !key) return;
    int k=keystr2num(key);
    if(G.ed.kcmd[k]) free(G.ed.kcmd[k]);
    G.ed.kcmd[k] = strdup(cmd);
}


static void c_editor_declevel()
{
    --G.ed.level;
    if (G.ed.level<0) G.ed.level=8;
}

static void c_editor_bgnforward()
{
    if (G.mode!=MODE_EDITOR) return;
    cam1_t *c=&G.ed.cam1;
    c->acc[1]=.1;
}

static void c_editor_endforward()
{
    if (G.mode!=MODE_EDITOR) return;
    G.ed.cam1.acc[1]=0;
}

static void c_editor_bgnbackward()
{
    if (G.mode!=MODE_EDITOR) return;
    G.ed.cam1.acc[1]=-.1;
}

static void c_editor_endbackward()
{
    if (G.mode!=MODE_EDITOR) return;
    G.ed.cam1.acc[1]=0;
}

static void c_editor_bgnstrafeleft()
{
    if (G.mode!=MODE_EDITOR) return;
    G.ed.cam1.acc[0]=-.1;
}

static void c_editor_endstrafeleft()
{
    if (G.mode!=MODE_EDITOR) return;
    G.ed.cam1.acc[0]=0;
}

static void c_editor_bgnstraferight()
{
    if (G.mode!=MODE_EDITOR) return;
    G.ed.cam1.acc[0]=.2;
}

static void c_editor_endstraferight()
{
    if (G.mode!=MODE_EDITOR) return;
    G.ed.cam1.acc[0]=0;
}

static void c_editor_bgnrotleft()
{
    if (G.mode!=MODE_EDITOR) return;
    G.ed.cam1.angacc[0]=.5;
}

static void c_editor_endrotleft()
{
    if (G.mode!=MODE_EDITOR) return;
    G.ed.cam1.angacc[0]=0;
}

static void c_editor_bgnrotright()
{
    if (G.mode!=MODE_EDITOR) return;
    G.ed.cam1.angacc[0]=-.5;
}

static void c_editor_endrotright()
{
    if (G.mode!=MODE_EDITOR) return;
    G.ed.cam1.angacc[0]=0;
}

static void c_editor_switchdim()
{
    if (G.mode==MODE_EDITOR) {
	if (G.ed.dim==3) G.ed.dim=2;
	else G.ed.dim=3;
    }
}

static void editor_init()
{
    G.ed.cam.zoom=1;
    G.ed.level=4;
    cmd_add("declevel", c_editor_declevel);
    cmd_add("pushnook", c_editor_pushnk);
    cmd_add("popnook", c_editor_popnk);
    cmd_add("+ed_forward", c_editor_bgnforward);
    cmd_add("-ed_forward", c_editor_endforward);
    cmd_add("+ed_backward", c_editor_bgnbackward);
    cmd_add("-ed_backward", c_editor_endbackward);
    cmd_add("+ed_strafeleft", c_editor_bgnstrafeleft);
    cmd_add("-ed_strafeleft", c_editor_endstrafeleft);
    cmd_add("+ed_straferight", c_editor_bgnstraferight);
    cmd_add("-ed_straferight", c_editor_endstraferight);
    cmd_add("+ed_rotleft", c_editor_bgnrotleft);
    cmd_add("-ed_rotleft", c_editor_endrotleft);
    cmd_add("+ed_rotright", c_editor_bgnrotright);
    cmd_add("-ed_rotright", c_editor_endrotright);
    cmd_add("ed_switchdim", c_editor_switchdim);
    editor_kcmd("g", "declevel");
    editor_kcmd(" ", "pushnook");
    editor_kcmd("backspace", "popnook");
    editor_kcmd("up", "ed_forward");
    editor_kcmd("down","ed_backward");
//    editor_kcmd("left", "ed_strafeleft");
//    editor_kcmd("right", "ed_straferight");
    editor_kcmd("left", "ed_rotleft");
    editor_kcmd("right","ed_rotright");
    editor_kcmd("kp_enter", "ed_switchdim");
//    G.mode=MODE_EDITOR;
    cam1_init(&G.ed.cam1,G.screenw,G.screenh,90);
    G.ed.dim=2;
}


static void drawdumb(uint tex)
{
    glBegin(GL_TRIANGLES);
    // floor
    glColor3f(0,0,.5);
    for (int i=0; i<G.nsects; i++) {
	float *t=G.sects[i].tess;
	for (int j=0; j<G.sects[i].nt; j++) {
	    glVertex3f(t[0],t[1],-10);
	    glVertex3f(t[2],t[3],-10);
	    glVertex3f(t[4],t[5],-10);

	    t+=6;
	}
    }

    // ceiling
    glColor3f(.5,0,0);
    for (int i=0; i<G.nsects; i++) {
	float *t=G.sects[i].tess;
	for (int j=0; j<G.sects[i].nt; j++) {
	    glVertex3f(t[0],t[1],10);
	    glVertex3f(t[2],t[3],10);
	    glVertex3f(t[4],t[5],10);

	    t+=6;
	}
    }
    glEnd();


    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,tex);
    // walls
//    glColor3f(0,.5,0);
    glColor3f(1,1,1);
    for (int i=0; i<G.nsects; i++) {
	int w=G.sects[i].firstw;
	int nc=0;
	glBegin(GL_QUAD_STRIP);
	for (int j=0;j<G.sects[i].nw;j++) {
	    glTexCoord2f(j*3, 0.6);
	    glVertex3f(G.walls[w].x, G.walls[w].y, 10);
	    glTexCoord2f(j*3, 0.4);
	    glVertex3f(G.walls[w].x, G.walls[w].y, -10);

	    if (w+1!=G.walls[w].nextco) {
		int x=G.walls[G.walls[w].nextco].x;
		int y=G.walls[G.walls[w].nextco].y;
		glTexCoord2f((j+1)*3,0.6);
		glVertex3f(x,y,10);
		glTexCoord2f((j+1)*3,0.4);
		glVertex3f(x,y,-10);

		glEnd();
		glBegin(GL_QUAD_STRIP);
	    }
	    ++w;
	}
	glEnd();
    }

    glDisable(GL_TEXTURE_2D);

}



static void drawclever(uint walltex, uint flortex)
{
    glColor3f(1,1,1);
    glEnable(GL_TEXTURE_2D);

    glBindTexture(GL_TEXTURE_2D, flortex);
    //floor/ceiling PVS
    int sz=arr_len(G.ed.pvs.flors);
    for (int i=0;i<sz;i++) {
	float x[10],y[10];

	pvs2d_flor_t *f=&G.ed.pvs.flors[i];
	for (int j=0;j<f->nv;j++) {
	    sect_t *s=&G.sects[f->sid];
	    if (f->eid[j]>=0) {
		float *v1=s->tess + f->eid[j]*2;
		float *v2=(f->eid[j]%3==2) ? v1-4 : v1+2;
		x[j]=v1[0]+(v2[0]-v1[0])*f->t[j];
		y[j]=v1[1]+(v2[1]-v1[1])*f->t[j];
		
	    }
	    else {
		x[j]=f->t[j];
		y[j]=f->y;
	    }
	}

	// floor
	glBegin(GL_POLYGON);
	for (int j=0;j<f->nv;j++) {
	    glTexCoord2f(x[j]*0.01, y[j]*0.01);
	    glVertex3f(x[j],y[j],-20);
	}
	glEnd();

	// ceiling
	glBegin(GL_POLYGON);
	for (int j=0;j<f->nv;j++) {
	    glTexCoord2f(x[j]*0.01, y[j]*0.01);
	    glVertex3f(x[j],y[j],20);
	}
	glEnd();
    }

    //wall PVS
    glBindTexture(GL_TEXTURE_2D, walltex);
    glColor3f(1,1,1);
    sz=arr_len(G.ed.pvs.walls);
    pvs2d_wall_t *pvsw=G.ed.pvs.walls;
    if (pvsw) {
	glBegin(GL_QUADS);
	for (int i=0;i<G.ed.pvs.wall_ns;i++) {
	    sect_t *s=&G.sects[pvsw->sid];
	    for (int j=0;j<G.ed.pvs.nwalls[i];j++,pvsw++) {
		if (pvsw->flag !=SOLID_WALL) continue;

		int id1=pvsw->eid;
		int id2=pvsw->eid%3==2 ? id1-2 : id1+1;
		float *ev1=&s->tess[id1*2];
		float *ev2=&s->tess[id2*2];
	    
		float dx=ev2[0]-ev1[0];
		float dy=ev2[1]-ev1[1];
	    
		float v1[2]={ev1[0]+pvsw->t1*dx, ev1[1]+pvsw->t1*dy};
		float v2[2]={ev1[0]+pvsw->t2*dx, ev1[1]+pvsw->t2*dy};

		// fake solid texture

		// 1m == 0.01 texcoord
		float umin,umax,vmin,vmax;
		if (fabsf(dx)<fabsf(dy)) {
		    umin=v1[1]*0.01;
		    umax=v2[1]*0.01;
		}
		else {
		    umin=v1[0]*0.01;
		    umax=v2[0]*0.01;
		}
		vmin=0.4;
		vmax=0.6;

	    
		if (IS_SOLID(s->e[pvsw->eid]) && 
		    G.ed.rayinfo.type == HIT_SOLID &&
		    G.ed.rayinfo.id == s->firstw + (s->e[pvsw->eid]&0xffff) ) {
		    glColor3f(1,0,0);
		}

		glTexCoord2f(umin,vmax); glVertex3f(v1[0],v1[1],20);
		glTexCoord2f(umin,vmin); glVertex3f(v1[0],v1[1],-20);
		glTexCoord2f(umax,vmin); glVertex3f(v2[0],v2[1],-20);
		glTexCoord2f(umax,vmax); glVertex3f(v2[0],v2[1],20);

		glColor3f(1,1,1);
	    }
	}
	glEnd();
    }

    glDisable(GL_TEXTURE_2D);

}


/*

 */


static void editor_draw3d()
{
    static uint walltex=0;
    static uint flortex=0;
    if (!walltex) {
	int w,h,bbp;
	uchar *data=readtga("wall.tga",&w,&h,&bbp);
	walltex=tex3b(w,h,011,data);
	free(data);

	data=readtga("flor.tga",&w,&h,&bbp);
	flortex=tex3b(w,h,011,data);
	free(data);
    }

    glEnable(GL_DEPTH_TEST);
    rnd_clearall();

    cam1_t *cam=&G.ed.cam1;
    cam1_use(cam);
    drawclever(walltex,flortex);
    glDisable(GL_DEPTH_TEST);

}	

static void editor_draw2d()
{
    rnd_tg(0);
    rnd_clearall();

    cam1_update(&G.ed.cam1);

    cam2_t *c=&G.ed.cam;
    cam2_use(c->x,c->y,c->zoom);
    rnd_simplegrid(c,G.ed.level);


    // walls while editing
    int sz=arr_len(G.nkstk);
    glColor3f(1,1,0);
    glBegin(GL_LINE_STRIP);
    for (int i=0;i<sz;i++) {
	glVertex2f(G.nkstk[i].x, G.nkstk[i].y);
    }
    glVertex2f(G.ed.mousex,G.ed.mousey);
    glEnd();


    // cursor
    glPointSize(5);
    glBegin(GL_POINTS);
    glVertex2f(G.ed.mousex,G.ed.mousey);

    // nooks
    glColor3f(1,0,0);
    for (int i=0;i<sz;i++) {
	glVertex2f(G.nkstk[i].x, G.nkstk[i].y);
    }

    glEnd();
    glPointSize(1);

    // sectors
    short undermouse[1000];
    int loc=chkloc(G.ed.mousex,G.ed.mousey,undermouse,1000);

    for (int i=0;i<G.nsects;i++) {
	sect_t *s=G.sects + i;

	// triangulation
	glColor3f(0,1,1);
	glBegin(GL_LINES);
	float *t=s->tess;
	for (int j=0;j<s->nt;j++) {
	    glVertex2f(t[0],t[1]);
	    glVertex2f(t[2],t[3]);
	    glVertex2f(t[2],t[3]);
	    glVertex2f(t[4],t[5]);
	    glVertex2f(t[4],t[5]);
	    glVertex2f(t[0],t[1]);
	    t+=6;
	}
	glEnd();

	// walls
	int w=s->firstw;
	glColor3f(.8,.8,.8);
	glBegin(GL_LINE_LOOP);
	for (int j=0;j<s->nw;j++) {
	    glVertex2f(G.walls[w].x,G.walls[w].y);
	    if (w+1 != G.walls[w].nextco) {
		glEnd();
		glBegin(GL_LINE_LOOP);
	    }
	    ++w;
	}
	glEnd();

	// shared edges
	w=s->firstw;
	glColor3f(1,0,0);
	glBegin(GL_LINES);
	for (int j=0;j<s->nw;j++) {
	    if (G.walls[w].adjw!=-1) {
		glVertex2f(G.walls[w].x,G.walls[w].y);
		glVertex2f(G.walls[G.walls[w].nextco].x, G.walls[G.walls[w].nextco].y);
	    }
	    ++w;
	}
	glEnd();
    }

    // show current sector
    if (loc==LOC_IN) {
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//	glColor4f(.5,.5,.5,.5);
	glBegin(GL_TRIANGLES);
	sect_t *s=G.sects + undermouse[0];
	float *t=s->tess;
	for (int j=0; j<s->nt; j++) {
	    glColor4f(1,0,0,.5);
	    glVertex2fv(t);
	    t+=2;
	    glColor4f(0,1,0,.5);
	    glVertex2fv(t);
	    t+=2;
	    glColor4f(0,0,1,.5);
	    glVertex2fv(t);
	    t+=2;
	}
	glEnd();
	glDisable(GL_BLEND);
    }

    // show pvs
    if (1) {
	// floor PVS
//	glPointSize(5);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glBegin(GL_TRIANGLES);
//	glBegin(GL_POINTS);
	int sz=arr_len(G.ed.pvs.flors);
	float r=1,g=0,b=0;
	for (int i=0;i<sz;i++) {
//	    glColor4f(r,g,b,.5);
	    glColor4f(1,1,1,.5);
	    float tmp=b;
	    b=g;
	    g=r;
	    r=tmp;

	    pvs2d_flor_t *f=G.ed.pvs.flors+i;
	    sect_t *s=G.sects + f->sid;
	    float verts[10];
	    for (int j=0;j<f->nv;j++) {
		if (f->eid[j]>=0) {
		    float *v1=s->tess + f->eid[j]*2;
		    float *v2=(f->eid[j]%3==2) ? v1-4 : v1+2;
		    verts[2*j]=v1[0]+(v2[0]-v1[0])*f->t[j];
		    verts[2*j+1]=v1[1]+(v2[1]-v1[1])*f->t[j];
		}
		else {
		    verts[2*j]=f->t[j];
		    verts[2*j+1]=f->y;
		}
//		glVertex2f(verts[2*j],verts[2*j+1]);
	    }

//	    printf("%iverts\n", f->nv);
#if 1
	    for (int j=1;j<f->nv-1;j++) {
		glVertex2f(verts[0],verts[1]);
		glVertex2f(verts[j*2],verts[j*2+1]);
		glVertex2f(verts[j*2+2],verts[j*2+3]);
	    }
#endif
	}
	glPointSize(1);
	glEnd();
	glDisable(GL_BLEND);

	// wall PVS
	pvs2d_wall_t *pvsw=G.ed.pvs.walls;
	if (pvsw) {
	    glColor3f(1,0,1);
	    glBegin(GL_LINES);
	    for (int i=0;i<G.ed.pvs.wall_ns;i++) {
		sect_t *s=&G.sects[pvsw->sid];
		for (int j=0;j<G.ed.pvs.nwalls[i];j++) {
		    int id1=pvsw->eid;
		    int id2=pvsw->eid%3==2 ? id1-2 : id1+1;
		    float *ev1=&s->tess[id1*2];
		    float *ev2=&s->tess[id2*2];

		    float dx=ev2[0]-ev1[0];
		    float dy=ev2[1]-ev1[1];

		    float v1[2]={ev1[0]+pvsw->t1*dx, ev1[1]+pvsw->t1*dy};
		    float v2[2]={ev1[0]+pvsw->t2*dx, ev1[1]+pvsw->t2*dy};
		    glVertex2fv(v1);
		    glVertex2fv(v2);

		    ++pvsw;
		}
	    }
	    glEnd();
	}
	
    }
    
    // cam1
    // pos
    float x=G.ed.cam1.pos.x;
    float y=G.ed.cam1.pos.y;

    glPointSize(5);
    glBegin(GL_POINTS);
    glVertex2f(x,y);
    glEnd();
    glPointSize(1);

    // dir
    float dirx=G.ed.cam1.rotaxis[1].x*20;
    float diry=G.ed.cam1.rotaxis[1].y*20;
    glBegin(GL_LINES);
    glVertex2f(x,y);
    glVertex2f(x+dirx,y+diry);
    glEnd();

}


void editor_draw()
{
    if (G.ed.dim==3) editor_draw3d();
    else if (G.ed.dim==2) editor_draw2d();
}

//////////////////////////////////////////// console ///////////////////////////////////
static void cs_draw()
{
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    vfs_use(-1);
    if(G.cs_winh<G.screenh/2){
	G.cs_winh += G.screenh/12;
    }

    // console background
    float bnd[4]={0, G.screenh-G.cs_winh, G.screenw, G.screenh};
    rnd_quad(G.cs_bg, bnd);

    // minibuffer background
    int h=fnt_getcur()->height;
    bnd[3]=bnd[1];
    bnd[1]-=h;
    rnd_quad(G.cs_mbbg, bnd);

    // console text
    fnt_color4fv(G.cs_fg);
    int sz=0;
    int nol=G.cs_nol;
    char **lines=G.cs_lines;
    int asc=fnt_getcur()->ascender;
    int y= 1.5*G.screenh-G.cs_winh;//-asc;

    for(int i=0; i<nol && i<G.cs_outh && i<G.cs_h; i++){
	fnt_printf(0, y-=h, 0, "%s", lines[(G.cs_out+i)%G.cs_h]); // DEBUGGED: (G.cs_out+i) =>  (G.cs_out+i)%G.cs_h
    }
	
    // minibuffer text
    fnt_printf(0, G.screenh-G.cs_winh-asc, 0, "%s", G.cs_mbcur->mb);
    fnt_flush();

    // blinking cursor:
    --G.cs_csrblink;
    if(G.cs_csrblink>0){
	vfs_use(-1);

	glEnable(GL_BLEND);
	glBlendFunc(GL_ONE_MINUS_DST_COLOR, GL_ZERO);  // inverse the color
	float rgba[4]={1,1,1,1};
		
	float bgn=0;
	char *mb=G.cs_mbcur->mb;
	fnt_t *fnt=fnt_getcur();

	// TODO: font scale factor
	int csr=G.cs_mbpos;			// cursor position
	for(int i=0; i<csr; i++){
	    bgn += fnt->ch[mb[i]].advance;
	}
		
	int cw = (!csr || !mb[csr]) ? fnt->avgw : fnt->ch[mb[csr]].advance;			// cursor width
	bnd[0]=bgn;
	bnd[2]=bgn+cw; //fnt->ch[*c].advance;
		
	rnd_quad(rgba, bnd);
	glDisable(GL_BLEND);
		
    }

    glDisable(GL_BLEND);

}


#define VU_CKSUM 0x58cffe26

void vu_use(const char *nam)
{
    rndvu_t *v = cvar_getp__(nam, VU_CKSUM);
    if(!v) return;

    // TESTING
    G.vu_cur=v;

    G.mode=MODE_VU;

    // ~TESTING
}


/* if the view doesn't exist, create a new view */
void vu_force_bgn(char *nam)
{ 
    if(G.vu_edit) return;

    rndvu_t *v=cvar_getp__(nam, VU_CKSUM);
    if(!v){
	v=malloc(sizeof(rndvu_t));

	memset(v->kcmd, 0, sizeof(v->kcmd));

	v->cmd=NULL;
	v->nam=strdup(nam);
	v->btnlist=NULL;
		
	v->g=4;
	v->x=v->y=0;
	v->w=G.screenw;
	v->h=G.screenh;

	v->next=G.vulist;
	if(G.vulist){
	    v->prev=G.vulist->prev;
	    G.vulist->prev=v;
	}
	else{
	    v->prev=v;
	}
	v->prev->next=v;
	G.vulist=v;


	// TESTING
	for (int i=0;i<4;i++) {
	    v->fg[i]=0;
	    v->bg[i]=0;
	}
	// ~TESTING

	cvar_setp__(nam, v, VU_CKSUM);

	G.guied_vuall[4]=v;
    }
    G.vu_edit = v;
}

void vu_end()
{
    G.vu_edit=NULL;
}

void vu_bgn(char *nam)
{
    if(G.vu_edit) return;

    rndvu_t *v=cvar_getp__(nam, VU_CKSUM);
    if(!v) return;
    G.vu_edit = v;
}


void vu_geterr()
{
    return !!G.vu_edit;
}

//#define vu_end() (G.vu_edit=NULL)
#define vu_del() (G.vu_edit ? cvar_del(G.vu_edit->nam) : 0)

static void vu_grid(int g)
{
    rndvu_t *v=G.vu_edit;
    if (v){
	if (g<0) g=0;
	else if (g>8) g=8;

	if (v->g!=g) {
	    rndvu_t *w=v->next;
	    while (w!=v) {
		if (w->g==v->g) {
		    G.guied_vuall[v->g]=w;
		    break;
		}
		w=w->next; // DEBUGGED: don't forget "next" -- using for() instead of while() is good for avoiding infinite loops
	    }
			
	    if (w==v) G.guied_vuall[v->g]=NULL;

	    rndvu_t *v2=G.guied_vuall[g];
	    if (v2) {
		v->prev->next=v->next;
		v->next->prev=v->prev;

		// put v in front of v2
		v2->prev->next=v;
		v->prev=v2->prev;
		v2->prev=v;
		v->next=v2;
	    }
	    G.guied_vuall[g]=v;
	    v->g=g;
	}
		
    }
}

static void vu_nam(char *newnam)
{
    rndvu_t *v=G.vu_edit;
    if(!v) return;

    if(v->nam) free(v->nam);
    v->nam = strdup(newnam);

    // ?????
}

void vu_kcmd(int k, char *cmd)
{
    if (!cmd) return;

    rndvu_t *v=G.vu_edit;
    if(!v) return;

    if(v->kcmd[k]) free(v->kcmd[k]);
    v->kcmd[k] = strdup(cmd);
}

void vu_cmd(char *cmd)
{
    rndvu_t *v=G.vu_edit;
    if(!v) return;

    if(v->cmd) free(v->cmd);
    v->cmd = strdup(cmd);
}


/*
  If the button doesn't exist, create a new button

*/
void btn_force_bgn(char *nam)
{
    if(!G.vu_edit) return;

    guibtn_t *b=NULL;
    char buf[32]="PushButton";
    if (nam && *nam) {
	for(b=G.vu_edit->btnlist; b; b=b->next){
	    if(!strcmp(b->nam, nam)) break;
	}
    }
    else {  // auto new name
	sprintf(buf+10, "%i", G.vu_edit->btnid++);
	nam=buf;
    }

    if(!b){ // new button
	b = malloc(sizeof(guibtn_t));

	b->nam=strdup(nam);
	b->txt=strdup("Untitiled\n");
	b->cmd=NULL;

	b->w=50;
	b->h=50;
	b->x=(G.screenw-b->w)/2;
	b->y=(G.screenh-b->h)/2;

	for(int i=0;i<4;i++){
	    b->fg[i]=1;
	    b->bg[i]=.5;
	}

	b->txtnol=1;
	b->txtscale=1;

	b->prev=NULL;
	b->next=G.vu_edit->btnlist;
	if(b->next){
	    b->next->prev=b;
	}
	G.vu_edit->btnlist=b;

	cs_printf("Button \"%s\" is created\n", nam);
    }

    G.btn_edit = b;
}


/*
  If the button doesn't exist, return
*/
static void btn_bgn(char *nam)
{
    if(!G.vu_edit) return;

    guibtn_t *b;
    for(b=G.vu_edit->btnlist; b; b=b->next){
	if(!strcmp(b->nam, nam)) break;
    }

    if(!b) return;

    G.btn_edit = b;
	
}

#define btn_end() (G.btn_edit=NULL)


void btn_del()
{
    rndvu_t *v=G.vu_edit;
    guibtn_t *b=G.btn_edit;
    if (!v || !b) return;

    if (b->prev) b->prev->next=b->next;
    else v->btnlist=b->next;

    if (b->next) b->next->prev=b->prev;

    free(b);
}


static void btn_cmd(char *cmd)
{
    guibtn_t *b=G.btn_edit;
    if(!b) return;

    if(b->cmd) free(b->cmd);
    b->cmd = strdup(cmd);
}

static void btn_txt(char *txt)
{
    guibtn_t *b=G.btn_edit;
    if(!b) return;

    char *c=txt;
    if(*c){
	b->txtnol=1;
	do{
	    if(*c=='\n'){
		if(c[1]) ++b->txtnol;
	    }
	    ++c;
	}while(*c);
    }
    else{
	b->txtnol=0;
    }

    if(b->txt) free(b->txt);
    b->txt = strdup(txt);
}

static void btn_nam(char *nam)
{
    guibtn_t *b=G.btn_edit;
    if(!b) return;

    if(b->nam) free(b->nam);
    b->nam = strdup(nam);
}

static void btn_fgcolor4f(float r,float g,float b,float a)
{
    guibtn_t *btn=G.btn_edit;
    if(!btn) return;

    btn->fg[0]=r;
    btn->fg[1]=g;
    btn->fg[2]=b;
    btn->fg[3]=a;
}

static void btn_bgcolor4f(float r,float g,float b,float a)
{
    guibtn_t *btn=G.btn_edit;
    if(!btn) return;

    btn->bg[0]=r;
    btn->bg[1]=g;
    btn->bg[2]=b;
    btn->bg[3]=a;
}

static int btn_intx(guibtn_t *b1, guibtn_t *b2)
{
    return b1->x < b2->x+b2->w  &&  b1->x+b1->w > b2->x  && b1->y < b2->y+b2->h  && b1->y+b1->h > b2->y;
}


static void btn_rankup()
{
    guibtn_t *btn=G.guied_btn;
    for (guibtn_t *b=btn->prev; b; b=b->prev) {
	if (btn_intx(btn, b)) {
	    if (btn->next) btn->next->prev = btn->prev;
	    if (btn->prev) btn->prev->next = btn->next;

	    btn->next = b;
	    btn->prev = b->prev;
	    if (b->prev) b->prev->next = btn;
	    else G.guied_vu->btnlist = btn;
	    b->prev = btn;
	    return;
	}
    }
}


static void btn_rankdown()
{
    guibtn_t *btn=G.guied_btn;
    for (guibtn_t *b=btn->next; b; b=b->next) {
	if (btn_intx(btn, b)) {
	    if (btn->next) btn->next->prev = btn->prev;
	    if (btn->prev) btn->prev->next = btn->next;
	    else G.guied_vu->btnlist=btn->next;

	    btn->prev = b;
	    btn->next = b->next;
	    if (b->next) b->next->prev = btn;
	    b->next = btn;
	    return;
	}
    }
}


/*
  the cmd will be executed everytime the key is hit
*/
static void vu_keycmd(char *v, int k, char *cmd)
{
    if(!v || !cmd) return;

    rndvu_t *vu = cvar_getp__(v, VU_CKSUM);
    if(!vu) return;
	
    if(vu->kcmd[k]) free(vu->kcmd[k]);
    vu->kcmd[k] = strdup(cmd);
}


/*
  the cmd will be executed in each frame.

  multiple commands can be separated by ';'
*/
static void vu_framecmd(char *v, char *cmd)
{
    rndvu_t *vu = cvar_getp__(v, VU_CKSUM);
    if(!vu) return;

    if(vu->cmd) free(vu->cmd);
    vu->cmd = strdup(cmd);
}


#define RESIZE_EDGE_WIDTH 6
#define BTN_MIN_WIDTH 5
#define BTN_MIN_HEIGHT 5
#define VU_MIN_WIDTH 20
#define VU_MIN_HEIGHT 20


//void btn_printf(int posx, int posy, int alignment, int bnd[4], const char *fmt, ...)

static void btn_prepvb(guibtn_t *btn, int border, float z)
{
    char *string=btn->txt;

    float scalefactor = 1; //btn->txtscale;
    char *line=string;
    int i=0;
    char *bgn = string;
    btn_vbfmt_t f[4];

    /*
      quad:

      3   2
	  
      0   1
    */

    for(int i=0; i<4; i++){
	memcpy(f[i].rgba, btn->fg, 16);
	f[i].z = z;
    }

    fnt_t *fnt=G.fnt_cur;
    float btncenter[2]={btn->x+btn->w*.5, btn->y+btn->h*.5};
	
    /*
      totalheight = nol * fntheight * scalefactor

      => top = btn->center.y + totalheight/2

      => baselinetop = top - ascender

      => baselinetop = btn->center.y + totalheight/2 - ascender


    */

    float totalh = btn->txtnol * fnt->height * scalefactor;
    float baseline = btncenter[1] + totalh/2 - fnt->ascender;
    float totalw=0;

    float left =btn->x;
    float right =btn->x+btn->w;
    float bottom =btn->y;
    float top =btn->y+btn->h;

    float y=baseline;
	
    // render each line of text
    while(*line){
	// compute line width
	char *c=line;
	while(*c && *c!='\n'){
	    totalw += fnt->ch[*c].advance;
	    ++c;
	}
	char *end=c;

	float ofs= -totalw * 0.5f * scalefactor;
	float x= btncenter[0] + ofs;
	fntch_t *ch;
	for (c=line; *c && *c!='\n'; ++c , x+=ch->advance*scalefactor) {
	    ch=&fnt->ch[*c];
	    float w=scalefactor * ch->w;
	    float h=scalefactor * ch->h;

	    float xmin = x + ch->x;
	    float xmax = xmin + w;

	    float ymin = y + ch->y;
	    float ymax = ymin + h;

	    if (xmax<left) continue;

	    if (xmin>right || ymin>top) {

		// go to next line
		if (!*end) {
		    goto done;
		}

		line=end+1;
		break;
	    }

	    if (ymax<bottom) {
		// done!
		goto done;
	    }

	    // clipped by the button
	    float texw = ch->texw;
	    float texh = ch->texh;

	    if (xmin<left) {
		f[0].u = f[3].u = ch->u + (left-xmin)/w*texw;
		f[0].x = f[3].x = left;
	    }
	    else {
		f[0].u = f[3].u = ch->u;
		f[0].x = f[3].x = xmin;
	    }

	    if (xmax>right) {
		f[1].u = f[2].u = ch->u+texw - (xmax-right)/w*texw;
		f[1].x = f[2].x = right;
	    }
	    else {
		f[1].u = f[2].u = ch->u+texw;
		f[1].x = f[2].x = xmax;
	    }

	    if (ymin<bottom) {
		f[0].v = f[1].v = ch->v + (bottom-ymin)/h*texh;
		f[0].y = f[1].y = bottom;
	    }
	    else {
		f[0].v = f[1].v = ch->v;
		f[0].y = f[1].y = ymin;				
	    }

	    if (ymax>top) {
		f[2].v = f[3].v = ch->v+texh - (ymax-top)/h*texh;
		f[2].y = f[3].y = top;
	    }
	    else {
		f[2].v = f[3].v = ch->v+texh;
		f[2].y = f[3].y = ymax;
	    }

	    arr_pushn(G.btn_vb, f, 4);
	}


	// DEBUGGED: arr_pushn(G.btn_vb, f, 4) => in the for loop above

	if(*end=='\n') ++end;
	line = end;
	y-=fnt->height*scalefactor;
    }


  done:
    ;

    // button polygon
    btn_vbfmt_t *b=f;
    for(int i=0; i<4; i++){
	memcpy(b[i].rgba, btn->bg, 16);
	b[i].z = z;
	b[i].u = b[i].v = 1;
    }

    b[0].x = b[3].x = btn->x;
    b[1].x = b[2].x = btn->x + btn->w;
    b[0].y = b[1].y = btn->y;
    b[2].y = b[3].y = btn->y + btn->h;
    arr_pushn(G.btn_vb, b, 4);


    // TODO: make this look prettier!
    if (border) {
	float xmin,xmax,ymin,ymax;
	xmin = b[0].x - RESIZE_EDGE_WIDTH;
	b[3].x = b[0].x = xmin+2;

	xmax = b[1].x + RESIZE_EDGE_WIDTH;
	b[2].x = b[1].x = xmax-2;

	ymin = b[0].y - RESIZE_EDGE_WIDTH;
	b[1].y = b[0].y = ymin+2;

	ymax = b[2].y += RESIZE_EDGE_WIDTH;
	b[3].y = b[2].y = ymax-2;

	for (int i=0; i<4; i++) {
	    b[i].z -= z*.5;
	    for (int j=0; j<4; j++){
		b[i].rgba[j] *=.5;
	    }
	}

	arr_pushn(G.btn_vb, b, 4);


	for (int i=0; i<4; i++){
	    b[i].z += z*.5;
	    for (int j=0; j<4; j++){
		b[i].rgba[j] =1;
	    }
	}

	// four corners (enlarge the corners a little)
	b[3].x = b[0].x = xmin;
	b[2].x = b[1].x = xmax;
	b[1].y = b[0].y = ymin;
	b[3].y = b[2].y = ymax;

	// 0
	b[1].x=b[2].x=btn->x;
	b[2].y=b[3].y=btn->y;
	arr_pushn(G.btn_vb, b, 4);

	// 1
	b[0].x=b[3].x=btn->x+btn->w;
	b[1].x=b[2].x=xmax;
	arr_pushn(G.btn_vb, b, 4);

	// 2
	b[0].y=b[1].y=btn->y+btn->h;
	b[2].y=b[3].y=ymax;
	arr_pushn(G.btn_vb, b, 4);

	// 3
	b[0].x=b[3].x=xmin;
	b[1].x=b[2].x=btn->x;
	arr_pushn(G.btn_vb, b, 4);

    }
}


enum {
    MT_VU_TO_GUIED_NORMAL,
    MT_GUIED1_TO_VU,
    MT_GUIED1_TO_GUIED9,
    MT_GUIED9_TO_VU,
    MT_GUIED9_TO_GUIED1,
    MT_VU_TO_VU,
    MT_CONSOLE_ON,
};



/*
  This function is called after the frame cmd is executed.

  Each button is just a box.
  Box text shouldn't be rendered using fnt_printf() & fnt_flush() => cannot handle transparent btns
  Characters should be clipped by the btn.
  
*/
static void vu_drawgui(rndvu_t *v)
{
    if(!v) return;
	
    /*
      font uses texture, but buttons don't (or uses different textures).

      fnt_color=(R,G,B,A)*tex
      btn_color=(R,G,B,A)*1  --> let btn's tex=1

      depthtest = LESS
      alphatest = enable
      alphafunc = !0
      blending = enable

      drawfnt
      drawbtn

      blend(fnt_color, bg)
      blend(btn_color, bg)
      dont_blend(fnt_color,btn_color)

    */
    rnd_scrncoord();

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    glEnable(GL_ALPHA_TEST);
    glAlphaFunc(GL_NOTEQUAL, 0);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    arr_popall(G.btn_vb);
    float z=.9;
    for(guibtn_t *b=v->btnlist; b; b=b->next){
	float dz=1e-3f;
	int border= G.mode==MODE_GUIED && b==G.guied_btn;
	btn_prepvb(b,border,z-=dz);
    }

    vfs_use(G.btn_shader);
	
    if (G.mt==MT_GUIED9_TO_GUIED1) {
	printf("hello");

    }

    samp2d("fnt_tex", G.fnt_tex);  // DEBUGGED: "font_tex" => "fnt_tex"

    int n=arr_len(G.btn_vb);
    int sz=sizeof(btn_vbfmt_t);
    int vbo=G.fnt_vbo;			// share VBO with font
    int cap=vb_getcap(vbo)/sz;

    // share the same VBO with fnt rendering
    btn_vbfmt_t *vb=G.btn_vb;
    while(n>0){
	if(n/cap){
	    vb_update(vbo,cap,sz,vb);
	    vb+=cap;
	}
	else{
	    vb_update(vbo,n,sz,vb);
	}

	n-=cap;
	rendq("{co texco color}", -1, vbo);
    }

    arr_popall(G.btn_vb);

    glDisable(GL_BLEND);
    glDisable(GL_ALPHA_TEST);
    glDisable(GL_DEPTH_TEST);

}


static void vu_mmfunc(SDL_MouseMotionEvent motion)
{
	
}


static void vu_mbfunc(SDL_MouseButtonEvent button, int down)
{
    rndvu_t *v=G.vu_cur;
    if (G.mode!=MODE_VU || !v) return;

    // Test if the mouse cursor is in any button

    int x=button.x;
    int y=button.y;

    if (button.button==SDL_BUTTON_LEFT) {
	for (guibtn_t *b=v->btnlist; b; b=b->next) {
	    if (x>b->x && x<b->x+b->w && y>b->y && y<b->y+b->h) {
		char *c=b->cmd;
		if (c&&*c) {
		    cmd_execnow(c);
		    return;
		}
		else break; // if the button has no cmd, ignore it.
	    }
	}
    }

    // If the button has no command, the mouse button event will be ignored by the GUI system

    char *c=v->kcmd[button.button];  // 1<=button.button<=5,   and SDLK_* doesn't use include this range. So it's ok to
    // to treat mouse button as ordinary key buttons.
    if (c) {
	char cmd[strlen(c)+2];
	cmd[0]=down?'+':'-';
	strcpy(cmd+1, c);
	cmd_execnow(cmd);
    }
}


static void vu_keyfunc(SDL_keysym keysym, int down)
{
    if (G.mode!=MODE_VU || !G.vu_cur) return;

    int k=keysym.sym;
    char *cmd = G.vu_cur->kcmd[k];
    if(cmd){
	char *buf = alloca(strlen(cmd) + 2);
	buf[0] = down ? '+' : '-';
	strcpy(buf+1, cmd);
	cmd_addtxt(buf);
    }
	
}


/*

  Transitions:

  view, guied_normal, guied9, console


  console
  __________________|____________
  |            |                 |
  \/           \/                \/
  view <-> guied_normal <-> guied9
  /\                              |
  |______________________________|
  

  console can be opened/closed during any transition



  view<->view -- anim
  view->guied_normal -- instant (in fact, nothing has changed except the status bar)
  guied_normal->view -- instant/anim (the view being editted and the view before entering the editting mode can be different)

  guied9->view -- anim
  guied_normal->guied9 -- anim
  guied9->guied_normal -- anim

*/




enum{BTN_LEFT=1, BTN_RIGHT=2, BTN_BOTTOM=4, BTN_TOP=8, BTN_MOVE=15,
     VU_LEFT=16, VU_RIGHT=32, VU_BOTTOM=64, VU_TOP=128, VU_MOVE=15<<4,
     VU_CHANGE_GRID=256};

char *guiop[257]={[1]="BTN_LEFT", [2]="BTN_RIGHT", [4]="BTN_BOTTOM", [8]="BTN_TOP", [15]="BTN_MOVE", 
		  [16]="VU_LEFT", [32]="VU_RIGHT", [64]="VU_BOTTOM", [128]="VU_TOP", [127]="VU_MOVE", [240]="VU_CHANGE_GRID"};

static void guied_mmfunc(SDL_MouseMotionEvent motion);

static void guied_refresh(int x, int y)
{
    SDL_MouseMotionEvent motion;
    motion.type = SDL_MOUSEMOTION;
    motion.state = 0;
    motion.x = x;
    motion.y = y;
    motion.xrel = motion.yrel = 0;
    guied_mmfunc(motion);
}


void vu_loadtn(rndvu_t *v)
{
    assert(v);
    assert(v->nam);

    FILE *f=fopen(v->nam, "r");

    if (!f) return;

    int w,h;
    fread(&w, sizeof(int), 1, f);
    fread(&h, sizeof(int), 1, f);

    int sz = G.vu_tnsz;
    assert(w==sz && h==sz);

    int nbytes=8*(w/4)*(h/4);
    uchar *img=alloca(nbytes+16);
    img=(uchar*)(((uint)img & ~15) +16);  // should be 16-byte aligned

    fread(img, 1, nbytes, f);
    fclose(f);

    int x,y;
    x=v->g%3 * sz;
    y=v->g/3 * sz;
	
    glBindTexture(GL_TEXTURE_2D, G.guied_tntex);
    glCompressedTexSubImage2D(GL_TEXTURE_2D, 0, x, y, w, h, GL_COMPRESSED_RGB_S3TC_DXT1_EXT, nbytes, img);

}


void vu_cleartn(int g)
{
    int sz = G.vu_tnsz;
    int w=sz;
    int h=sz;
    int nbytes=8*(w/4)*(h/4);
    uchar *img=alloca(nbytes);
	
    memset(img,0,nbytes);

    int x,y;
    x=g%3 * sz;
    y=g/3 * sz;
	
    glBindTexture(GL_TEXTURE_2D, G.guied_tntex);
    glCompressedTexSubImage2D(GL_TEXTURE_2D, 0, x, y, w, h, GL_COMPRESSED_RGB_S3TC_DXT1_EXT, nbytes, img);
}



void rnd_fb2tex(uint tex)
{
    glBindTexture(GL_TEXTURE_2D, tex);
    glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, G.screenw, G.screenh);
}

/*
  Shouldn't save console window in any case -> called before rendering the console
*/
void vu_savetn(rndvu_t *v)
{
    assert(v);
    assert(v->nam);

    vfs_use(-1); // DEBUGGED: have to disable shaders

    // render to a temp texture
    rnd_fb2tex(G.tex_screensz);

    // render a small quad to some RT
    int sz=G.vu_tnsz;

    rnd_tg("td", G.tex_128x128, G.dtex_128x128);
    rnd_clearall();

    float tcbnd[4]={0,0,1,1};
    float qbnd[4]={0,0,sz,sz};
    rnd_texquadsz(G.tex_screensz, sz,sz,tcbnd, qbnd); //TODO

    // read back pixels
    int nbytes;
    nbytes=4*sz*sz;
    uchar *rgba=alloca(16 + 8 + nbytes); // padding, w,h,pixels[]
    rgba=(uchar*)(((uint)rgba & ~15) +16);  // should be 16-byte aligned
    glReadPixels(0,0,sz,sz,GL_RGBA,GL_UNSIGNED_BYTE,rgba);

    // compress to dxt1
    nbytes=8*(sz/4)*(sz/4);
    uchar *dxt=alloca(16 + 8 + nbytes); // padding, w,h,pixels[]
    dxt=(uchar*)(((uint)dxt & ~15) +16);  // should be 16-byte aligned
    img_todxt1(rgba, dxt, sz, sz, &nbytes);

    // store on disk
    FILE *f=fopen(v->nam, "w");
    assert(f);

    fwrite(&sz,sizeof(int),1,f); // DEBUGGED: img_todxt1() doesn't wrihte the size of the picture at the beginning
    fwrite(&sz,sizeof(int),1,f); // write it yourself
    fwrite(dxt,nbytes,1,f);
    fclose(f);

    rnd_tg(0);
	
}


typedef struct {
    char *str;
    int k;
}k__t;

static int cmpk(k__t *k1, k__t *k2)
{
    return strcmp(k1->str,k2->str);
}


static k__t keys[]={
    " ", 32,
    "0", 48,
    "1", 49,
    "2", 50,
    "3", 51,
    "4", 52,
    "5", 53,
    "6", 54,
    "7", 55,
    "8", 56,
    "9", 57,
    "A", 97,
    "AMPERSAND", 38,
    "ASTERISK", 42,
    "AT", 64,
    "B", 98,
    "BACKQUOTE", 96,
    "BACKSLASH", 92,
    "BACKSPACE", 8,
    "BREAK", 318,
    "C", 99,
    "CAPSLOCK", 301,
    "CARET", 94,
    "CLEAR", 12,
    "COLON", 58,
    "COMMA", 44,
    "COMPOSE", 314,
    "D", 100,
    "DELETE", 127,
    "DOLLAR", 36,
    "DOWN", 274,
    "E", 101,
    "END", 279,
    "EQUALS", 61,
    "ESCAPE", 27,
    "EURO", 321,
    "EXCLAIM", 33,
    "F", 102,
    "F1", 282,
    "F10", 291,
    "F11", 292,
    "F12", 293,
    "F13", 294,
    "F14", 295,
    "F15", 296,
    "F2", 283,
    "F3", 284,
    "F4", 285,
    "F5", 286,
    "F6", 287,
    "F7", 288,
    "F8", 289,
    "F9", 290,
    "FIRST", 0,
    "G", 103,
    "GREATER", 62,
    "H", 104,
    "HASH", 35,
    "HELP", 315,
    "HOME", 278,
    "I", 105,
    "INSERT", 277,
    "J", 106,
    "K", 107,
    "KP0", 256,
    "KP1", 257,
    "KP2", 258,
    "KP3", 259,
    "KP4", 260,
    "KP5", 261,
    "KP6", 262,
    "KP7", 263,
    "KP8", 264,
    "KP9", 265,
    "KP_DIVIDE", 267,
    "KP_ENTER", 271,
    "KP_EQUALS", 272,
    "KP_MINUS", 269,
    "KP_MULTIPLY", 268,
    "KP_PERIOD", 266,
    "KP_PLUS", 270,
    "L", 108,
    "LALT", 308,
    "LCTRL", 306,
    "LEFT", 276,
    "LEFTBRACKET", 91,
    "LEFTPAREN", 40,
    "LESS", 60,
    "LMETA", 310,
    "LSHIFT", 304,
    "LSUPER", 311,
    "M", 109,
    "MENU", 319,
    "MINUS", 45,
    "MODE", 313,
    "MOUSE1", 1,
    "MOUSE2", 2,
    "MOUSE3", 3,
    "MOUSEWHEELUP", 4,
    "MOUSEWHEELDOWN", 5,
    "N", 110,
    "NUMLOCK", 300,
    "O", 111,
    "P", 112,
    "PAGEDOWN", 281,
    "PAGEUP", 280,
    "PAUSE", 19,
    "PERIOD", 46,
    "PLUS", 43,
    "POWER", 320,
    "PRINT", 316,
    "Q", 113,
    "QUESTION", 63,
    "QUOTE", 39,
    "QUOTEDBL", 34,
    "R", 114,
    "RALT", 307,
    "RCTRL", 305,
    "RETURN", 13,
    "RIGHT", 275,
    "RIGHTBRACKET", 93,
    "RIGHTPAREN", 41,
    "RMETA", 309,
    "RSHIFT", 303,
    "RSUPER", 312,
    "S", 115,
    "SCROLLOCK", 302,
    "SEMICOLON", 59,
    "SLASH", 47,
    "SPACE", 32,
    "SYSREQ", 317,
    "T", 116,
    "TAB", 9,
    "U", 117,
    "U", 121,
    "UNDERSCORE", 95,
    "UNDO", 322,
    "UNKNOWN", 0,
    "UP", 273,
    "V", 118,
    "W", 119,
    "WORLD_0", 160,
    "WORLD_1", 161,
    "WORLD_10", 170,
    "WORLD_11", 171,
    "WORLD_12", 172,
    "WORLD_13", 173,
    "WORLD_14", 174,
    "WORLD_15", 175,
    "WORLD_16", 176,
    "WORLD_17", 177,
    "WORLD_18", 178,
    "WORLD_19", 179,
    "WORLD_2", 162,
    "WORLD_20", 180,
    "WORLD_21", 181,
    "WORLD_22", 182,
    "WORLD_23", 183,
    "WORLD_24", 184,
    "WORLD_25", 185,
    "WORLD_26", 186,
    "WORLD_27", 187,
    "WORLD_28", 188,
    "WORLD_29", 189,
    "WORLD_3", 163,
    "WORLD_30", 190,
    "WORLD_31", 191,
    "WORLD_32", 192,
    "WORLD_33", 193,
    "WORLD_34", 194,
    "WORLD_35", 195,
    "WORLD_36", 196,
    "WORLD_37", 197,
    "WORLD_38", 198,
    "WORLD_39", 199,
    "WORLD_4", 164,
    "WORLD_40", 200,
    "WORLD_41", 201,
    "WORLD_42", 202,
    "WORLD_43", 203,
    "WORLD_44", 204,
    "WORLD_45", 205,
    "WORLD_46", 206,
    "WORLD_47", 207,
    "WORLD_48", 208,
    "WORLD_49", 209,
    "WORLD_5", 165,
    "WORLD_50", 210,
    "WORLD_51", 211,
    "WORLD_52", 212,
    "WORLD_53", 213,
    "WORLD_54", 214,
    "WORLD_55", 215,
    "WORLD_56", 216,
    "WORLD_57", 217,
    "WORLD_58", 218,
    "WORLD_59", 219,
    "WORLD_6", 166,
    "WORLD_60", 220,
    "WORLD_61", 221,
    "WORLD_62", 222,
    "WORLD_63", 223,
    "WORLD_64", 224,
    "WORLD_65", 225,
    "WORLD_66", 226,
    "WORLD_67", 227,
    "WORLD_68", 228,
    "WORLD_69", 229,
    "WORLD_7", 167,
    "WORLD_70", 230,
    "WORLD_71", 231,
    "WORLD_72", 232,
    "WORLD_73", 233,
    "WORLD_74", 234,
    "WORLD_75", 235,
    "WORLD_76", 236,
    "WORLD_77", 237,
    "WORLD_78", 238,
    "WORLD_79", 239,
    "WORLD_8", 168,
    "WORLD_80", 240,
    "WORLD_81", 241,
    "WORLD_82", 242,
    "WORLD_83", 243,
    "WORLD_84", 244,
    "WORLD_85", 245,
    "WORLD_86", 246,
    "WORLD_87", 247,
    "WORLD_88", 248,
    "WORLD_89", 249,
    "WORLD_9", 169,
    "WORLD_90", 250,
    "WORLD_91", 251,
    "WORLD_92", 252,
    "WORLD_93", 253,
    "WORLD_94", 254,
    "WORLD_95", 255,
    "X", 120,
    "Z", 122,
};


int keystr2num(char *k)
{
    int sz=strlen(k);
    char s[sz+1];
    for (int i=0;i<sz;i++) s[i]=toupper(k[i]);
    s[sz]=0;
    k__t key={s,0};
    int n=sizeof(keys)/sizeof(k__t);
    k__t *found=bsearch(&key, keys, n, sizeof(k__t), cmpk);

    return found?found->k:0;
}

static char *keynum2str(int k)
{
    k%=SDLK_LAST;
    static int init=0;
    static char *str[SDLK_LAST];
    if (!init){
	init=1;
	int n=sizeof(keys)/sizeof(k__t);
	for (int i=0;i<n;i++) {
	    str[keys[i].k]=keys[i].str;
	}
    }

    return str[k];
}


/*
  vu {
  nam string;
  cmd string;
  grid integer;
  fg vec4;
  bg vec4;
  rect x y w h;
  btnid integer;

  bind key cmd;
  bind key cmd;
  btn {
  nam string;
  txt string;
  cmd string;
  rect x y w h;
  scale float;
  fg vec4;
  bg vec4;
  }

  btn {
  }

  }

  guied {
  bind key cmd;
  }

*/

static void outputbtn(guibtn_t *b, FILE *f)
{
    fprintf(f, 
	    "	btn {\n\
		nam %s;\n\
		txt %s;\n\
		cmd %s;\n\
		rect %i %i %i %i;\n\
		scale %f;\n\
		fg %f %f %f %f;\n\
		bg %f %f %f %f;\n\
	}\n",
	    b->nam,
	    b->txt,
	    b->cmd,
	    b->x,b->y,b->w,b->h,
	    b->txtscale,
	    b->fg[0],b->fg[1],b->fg[2],b->fg[3],
	    b->bg[0],b->bg[1],b->bg[2],b->bg[3]);
}


static void outputvu(rndvu_t *v, FILE *f)
{
    fprintf(f, 
	    "vu {\n\
	nam %s;\n\
	cmd %s;\n\
	grid %i;\n\
	fg %f %f %f %f;\n\
	bg %f %f %f %f;\n\
	rect %i %i %i %i;\n\
	btnid %i;\n",
	    v->nam,
	    v->cmd, 
	    v->g, 
	    v->fg[0],v->fg[1],v->fg[2],v->fg[3],
	    v->bg[0],v->bg[1],v->bg[2],v->bg[3],
	    v->x,v->y,v->w,v->h,
	    v->btnid);

    guibtn_t *b=v->btnlist;
    while(b){
	outputbtn(b,f);
	b=b->next;
    }

    // key bindings
    for (int i=0; i<SDLK_LAST; i++) {
	if (v->kcmd[i]) {
	    char *k=keynum2str(i);
	    if (k) {
		fprintf(f, 
			"	bind %s %s;\n", k, v->kcmd[i]);
	    }
	}
    }

    fprintf(f, 
	    "}\n");
	
}


void vu_saveall(char *nam)
{
    FILE *f=NULL;
    if (!nam || !*nam) nam="default.view";
    f=fopen(nam, "w+");

    if (!f) {
	cs_printf("Cannot save view");
	return;
    }

    rndvu_t *v=G.vulist;
    if (!v) return;

    do {
	outputvu(v,f);
	v=v->next;
    } while(v!=G.vulist);

    fclose(f);
}


void vu_loadall(char *nam)
{
    char *source=loadtext_alloca(nam);
    if (!source) return;

	

}


static void guied_mbfunc(SDL_MouseButtonEvent button, int down)
{
    // ignore all mouse events during mode transitions
    //	if(G.mt) return;
    if (!G.guied_vu) return;

    if (button.button==SDL_BUTTON_LEFT) {  // DEBUGGED: SDL_BUTTON(1) =>  SDL_BUTTON_LEFT 
	if (down) {
	    G.guied_LMBdown = 1;

	    if (G.guied9) {
		int x=button.x*3/G.screenw;
		int y=button.y*3/G.screenh;

		G.guied9_dragpnt[0]=button.x; // - x*G.screenw/3;
		G.guied9_dragpnt[1]=button.y; // - y*G.screenw/3;				
	    }
	}
	else{
	    G.guied_LMBdown = 0;

	    if (G.guied9) {
		if (G.guied_op == VU_CHANGE_GRID) { // drag the selected view to another grid
		    // put the view into the grid the cursor is currently in
		    int x = button.x*3/G.screenw;
		    int y = button.y*3/G.screenh;
		    int g = x + 3*y;

		    rndvu_t *vu=G.guied_vu;
		    int oldg = vu->g;
		    if (g != oldg) {

#if 0
			// remove the view from its original grid
			rndvu_t *v;
			for (v=vu->next; v!=vu; v=v->next) {
			    if(v->g == g) break;
			}

			if(v!=vu){
			    G.guied_vuall[oldg] = v;
			    vu_loadtn(v);
			}
			else{
			    G.guied_vuall[oldg] = NULL;
			    vu_cleartn(oldg);
			}

			// add the view to grid g

			if (G.guied_vuall[g]) {
			    if (vu->next != vu) {
				vu->next->prev=vu->prev;
				vu->prev->next=vu->next;

				rndvu_t *v2=G.guied_vuall[g];
				v2->prev->next=vu;

				vu->prev=v2->prev;
				vu->next=v2;

				v2->prev=vu;
							
			    }
			}
			G.guied_vuall[g] = vu;
			vu->g = g;
#endif

			vu_bgn(vu->nam);
			{
			    vu_grid(g);
			}
			vu_end();

			if (!G.guied_vuall[oldg]) vu_cleartn(oldg);
			else vu_loadtn(G.guied_vuall[oldg]);
		    }
		    G.guied_op = 0;
		}
	    }

	    guied_refresh(button.x, button.y);
	}
	return;
    }

    if (G.guied_LMBdown) return; // when LMB is down, ignore all other mouse button events.

    if (button.button==SDL_BUTTON_WHEELDOWN) { // scroll down
	if (down) {
	    if(G.keydown[SDLK_LCTRL] || G.keydown[SDLK_RCTRL]){
		if(G.guied_op == BTN_MOVE){ // ctrl + mouse wheel down == button rank down
		    btn_rankdown();
		}
	    }
	    else{ // next view
		// save the view's thumbnail to disk
			
		if(!G.guied9){ // 1-view

#if 0
		    if(G.guied_changed){
			vu_savetn(G.guied_vu); // TODO
		    }
#endif

		}
			
		int g = G.guied_vu->g;
		for (rndvu_t *v=G.guied_vu->next; v!=G.guied_vu; v=v->next) {
		    if (v->g == g) {
			G.guied_vuall[g] = G.guied_vu = v;
			//						vu_loadtn(v); // TODO
			guied_refresh(button.x, button.y);
			break;
		    }
		}
	    }
	}
    }
    else if(button.button==SDL_BUTTON_WHEELUP){ // scroll up
	if (down) {
	    if (G.keydown[SDLK_LCTRL] || G.keydown[SDLK_RCTRL]) { // ctrl + mouse wheel up == button rank up
		if (G.guied_op == 15) {
		    btn_rankup();
		}
	    }
	    else { // previous view
		if (!G.guied9) {
#if 0
		    if (G.guied_changed) {
			vu_savetn(G.guied_vu); // TODO
		    }
#endif
		}

		int g = G.guied_vu->g;
		for (rndvu_t *v=G.guied_vu->prev; v!=G.guied_vu; v=v->prev) {
		    if (v->g == g) {
			G.guied_vuall[g] = G.guied_vu = v;
			//						vu_loadtn(v); // TODO
			guied_refresh(button.x, button.y);
			break;
		    }
		}
	    }
	}
    }
    else if(button.button==SDL_BUTTON_RIGHT){ // rmb -- switch between 1 and 9 views
	if(down){
	    if(G.guied9){
		G.mt = MT_GUIED9_TO_GUIED1;
	    }
	    else{
		G.mt = MT_GUIED1_TO_GUIED9;
	    }

	}
    }

}


static void guied_mmfunc(SDL_MouseMotionEvent motion)
{
    // SDL mouse coord's origin is the upperleft corner.

    if(!G.guied_vu && !G.guied9){
	return;
    }

    if(!G.guied_LMBdown){ // not dragging/resizing anything => find what's under the cursor
	if(!G.guied9){
	    // find which button is in focus
	    // brute force
	    int x = motion.x;
	    int y = motion.y; //G.screenh-motion.y;
	    G.guied_op=0;
	    G.guied_btn=NULL;
	    for(guibtn_t *b=G.guied_vu->btnlist; b; b=b->next)
		{
		    int left  =b->x-RESIZE_EDGE_WIDTH;
		    int right =b->x+b->w + RESIZE_EDGE_WIDTH;
		    int bottom=b->y-RESIZE_EDGE_WIDTH;
		    int top   =b->y+b->h + RESIZE_EDGE_WIDTH;

		    if(x>= left && x<= right && y>= bottom && y<= top)
			{
			    if(x>=b->x && x<=b->x+b->w && y>=b->y && y<=b->y+b->h){
				G.guied_op=BTN_MOVE;
				G.guied_btn = b;
				return;
			    }
			    else{
				int op=0;
				int cx=b->x+b->w/2;
				int cy=b->y+b->h/2;
				int ew=RESIZE_EDGE_WIDTH/2; // edge width
				if(x<b->x || x>b->x+b->w){
				    if(x<cx) op |= BTN_LEFT;
				    else op |= BTN_RIGHT;
				}
				if(y<b->y || y>b->y+b->h){
				    if(y<cy) op |= BTN_BOTTOM;
				    else op |= BTN_TOP;
				}
#if 0
						
				if(!G.guied_op){
				    G.guied_btn = b;
				    return;
				}
#endif


				G.guied_op=op;
				G.guied_btn=b;
			    }
			}
		}

	    // no button is under the cursor => view
	    rndvu_t *v = G.guied_vu;
	    if (x>=v->x && x<=v->x + v->w && y>=v->y && y<=v->y + v->h)
		{
		    int cx=v->x+v->w/2;
		    int cy=v->y+v->h/2;
		    int ew = RESIZE_EDGE_WIDTH/2;
		    int op=G.guied_op;
		    if (y<v->y+RESIZE_EDGE_WIDTH || y>v->y+v->h-RESIZE_EDGE_WIDTH) {
			if(y - v->y < RESIZE_EDGE_WIDTH) op |= 64;
			else if(v->y+v->h - y < RESIZE_EDGE_WIDTH) op |= 128;
		    }

		    if (x<v->x+RESIZE_EDGE_WIDTH || x>v->x+v->w-RESIZE_EDGE_WIDTH) {
			if (x - v->x < RESIZE_EDGE_WIDTH) op |= 16;
			else if (v->x+v->w - x < RESIZE_EDGE_WIDTH) op |= 32;
		    }

		    //				if(!G.guied_op) G.guied_op = VU_CHANGE_GRID;
		}
			   
	}
	else {
	    // find which view is in focus
	    int gx = motion.x*3/G.screenw;
	    int gy = motion.y*3/G.screenh;

	    G.guied_vu = G.guied_vuall[gy*3+gx];
	    G.guied_op = VU_CHANGE_GRID;
	}
	return;
    }


    int dx = motion.xrel;
    int dy = motion.yrel;
    int x = motion.x;
    int y = motion.y; //G.screenh-motion.y;
    guibtn_t *btn = G.guied_btn;
    rndvu_t *v = G.guied_vu;
    if(!G.guied9){ // edit buttons/view
	int op=G.guied_op;
	if (op==BTN_MOVE) {
	    btn->x += dx;
	    btn->y += dy;
	}
	else if (op==VU_MOVE) {
	    v->x += dx;
	    v->y += dy;
	    for(guibtn_t *b=v->btnlist; b; b=b->next){
		b->x += dx;
		b->y += dy;
	    }
	}
	else{
	    if (op & BTN_LEFT) {
		if (x <= btn->x+btn->w-BTN_MIN_WIDTH) {
		    btn->x += dx;
		    btn->w -= dx;
		    if(btn->w < BTN_MIN_WIDTH){
			btn->x += btn->w-BTN_MIN_WIDTH;
			btn->w = BTN_MIN_WIDTH;
		    }
		}
		else {
		    btn->x = btn->x+btn->w-BTN_MIN_WIDTH;
		    btn->w = BTN_MIN_WIDTH;
		}
	    }

	    if (op & BTN_RIGHT) {
		if(x >= btn->x + BTN_MIN_WIDTH){
				
		    btn->w += dx;
		    if(btn->w < BTN_MIN_WIDTH){
			btn->w = BTN_MIN_WIDTH;
		    }
		}
		else {
		    btn->w = BTN_MIN_WIDTH;
		}
	    }

	    if (op & BTN_BOTTOM) {
		if(y <= btn->y + btn->h - BTN_MIN_HEIGHT){
				
		    btn->y += dy;
		    btn->h -= dy;
		    if(btn->h < BTN_MIN_HEIGHT){
			btn->y += btn->h - BTN_MIN_HEIGHT;
			btn->h = BTN_MIN_HEIGHT;
		    }
		}
		else {
		    btn->y = btn->y+btn->h-BTN_MIN_HEIGHT;
		    btn->h = BTN_MIN_HEIGHT;
		}
	    }

	    if (op & BTN_TOP) {
		if(y >= btn->y + BTN_MIN_HEIGHT){
					
		    btn->h += dy;
		    if(btn->h < BTN_MIN_HEIGHT){
			btn->h = BTN_MIN_HEIGHT;
		    }
		}
		else {
		    btn->h = BTN_MIN_HEIGHT;
		}
	    }


	    if (op & VU_LEFT) {
		if(x <= v->x+v->w - VU_MIN_WIDTH) {
				
		    v->x += dx;
		    v->w -= dx;
		    if(v->w < VU_MIN_WIDTH){
			v->x += v->w-VU_MIN_WIDTH;
			v->w = VU_MIN_WIDTH;
		    }
		}
	    }

	    if (op & VU_RIGHT) {
		if(x >= v->x + VU_MIN_WIDTH) {
				
		    v->w += dx;
		    if(v->w < VU_MIN_WIDTH){
			v->w = VU_MIN_WIDTH;
		    }
		}
	    }

	    if (op & VU_BOTTOM) {
		if(y <= v->y + v->h - VU_MIN_HEIGHT) {
				
		    v->y += dy;
		    v->h -= dy;
		    if(v->h < VU_MIN_HEIGHT){
			v->y += v->h - VU_MIN_HEIGHT;
			v->h = VU_MIN_HEIGHT;
		    }
		}
	    }

	    if (op & VU_TOP) {
		if (y >= btn->y + BTN_MIN_HEIGHT) {
				
		    btn->y += dy;
		    if(btn->h < BTN_MIN_HEIGHT){
			btn->h = BTN_MIN_HEIGHT;
		    }
		}
	    }			
	}
    }
    else{
	// move view to other grid
	//		G.guied_vucenter[0] += dx;
	//		G.guied_vucenter[1] += dy;
    }
}


void guied_keyfunc(SDL_keysym keysym, int down)
{
    if (G.mode!=MODE_GUIED ) return;

    int k=keysym.sym;
    if (k==SDLK_ESCAPE) {
	if (down) {
	    cmd_execnow("guied\n");
	}
    }
    else {
	char *c=G.guied_kcmd[k];
	if (c&&*c) {
	    char cmd[strlen(c)+2];
	    cmd[0]=down?'+':'-';
	    strcpy(cmd+1,c);
	    cmd_execnow(cmd);
	}
    }
}


/*
  in gui editting mode, changing views is different from that in view mode
  -- the background of the view being edited is only rendered once at the time it's loaded,
  and stored as a texture.

*/


void rnd_clearall()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
}

static void guied_drawstatusbar()
{
    int h=fnt_getcur()->height;
    int y=0;
    float qbnd[4]={0,y,G.screenw,y+h};

    rnd_scrncoord();

    glEnable(GL_BLEND);
    rnd_quad(G.guied_statusbarbg, qbnd);
    fnt_color4fv(G.guied_statusbarfg);
    fnt_printf(0,y+h-fnt_getcur()->ascender, 0, "VIEW: %s | BUTTON: %s | XY: (%i,%i) | OP: %s | LMB:%i", 
	       G.guied_vu?G.guied_vu->nam:0, 
	       G.guied_btn?G.guied_btn->nam:0,
	       G.mousex, G.mousey,
	       guiop[G.guied_op],
	       G.guied_LMBdown); //TODO
    fnt_flush();
    glDisable(GL_BLEND);
}


/*---------------------------------
  draw 1-view in gui editing mode

  There're currently three modes: console, gui-editting and normal mode.

  Console controls everything and has the highest priority;
  The program consists of GUI and non-GUI stuff. 


  --------------------------------*/
void guied_draw()
{
    rnd_clearall();

    // background is just a texture; it should be already available now
    // the problem is, the view may be resized/moved.

#if 0
    rndvu_t *v=G.guied_vu;

    float tcbnd[4]={0,0,1,1};
    float qbnd[4]={v->x, v->y, v->x+v->w, v->y+v->h};

    rnd_texquad(G.tex_screensz, tcbnd, qbnd);
#endif

    // gui
    vu_drawgui(G.guied_vu);
	
    // status bar
    if (!G.mt) guied_drawstatusbar();
}


/*
  draw overview
*/
void guied_draw9()
{
    vfs_use(-1);
    rnd_clearall();

    // thumbnails
    // each thumbnail is a 128x128 dxt texture. Everytime switching view, it's got updated.
	
    float tcbnd[4]={0,0,1,1};
    float qbnd[4]={0,0,G.screenw, G.screenh};
    uint tn=G.guied_tntex;

    rnd_texquad(tn, tcbnd, qbnd);

    if (!G.mt){
	guied_drawstatusbar();
	vfs_use(-1);
    }

    if(G.guied_LMBdown && G.guied_op==VU_CHANGE_GRID){ // dragging
	glEnable(GL_BLEND);
	float rgba[4]={.5,.5,.5,.5};
	int x=G.guied_vu->g%3;
	int y=G.guied_vu->g/3;
	int w=G.screenw/3;
	int h=G.screenh/3;

	qbnd[0]=x*w;
	qbnd[1]=y*h;
	qbnd[2]=(x+1)*w;
	qbnd[3]=(y+1)*h;
	rnd_quad(rgba, qbnd);

	qbnd[0]+=G.mousex-G.guied9_dragpnt[0];
	qbnd[1]+=G.mousey-G.guied9_dragpnt[1];
	qbnd[2]=qbnd[0]+w;
	qbnd[3]=qbnd[1]+h;
	tcbnd[0]=x/3.0f;
	tcbnd[1]=y/3.0f;
	tcbnd[2]=(x+1)/3.0f;
	tcbnd[3]=(y+1)/3.0f;
	rnd_texquad(tn, tcbnd, qbnd);
	glDisable(GL_BLEND);
    }
}


/*
  draw overview->view transition


  0<= tween <= 1
*/
void mt_guied9vu(rndvu_t *v, float t)
{
    vfs_use(-1);

    // draw background
    // the background zooms in -- only texcoord changes
    float qbnd[4]={0,0,G.screenw,G.screenh};
    float tcbnd[4];
	
    float x=(v->g%3)/3.0f;
    float y=(v->g/3)/3.0f;
    tcbnd[0]=x*t;
    tcbnd[1]=y*t;
    tcbnd[2]=1+(x-2.0/3.0)*t;
    tcbnd[3]=1+(y-2.0/3.0)*t;
	
    rnd_texquad(G.guied_tntex, tcbnd, qbnd);
	
    tcbnd[0]=tcbnd[1]=0;
    tcbnd[2]=tcbnd[3]=1;
    qbnd[0]=x*(1-t)*G.screenw;
    qbnd[1]=y*(1-t)*G.screenh;
    qbnd[2]=((x+1.0/3.0)+t*(2.0/3.0-x))*G.screenw;
    qbnd[3]=((y+1.0/3.0)+t*(2.0/3.0-y))*G.screenh;

    rnd_texquad(G.tex_screensz, tcbnd, qbnd);
}


// exactly the inverse of mt_guied9vu()
void mt_vuguied9(rndvu_t *v, float t)
{
    vfs_use(-1);

    // draw background
    // the background zooms in -- only texcoord changes
    float qbnd[4]={0,0,G.screenw,G.screenh};
    float tcbnd[4];
	
    float x=(v->g%3)/3.0f;
    float y=(v->g/3)/3.0f;
    tcbnd[0]=x*(1-t);
    tcbnd[1]=y*(1-t);
    tcbnd[2]=(x+1.0/3.0)+t*(2.0/3.0-x);
    tcbnd[3]=(y+1.0/3.0)+t*(2.0/3.0-y);

    glColor4f(1,1,1,1);
	
    rnd_texquad(G.guied_tntex, tcbnd, qbnd);

    tcbnd[0]=tcbnd[1]=0;
    tcbnd[2]=tcbnd[3]=1;

    qbnd[0]=x*t*G.screenw;
    qbnd[1]=y*t*G.screenh;
    qbnd[2]=(1+(x-2.0/3.0)*t)*G.screenw;
    qbnd[3]=(1+(y-2.0/3.0)*t)*G.screenh;

    rnd_texquad(G.tex_screensz, tcbnd, qbnd);
}


static void lerp4f(float dest[4], float a[4], float b[4], float t)
{
    for(int i=0; i<4; i++){
	dest[i] = a[i]+(b[i]-a[i])*t;
    }
}


/*
  draw view->view transition
*/
void mt_vuvu(rndvu_t *v1, rndvu_t *v2, float t)
{
    // the rendering uses screen coordinate (i.e. [0,w]x[0,h]

    rnd_clearall();
    vfs_use(-1);

    if(v1->g==v2->g){
	float rgba[4];

	lerp4f(rgba, v1->bg, v2->bg, t);

	float qbnd[4];
	qbnd[0]=v1->x + (v2->x-v1->x)*t;
	qbnd[1]=v1->y + (v2->y-v1->y)*t;
	qbnd[2]=qbnd[0] + v1->w+(v2->w - v1->w)*t;
	qbnd[3]=qbnd[1] + v1->h+(v2->h - v1->h)*t;

	float tcbnd[4]={0,0,1,1};


	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	t=0;
	// new
	glColor4f(1,1,1,1);
	rnd_texquad(G.tex_screensz1, tcbnd, qbnd);

	glColor4f(1,1,1,1);
	rnd_texquad(G.tex_screensz, tcbnd, qbnd);

	glDisable(GL_BLEND);
    }
    else{ // jump from one grid to another

	// just draw two quads
	int x1=v1->g%3;
	int x2=v2->g%3;
	int y1=v1->g/3;
	int y2=v2->g/3;

	int d;
	d=x2-x1;
	int l2r = -.5*d + 1.5*(d%2);
		
	d=y2-y1;
	int b2t = -.5*d + 1.5*(d%2);

	float tx=0,ty=0;
		
	if (l2r) tx=t;
	if (b2t) ty=t;

	float tcbnd[4]={0,0,1,1};
	float qbnd1[4], qbnd2[4];
	int scrw=G.screenw;
	int scrh=G.screenh;

	if (l2r==1) {
	    qbnd1[2]=qbnd2[0]=scrw*(1-tx);
	    qbnd1[0]=qbnd1[2]-scrw;
	    qbnd2[2]=qbnd2[0]+scrw;
	}
	else if (l2r==-1) {
	    qbnd1[0]=qbnd2[2]=scrw*tx;
	    qbnd1[2]=qbnd1[0]+scrw;
	    qbnd2[0]=qbnd2[2]-scrw;
	}
	else{
	    qbnd1[0]=qbnd2[0]=0;
	    qbnd1[2]=qbnd2[2]=scrw;
	}

	if(b2t==1){
	    qbnd1[3]=qbnd2[1]=scrh*(1-ty);
	    qbnd1[1]=qbnd1[3]-scrh;
	    qbnd2[3]=qbnd2[1]+scrh;
	}
	else if(b2t==-1){
	    qbnd1[1]=qbnd2[3]=scrh*ty;
	    qbnd1[3]=qbnd1[1]+scrh;
	    qbnd2[1]=qbnd2[3]-scrh;
	}
	else{
	    qbnd1[1]=qbnd2[1]=0;
	    qbnd1[3]=qbnd2[3]=scrh;
	}


	rnd_texquad(G.tex_screensz, tcbnd, qbnd1);

	// new
	rnd_texquad(G.tex_screensz1, tcbnd, qbnd2);

    }
}


static void cs_keyfunc(SDL_keysym keysym)
{
    switch(keysym.sym){
    case SDLK_LEFT:
	mb_left();
	break;
    case SDLK_RIGHT:
	mb_right();
	break;
    case SDLK_UP:
	mb_prev();
	break;
    case SDLK_DOWN:
	mb_next();
	break;
    case SDLK_DELETE:
	mb_del();
	break;
    case SDLK_BACKSPACE:
	mb_backspace();
	break;
    case SDLK_RETURN:
	mb_enter();
	break;
    case SDLK_PAGEUP:
	cs_scrollup();
	break;
    case SDLK_PAGEDOWN:
	cs_scrolldown();
	break;
    case SDLK_END:
	cs_end();
	break;
    default:
	if(keysym.unicode){
	    mb_ins(keysym.unicode);
	}
	break;
    }	


}


static void cs_mbfunc(SDL_MouseButtonEvent button, int down)
{
    if(!down) return;

    // scroll up
    if(button.button==4){
	cs_scrollup();
    }
    // scroll down
    else if(button.button==5){
	cs_scrolldown();
    }
}



/*
  If key==`:
  If in console mode: quit console mode
  Else: enter console mode

  Else:
  	

*/

static void keyfunc(SDL_keysym keysym, int down)
{
    //	if(G.mt) return;

    // console?
    // ignore text input
    int k=keysym.sym;
    if(k == SDLK_BACKQUOTE){ // highest priority
	if(down){
	    if(G.cs_enabled){
		G.cs_enabled=0;
		G.cs_winh=0;
		SDL_EnableKeyRepeat(0,0);
		SDL_EnableUNICODE(0);
		//				SDL_WM_GrabInput(SDL_GRAB_ON);
	    }
	    else{
		G.cs_enabled=1;
		G.cs_csrblink=31;
		// TODO: console transition

		SDL_EnableKeyRepeat(SDL_DEFAULT_REPEAT_DELAY,
				    SDL_DEFAULT_REPEAT_INTERVAL);
		SDL_EnableUNICODE(1);
		//				SDL_WM_GrabInput(SDL_GRAB_OFF); // turn off grab input, so the cursor can move outside the window
	    }
	}
    }
    else{
	if(G.cs_enabled){
	    if(down){
		cs_keyfunc(keysym);

		// don't forget the blinking cursor
		G.cs_csrblink=31;
	    }
	}
	else{
	    switch(G.mode){
	    case MODE_GUIED:
		guied_keyfunc(keysym, down);
		break;
	    case MODE_VU:
		vu_keyfunc(keysym, down);
		break;
	    case MODE_FILTER:
		filter_keyfunc(keysym, down);
		break;
	    case MODE_EDITOR:
		editor_keyfunc(keysym, down);
		break;
	    }
	}
    }      

    G.keydown[k] = down;
}


/*
  If in console mode, all button events are only processed by console.

  Otherwise:
  If in transition mode, all button events are ignored.
  If in gui editting mode or normal modes, let the respected mode to process the events

*/
static void mousebuttonfunc(SDL_MouseButtonEvent button, int down)
{
    button.y=G.screenh-button.y;
    if(G.cs_enabled){
	cs_mbfunc(button, down);
	return;
    }

    //	if(G.mt) return;

    if(G.mode == MODE_GUIED){
	guied_mbfunc(button, down);
    }
    else if(G.mode == MODE_VU){
	vu_mbfunc(button, down);
    }
    else if (G.mode==MODE_FILTER) {
	filter_mbfunc(button, down);
    }
    else if (G.mode==MODE_EDITOR) {
	editor_mbfunc(button, down);
    }

}

/*
  If in console mode, all mouse events are ignored.

  Otherwise:
  If in transition mode, all mouse events are ignored.
  If in gui editting mode or normal modes, let the respected mode to process the events  

*/
static void mousemotionfunc(SDL_MouseMotionEvent motion)
{
    motion.y=G.screenh-motion.y;
    motion.yrel*=-1;

    if(G.cs_enabled){
	goto donothing;
    }

    //	if(G.mt) goto donothing;

    if(G.mode == MODE_GUIED){
	guied_mmfunc(motion);
    }
    else if(G.mode == MODE_VU){ // mouse binding / GUI
	vu_mmfunc(motion);
    }
    else if (G.mode == MODE_FILTER) {
	filter_mmfunc(motion);
    }
    else if (G.mode == MODE_EDITOR) {
	editor_mmfunc(motion);
    }

  donothing:
    G.mousex=motion.x;
    G.mousey=motion.y; //G.screenh-motion.y;
}


/*
  eventloop() may be called multiple times each frame.
*/
static int eventloop()
{
    SDL_Event e;
    while(SDL_PollEvent(&e)){
	switch(e.type){
	case SDL_KEYDOWN:
	case SDL_KEYUP:
	{
	    SDL_KeyboardEvent k=e.key;
	    keyfunc(k.keysym, e.type==SDL_KEYDOWN);
	    break;
	}
	case SDL_MOUSEMOTION:
	{
	    SDL_MouseMotionEvent motion=e.motion;
	    mousemotionfunc(motion);
	    break;
	}
	case SDL_MOUSEBUTTONDOWN:
	case SDL_MOUSEBUTTONUP:
	{
	    SDL_MouseButtonEvent button=e.button;
	    mousebuttonfunc(button, e.type==SDL_MOUSEBUTTONDOWN);
	    break;
	}

	case SDL_QUIT:
	    break;

	}
    }
}


/*
  
 */


void mt_init()
{
	
}

static int curmsec()
{
    struct timeval	tv;
	
    gettimeofday(&tv, NULL);
	
    if(!G.timebase){
	G.timebase = tv.tv_sec;
	return tv.tv_usec/1000;
    }

    G.timecur = (tv.tv_sec - G.timebase)*1000 + tv.tv_usec/1000;
    return G.timecur;
}



static void frame()
{	
    int		ms=0;
    static int	lastms=0;
    int		minms=17; // TODO: svar

    int test1=lastms;
    do {
	eventloop(); // cmds bound to key are added to cbuf_
	ms = curmsec()-lastms;
    } while(ms < minms);
    lastms += ms;
    int test2=lastms;

    //	cs_printf("frame time=%ims\n", test2-test1);

    cmd_exec(); // Stuff like loading maps are done here. May also contain view switch/editing commands.

    // eventloop() may call for mode transition (directly if in GUI editting mode, or via command normal mode)

    rnd_clearall();

    if(G.mt && G.mt_tween){ // mode transition. just play some fancy animation
	float t = G.mt_tween;
	switch(G.mt){

	case MT_VU_TO_VU:
	    mt_vuvu(G.mt_vuvu[0], G.mt_vuvu[1], t);
	    if (t>=1.0f) {
		G.vu_cur=G.mt_vuvu[1];
	    }
	    break;

	case MT_GUIED1_TO_VU:
	    //			mt_vuvu(G.guied_vu, G.vu_cur, t);
	    break;

	case MT_GUIED9_TO_VU:
	    //			mt_guied9vu(G.vu_cur, t);
	    break;

	case MT_GUIED9_TO_GUIED1:
	    mt_guied9vu(G.guied_vu, t);
	    if (t>=1.0f) {
		G.guied9=0;
	    }
	    break;

	case MT_GUIED1_TO_GUIED9:
	    mt_vuguied9(G.guied_vu, t);
	    if (t>=1.0f) {
		G.guied9=1;
	    }
	    break;

	}

	if(t>=1.0f){
	    G.mt_tween=0;
	    G.mt=0;
	}
	else {
	    G.mt_tween += G.mt_dt;
	}

	//		SDL_GL_SwapBuffers(); // DEBUGGED: should swap framebuffer here

	goto console;
    }
    if(G.mode==MODE_GUIED){
	if (G.guied9) guied_draw9();
	else guied_draw();

	if(G.mt && !G.mt_tween){
	    if(G.mt==MT_GUIED1_TO_VU){
		vu_savetn(G.guied_vu);
		rnd_tg("t", G.tex_screensz1);
		cmd_execnow(G.vu_cur->cmd);
		vu_drawgui(G.vu_cur);
	    }
	    else if(G.mt==MT_GUIED9_TO_VU){
		rnd_tg("t", G.tex_screensz);
		cmd_execnow(G.vu_cur->cmd);
		vu_drawgui(G.vu_cur);
	    }


	    else if(G.mt==MT_GUIED9_TO_GUIED1){
		rnd_tg("td", G.tex_screensz, G.dtex_screensz);
		rnd_clearall();
		cmd_execnow(G.guied_vu->cmd);
		vu_drawgui(G.guied_vu);
		rnd_tg(0);
	    }
	    else if(G.mt==MT_GUIED1_TO_GUIED9){
		rndvu_t *v=G.guied_vu;
		vu_savetn(v);
		G.guied_vuall[v->g]=v;

		for (int i=0; i<9; i++) {
		    if (G.guied_vuall[i]) {
			vu_loadtn(G.guied_vuall[i]);
		    }
		}
	    }

	    G.mt_tween=G.mt_dt;
	}
    }
    else if (G.mode==MODE_VU) {
	cmd_execnow(G.vu_cur->cmd);
	vu_drawgui(G.vu_cur);

	if (G.mt && !G.mt_tween) {
	    vu_savetn(G.vu_cur);
	    if (G.mt==MT_VU_TO_VU) {
		rnd_tg("td", G.tex_screensz1, G.dtex_screensz);
		rnd_clearall();
		cmd_execnow(G.mt_vuvu[1]->cmd);
		vu_drawgui(G.mt_vuvu[1]);
		rnd_tg(0);
	    }

	    G.mt_tween=.1;
	}
    }
    else if (G.mode==MODE_FILTER) {
	// you can only edit filters in view mode (not GUIed mode)
	cmd_execnow(G.vu_cur->cmd);
    }
    else if (G.mode==MODE_EDITOR) {
	cmd_execnow(G.ed.cmd);
	editor_draw();
    }

  console:
    if(G.cs_enabled){
	cs_draw();
    }

    SDL_GL_SwapBuffers();
}		 


void guied_addvu()
{
    if (G.mode==MODE_GUIED) {
		
    }
}



/////////////////////////////////////////////////////////////////
/// Built-in editors
void c_editor_toggle()
{
    if (G.mode!=MODE_EDITOR) {
	G.prevmode=G.mode;
	G.mode=MODE_EDITOR;
    }
    else {
	G.mode=G.prevmode;
    }
}

void c_guied_toggle()
{
    // TODO: 
    // TESTBGN

    if (G.mode!=MODE_GUIED) {
	G.prevmode=G.mode;
	G.mode=MODE_GUIED;
		
	// TESTEND
		
	G.guied_vu=G.vu_cur;
	G.mt=MT_VU_TO_GUIED_NORMAL;

	cs_printf("Enters GUI editing mode\n");
    }
    else {
	G.mode=G.prevmode;
	G.guied_vu=NULL;

	cs_printf("Quits GUI editing mode\n");
    }
}



/*
  c_* means command functions
*/
void c_print()
{
    char buf[256];
    if (cmd_gets(buf)) {
	cs_printf("%s\n", buf);
    }
}


void c_guied_btnadd()
{
    if (G.mode!=MODE_GUIED || !G.guied_vu) return;

    char nam[256]={0};
    if (cmd_gets(nam)) {
	char *c=nam;
	while (isalnum(*c)) ++c;
	*c=0;
    }

    if (G.mode==MODE_GUIED) {
	vu_bgn(G.guied_vu->nam);
	{
	    btn_force_bgn(nam);
	    btn_end();
	}
	vu_end();
    }
}


void c_guied_btndel()
{
    if (G.mode!=MODE_GUIED || !G.guied_vu) return;

    char nam[256]={0};
    if (cmd_gets(nam)) { // delete the named button
	vu_bgn(G.guied_vu->nam);
	{
	    btn_bgn(nam);
	    {
		btn_del();
	    }
	    btn_end();
	}
	vu_end();
    }
    else if (G.guied_btn) { // delete the button currently in focus
	vu_bgn(G.guied_vu->nam);
	{
	    btn_bgn(G.guied_btn->nam); 
	    {
		btn_del();
	    }
	    btn_end();
	}
	vu_end();		
    }

}


void c_guied_btntxt()
{
    if (G.mode!=MODE_GUIED || !G.guied_vu || !G.guied_btn) return;

    char txt[256];
    if (cmd_gets(txt)) {
	vu_bgn(G.guied_vu->nam);
	{
	    btn_bgn(G.guied_btn->nam);
	    {
		btn_txt(txt);
	    }
	    btn_end();
	}
	vu_end();
    }

}


// TODO:
void c_guied_btncmd()
{
    if (G.mode!=MODE_GUIED || !G.guied_vu || !G.guied_btn) return;

    char cmd[256];
    if (cmd_gets(cmd)) {
	vu_bgn(G.guied_vu->nam);
	{
	    btn_bgn(G.guied_btn->nam);
	    {
		btn_cmd(cmd);
	    }
	    btn_end();
	}
	vu_end();
    }
    else {
	cs_printf(G.guied_btn->cmd);
    }
}


void c_guied_btnbg()
{
    if (G.mode!=MODE_GUIED || !G.guied_vu || !G.guied_btn) return;

    for (int i=0; i<4; i++) {
	float color;
	if (cmd_getf(&color)){
	    G.guied_btn->bg[i]=color;
	}
	else break;
    }
}


void c_guied_addview()
{
    if (G.mode!=MODE_GUIED) return;

    char nam[256];
    if (!cmd_gets(nam)) {
	strcpy(nam, "View");
	sprintf(nam+4, "%i", G.vu_id++);
    }

    if (G.guied9) {
	int x=G.mousex*3/G.screenw;
	int y=G.mousey*3/G.screenh;
	int g=3*y+x;

	if (g<0) g=0;
	if (g>8) g=8;
		
	vu_force_bgn(nam);
	vu_grid(g);
	G.guied_vuall[g]=G.vu_edit;
	vu_end();
    }
    else {
	vu_force_bgn(nam);
	if (G.guied_vu) vu_grid(G.guied_vu->g);
	G.guied_vuall[G.guied_vu->g]=G.vu_edit;
	vu_end();
    }
}
	


void c_vu_gotoview()
{
    if (G.mode!=MODE_VU || !G.vu_cur) return;

    char nam[256];
    if (cmd_gets(nam)) {
	rndvu_t *v=cvar_getp__(nam, VU_CKSUM);
	if (!v) {
	    cs_printf("Cannot find view \"%s\"\n", nam);
	    return;
	}

	if (v==G.vu_cur) return;

	G.mt_vuvu[0]=G.vu_cur;
	G.mt_vuvu[1]=v;

	G.mt=MT_VU_TO_VU;
    }

}



void c_editor_addnook()
{
    
}



/*
  TODO:

  bind key cmd

  Current view will be affected
*/
void c_bind()
{
    if (G.mode==MODE_GUIED) {
    }
    else if (G.mode==MODE_VU) {
    }
}

/*
  binda view key cmd
*/
void c_binda()
{
    char nam[256];
    char *vnam, *kstr, *cmd;
    int k;
    if (cmd_gets(nam)) {
	char *c=nam;

	// view
	vnam=c;
	while (isalnum(*c)||*c=='_') ++c;

	if (!*c){
	    cs_printf("Usage: binda ViewName KeyString Cmd\n");
	    return;
	}

	*c++=0;

	// keystring
	kstr=c;
	while (isalnum(*c)||*c=='_') {
	    *c=toupper(*c);
	    ++c;
	}

	if (!*c){ // error
	    cs_printf("Usage: binda ViewName KeyString Cmd\n");
	    return;
	}

	*c++=0;
		
	k=keystr2num(kstr);
	if (!k){
	    cs_printf("Unknown key string %s\n", kstr);
	    return;
	}

	// command
	cmd=c;
	vu_bgn(nam);
	{
	    vu_kcmd(k,cmd);
	}
	vu_end();
    }
}


void c_vu_saveview()
{
    vu_saveall(NULL);
}


void c_guiedbind()
{
    if (G.mode!=MODE_GUIED) return;

    char kc[256];
    if (cmd_gets(kc)) {
	char *c=kc;

	while (isalnum(*c) || *c=='_') {
	    *c=toupper(*c);
	    ++c;  // key
	}
	if (!*c) return;  // no command
	*c++=0;

	int k=keystr2num(kc);
	if (k && *c) {
	    if (G.guied_kcmd[k]) free(G.guied_kcmd[k]);
	    G.guied_kcmd[k]=strdup(c);
	}
    }
}


void c_quit()
{
    exit(0);
}


void vu_init()
{
    G.btn_shader=vfs_load("shaders/guibuttons.c");
    vu_force_bgn("hello");
    {
	btn_force_bgn("world");
	btn_end();
		
	btn_force_bgn("world2");
	btn_end();

	vu_grid(6);
    }
    vu_end();


    vu_force_bgn("yeah");
    {
	btn_force_bgn("oh yeah");
	btn_end();
		
	btn_force_bgn("oh yeah2");
	btn_end();

	btn_force_bgn("oh yeah3");
	btn_end();
    }
    vu_end();

    vu_use("hello");

    for(int i=0;i<4;i++){
	G.guied_statusbarfg[i]=1;
	G.guied_statusbarbg[i]=.8;
    }

    G.vu_tnsz=128;

    cmd_add("btnadd", c_guied_btnadd);
    cmd_add("btndel", c_guied_btndel);
    cmd_add("btntxt", c_guied_btntxt);
    cmd_add("btnbg",  c_guied_btnbg);
    cmd_add("btncmd", c_guied_btncmd);
    cmd_add("gotoview", c_vu_gotoview);
    cmd_add("guiedbind", c_guiedbind);
    cmd_add("addview", c_guied_addview);
    cmd_add("saveview", c_vu_saveview);
    cmd_add("binda", c_binda);

    G.mt_dt=0.1;
}


static void dll_init()
{
    G.dll=dlopen("./game.so", RTLD_NOW);
    if (!G.dll) {
	cs_printf("Cannot load game.so\n");
	cs_printf("%s\n", dlerror());
	return;
    }

    cvec_init();

    void (*init)()=dlsym(G.dll, "dllInit");
    if (!dlerror()) {
	if (init) {
	    init();
	}
    }
    else {
	cs_printf("Cannot find dllInit()!\n");
    }
}

static void dll_shutdown()
{
    if (G.dll) dlclose(G.dll);
}

static void cmd_init()
{
    cmd_add("guied", c_guied_toggle);
    cmd_add("editor", c_editor_toggle);
    cmd_add("print", c_print);
    cmd_add("bind", c_bind);
    cmd_add("quit", c_quit);
}

int main(int argc, char *argv[])
{
    vid_init();
    cvar_init();
    vfs_init();
    vb_init();
    cs_init();
    cmd_init();
    fnt_init();
    rt_init();
    filter_init();
    vu_init();
    editor_init();

    dll_init();

    while(1){
	frame();
    }

    dll_shutdown();
    fnt_shutdown();
    cvec_shutdown();
    cvar_shutdown();
    vb_shutdown();
    cs_shutdown();
    return 0;
}
