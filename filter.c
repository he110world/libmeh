static void filter_save(filter_t *ft)
{
    if (!ft) return;

    FILE *f=fopen(ft->nam, "w");
    if (!f) return;

    int n=arr_len(ft->nodes);

    fprintf(f, "/* FILTER V5 */\n");
	
    // remove unused IDs
    int idx[n];
    int k=0;
    for (int i=0; i<n; ++i) {
	filternode_t *node=ft->nodes+i;
	if (node->indeg!=-1) {
	    idx[i]=k++;
	}
    }

    // input nodes
    for (int i=0; i<arr_len(ft->in); ++i) {
	fprintf(f, "in {\n");
	fprintf(f, "	id %i\n", idx[ft->in[i].id]);
	fprintf(f, "	name %s\n", ft->in[i].nam);
	fprintf(f, "}\n");
    }

    // output nodes
    for (int i=0; i<arr_len(ft->out); ++i) {
	fprintf(f, "out {\n");
	fprintf(f, "	id %i\n", idx[ft->out[i].id]);
	fprintf(f, "	name %s\n", ft->out[i].nam);
	fprintf(f, "}\n");
    }

    // other nodes
    for (int i=0; i<n; ++i) {
	filternode_t *node=ft->nodes+i;
	if (node->indeg==-1) continue;

	if (node->type==FILTER_SHADER) {
	    fprintf(f, "shader {\n");
	    fprintf(f, "	id %i\n", idx[i]);
	    fprintf(f, "	name \"%s\"\n", node->shader.nam);
	    fprintf(f, "	mask %i %i %i %i\n", 
		    node->shader.colormask[0],
		    node->shader.colormask[0], 
		    node->shader.colormask[0], 
		    node->shader.colormask[0]);
	    fprintf(f, "}\n");

	    // edge in
	    for (int j=0; j<node->shader.ni; ++j) {
		int bid=node->shader.in[j].b;
		if (bid!=-1) {
		    fprintf(f, "edge {\n");
		    fprintf(f, "	source %i[%i]\n", idx[bid], 0);
		    fprintf(f, "	target %i[%i]\n", idx[i], j);
		    fprintf(f, "}\n");
		}
	    }

	    // edge out
	    for (int j=0; j<node->shader.no; ++j) {
		int bid=node->shader.out[j].b;
		if (bid!=-1) {
		    fprintf(f, "edge {\n");
		    fprintf(f, "	target %i[%i]\n", idx[i], j);
		    fprintf(f, "	source %i[%i]\n", idx[bid], 0);
		    fprintf(f, "}\n");
		}
	    }

	}
	else if (node->type==FILTER_BUFFER) {
	    fprintf(f, "buffer {\n");
	    fprintf(f, "	id %i\n", idx[i]);
	    fprintf(f, "}\n");
	}

    }

    fclose(f);

}

static int isdag(filter_t *ft)
{
    int n=arr_len(ft->nodes);
    int indeg[n];
    short stack[n];
    int ns=0;

    // init indeg[]
    for (int i=0; i<n; ++i) {
	indeg[i]=ft->nodes[i].indeg;
	if (indeg[i]==0) { 	// find start nodes (indeg==0)	-- not only input nodes, but also shaders with no input
	    stack[ns++]=i;
	}
    }

    int cnt=0;
    while (ns) {
	short id=stack[--ns];
	filternode_t *node=ft->nodes+id;
	++cnt;
	for (int i=0; i<node->outdeg; ++i) {
	    int k=node->inout[63-i];
	    if (!--indeg[k]) stack[ns++]=k;
	}
    }

    return cnt==n;
}

static int filter_update(filter_t *ft)
{
    if (!isdag(ft)) {
	return 1;
    }

    // find valid shaders
    int n=arr_len(ft->nodes);
    int indeg[n];
    short stack[n];
    int ns=0;

    // init indeg[]
    for (int i=0; i<n; ++i) {
	indeg[i]=ft->nodes[i].indeg;
	if (ft->nodes[i].type==FILTER_IN) { // push input nodes into stack -- traverse from input nodes
	    if (ft->nodes[i].indeg!=0) {
		printf("node%i=%i\n", i,ft->nodes[i].indeg);
	    }
	    assert(ft->nodes[i].indeg==0);
	    stack[ns++]=i;
	}
    }

    // (incomplete) topological sort
    short valid[n];
    int nv=0;
    while (ns) {
	short id=stack[--ns];
	filternode_t *node=ft->nodes+id;
	if (node->type==FILTER_SHADER || node->type==FILTER_OUT) {
	    valid[nv++]=id;
	}
	for (int i=0; i<node->outdeg; ++i) {
	    int k=node->inout[63-i];
	    if (!--indeg[k]) stack[ns++]=k;
	}
    }

    arr_popall(ft->validshaders);
    arr_pushn(ft->validshaders, valid, nv);

    return 0;
}



static void filter_allocnode(filter_t *ft, int id)
{
    int sz=arr_len(ft->nodes);
    if (sz>id) {
	// TODO: node id conflict

	filternode_t *node=ft->nodes+id;
	arr_popall(node->shader.in);
	arr_popall(node->shader.out);
	if (node->shader.nam) {
	    free(node->shader.nam);
	    node->shader.nam=NULL;
	}
	return;
    }

    int n=id-sz+1;
    arr_pushn(ft->nodes,NULL,n);
    for (int i=sz; i<sz+n; ++i) { // invalidate all newly allocated nodes
	filternode_t *node=ft->nodes + i;
	memset(node,0,sizeof(filternode_t));
	node->indeg=-1;
    }
}

static int filter_addnode(filter_t *ft)
{
    int sz=arr_len(ft->nodes);
    if (ft->hasfreenode) {
	for (int i=0; i<sz; ++i) {
	    filternode_t *node=ft->nodes+i;
	    if (node->indeg==-1) {
		return i;
	    }
	}
	assert(0); // shouldn't arrive here
    }
    else {
	arr_pushn(ft->nodes,NULL,1);
	return sz;
    }
}

/*
  only called by filter_load()
*/
static int filter_addshader_simple(filter_t *ft, char *nam, char mask[4], int id)
{
    // alloc node
    filter_allocnode(ft, id);
    filternode_t *node=ft->nodes+id;
    node->type=FILTER_SHADER;
    node->indeg=node->outdeg=0;
    for (int i=0;i<4;++i){
	node->shader.colormask[i]=mask[i];
    }

    int type=0;
    if (!strcmp(nam, "ds")) {
	type=SHADER_DS;
    }
    else if (!strcmp(nam, "pack")) {
	type=SHADER_PACK;
    }

    node->shader.type=type;

    // load shader
    if (type==SHADER_DS) {
	node->shader.ni=node->shader.no=1;
	node->shader.p_dscolor=vfs_load("shaders/colordownsampler.c");
	assert(!glGetError()+888*0);
	node->shader.p_dsdepth=vfs_load("shaders/depthdownsampler.c");	
	assert(!glGetError()+999*0);
	
    }
    else if (type==SHADER_PACK) {
	node->shader.ni=0;
	node->shader.no=1;
    }
    else {
	int p=vfs_load(nam);
	if (p<0){
	    return 1;
	}
	vfs_t *v=VFS_PTR(p);
	node->shader.p=p;
	node->shader.ni=v->ntex;
	node->shader.no=v->nrts;
    }

    node->disp.x=rand()%300-150;
    node->disp.y=rand()%300-150;
    node->disp.rank=++ft->disp.maxrank;

    //	assert(!node->shader.in);
    arr_pushn(node->shader.in,NULL,node->shader.ni);
    arr_pushn(node->shader.out,NULL,node->shader.no);
	
    // init I/O buffers
    for (int i=0; i<node->shader.ni; ++i) {
	node->shader.in[i].b=-1;
    }

    for (int i=0; i<node->shader.no; ++i) {
	node->shader.out[i].b=-1;
    }
    return 0;
}


static void filter_addio(filter_t *ft, char *nam, int id, int type)
{
    if (!ft || !nam || !*nam) return;

    // check node id. If the node already exists, return.
    int sz=arr_len(ft->nodes);
    if (id>=0) {
	if (id<sz && ft->nodes[id].indeg!=-1) {
	    return;
	}
    }

    float *f=cvec_get(nam);
    if (!f) return;

    filterio_t *io = type==FILTER_INPUT? ft->in : ft->out;
    int n = arr_len(io);
    int j=-1;

    // quit if the node already exist
    for (int i=0; i<n; ++i) {
	if (io[i].nam) {
	    if (!strcmp(io[i].nam, nam)) {
		io[i].tex=f;
		return;
	    }
	}
	else {
	    j=i;
	}
    }
    if (j==-1) {
	if (type==FILTER_INPUT) {
	    arr_pushn(ft->in,NULL,1);
	    io=ft->in;
	}
	else {
	    arr_pushn(ft->out,NULL,1);
	    io=ft->out;
	}
	j=n;
    }

    io+=j;
    io->nam=strdup(nam);
    io->tex=f;

    // allocate space for the node if necessary
    if (id!=-1) {
	filter_allocnode(ft,id);
	io->id=id;
    }
    else {
	if (ft->hasfreenode) {
	    for (int i=0; i<sz; ++i) {
		filternode_t *node=ft->nodes + i;
		if (node->indeg==-1) {
		    io->id=i;
		    --ft->hasfreenode;
		    break;
		}
	    }
	}
	else {
	    arr_pushn(ft->nodes,NULL,1);
	    io->id=sz;
	}
    }
	
    // init the node
    if (type==FILTER_INPUT) {
	int tex=f[0];
	rtinfo_t info;
	int err=rt_getinfo(tex, &info);
	filternode_t *b=ft->nodes + io->id;
	b->type=FILTER_IN;
	b->indeg=b->outdeg=0;
	b->buffer.tex=tex;
	if (!err) {
	    b->buffer.level=info.level;
	    b->buffer.type=info.type;
	}
	else { // not a render target (e.g. noise texture)
	    b->buffer.type=RT_MAX;
	}

	b->disp.x=rand()%300;
	b->disp.y=rand()%300;
	b->disp.w = G.screenw>>b->buffer.level;
	b->disp.h = (G.screenh>>b->buffer.level)+40; // TODO:
	b->disp.rank = ++ft->disp.maxrank;
    }
    else {
	filternode_t *s=ft->nodes + io->id;
	s->type=FILTER_OUT;
	s->indeg=s->outdeg=0;
	arr_popall(s->shader.out);
	s->shader.no=0;
	s->disp.x=0;
	s->disp.y=0;
	s->disp.rank= ++ft->disp.maxrank;
    }
	
    // no need to call filter_update() -- nothing changed in the graph's connectivity
}


/*
  Number of render targets and multitextures may change if the shader changes.
  So you should check the shader at loading time -- some edges may be invalid then

  TODO: this function isn't done!
*/

static filter_t *filter_load(char *nam, filter_t *ft)
{
    char *code=loadtext_alloca(nam);
    if (!code) return NULL;

    // alloc mem
    if (!ft) {
	ft=calloc(1, sizeof(filter_t));
	ft->nam=strdup(nam);
	ft->disp.zoom=.5; // TODO
    }
    else {
	arr_popall(ft->nodes);
	arr_popall(ft->in);
	arr_popall(ft->out);
    }

    /*
      arr_pushn(ft->nodes,NULL,n);
      for (int i=0; i<n; ++i) { 	// invalidate all nodes first
      ft->nodes[n].indeg=-1;
      }
    */

    // parse other stuff
    lex_t lex;

    struct {
	struct {
	    short id,channel;;
	} source,target;
    } edges[10000];
    int ne=0;

    enum {TOKSHADER=TOK_MAX, TOKBUFFER, TOKID, TOKNAME, TOKMASK, TOKEDGE, TOKSOURCE, TOKTARGET, TOKIN, TOKOUT};
    kw_t kw[]={{"shader", TOKSHADER}, 
	       {"buffer", TOKBUFFER}, 
	       {"id", TOKID},
	       {"name", TOKNAME},
	       {"mask", TOKMASK},
	       {"edge", TOKEDGE},
	       {"source", TOKSOURCE},
	       {"target", TOKTARGET},
	       {"in", TOKIN},
	       {"out", TOKOUT}};
    lex_bgn(&lex,code,kw,sizeof(kw)/sizeof(kw_t));
    lex_gettok();
    int jjj=0;
    while (lex.tok!=TOK_END) {
	switch (lex.tok) {
	case TOKIN:
	case TOKOUT:
	{
	    int type=lex.tok==TOKIN? FILTER_INPUT : FILTER_OUTPUT;
	    lex_gettok();
	    lex_expect('{');

	    int id=-1;
	    char nam[128]="";
	    while (lex.tok!='}') {
		switch (lex.tok) {
		case TOKID:
		    lex_gettok();
		    id=lex_getint();
		    break;
		case TOKNAME:
		    lex_gettok();
		    lex_getstr(nam, 127);
		    break;
		}
		//				lex_gettok();
	    }

	    if (id==-1 && !*nam) // ignore empty input/output nodes
		continue;

	    filter_addio(ft,nam,id,type);
	    lex_gettok(); // '}'
	    break;
	}
	case TOKBUFFER:
	    lex_gettok();
	    lex_expect('{');
	    while (lex.tok!='}') {
		switch (lex.tok) {
		case TOKID: 
		{
		    lex_gettok();
		    int i=lex_getint();
		    filter_allocnode(ft, i);
		    ft->nodes[i].indeg=0;
		    ft->nodes[i].type=FILTER_BUFFER;
		    break;
		}

		default:
		    lex_error();
		    break;
		}
		//				lex_gettok();
	    }
	    lex_gettok(); // '}'
	    break;

	case TOKSHADER:
	{
	    lex_gettok();
	    lex_expect('{');
	    int id=-1;
	    char nam[128]={[127]=0};
	    char mask[4]={1,1,1,1};
	    while (lex.tok!='}') {
		switch (lex.tok) {
		case TOKID: 
		    lex_gettok();
		    id=lex_getint();
		    break;

		case TOKNAME:
		    lex_gettok();
		    lex_getstr(nam, 127);
		    break;

		case TOKMASK:
		    lex_gettok();
		    for (int i=0; i<4; ++i) {
			mask[i]=lex_getint();
			if (mask[i]!=1) mask[i]=1;
		    }
		    break;
		default:
		    lex_error();
		    break;
		}
		//				lex_gettok();
	    }
	    lex_gettok(); // '}'

	    if (*nam && id!=-1) {
		int err=filter_addshader_simple(ft, nam, mask, id);
		if (err) {
		    lex_error();
		}
	    }
	    break;
	}
	case TOKEDGE:
	    lex_gettok();
	    lex_expect('{');
	    while (lex.tok!='}') {
		if (lex.tok==TOKSOURCE || lex.tok==TOKTARGET) {
		    int type=lex.tok;
		    lex_gettok();
		    int id=lex_getint();
		    lex_expect('[');
		    int channel=lex_getint();
		    lex_expect(']');

		    if (type==TOKSOURCE) {
			edges[ne].source.id=id;
			edges[ne].source.channel=channel;
		    }
		    else {
			edges[ne].target.id=id;
			edges[ne].target.channel=channel;
		    }
		}
		else {
		    lex_error();
		}
	    }
	    lex_gettok(); // '}'
	    ++ne;
	    break;
			
	case TOK_END:
	    break;
			
	default:
	    lex_error();
	    break;
	}
    }

    // link edges
    for (int i=0; i<ne; ++i) {
	int sc=edges[i].source.channel;
	int tc=edges[i].target.channel;			
	int sid=edges[i].source.id;
	int tid=edges[i].target.id;
	filternode_t *source=ft->nodes + sid;
	filternode_t *target=ft->nodes + tid;
		
	if (source->type==FILTER_SHADER) { // shader to buffer
	    if (sc>=arr_len(source->shader.out)) continue;
	    source->shader.out[sc].b=tid;
	}
	else {
	    assert(target->type==FILTER_SHADER); // buffer to shader
	    if (tc>=arr_len(target->shader.in)) continue;
	    target->shader.in[tc].b=sid;
	}
		
	source->inout[63-source->outdeg++]=tid;
	target->inout[target->indeg++]=sid;
    }

    // done!
    lex_end();

    // for shader with no output buffers, generate (invalid) output buffers for them
    int sz=arr_len(ft->nodes);
    for (int i=0; i<sz; ++i) {
	filternode_t *node=ft->nodes + i;
	if (node->indeg==-1 || node->type!=FILTER_SHADER ) continue;
	if (!node->shader.no) continue;

	for (int j=0; j<node->shader.no; ++j ){
	    if (node->shader.out[j].b!=-1) continue;

	    int id=filter_addnode(ft);
	    filternode_t *b=ft->nodes+id;
	    b->type=FILTER_BUFFER;
	    b->indeg=1;
	    b->outdeg=0;
	    b->inout[0]=i;
	    b->disp.x=0;
	    b->disp.y=0;
	    b->disp.rank=++ft->disp.maxrank;
			
	    filternode_t *s=ft->nodes + i;  // DEBUGGED: s=node => s=ft->nodes+i:  
	    // realloc is called, previous "filternode_t *node" may be invalid now
	    s->inout[63-s->outdeg++]=id;
	    //			arr_pushn(s->shader.out, NULL, 1);
	    s->shader.out[j].b=id;
	}
    }

    assert(!filter_update(ft));
    return ft;
}


static void filter_addshader(filter_t *ft, char *nam)
{
    if (!ft || !nam || !*nam) return;

    vfs_t *p=NULL;
    int ntex,nrts;
    int ds=!strcmp(nam, "ds"); // is the shader downsampler?
    int pack;
    if (!ds) {
	pack=!strcmp(nam, "pack");
    }

    // if the shader is downsampler, the shader and output buffer type are decided at runtime.
    if (ds) {
	// there're two versions of downsampler: color and depth -- they're two different shaders
	ntex=1;
	nrts=1;
    }
    else if (pack){ // can have any number of input buffers; only one output buffer
	nrts=1;
    }
    else {
	p=cvar_getp__(nam, VFS_CKSUM);
	if (!p) {
	    p=vfs_load(nam);
	    if (!p) { // error
		return;
	    }
	}
	ntex=p->ntex;
	nrts=p->nrts;
    }
	
    // let's find enough nodes to store them (shader node + children buffers)
    short id[nrts+1];
    int n=0;
    int sz=arr_len(ft->nodes);
    for (int i=0; i<sz && n<=nrts; ++i) {
	if (ft->nodes[i].indeg==-1) {
	    id[n++]=i;
	}
    }
    if (n<=nrts) {
	int k=nrts+1-n;
	arr_pushn(ft->nodes, NULL, k); // TODO: init ft->nodes
	for (int i=0; i<k; i++) {
	    id[nrts+i]=sz+i;
	}
    }
    n=nrts+1;
	
    // input (multitexture)
    filternode_t *s=ft->nodes + id[nrts];
    s->type=FILTER_SHADER;
    arr_popall(s->shader.in);
    arr_pushn(s->shader.in, NULL, ntex);

    if (s->shader.nam) free(s->shader.nam);
    s->shader.nam=strdup(nam);
    if (ds) {
	s->shader.in->b=-1;
	s->shader.type=SHADER_DS;
	s->shader.p_dscolor=vfs_load("shaders/colordownsampler.c");
	s->shader.p_dsdepth=vfs_load("shaders/depthdownsampler.c");
    }
    else if (pack) {
		
    }
    else {
	int ni=0;
	for (int i=0; i<p->nu; ++i) {
	    if (p->multitex[i]>=0) {
		// name is used purely for display pleasure (not necessary for the algorithms at all)
		// name, uniform, buffernode_id
		glGetActiveUniform(p->p, i, 64, 0,0,0, s->shader.in->nam);
		//				s->shader.in[ni].u=i;
		s->shader.in[ni++].b=-1; // WATCHOUT!
	    }
	}
    }
    s->shader.ni=ntex;
    s->shader.no=nrts;
    s->indeg=0;
    s->outdeg=nrts;

    // add output buffer nodes
    arr_popall(s->shader.out);
    arr_pushn(s->shader.out,NULL,nrts);
    for (int i=0; i<s->shader.no; ++i) {
	filternode_t *b=ft->nodes + id[i];
	b->type=FILTER_BUFFER;
	b->indeg=1;
	b->outdeg=0;
	b->buffer.type=RT_4B; // default is RGBA8
	b->buffer.tex=-1; // TODO
	// buffer texture etc. will be allocated during filter_uses() calls
    }
	
    if (!ds) {
	if (p->depthrt) {
	    filternode_t *b=ft->nodes + id[s->shader.no-1];
	    b->buffer.type=RT_DEPTH;
	}
    }
      
}


/*
  should call filter_update()
*/
static void filter_delshader(filter_t *ft, int id)
{
    filternode_t *s=ft->nodes+id;
    assert(s->type==FILTER_SHADER);

    // remove prev's next
    for (int i=0; i<s->indeg; ++i) {
	filternode_t *b=ft->nodes + s->inout[i];

	int j;
	for (j=0; j<b->outdeg; ++j) {
	    if (b->inout[63-j]==id) break;
	}

	for (;j<b->outdeg-1; ++j) {
	    b->inout[63-j]=b->inout[62-j];
	}

	--b->outdeg;
    }


    // remove next's prev
    for (int i=0; i<s->outdeg; ++i) {
	int c=s->inout[63-i];
	filternode_t *child=ft->nodes+c;
	for (int j=0; j<child->outdeg; ++j) {
	    filternode_t *grandchild=ft->nodes + child->inout[63-j];
	    int k;
	    for (k=0; k<grandchild->indeg; ++k) {
		if (grandchild->inout[k]==c) break;
	    }

	    for (;k<grandchild->indeg-1; ++k) {
		grandchild->inout[k]=grandchild->inout[k+1];
	    }

	    --grandchild->indeg;
	}

	// remove the shader node and its children -- just set their indeg to -1		
	child->indeg=-1;
    }

    s->indeg=-1;

    assert(!filter_update(ft));
}


// buffer->shader
static int filter_addedge(filter_t *ft, int bid, int sid, int channel)
{
    // if the edge already exists
    filternode_t *b=ft->nodes+bid;
    filternode_t *s=ft->nodes+sid;
    int i;
    for (i=0; i<b->outdeg; ++i) {
	if (b->inout[63-i]==sid) break;

	filternode_t *nexts=ft->nodes+b->inout[63-i];
	if (nexts->shader.type==SHADER_PACK) {
	    return 1;
	}
    }
    if (i!=b->outdeg) return 1;

    if (s->shader.type==SHADER_PACK) {
	if (b->outdeg) {
	    return 1;
	}
    }

    // TODO: test if there're too many in/out
	
    b->inout[63-b->outdeg++]=sid;
    s->inout[s->indeg++]=bid;

    // if cyclic graph
    if (filter_update(ft)) {
	--b->outdeg;
	--s->indeg;
	return 1;
    }

    s->shader.in[channel].b=bid;
    return 0;
}



/*
  should call filter_update()
*/

// buffer -X-> shader
static int filter_deledge(filter_t *ft, int bid, int sid)
{
    filternode_t *b=ft->nodes+bid;
    filternode_t *s=ft->nodes+sid;

    if (!(b->type==FILTER_IN || b->type==FILTER_BUFFER)) {
	assert(b->type==FILTER_IN || b->type==FILTER_BUFFER);
    }
    if (!(s->type==FILTER_SHADER || s->type==FILTER_OUT)) {
	assert(s->type==FILTER_SHADER || s->type==FILTER_OUT);
    }

    // remove buffer's next
    {
	int i;
	int found=0;
	for (i=0; i<b->outdeg; ++i) {
	    if (b->inout[63-i]==sid){ 
		found=1; break; }
	}
		
	if (!found) {
	    return 1;
	}
		
	for (;i<b->outdeg-1; ++i) {
	    b->inout[63-i]=b->inout[62-i];
	}
	--b->outdeg;
    }
	

    // remove shader's prev
    {
	int i;
	for (i=0; i<s->indeg; ++i) {
	    if (s->inout[i]==bid) break;
	}

	assert(i!=s->indeg);

	for (;i<s->indeg-1; ++i) {
	    s->inout[i]=b->inout[i+1];
	}
	--s->indeg;

	// remove the buffer from the shader's channel
	for (i=0; i<s->shader.ni; ++i) {
	    if (s->shader.in[i].b==bid){
		s->shader.in[i].b=-1;
		break;
	    }
	}
    }

    filter_update(ft);
    return 0;
}


static void filter_delio(filter_t *ft, int id)
{
    filternode_t *node=ft->nodes + id;
    if (node->type==FILTER_IN) {
	// remove next's prev
	for (int i=0; i<node->outdeg; ++i) {
	    int c=node->inout[63-i];
	    filternode_t *child=ft->nodes+c;
	    for (int j=0; j<child->outdeg; ++j) {
		filternode_t *grandchild=ft->nodes + child->inout[63-j];
		int k;
		for (k=0; k<grandchild->indeg; ++k) {
		    if (grandchild->inout[k]==c) break;
		}
				
		for (;k<grandchild->indeg-1; ++k) {
		    grandchild->inout[k]=grandchild->inout[k+1];
		}
				
		--grandchild->indeg;
	    }
			
	    // remove the shader node and its children -- just set their indeg to -1		
	    child->indeg=-1;
	}
    }
    else if (node->type==FILTER_OUT) {
	// remove prev's next
	for (int i=0; i<node->indeg; ++i) {
	    filternode_t *b=ft->nodes + node->inout[i];
			
	    int j;
	    for (j=0; j<b->outdeg; ++j) {
		if (b->inout[63-j]==id) break;
	    }
			
	    for (;j<b->outdeg-1; ++j) {
		b->inout[63-j]=b->inout[62-j];
	    }
			
	    --b->outdeg;
	}
    }
		
    node->indeg=-1;

    filter_update(ft);
	
}


static void filter_reload(char *nam)
{
    if (!nam || !*nam) return;

    filter_t *ft=cvar_getp__(nam, FILTER_CKSUM);
    if (!ft) return;

    time_t t=fmodtime(nam);
    if (t==ft->time) return;

    filter_load(nam, ft);
}
