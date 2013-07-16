#include "geometry.h"
#include "compiler.h"

// Articulate bodies

#include "math2.h"
#include "util.h"
#include "anim.h"
#include <stdlib.h>
#include <string.h>
#include <GL/gl.h>
#include <alloca.h>
#include "image.h"


typedef struct{
    char *name;
    int parent;
    v3_t pos;
    q_t rotq;
    float rotmat[3][3];
} joint_t;

typedef struct{
    unsigned widx;
    unsigned nw;
    float uv[2];
} vert_t;

typedef struct{
    unsigned v[3];
} tri_t;

typedef struct{
    unsigned bone;
    float w;
    v3_t pos;
    v3_t normal;
    v3_t tangent;
} weight_t;

typedef struct{
    joint_t *j;
    vert_t *v;
    tri_t *t;
    weight_t *w;
    unsigned nj, nv, nt, nw;

    unsigned vbo, ebo;
    unsigned colortex, normaltex;
} ab_t;



/*
  file format for articulate body:

  suffix: .ab

  abversion1

  j <num_joints(int)>
  <name(str)> <parent(int)> <pos(float*3)> <rot_quat(float*3)>

  v <num_verts(int)>
  <texcoord(float*2)> <weight_offset(int)> <num_weights(int)>

  t <num_tris(int)>
  <index_to_verts(int*3)>

  w <num_weights(int)>
  <bone_idx(int)> <weight(float)> <pos(float*3)>
  

  j,v,t,w can be in any order.

*/

enum{KW_ABVERSION, KW_J, KW_V, KW_T, KW_W};


static char *kwtab[]={"abversion", KW_ABVERSION, "j", KW_J, "v", KW_V, "t", KW_T, "w", KW_W};
static int nkws = sizeof(kwtab)/sizeof(char*)/2;

/* load joints */
static int loadj(text_t *t, ab_t *ab, int cnt)
{
    joint_t *j = ab->j = malloc(sizeof(joint_t)*cnt);
    ab->nj = cnt;
    char name[64];
    for(int i=0; i<cnt; i++){
	int err;
	err = readstr(t, name, 64); 
	if(err) return err;		

	j->name = strdup(name);
	err = readint(t, &j->parent);
	if(err) return err;

	err = readv3(t, &j->pos);
	if(err) return err;

	q_t q;
	err = readv3(t, &q.v);
	if(err) return err;

	float t = 1.0f - dot(q.v,q.v);
	q.w = (t<0.0f) ? 0.0f : -sqrtf(t);
	j->rotq = q;
	j++;
    }
    return 0;
}

/* load vertices */
static int loadv(text_t *t, ab_t *ab, int cnt)
{
    vert_t *v = ab->v = malloc(sizeof(joint_t)*cnt);
    ab->nv = cnt;
    for(int i=0; i<cnt; i++){
	int err;
	err = readfloat(t, v->uv);
	if(err) return err;

	err = readfloat(t, v->uv+1);
	if(err) return err;

	err = readint(t, &v->widx);
	if(err) return err;

	err = readint(t, &v->nw);
	if(err) return err;
	v++;
    }
    return 0;
}


/* load triangles */
static int loadt(text_t *t, ab_t *ab, int cnt)
{
    tri_t *tri = ab->t = malloc(sizeof(tri_t)*cnt);
    ab->nt = cnt;
    for(int i=0; i<cnt; i++){
	int err;
	err = readint(t, tri->v);
	if(err) return err;

	err = readint(t, tri->v+1);
	if(err) return err;

	err = readint(t, tri->v+2);
	if(err) return err;
	tri++;
    }
    return 0;
}

/* load weight vertices */
static int loadw(text_t *t, ab_t *ab, int cnt)
{
    weight_t *w = ab->w = malloc(sizeof(weight_t)*cnt);
    ab->nw = cnt;
    for(int i=0; i<cnt; i++){
	int err;
	err = readint(t, &w->bone);
	if(err) return err;

	err = readfloat(t, &w->w);
	if(err) return err;

	err = readv3(t, &w->pos);
	if(err) return err;
	w++;
    }
    return 0;
}

#define MAXNJ 1000
#define MAXNV 3000000
#define MAXNT 1000000
#define MAXNW 10000000

static int chkab(ab_t *ab)
{
    /* unreasonable num of stuff (at year 2008) */
    if(!ab->nj || !ab->nv || !ab->nt || !ab->nw ||
       ab->nj > MAXNJ || ab->nv > MAXNV || 
       ab->nt > MAXNT || ab->nw > MAXNW) return 1;


    /* check idx bounds */
    /* joint */
    int lastj = ab->nj-1;
    for(int i=0; i<ab->nj; i++){
	int p=ab->j[i].parent;
	if(p < -1 || p > lastj) return 1;
    }

    /* vertex */
    for(int i=0; i<ab->nv; i++)
	if(ab->v[i].widx + ab->v[i].nw > ab->nw) 
	    return 1;

    /* triangle */
    unsigned lastv = ab->nv-1;
    for(int i=0; i<ab->nt; i++){
	unsigned *v = ab->t[i].v;
	if(v[0] > lastv || v[1]>lastv || v[2]>lastv) 
	    return 1;
    }

    /* weight */
    for(int i=0; i<ab->nw; i++){
	int b = ab->w[i].bone;
	if(b >= ab->nj) return 1;
    }

    return 0;
}



void ab_free(ab_t *ab)
{
    if(!ab) return;

    if(ab->nj) free(ab->j);
    if(ab->nv) free(ab->v);
    if(ab->nt) free(ab->t);
    if(ab->nw) free(ab->w);
    free(ab);
}


static int readver(text_t *t, int *v)
{
    int kw;
    int err = readkw(t, kwtab, nkws, &kw);
    if(err) return err;

    if(kw != KW_ABVERSION) return 1;
    return readint(t, v);
}


static ab_t *loadab1(text_t *t)
{
    ab_t *ab = calloc(1, sizeof(ab_t));
    int err=0;

    while(!err){
	int type;
	err = readkw(t, kwtab, nkws, &type);
	if(err){
	    if(err = ERROR_EOF){
		err = 0;
	    }
	    break;
	}

	int cnt;
	err = readint(t, &cnt);
	if(err) break;

	switch(type){
	case KW_J:
	    err = loadj(t, ab, cnt);
	    break;

	case KW_V:
	    err = loadv(t, ab, cnt);
	    break;

	case KW_T:
	    err = loadt(t, ab, cnt);
	    break;

	case KW_W:
	    err = loadw(t, ab, cnt);
	    break;

	default:
	    err = 1;
	    break;
	}
    }

    if(err){
	ab_free(ab);
	return NULL;
    }
    else {
	return ab;
    }
}




/*
  compute weight point's normal & tangent

*/
static void compnt(ab_t *ab)
{
    /* compute object space vertex pos */
    v3_t *objv, *objn, *objs, *objt;
    void *data;
    {
	int size = 4*ab->nv*sizeof(v3_t);
	data = malloc(size);
	memset(data, 0, size);
	objv = data;
	objn = objv+ab->nv;
	objs = objn+ab->nv;
	objt = objs+ab->nv;
    }

    /* compute object space weight pos => add up to object vert pos */
    for(int i=0; i<ab->nv; i++){
	vert_t *v = ab->v + i;
	weight_t *w = ab->w + v->widx;
	for(int j=0; j<v->nw; j++){
	    joint_t *bone = ab->j+w->bone;
	    v3_t objw = add(bone->pos, rotbyq(bone->rotq, w->pos));
	    objv[i] = add(objv[i], scale(objw, w->w));
	    w++;
	}
    }

    /* compute per-triangle object space normal & tangents */
    for(int i=0; i<ab->nt; i++){
	/* normal */
	int *v = ab->t[i].v;
	v3_t e1 = sub(objv[v[1]], objv[v[0]]);
	v3_t e2 = sub(objv[v[2]], objv[v[0]]);
	v3_t n = normalize(cross(e1, e2));

	/* tangent */
	float *st0 = ab->v[v[0]].uv;
	float *st1 = ab->v[v[1]].uv;
	float *st2 = ab->v[v[2]].uv;
	float s1 = st1[0] - st0[0];
	float s2 = st2[0] - st0[0];
	float t1 = st1[1] - st0[1];
	float t2 = st2[1] - st0[1];
	float det = s1*t2 - s2*t1;
	float invdet = det ? 1.0f/det : 1.0f;
	v3_t s = normalize(sub(scale(e1, t2*invdet), scale(e2, t1*invdet)));
	v3_t t = normalize(sub(scale(e2, s1*invdet), scale(e1, s2*invdet)));
		
	for(int j=0; j<3; j++){
	    int k=v[j];
	    objn[k] = add(objn[k], n);
	    objs[k] = add(objs[k], s);
	    objt[k] = add(objt[k], t);
	}

    }


    /* per-vertex object space normal => per-weight-vertex bone space normal */
    for(int i=0; i<ab->nv; i++){
	v3_t normal = normalize(objn[i]);
	v3_t tangent = normalize(sub(objs[i], scale(normal, proj(objs[i], normal))));
	if(proj(cross(normal, tangent), objt[i]) < 0.0f){
	    tangent = scale(tangent, -1.0f);
	}

	vert_t *v = ab->v + i;
	weight_t *w = ab->w + v->widx;
	for(int j=0; j<v->nw; j++){
	    q_t invrotq = conjq(ab->j[w->bone].rotq);
	    w->normal = rotbyq(invrotq, normal); /* no need to normalize -- rotation doesn't change vector length */
	    w->tangent = rotbyq(invrotq, tangent);
	    w++;
	}
    }

    free(data);

}



/*
  check (get) the version number of ab.

  correct version header should be:
  abversionN

  N is a int.
*/


ab_t *ab_loadfile(char *filename)
{
    text_t *t = loadtext(filename);
    if(!t) return NULL;

    int ver;
    int err = readver(t, &ver);
    if(err){
	freetext(t);
	return NULL;
    }

    ab_t *ab=NULL;
    switch(ver){
    case 1:
	ab = loadab1(t);
	if(ab){
	    err = chkab(ab);
	    if(err){
		ab_free(ab);
		ab = NULL;
	    }			
	}
	break;
    }

    freetext(t);
    return ab;
}


/*
  For each AB there is a vbo.
  pos, texcoord, normal, tangent, index

  data is per-vertex so they can be interleaved.

*/
typedef struct{
    v3_t pos;
    v3_t normal;
    v3_t tangent;
    float uv[2];
    float pad[2];
} rendvert_t;


void ab_updvbo(ab_t *ab)
{	

    int size=ab->nv * sizeof(rendvert_t);
    if(!ab->vbo){
	glGenBuffers(1, &ab->vbo);
	glGenBuffers(1, &ab->ebo);
	/* index buffer -- static */
	glBindBuffer(GL_ARRAY_BUFFER, ab->vbo);
	glBufferData(GL_ARRAY_BUFFER, size, NULL, GL_STREAM_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ab->ebo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, ab->nt*sizeof(tri_t), ab->t, GL_STATIC_DRAW);
    }

    rendvert_t *rendv = alloca(size);
    memset(rendv, 0, size);
    for(int i=0; i<ab->nv; i++){
	vert_t *v = ab->v+i;
	rendvert_t *rv = rendv+i;
	weight_t *w = ab->w + v->widx;
	for(int j=0; j<v->nw; j++){
	    joint_t *bone=ab->j + w->bone;
	    {
		v3_t wpos = scale(add(bone->pos, rotbyq(bone->rotq, w->pos)), w->w);
		rv->pos = add(rv->pos, wpos);
	    }

	    {
		v3_t wn = scale(rotbyq(bone->rotq, w->normal), w->w);
		rv->normal = add(rv->normal, wn);
	    }

	    {
		v3_t wt = scale(rotbyq(bone->rotq, w->tangent), w->w);
		rv->tangent = add(rv->tangent, wt);
	    }
	    w++;
	}
		
	rv->uv[0] = v->uv[0]; rv->uv[1] = v->uv[1];
    }

    glBindBuffer(GL_ARRAY_BUFFER, ab->vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, size, rendv, GL_STREAM_DRAW);


}


void ab_render(ab_t *ab, unsigned loc)
{
    glBindBuffer(GL_ARRAY_BUFFER, ab->vbo);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    glVertexPointer(3, GL_FLOAT, sizeof(rendvert_t), 0);

    printf("normalptr=%i\n", &(((rendvert_t*)NULL)->normal));
    glNormalPointer(GL_FLOAT, sizeof(rendvert_t), &(((rendvert_t*)NULL)->normal));

    glEnableVertexAttribArray(loc);
    printf("tangentptr=%i\n", &(((rendvert_t*)NULL)->tangent));
    glVertexAttribPointer(loc, 3, GL_FLOAT, 0, sizeof(rendvert_t), &(((rendvert_t*)NULL)->tangent));

    printf("uvptr=%i\n", ((rendvert_t*)NULL)->uv);
    glTexCoordPointer(2, GL_FLOAT, sizeof(rendvert_t), ((rendvert_t*)NULL)->uv);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ab->ebo);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, ab->colortex);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, ab->normaltex);

    glDrawRangeElements(GL_TRIANGLES, 0, ab->nv-1, ab->nt*3, GL_UNSIGNED_INT, NULL);

    printf("render!\n");

}


ab_t *g_ab;
void test_initab()
{
    g_ab = ab_loadfile("hello2.ab");

    if(g_ab){
//		printf("nj=%i, nv=%i, nt=%i, nw=%i\n", ab->nj, ab->nv, ab->nt, ab->nw);
    }

    compnt(g_ab);
    ab_updvbo(g_ab);
    g_ab->colortex = loadtexture("Bandersha_512.tga", TEX_MIPMAP, 1);
    g_ab->normaltex = loadtexture("Bandersha_512_NRM.tga", TEX_MIPMAP, 1);
    return 0;
}

void test_ab_render(unsigned loc)
{
    if(g_ab) ab_render(g_ab, loc);
}
