#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <GL/gl.h>

#include "math2.h"
#include "meh.h"
#include "camera.h"
#include "cd.h"

static struct{
	sp_t			*sp;
	sctr_t			*sctrlist;
	int			nsctrs;
	obj_t			*objlist;
} map;

/* engine variables are stored here, and can be accessed by external functions by including the meh.h file. */

static unsigned short meh_touchers[1024];
float meh_wallthickness = 0.2f;


/* a bound is invalid if it collides with any sector */
#define BND_TOO_TINY 0.1f
static inline int bndok(float bnd[2][3])
{
	return (bnd[1][0]-bnd[0][0]>BND_TOO_TINY ||
		bnd[1][1]-bnd[0][1]>BND_TOO_TINY ||
		bnd[1][2]-bnd[0][2]>BND_TOO_TINY);
}



static inline int bndbnd(const float a[2][3], const float b[2][3])
{
	return (a[0][0]<b[1][0] && a[1][0]>b[0][0] && 
		a[0][1]<b[1][1] && a[1][1]>b[0][1] && 
		a[0][2]<b[1][2] && a[1][2]>b[0][2]);
}

#if 0
/* tests whether a bound touches a sector */
static int bndsctr(const float bnd[2][3], const sctr_t*sctr)
{
	if(!bndbnd(bnd,sctr->bnd)) return 0;

}
#endif


/*
  tries to add a new sector to the current map.
 */
#if 0
int bnd2sctr(float bnd[2][3])
{
	if(!bndok(bnd)) return 0;
	int n = rngqry(map.sp, bnd, meh_touchers);
	for(int i=0;i<n;i++){
		sctr_t *sctr = getowner(map.sp, meh_touchers[i]);
		if(bndsctr(bnd,sctr)) return 0;
	}

	sctr_t *sctr = malloc(sizeof(sctr_t));

	if(map.sctrlist){
		map.sctrlist->prev = sctr;
	}
	sctr->next = map.sctrlist;
	sctr->prev = NULL;
	map.sctrlist = sctr;

	/* find portals */
}
#endif

/* are p and p2 separate? */
static inline int plnsep(plane_t *p, plane_t *p2)
{
	if(dot(p->n,p2->n)<-0.99){
		float d = p->d + p2->d;
		if(d<0 && d>-0.01){
			return 0;
		}
	}
	return 1;

#if 0
	if(fabs(p->n.x-p2->n.x)<0.01 &&
	   fabs(p->n.y-p2->n.y)<0.01 &&
	   fabs(p->n.z-p2->n.z)<0.01){
		float d = p->d + p2->d;
		if(d<0 && d>-0.01){
			return 0;
		}
	}
	return 1;
#endif
}


/*
  Editing sectors

  constrained/unconstrained add/move/rotate sectors;
  constrained/unconstrained move/rotate walls

 */

/*
  whether two __coplannar__ polygons intersct
 */
static int polyintx(phd_t *phd, int pid, phd_t *phd2, int pid2, v3_t intx[])
{
	poly_t *p,*p2;
	p = phd->polys + pid;
	p2 = phd2->polys + pid2;
	for(int i=0;i<p->nv;i++){
		intx[i] = phd->verts[p->v[i]];
	}
	int nv = p->nv;
	for(int i=0;i<p2->nv;i++){
		plane_t pln;
		v3_t *v,*v2;
		int i2 = (i+1)%p2->nv;
		v = phd2->verts + p2->v[i];
		v2 = phd2->verts + p2->v[i2];
		pln.n = cross(p2->plane.n, sub(*v2,*v));
		pln.d = dot(pln.n, *v);
		nv = splitpoly(intx, nv, &pln);
		if(!nv) return 0;
	}
}


/*
  portal generation
  
  adding/removing brushes may affect portals

  how to store portals, and efficiently update them?



 */


static void add2mptllist(metaptl_t *mp, metaptl_t *twin, sctr_t *sctr, int wallid, v3_t *verts, int ccw, int nverts)
{
	mp->prev = NULL;
	mp->next = sctr->walls[wallid].mptllist;
	if(mp->next){
		mp->next->prev = mp;
	}
	sctr->walls[wallid].mptllist = mp;
	mp->twin = twin;
	mp->verts = verts;
	mp->nverts = nverts;
	mp->sctr = sctr;
	mp->wallid = wallid;
	mp->ccw = ccw;
}

static void addsctr2map(sctr_t *sctr)
{
	sctr->prev = map.sctrlist;
	if(sctr->prev){
		sctr->next = sctr->prev->next;
		sctr->prev->next = sctr;
	}
	else{
		map.sctrlist = sctr;
		sctr->next = NULL;
	}
	map.nsctrs++;
}


static int genmptls(sctr_t *sctr)
{
	phd_t *phd = sctr->phd;
	float loosebnd[2][3];
	for(int a=0;a<3;a++){
		loosebnd[0][a] = phd->bnd[0][a]-0.01;
		loosebnd[1][a] = phd->bnd[1][a]+0.01;
	}
	int ntouchers = rngqry(map.sp, loosebnd, meh_touchers);

	int seppoly[32];
	for(int i=0;i<ntouchers;i++){
		sctr_t *s=SP_GETOWNER(map.sp,meh_touchers[i]);
		if(s==sctr) continue;
		seppoly[i] = phdqry(phd, s->phd);
		if(seppoly[i]<0) return -1;
	}

	/* find metaportals */
	v3_t ptlverts[32];
	int nmptls=0;
	for(int t=0;t<ntouchers;t++){
		sctr_t *s=SP_GETOWNER(map.sp,meh_touchers[t]);
		if(s==sctr) continue;
		if(seppoly[t]>=phd->npolys) continue;

		int pid = -1;
		for(int j=seppoly[t];j<phd->npolys;j++){
			if(phdabovepln(s->phd, &phd->polys[j].plane)){
				pid = j;
				break;
			}
		}
		assert(pid>=0);

		plane_t *p = &phd->polys[pid].plane;
		for(int j=0; j<s->phd->npolys; j++){
			plane_t *p2 = &s->phd->polys[j].plane;
			if(plnsep(p,p2)) continue;

			int nptlverts = polyintx(phd, pid, s->phd, j, ptlverts);
			if(nptlverts){
				printf("sctr=%x\n");
				void *mem = malloc(2*sizeof(metaptl_t) + 32*sizeof(v3_t)); /* mp,mp2 allocated using a single malloc */
				metaptl_t *mp = mem;
				metaptl_t *mp2 = mp+1;
				v3_t *verts = mp2 + 1;
				memcpy(verts, ptlverts, nptlverts*sizeof(v3_t));
				add2mptllist(mp,mp2,sctr,pid,verts,0,nptlverts);
				add2mptllist(mp2,mp,s,j,verts,1,nptlverts);
				nmptls++;
				break;
			}
		}
	}
	return nmptls;
}


int phd2sctr(phd_t *phd)
{
	sctr_t *sctr;
	int nptls;
	sctr = malloc(sizeof(sctr_t));
	sctr->phd = phd;
	memset(sctr->walls, 0, phd->npolys*sizeof(wall_t));
	nptls = genmptls(sctr);
	if(nptls<0){
		free(sctr);
		return 1;
	}

	printf("hello\n");

	/* the phd is valid */
	sctr->sp = allocsp(128);
	sctr->brushlist = NULL;
	addbnd(map.sp, phd->bnd, 1, sctr, &sctr->spid);
	flushsp(map.sp);
	addsctr2map(sctr);

	for(int j=0;j<phd->nedges;j++){
		printf("edge->p[0]=%i, edge->p[1]=%i\n", phd->edges[j].p[0], phd->edges[j].p[1]);
	}


	/* generate wall brushes */
	sctr->nwalls = phd->npolys;
	plane_t wallplns[32];
	for(int i=0;i<sctr->nwalls;i++){
		/* extrude polys inward */
		brush_t *b = malloc(sizeof(brush_t));
		b->type = BRUSH_TYPE;
		b->prev = b->next = NULL;
		b->wallid = i;
		b->phd = poly2phd(phd, i, meh_wallthickness);
		sctr->walls[i].brushlist = b;
		addbnd(sctr->sp, b->phd->bnd, 0, b, &b->spid);
	}
	flushsp(sctr->sp);
	return 0;
}


int bnd2sctr(float bnd[2][3])
{
	phd_t *phd = allocphd();
	plane_t planes[6];
	bnd2planes(bnd,planes);
	addplns(phd,planes,6);
	phd2sctr(phd);
}

/* adding brush, with collision detection */
int addbrush(sctr_t *sctr, plane_t plns[], int np, int cd)
{
	int stuck=0;
	phd_t *phd = allocphd();
	addplns(phd, plns, np);

	/* is the brush inside the sctr? */
	{
		int a=phdqry(phd, sctr->phd);
		if(a>=0){
			freephd(phd);
			return 1;
		}
	}

	/* detection collisions.
	   if collide with other stuff, add nothing and return "stuck"  */
	int n = rngqry(sctr->sp, phd->bnd, meh_touchers);
	if(cd){
		for(int i=0; i<n; i++){
			stuff_t *stuff = SP_GETOWNER(sctr->sp, meh_touchers[i]);
			if(stuff->type != OBJ_TYPE) continue;
			
			int a=phdqry(phd, stuff->common.phd);
			if(a>=0){
				stuck=1;
				break;
			}
		}
		
		if(stuck){
			freephd(phd);
			return 1;
		}
	}
	else{
		for(int i=0; i<n; i++){
			stuff_t *stuff = SP_GETOWNER(sctr->sp, meh_touchers[i]);
			if(stuff->type != OBJ_TYPE) continue;

			int a=phdqry(phd, stuff->common.phd);
			if(a>=0){
				obj_t *obj = stuff;
				obj->flag |= OBJ_PHANTOM;
			}
		}
	}

	/* TODO: clip brush against sctr & portal stuff */

	brush_t *b = malloc(sizeof(brush_t));
	b->prev = NULL;
	b->next = sctr->brushlist;
	if(b->next) { b->next->prev = b;}
	sctr->brushlist = b;

	b->type = BRUSH_TYPE;
	b->phd = phd;
	b->wallid = -1;
	addbnd(sctr->sp, b->phd->bnd, 0, b, &b->spid);
	flushsp(sctr->sp);

	return 0;
}

void rmbrush(sctr_t *sctr, brush_t *b)
{
	if(b->prev){
		b->prev->next = b->next;
	}
	else{
		if(b->wallid>=0){
			sctr->walls[b->wallid].brushlist = b->next;
		}
		else{
			sctr->brushlist = b->next;
		}
	}

	if(b->next){
		b->next->prev = b->prev;
	}
	
	freephd(b->phd);
	rmbnd(sctr->sp, b->spid);
	flushsp(sctr->sp);
	free(b);
}

static void rmbrushlist_bat(sctr_t *sctr, brush_t *blist)
{
	brush_t *next;
	for(brush_t *b=blist; b; b=next){
		next = b->next;
		freephd(b->phd);
		rmbnd(sctr->sp, b->spid);
		free(b);
	}
}


/* should be used only by rmsctr(): no rmbnd() stuff */
static void rmsctr_rmbrushlist(brush_t *blist)
{
	brush_t *next;
	for(brush_t *b=blist; b; b=next){
		next = b->next;
		freephd(b->phd);
		free(b);
	}
}

static void rmmptllist(metaptl_t *mlist)
{
	printf("mlist=0x%x\n",mlist);
	metaptl_t *next;
	if(mlist && !mlist->prev){
		mlist->sctr->walls[mlist->wallid].mptllist = NULL;
	}

	for(metaptl_t *m=mlist; m; m=next){
		metaptl_t *t = m->twin;
		next = m->next;
		if(t->prev){
			t->prev->next = t->next;
		}
		else{
			printf("t->wallid=%i, t->next=%i\n",t->wallid,t->next);
			t->sctr->walls[t->wallid].mptllist = t->next;
		}
		if(t->next){
			t->next->prev = t->prev;
		}
		if(m->ccw)free(t);
		else free(m);
	}
}

void rmsctr(sctr_t *s)
{
	if(s->prev){
		s->prev->next = s->next;
	}
	else{
		map.sctrlist = s->next;
	}
	if(s->next){
		s->next->prev = s->prev;
	}
	freephd(s->phd);
	rmbnd(map.sp, s->spid);	
	flushsp(map.sp);
	freesp(s->sp);

	for(int i=0;i<s->nwalls;i++){
		rmsctr_rmbrushlist(s->walls[i].brushlist);
		rmmptllist(s->walls[i].mptllist);
	}
	rmsctr_rmbrushlist(s->brushlist);
	free(s);
}

static void rmallsctrs()
{
	sctr_t *next;
	for(sctr_t *s=map.sctrlist; s; s=next){
		next = s->next;
		freephd(s->phd);
		freesp(s->sp);

		for(int i=0;i<s->nwalls;i++){
			rmsctr_rmbrushlist(s->walls[i].brushlist);
			rmmptllist(s->walls[i].mptllist);
		}
		rmsctr_rmbrushlist(s->brushlist);
		free(s);
	}
}


/* 
   set the wall's normal.
   returnsfalse and keeps the sector unchanged if the sector with the new normal collides with other sectors (static cd)
*/
int setwalln(sctr_t *s, wall_t *w, const v3_t *n)
{
}

void poly2bnd(phd_t *phd, int pid, float bnd[2][3])
{
	bnd[0][0]=bnd[0][1]=bnd[0][2]=1.0e20f;
	bnd[1][0]=bnd[1][1]=bnd[1][2]=1.0e-20f;
	poly_t *p=phd->polys+pid;
	for(int i=0; i<p->nv; i++){
		v3_t *v = phd->verts + p->v[i];
		if(v->x<bnd[0][0]) bnd[0][0]=v->x;
		else if(v->x>bnd[1][0]) bnd[1][0]=v->x;
		if(v->y<bnd[0][1]) bnd[0][1]=v->y;
		else if(v->y>bnd[1][1]) bnd[1][1]=v->y;
		if(v->z<bnd[0][2]) bnd[0][2]=v->z;
		else if(v->z>bnd[1][2]) bnd[1][2]=v->z;
	}
}


/*
  1. find affected brushes
  2. extrude them to infinite
  3. clip them against borders of the wall
 */


int rotwall2d(sctr_t *sctr, int wallid, float deg)
{

}

int xlatwall2d(sctr_t *sctr, int wallid, v3_t xlat)
{

}


v3_t polycntr(phd_t *phd, int pid)
{
	v3_t cntr = {0,0,0};
	poly_t *poly = phd->polys + pid;
	for(int i=0; i<poly->nv; i++){ 
		cntr = add(cntr, phd->verts[poly->v[i]]); 
	}
	return scale(cntr, 1.0f/poly->nv);
}


#define SAME_VEC(v,v2) (dot((v),(v2))>0.99)
#define Y_VEC(p) (fabs((p)->y>0.999))

#define PLN_XLAT 1
#define PLN_ROT 2
int getplntransf(phd_t *oldphd, phd_t *newphd, int pid, v3_t rightvec[2], v3_t *xlat, float rot[3][3])
{
	int transf=0;
	v3_t oldcenter={0,0,0}, newcenter={0,0,0};
	plane_t *oldpln, *newpln;

	oldpln = &oldphd->polys[pid].plane;
	newpln = &newphd->polys[pid].plane;
	/* extract the transformation of the plane */

	/* 1. translation */
	if(SAME_PLN(oldpln, newpln)){
		return 0;
	}
	else{
		newcenter = polycntr(newphd, pid);
		oldcenter = polycntr(oldphd, pid);
		*xlat = sub(newcenter, oldcenter);
		if(dot(*xlat, *xlat)>0.0001){transf |= PLN_XLAT; };
	}

	/* 2. rotation	*/
	if(!SAME_VEC(oldpln->n, newpln->n)){
		v3_t y={0,1,0};
		v3_t right[2], up[2], depth[2];
		v3_t x;

		transf |= PLN_ROT;
		depth[0] = oldpln->n;
		depth[1] = newpln->n;
		x = Y_VEC(depth) ? rightvec[depth[0].y>0] : normalize(cross(depth[0],y));
		up[0] = cross(x, depth[0]);
		right[0] = cross(depth[0], up[0]);

		x = Y_VEC(depth+1) ? rightvec[depth[1].y>0] = right[0] : normalize(cross(depth[1], y));
		up[1] = cross(x, depth[1]);
		right[1] = cross(depth[1], up[1]);

		if(Y_VEC(depth)){
			float proj[2];
			proj[0] = dot(right[0], right[1]);
			proj[1] = dot(up[0], right[1]);
			if(fabs(proj[1]) > fabs(proj[0])){ /* swap right and up vectors */
				if(proj[0]<0){ up[0] = scale(right[0], -1.0f); }
			}
			else{
				if(proj[1]<0){ up[0] = scale(up[0], -1.0f);}
			}
			right[0] = cross(depth[0], up[0]);
		}

		/* rot = R_new * R_old^-1 = R_new * R_old^T */
		rot[0][0] = right[1].x*right[0].x + up[1].x*up[0].x + depth[1].x*depth[0].x;
		rot[0][1] = right[1].x*right[0].y + up[1].x*up[0].y + depth[1].x*depth[0].y;
		rot[0][2] = right[1].x*right[0].z + up[1].x*up[0].z + depth[1].x*depth[0].z;
		rot[1][0] = right[1].y*right[0].x + up[1].y*up[0].x + depth[1].y*depth[0].x;
		rot[1][1] = right[1].y*right[0].y + up[1].y*up[0].y + depth[1].y*depth[0].y;
		rot[1][2] = right[1].y*right[0].z + up[1].y*up[0].z + depth[1].y*depth[0].z;
		rot[2][0] = right[1].z*right[0].x + up[1].z*up[0].x + depth[1].z*depth[0].x;
		rot[2][1] = right[1].z*right[0].y + up[1].z*up[0].y + depth[1].z*depth[0].y;
		rot[2][2] = right[1].z*right[0].z + up[1].z*up[0].z + depth[1].z*depth[0].z;
		*xlat = sub(newcenter, matv3(rot, oldcenter));
	}
	return transf;
}


static phd_t *replacepln(phd_t *phd, plane_t *oldpln, int pid)
{
	phd_t *rphd;
	plane_t *newpln = &phd->polys[pid].plane;
	if(SAME_PLN(oldpln, newpln)) return NULL;

	/* reconstruct the old sector */
	plane_t plns[32];
	for(int i=0; i<phd->npolys; i++){
		plns[i] = phd->polys[i].plane;
	}
	plns[pid] = *oldpln;
	rphd = allocphd();
	addplns(rphd, plns, phd->npolys);
	return rphd;
}

static int poly2borderplns(phd_t *phd, int pid, plane_t borderplns[])
{
	poly_t *poly = phd->polys + pid;
	for(int i=0; i<poly->nv; i++){
		edge_t *e = phd->edges + poly->e[i];
		borderplns[i] = phd->polys[e->p[e->p[0]==pid]].plane;
	}
	return poly->nv;
}


static void mvborderplns_updbrush(sctr_t *sctr, brush_t *b, int intx)
{
	if(intx == PHD_ABOVE_PLN){
		b->type = NULL_BRUSH_TYPE;
		rmbnd(sctr->sp, b->spid);
	}
	else{
		if(b->type == NULL_BRUSH_TYPE){
			b->type = BRUSH_TYPE;
			addbnd(sctr->sp, b->phd->bnd, 0, b, &b->spid);
		}
		else{
			updbnd(sctr->sp, b->spid, b->phd->bnd);
		}
	}
}

/* don't expect this function to run in real-time when the sctr is very complex */
/* can be used only when no walls have been added/removed */
int updsctr(sctr_t *sctr, phd_t *oldsctrphd, int wallid)
{
	plane_t *oldpln = &oldsctrphd->polys[wallid].plane;
	plane_t *newpln = &sctr->phd->polys[wallid].plane;

	/* for each brush in the neighbouring walls:
	     1. find & remove the old border pln
	     2. clipped by the new border pln
	 */
	/* first we'd find which walls are going to be discarded/clipped/kept untouched */
	char above[64];
	memset(above,0,oldsctrphd->nverts);
	for(int i=0; i<oldsctrphd->nverts; i++){
		above[i] = dot(oldsctrphd->verts[i], newpln->n) > newpln->d;
	}

	for(int i=0; i<oldsctrphd->npolys; i++){
		if(i==wallid) continue;
		poly_t *wallpoly = oldsctrphd->polys + i;
		int nf=0;
		for(int j=0; j<wallpoly->nv; j++){ nf += above[wallpoly->v[j]];}
		if(nf>0 && nf<wallpoly->nv){ /* clip */
			for(brush_t *b=sctr->walls[i].brushlist; b; b=b->next){
				int intx = mvborderpln(b->phd, oldpln, newpln);
				mvborderplns_updbrush(sctr, b, intx);
			}
		}
		else if(nf == wallpoly->nv){/* the wall is totally behind some other walls, hence should be discarded */
			rmbrushlist_bat(sctr, sctr->walls[i].brushlist);
		}
		/* else: no clipping is needed */
	}

	/* extrude the wall's brushes */
	if(sctr->phd->polys[wallid].nv < 3){
		rmbrushlist_bat(sctr, sctr->walls[wallid].brushlist);
	}
	else{
		if(!sctr->walls[wallid].brushlist){
			brush_t *b = malloc(sizeof(brush_t));
			b->type = BRUSH_TYPE;
			b->phd = poly2phd(sctr->phd, wallid, meh_wallthickness);
			b->wallid = wallid;
			addbnd(sctr->sp, b->phd->bnd, 0, b, &b->spid);
		}
		else{

			/* for each brush in the moved wall:
			   1. find & remove old border plns 
			   2. translate/rotate
			   3. clipped by new border plns
			*/
			plane_t borderplns[16], newborderplns[16];
			int nb,nnb;
			nb = poly2borderplns(oldsctrphd, wallid, borderplns);
			nnb = poly2borderplns(sctr->phd, wallid, newborderplns);

			v3_t xlat;
			float rotmat[3][3];
			int transf;
			transf = getplntransf(oldsctrphd, sctr->phd, wallid, sctr->rightvec, &xlat, rotmat);
			void *rot = (transf & PLN_ROT) ? rotmat : NULL;
			for(brush_t *b=sctr->walls[wallid].brushlist; b; b=b->next){
				int intx=mvborderplns_n_phd(b->phd, rot, xlat, borderplns, newborderplns, nb, nnb);
				mvborderplns_updbrush(sctr, b, intx);
			}
		}
	}

	/* clean up */
	flushsp(sctr->sp);
}

int updsctr_lock(sctr_t *sctr, phd_t *oldsctrphd, float rot[3][3], v3_t xlat, int lockids[], int nlocks)
{
	char locked[32];
	memset(locked, 0, sctr->phd->npolys);
	for(int i=0; i<nlocks; i++){ locked[i]=1;}

	for(int i=0; i<sctr->phd->npolys; i++){
		if(locked[i]){
			v3_t oldcntr, newcntr;
			plane_t borderplns[16], newborderplns[16];
			int nb,nnb;
			oldcntr = polycntr(oldsctrphd, i);
			newcntr = polycntr(sctr->phd, i);

			nb = poly2borderplns(oldsctrphd, i, borderplns);
			nnb = poly2borderplns(sctr->phd, i, newborderplns); /* the new phd is already transformed */
			for(brush_t *b=sctr->walls[i].brushlist; b; b=b->next){
				int intx = mvborderplns_n_phd(sctr->phd, NULL, sub(newcntr, oldcntr), 
							      borderplns, newborderplns, nb, nnb);
				mvborderplns_updbrush(sctr, b, intx);
			}
		}
		else{
			poly_t *oldpoly, *newpoly;
			plane_t borderplns[16], newborderplns[16];
			int nb=0,nnb=0;
			oldpoly = oldsctrphd->polys + i;
			newpoly = sctr->phd->polys + i;

			for(int j=0; j<oldpoly->nv; j++){
				edge_t *e = oldsctrphd->edges + oldpoly->e[j];
				int adj = e->p[e->p[0]==i];
				if(locked[adj]){
					borderplns[nb++] = oldsctrphd->polys[adj].plane;
				}
			}

			for(int j=0; j<newpoly->nv; j++){
				edge_t *e = sctr->phd->edges + newpoly->e[j];
				int adj = e->p[e->p[0]==i];
				if(locked[adj]){
					newborderplns[nnb++] = sctr->phd->polys[adj].plane;
				}
			}

			for(brush_t *b=sctr->walls[i].brushlist; b; b=b->next){
				int intx = mvborderplns_n_phd(sctr->phd, rot, xlat,borderplns, newborderplns, nb, nnb);
				mvborderplns_updbrush(sctr, b, intx);
			}
		}
	}
	flushsp(sctr->sp);
}


int addwalls(sctr_t *sctr, plane_t wallplns[], int nwalls)
{
	int np = 0;
	plane_t plns[32];
	int oldnp = sctr->phd->npolys;

	/* check whether these added walls already exist */
	for(int i=0; i<sctr->nwalls; i++){
		for(int j=0; j<nwalls; j++){
			if(!SAME_PLN(&sctr->phd->polys[i].plane, wallplns+j)){
				plns[np++] = wallplns[j];
			}
		}
	}

	if(!np) return 0;

	addplns(sctr->phd, plns, np);
	for(int i=0; i<sctr->nwalls; i++){
		if(sctr->phd->polys[i].nv<3){
			sctr->phd->polys[i].nv = 0;
			rmbrushlist_bat(sctr, sctr->walls[i].brushlist);
		}
	}

	sctr->nwalls = sctr->phd->npolys;
	for(int i=0; i<oldnp; i++){
		poly_t *poly = sctr->phd->polys + i;
		plane_t newplns[32];
		int newnp=0;
		int nf=0;
		for(int j=0; j<poly->nv; j++){
			edge_t *e = sctr->phd->edges + poly->e[j];
			int adj = e->p[e->p[0]==i];
			if(adj>=oldnp){
				newplns[newnp++] = sctr->phd->polys[adj].plane;
			}
		}

		for(brush_t *b=sctr->walls[i].brushlist; b; b=b->next){
			int intx = addborderplns(b->phd, newplns, newnp);
			mvborderplns_updbrush(sctr, b, intx);
		}
	}

	for(int i=oldnp; i<oldnp+np; i++){
		brush_t *b = malloc(sizeof(brush_t));
		b->type = BRUSH_TYPE;
		b->prev = b->next = NULL;
		b->wallid = i;
		b->phd = poly2phd(sctr->phd, i, meh_wallthickness);
		sctr->walls[i].brushlist = b;
		addbnd(sctr->sp, b->phd->bnd, 0, b, &b->spid);		
	}

	flushsp(sctr->sp);
}


int rmdeadwalls(sctr_t *sctr)
{
}


/*
  1. xlat sctr (only xlat)
  2. rot sctr (only rot)
  3. xlat walls 
  ((1)xlat + change border plns -- xlat wall; 
  (2)only change border plns -- locked walls; 
  (3)change nothing at all -- locked walls)
  4. rot walls ((1)rot + change border plns; (2)only change border plns; (3)change nothing at all)
  5. xlat sctr while some walls are locked 
  ((1)xlat + change border plns -- xlat walls; 
  (2)only change border plns -- locked walls; 
  (3)only xlat -- xlat walls)
  6. rot sctr while some walls are locked ((1)xlat + change border plns; (2)only change border plns; (3)only rot)

  3 == xlat walls relative to sctr (lock sctr (==some walls), xlat walls) == lock some walls, xlat other walls
  5 == xlat sctr relative to walls (xlat sctr (==some walls), lock others)== lock some walls, xlat other walls
  4 == rot walls relative to sctr (lock walls, xlat sctr (==some walls)) == lock some walls, rot other walls
  6 == rot sctr relative to walls ========================================= lock some walls, rot other walls

  ==> xlat/rot some walls


  7. add new walls ((1)only change border plns (2)change nothing)
  8. remove walls ((1)only change border plns (2)change nothing)
  9. split sctr ((1)only change border plns (2)change nothing)
  10. xlat wall 2d ((1)only change border plns (2)change nothing)
  11. rot wall 2d ((1)only change border plns (2)change nothing)

  all operations above use ccd to avoid penetration
 */
static int updwall(sctr_t *sctr, phd_t *oldphd, int wallid)
{
}

/* removing a wall == move it far away, instead of really deleting it
   adding a wall == really add it

   but after the updating, those quite far away walls (not only those "deleted" walls, but also
   walls which will unlikely be useful again) will be actually deleted, and wallid of the affected brushes are 
   updated
 */


int xlatwall_cd()
{
}

int xlatwall_rearrange()
{
}

int xlatwall(sctr_t *sctr, int wallid, v3_t x, v3_t *xdone)
{
}

void csgsub()
{
}

static void updmptls(sctr_t *sctr)
{
	/* remove all old metaptls */
	for(int i=0;i<sctr->nwalls;i++){
		wall_t *w = sctr->walls+i;
		rmmptllist(w->mptllist);
/* 		w->mptllist = NULL; */
	}
  	genmptls(sctr); 
}

static inline int bndcol(float b[2][3], float b2[2][3])
{

}

#define BNDINTX(b,b2) ((b)[0][0]<=(b2)[1][0] && (b)[1][0]>=(b2)[0][0] && \
		       (b)[0][1]<=(b2)[1][1] && (b)[1][1]>=(b2)[0][1] && \
		       (b)[0][2]<=(b2)[1][2] && (b)[1][2]>=(b2)[0][2])


static float findtoi(phd_t *phd, void (*phdfunc)(phd_t *phd, void *data), void (*step)(void *data, float dt), void *data, 
		     phd_t *touchers[], int n, float len, float threshold)
{
	float dt=1.0f;
	int col;
	float t = 0.0f;
	float validt = 0.0f;
/* 	v3_t d=dx; */
	float dlen=len;
	float curlen=0;
	float validlen=0;
	float invalidlen=1;
	v3_t validx={0,0,0};
	int iii=0;
	do{
		iii++;
/* 		xlatphdpln(phd, wallid, d); */
		phdfunc(phd, data);

		t+=dt;
/* 		x = add(x,d); */
		curlen += dlen;
		col = 0;
		for(int i=0;i<n;i++){
			
			if(BNDINTX(phd->bnd, touchers[i]->bnd)){
				printf("jj\n");
				col=1;
				break;
			}
			int a=phdqry(phd, touchers[i]);
			if(a<0){
				printf("kk\n");
				col=1;
				break;
			}
		}
		if(col){
			dt = -.5*fabs(dt);
			if(curlen<invalidlen) invalidlen = curlen;
		}
		else{
			if(iii==1) break; /* move from start to the end, without any collision */
			dt = .5*fabs(dt);
			if(curlen>validlen) validlen = curlen;
/* 			validx=x; */

			validt = t;
		}
/* 		d = scale(dx,dt); */
		step(data, dt);
		dlen = dt*len;
	}while(invalidlen-validlen>threshold);
	printf("invalid=%f, valid=%f, threshold=%f, t=%f\n",invalidlen,validlen,threshold,t);
	return validt;
}


static struct xlatphdpln_data_s{
	int wallid;
	v3_t d;
	v3_t dx;
};


static void xlatphdpln_wrap(phd_t *phd, struct xlatphdpln_data_s *data)
{
	xlatphdpln(phd, data->wallid, data->d);
}

static void xlatphdpln_step(struct xlatphdpln_data_s *data, float dt)
{
	data->d = scale(data->dx, dt);
}


static struct xlatphd_data_s{
	v3_t d;
	v3_t dx;
};


static void xlatphd_wrap(phd_t *phd, struct xlatphd_data_s *data)
{
	xlatphd(phd, data->d);
}

static void xlatphd_step(struct xlatphd_data_s *data, float dt)
{
	data->d = scale(data->dx, dt);
}

static void rotphd_wrap(phd_t *phd, rot_t *r)
{
	rotphd(phd, r);
}

static void rotphd_step(rot_t *r, float dt)
{
	r->deg *= dt;
}


/* deferred */
int xlatwall_def(sctr_t *sctr, int wallid, v3_t dx, v3_t *xdone)
{
	/* adjust the offset vector */
	int block=0;
	v3_t x={0,0,0};
	wall_t *w = sctr->walls + wallid;
	float proj = dot(dx, sctr->phd->polys[wallid].plane.n);
	dx = scale(dx, proj);
	if(fabs(proj)<0.01){
		*xdone = x;
		return 1;
	}
	else if(proj>0){
		if(w->mptllist){
			*xdone = x;
			return 1;
		}
	}
	else{
		xlatphdpln(sctr->phd, wallid, dx);
		updbnd(map.sp, sctr->spid, sctr->phd->bnd);
		updmptls(sctr);
		*xdone = dx;
		return 0;
	}

	float len=length(dx);
	float bnd4d[2][3];
	int n;
	xlatbnd4d(bnd4d, sctr->phd->bnd, &dx);
	n = bndqry4(map.sp, sctr->spid, bnd4d, meh_touchers);
	if(!n){
		xlatphdpln(sctr->phd, wallid, dx);
		updbnd(map.sp, sctr->spid, sctr->phd->bnd);
		updmptls(sctr);
		*xdone = dx;
		return 0;	/* n=0 mean now there's no contacting sctrs, and there won't be any after the tanslation, either. */
	}

	phd_t *touchers[64];
	for(int i=0; i<n; i++){
		touchers[i] = ((sctr_t *)SP_GETOWNER(map.sp, meh_touchers[i]))->phd;
	}
	struct xlatphdpln_data_s data = {wallid, dx, dx};
	phd_t *phd = dupphd(sctr->phd);
	float t =  findtoi(phd, xlatphdpln_wrap, xlatphdpln_step, &data, touchers, n, length(dx), 0.01f);
	v3_t validx = scale(dx,t);
	xlatphdpln(sctr->phd, wallid, validx);
	updbnd(map.sp, sctr->spid, sctr->phd->bnd);
	updmptls(sctr);
	*xdone = validx;
	freephd(phd);
	return 0;
}



/* deferred */
int xlatsctr_def(sctr_t *sctr, v3_t dx, v3_t *xdone)
{
	/* adjust the offset vector */
	int block=0;
	v3_t x={0,0,0};
	for(int i=0;i<sctr->nwalls;i++){
		wall_t *w = sctr->walls + i;
		if(!w->mptllist){
			continue;
		}
		plane_t *pln = &sctr->phd->polys[i].plane;
		if(dot(dx, pln->n)>0){
			if(block++){
				*xdone = x;
				return 1;
			}
			v3_t v = normalize(cross(pln->n, cross(dx, pln->n)));
			dx = scale(v, dot(dx,v));
		}
	}

	float len=length(dx);
	if(len<0.01){
		*xdone=x;
		return 1;
	}

	float bnd4d[2][3];
	int n;
	xlatbnd4d(bnd4d, sctr->phd->bnd, &dx);
	n = bndqry4(map.sp, sctr->spid, bnd4d, meh_touchers);
	if(!n){
		xlatphd(sctr->phd, dx);
		updbnd(map.sp, sctr->spid, sctr->phd->bnd);
		updmptls(sctr);
		*xdone = dx;
		return 0;	/* n=0 mean now there's no contacting sctrs, and there won't be any after the tanslation, either. */
	}

	phd_t *touchers[64];
	for(int i=0; i<n; i++){
		touchers[i] = ((sctr_t *)SP_GETOWNER(map.sp, meh_touchers[i]))->phd;
	}
	struct xlatphd_data_s data = {dx, dx};
	phd_t *phd = dupphd(sctr->phd);
	float t =  findtoi(phd, xlatphd_wrap, xlatphd_step, &data, touchers, n, length(dx), 0.01f);
	v3_t validx = scale(dx,t);
	xlatphd(sctr->phd, validx);
	updbnd(map.sp, sctr->spid, sctr->phd->bnd);
	updmptls(sctr);
	*xdone = validx;
	freephd(phd);
	return 0;
}



void xlatsctr_upd(sctr_t *s, v3_t x)
{
	for(int i=0;i<s->nwalls;i++){
		wall_t *w = s->walls+i;
		for(brush_t *b=w->brushlist; b; b=b->next){
			phd_t *p = b->phd;
			for(int j=0; j<p->nverts; j++){
				p->verts[j] = add(p->verts[j],x);
			}
		}
	}
	for(brush_t *b=s->brushlist; b; b=b->next){
		phd_t *p = b->phd;
		for(int j=0; j<p->nverts; j++){
			p->verts[j] = add(p->verts[j],x);
		}
	}
	xlatsp(s->sp,x);
}


int xlatsctr(sctr_t *sctr, v3_t dx)
{
	v3_t xdone;
	int stuck = xlatsctr_def(sctr,dx,&xdone);
	if(!stuck) xlatsctr_upd(sctr,xdone);
	return stuck;
}


int xlatsctr_locked()
{
}


int rotwall()
{
}

int rotwall_cd()
{
}

int rotwall_rearrange()
{
}

int rotsctr_locked()
{
}

static v3_t bndcntr(float bnd[2][3])
{
	v3_t c;
	c.x = (bnd[0][0] + bnd[1][0])*.5f;
	c.y = (bnd[0][1] + bnd[1][1])*.5f;
	c.z = (bnd[0][2] + bnd[1][2])*.5f;
	return c;
}

static float bndradius(float bnd[2][3])
{
	v3_t diag;
	diag.x = (bnd[1][0]-bnd[0][0])*.5f;
	diag.y = (bnd[1][1]-bnd[0][1])*.5f;
	diag.z = (bnd[1][2]-bnd[0][2])*.5f;
	return length(diag);
}

/* bounding box of a rotating point */
void bnd_of_rotpnt(float bnd[2][3], rot_t *rot, v3_t v)
{

}


void bnd_of_rotphd(float bnd[2][3], rot_t *rot, phd_t *phd)
{

}


#if 0
void sumbnds(float dest[2][3], float a[2][3], float b[2][3])
{
	enum{A,B};
	enum{XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX};
	float *bnd[2]={a,b};

	/* [0][0], [0][1], [0][2], [1][0], [1][1], [1][2] */

	dest[0][0] = bnd[bnd[A][XMIN] > bnd[B][XMIN]][XMIN];
	dest[1][0] = bnd[bnd[A][XMAX] < bnd[B][XMAX]][XMAX];
	dest[0][1] = bnd[bnd[A][YMIN] > bnd[B][YMIN]][YMIN];
	dest[1][1] = bnd[bnd[A][YMAX] < bnd[B][YMAX]][YMAX];
	dest[0][2] = bnd[bnd[A][ZMIN] > bnd[B][ZMIN]][ZMIN];
	dest[1][2] = bnd[bnd[A][ZMAX] < bnd[B][ZMAX]][ZMAX];
}
#endif

void sumbnds(float a[6], float b[6])
{
	if(a[0]>b[0]) a[0]=b[0];
	if(a[1]>b[1]) a[1]=b[1];
	if(a[2]>b[2]) a[2]=b[2];
	if(a[3]<b[3]) a[3]=b[3];
	if(a[4]<b[4]) a[4]=b[4];
	if(a[5]<b[5]) a[5]=b[5];
}


int rotsctr_def(sctr_t *sctr, rot_t *rot, float *t1)
{
	/* radius of rotation */
	v3_t cntr = bndcntr(sctr->phd->bnd);
	float size = bndradius(sctr->phd->bnd);
	v3_t rr = sub(cntr, rot->org);
	rr = sub(rr, scale(rot->axis, dot(rr, rot->axis)));
	
		
	/* adjust the offset vector */
	int block=0;
	float arc = DEG2RAD(rot->deg) * (length(rr) + size);
	if(arc<0.01f){
		*t1 = 1.0f;
		return 1;
	}

#if 0
	float bnd4d[2][3];
	int n;
	n = bndqry4(map.sp, sctr->spid, bnd4d, meh_touchers);
	if(!n){
		xlatphd(sctr->phd, dx);
		updbnd(map.sp, sctr->spid, sctr->phd->bnd);
		updmptls(sctr);
		*t1 = 1.0f;
		return 0;	/* n=0 mean now there's no contacting sctrs, and there won't be any after the tanslation, either. */
	}

	phd_t *touchers[64];
	for(int i=0; i<n; i++){
		touchers[i] = SP_GETOWNER(map.sp, meh_touchers[i]);
	}
	rot_t rot =*rot;
	phd_t *phd = dupphd(sctr->phd);
	float t =  findtoi(phd, rotphd_wrap, rotphd_step, &rot, touchers, n, arc, 0.01f);

	xlatphd(sctr->phd, sub(validx,x));
	updbnd(map.sp, sctr->spid, sctr->phd->bnd);
	updmptls(sctr);
	*xdone = validx;
	freephd(phd);
	return 0;
#endif
}


void rotsctr_upd(sctr_t *s, rot_t *rot)
{
	float rotmat[3][3];
	rot2mat(rot, rotmat);
	for(int i=0; i<s->nwalls; i++){
		for(brush_t *b=s->walls[i].brushlist; b; b=b->next){
			rotphd_bat(b->phd, rotmat, rot->org);
			updphdbnd(b->phd);
			updbnd_def(s->sp, b->spid, b->phd->bnd);
		}
	}
	for(brush_t *b=s->brushlist; b; b=b->next){
		rotphd_bat(b->phd, rotmat, rot->org);
		updphdbnd(b->phd);
		updbnd_def(s->sp, b->spid, b->phd->bnd);
	}
	updallbnds(s->sp);
}

int rotsctr()
{
}



int xlatpnt(pnt_t *p, float *frac)
{
}


/* translate a point without collision detection */
int xlatpnt_nocd(pnt_t *p)
{
}


void unlinkbrush(sctr_t *sctr, brush_t *b)
{
	if(b->prev){ 
		b->prev->next = b->next; 
	}
	else{
		if(b->wallid>=0){
			sctr->walls[b->wallid].brushlist = b->next; 
		}
		else{
			sctr->brushlist = b->next; 
		}	
	}
	if(b->next){ b->next->prev = b->prev; }
}




void delbrush(sctr_t *sctr, brush_t *b)
{
}


int brush2obj(sctr_t *sctr, brush_t *b, int cd)
{
	int stuck=0;
	int n=rngqry(sctr->sp, b->phd->bnd, meh_touchers);
	for(int i=0; i<n; i++){
		stuff_t *stuff = SP_GETOWNER(sctr->sp, meh_touchers[i]);
		if(stuff->type == NULL_BRUSH_TYPE){
			int a=phdqry(b->phd, stuff->common.phd);
			if(a>=0){ 
				if(cd){	return 1; }
				stuck=1;
				break;
			}
		}
		/* else: toucher->type == OBJ_TYPE. since obj won't intersect with anything including brush, we omit the test */
	}
	obj_t *obj = malloc(sizeof(obj_t));
	obj->type = OBJ_TYPE;
	obj->spid = b->spid;
	if(stuck){ obj->flag |= OBJ_PHANTOM; }
	SP_SETOWNER(sctr->sp, obj->spid, obj);
	setdbnd(sctr->sp, obj->spid);
	obj->phd = b->phd;
	unlinkbrush(sctr, b);
	free(b);

	obj->prev = NULL;
	obj->next = map.objlist;
	if(obj->next) obj->next->prev = obj;
	map.objlist = obj;

	/* TODO: update portals */
	updmass(obj);

	return stuck;
}





typedef struct{
	v3_t intx;
	float toi;
} ccdinfo_t;

/* transform objects with continuous collision detection */
int xlatobj(v3_t *xlat, obj_t *obj, sctr_t *sctr, ccdinfo_t *ccdinfo)
{
}


static inline int bbintx(float b[6], float b2[6])
{
	return (b[0]<b2[3] && b[3]>b2[0] && 
		b[1]<b2[4] && b[4]>b2[1] &&
		b[2]<b2[5] && b[5]>b2[2]);
}

static inline int pntinbnd(float p[3], float bnd[6])
{
	return (p[0]>=bnd[0] && p[0]<=bnd[3] &&
		p[1]>=bnd[1] && p[1]<=bnd[4] &&
		p[2]>=bnd[2] && p[2]<=bnd[5]);
}

static inline void mergebnds(float dst[6], float a[6], float b[6])
{
	dst[0] = a[0]<b[0] ? a[0] : b[0];
	dst[1] = a[1]<b[1] ? a[1] : b[1];
	dst[2] = a[2]<b[2] ? a[2] : b[2];
	dst[3] = a[3]>b[3] ? a[3] : b[3];
	dst[4] = a[4]>b[4] ? a[4] : b[4];
	dst[5] = a[5]>b[5] ? a[5] : b[5];
}


#define neg3(c,a) {(c)[0]=-(a)[0]; (c)[1]=-(a)[1]; (c)[2]=-(a)[2];}
#define add3(c,a,b) {(c)[0]=(a)[0]+(b)[0]; (c)[1]=(a)[1]+(b)[1]; (c)[2]=(a)[2]+(b)[2];}
#define sub3(c,a,b) {(c)[0]=(a)[0]-(b)[0]; (c)[1]=(a)[1]-(b)[1]; (c)[2]=(a)[2]-(b)[2];}
#define dot3(a,b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])
#define cross3(c,a,b) {(c)[0]=((a)[1])*((b)[2]) - ((a)[2])*((b)[1]);	\
		(c)[1]=((a)[2])*((b)[0]) - ((a)[0])*((b)[2]);		\
		(c)[2]=((a)[0])*((b)[1]) - ((a)[1])*((b)[0]);}
#define scale3(c,a,k) {(c)[0]=(a)[0]*(k);	\
		(c)[1]=(a)[1]*(k);		\
		(c)[2]=(a)[2]*(k);}
#define mov3(c,a) {(c)[0]=(a)[0]; (c)[1]=(a)[1]; (c)[2]=(a)[2];}
#define mad3(c,a,k,b) {(c)[0]=(a)[0]*(k)+(b)[0];			\
		(c)[1]=(a)[1]*(k)+(b)[1];				\
		(c)[2]=(a)[2]*(k)+(b)[2];}
static int xlatvpln(const float x[3], const float v[3], const plane_t *pln, float *tmin, float itx[3])
{
	float v2[3];
	float prjx = dot3(x, (float*)&pln->n);
	float prjv;
	float t;
	if(prjx>-0.00001) return 0; /* ignore backface */

	add3(v2,v,x);
	prjv = dot3(v, (float*)&pln->n);
	if(prjv<pln->d || prjv+prjx>pln->d) return 0; /* the trace is totally on one side of the pln */

	t = (pln->d - prjv) / prjx;
	if(t>=0 && t<*tmin){
		*tmin = t;
		scale3(itx, x, t);
		add3(itx, itx, v);
		return 1;
	}
	return 0;
}

static int xlatee(float x[3], float d1[3], float d2[3], float e1[3], float e2[3], float *tmin, float itx[3])
{
	float d[3], e[3];
	float de[3];
	float x_x_e[3];
	float de_x_x[3];
	float denom;
	float s,u,v;
	sub3(d, d2, d1);
	sub3(e, e2, e1);
	cross3(x_x_e, x, e); 
	denom = dot3(d, x_x_e);
	if(fabs(denom) < 0.00001){
		return 0;
	}

	sub3(de, e1, d1);
	cross3(de_x_x, de, x);
	
	s = dot3(d, de_x_x)/denom;
	if(s<0 || s>1) return 0;

	u = dot3(e, de_x_x)/denom;
	if(u<0 || u>1) return 0;
	
	cross3(de_x_x, de, e);	/* de_x_x is used to store de_x_e */
	v = dot3(d, de_x_x)/denom;
	if(v>=0 && v<*tmin){
		*tmin = v;
		mad3(itx, e, s, e1);
		return 1;
	}
	return 0;
}


void testxlatee()
{
	v3_t d1={0,0,0},d2={2,0,0};
	v3_t e1={1.5,1,-1},e2={1,2,1};
	v3_t x={0,5,0};
	v3_t itx;
	float tmin = 1;
	int stuck=xlatee(&x,&d1,&d2,&e1,&e2,&tmin,&itx);
	if(stuck){
		printf("itx=(%f %f %f)\n", itx.x, itx.y, itx.z);
	}
	glLineWidth(3);
	glBegin(GL_LINES);
	glColor3f(0,0,0);
	glVertex3fv(&d1);
	glVertex3fv(&d2);
	glColor3f(1,1,1);
	glVertex3fv(&e1);
	glVertex3fv(&e2);
	glEnd();
	glLineWidth(1);

	glPointSize(3);
	glColor3f(1,.5,.2);
	glBegin(GL_POINTS);
	glVertex3fv((float*)&itx);
	glEnd();
	glPointSize(1);

}


#define EDGECACHE 1

static struct {
	v3_t *x0;
	v3_t xlat;
	phd_t *phd;
	v3_t intx;
	float t;
	float pbnd[2][3];
	float ebnd[2][3];	/* epsilon bnd -- used for contact detection */
	ccdinfo_t *ccd;
	int stuck;
} meh_xp;

#if EDGECACHE
static inline int pntinpoly(v3_t *pnt, phd_t *phd, int pid)
{
	poly_t *poly = phd->polys + pid;
	for(int i=0; i<poly->nv; i++){
		edge_t *e = phd->edges + poly->e[i];
		float d=dot(*pnt, e->pln.n);
		if((pid==e->p[0]) ^ (d>e->pln.d)) return 0;
	}
	return 1;
}
#else

static inline int pntinpoly(v3_t *pnt, plane_t *edgeplns, int ne)
{
	for(int i=0; i<ne; i++){
		float d=dot(*pnt, edgeplns[i].n);
		printf("pnt=(%f %f %f)\nd=%f, edgeplns[i].d=%f\n", pnt->x,pnt->y,pnt->z,d, edgeplns[i].d);
		if(d >= edgeplns[i].d){ return 0; }
	}
	return 1;
}
#endif


static void beginxlatphd(v3_t *x, phd_t *phd, ccdinfo_t *ccd)
{
	meh_xp.x0 = x;
	meh_xp.xlat = *x;
	meh_xp.phd = phd;
	meh_xp.t = 1;
	xlatbnd4d(meh_xp.pbnd, phd->bnd, x);
	meh_xp.stuck = 0;
	meh_xp.ccd = ccd;
}


static int xlatagstphd(phd_t *phd)
{
	/*
	  phd1->v vs. phd2->p
	  phd1->p vs. phd2->v
	  phd1->e vs. phd2->e
	*/
	float bnd[2][3];
	char v[32];
	int nv=0;
	float negx[3];
	float *x = &meh_xp.xlat;
	int coll=0;
	neg3(negx, ARRAY(meh_xp.xlat));
	xlatbnd4d(bnd, phd->bnd, negx);

	if(!bbintx(bnd, meh_xp.pbnd)) return 0;

	/* 1. */
	for(int i=0; i<meh_xp.phd->nverts; i++){
		if(pntinbnd(meh_xp.phd->verts+i, bnd)) v[nv++] = i;
	}

	for(int i=0; i<phd->npolys; i++){
		float polybnd[2][3];
		poly_t *poly = phd->polys + i;

#if EDGECACHE
#else
		plane_t edgeplns[16];
		for(int j=0; j<poly->nv; j++){ /* TODO: cache this! */
			v3_t *v = phd->verts + poly->v[j];
			v3_t *v2 = phd->verts + poly->v[(j+1)%poly->nv];
			edgeplns[j].n = cross(poly->plane.n, sub(*v2, *v));
			edgeplns[j].d = dot(*v, edgeplns[j].n);
		}
#endif

		poly2bnd(phd, i, polybnd);
		xlatbnd4d(polybnd, polybnd, negx);
		for(int j=0; j<nv; j++){
			v3_t *pnt = meh_xp.phd->verts + v[j];
			float tmin=meh_xp.t;
			v3_t intx;
			int stuck;
			if(!pntinbnd(pnt, polybnd)) continue;
			stuck = xlatvpln(x, pnt, &poly->plane, &tmin, &intx);
			if(stuck){
#if EDGECACHE
 				if(pntinpoly(&intx, phd, i)){ 
#else
 				if(pntinpoly(&intx, edgeplns, poly->nv)){ 
#endif
					meh_xp.t = tmin;
					meh_xp.intx = intx;
					coll=1;
				}
			}
		}		
	}
	
	/* 2. */
	nv = 0;
	for(int i=0; i<phd->nverts; i++){
		if(pntinbnd(phd->verts+i, meh_xp.pbnd)) v[nv++] = i;
	}

	for(int i=0; i<meh_xp.phd->npolys; i++){
		float polybnd[2][3];
		poly_t *poly = meh_xp.phd->polys + i;
#if EDGECACHE
#else
		plane_t edgeplns[16];

		for(int j=0; j<poly->nv; j++){ /* TODO: cache this! */
			v3_t *v = meh_xp.phd->verts + poly->v[j];
			v3_t *v2 = meh_xp.phd->verts + poly->v[(j+1)%poly->nv];
			edgeplns[j].n = cross(poly->plane.n, sub(*v2, *v));
			edgeplns[j].d = dot(*v, edgeplns[j].n);
		}
#endif

		poly2bnd(meh_xp.phd, i, polybnd);
		xlatbnd4d(polybnd, polybnd, x);
		for(int j=0; j<nv; j++){
			v3_t *pnt = phd->verts + v[j];
			float tmin = meh_xp.t;
			v3_t intx;
			int stuck;
			if(!pntinbnd(pnt, polybnd)) continue;
			stuck = xlatvpln(negx, pnt, &poly->plane, &tmin, &intx);
			if(stuck){
				intx = add(intx, meh_xp.xlat);
#if EDGECACHE
 				if(pntinpoly(&intx, meh_xp.phd, i)){ 
#else
 				if(pntinpoly(&intx, edgeplns, poly->nv)){ 
#endif
					meh_xp.t = tmin;
					meh_xp.intx = intx;
					coll=1;
				}
			}
		}		
	}

	/* 3. */
	for(int i=0; i<meh_xp.phd->nedges; i++){
		float *d1,*d2,*e1,*e2;
		int *v;
		v = meh_xp.phd->edges[i].v;
		d1 = meh_xp.phd->verts + v[0];
		d2 = meh_xp.phd->verts + v[1];
		for(int j=0; j<phd->nedges; j++){
			v = phd->edges[j].v;
			e1 = phd->verts + v[0];
			e2 = phd->verts + v[1];
			float tmin = meh_xp.t;
			v3_t intx;
			int stuck = xlatee(x, d1, d2, e1, e2, &tmin, &intx);
			if(stuck){
/* 				if(tmin < meh_xp.t){ */
					meh_xp.t = tmin;
					meh_xp.intx = intx;
					coll=1;
/* 				} */
			}
		}
	}
	meh_xp.stuck += coll;
	return coll;
}

#define CONTACT_EPSILON 0.01

#define MAX_NCONTACTS 32
static v3_t meh_contacts[MAX_NCONTACTS];
static int meh_ncontacts;


static void begincontact(phd_t *phd, v3_t **contacts)
{
	ccdinfo_t ccd;
	v3_t x = {0,0,0};
	beginxlatphd(&x, phd, &ccd);
	meh_ncontacts = 0;
	*contacts = meh_contacts;
#if 0
	for(int i=0; i<3; i++){
		meh_xp.pbnd[0][i] -= CONTACT_EPSILON;
		meh_xp.pbnd[1][i] += CONTACT_EPSILON;
	}
#endif
}

static int endcontact()
{
	return meh_ncontacts;
}


static inline int segbndintx(float v[3], float v2[3], float bnd[2][3])
{
	/* v->v2 */
	float d[3];
	sub3(d,v2,v);

	/* x,y,z */
	for(int i=0; i<3; i++){
		if(v[i]<v2[i]){
			if(v[i]>bnd[1][i] || v2[i]<bnd[0][i]) return 0;
		}
		else{
			if(v[i]<bnd[0][i] || v2[i]<bnd[1][i]) return 0;
		}
	}

	/*  */
}


static inline int bbintx2(float dest[6], float a[6], float b[6])
{
	for(int i=0; i<3; i++){
		int j=i+3;
		dest[i] = (a[i]>b[i] ? a[i] : b[i]) - CONTACT_EPSILON;
		dest[j] = (a[j]<b[j] ? a[j] : b[j]) + CONTACT_EPSILON;
		if(dest[i]>dest[j]) return 0;
	}
	return 1;
}



static inline int eeintx(phd_t *phd, int eid, phd_t *phd2, int eid2, float intx[3])
{
	edge_t *e = phd->edges + eid;
	edge_t *e2 = phd2->edges + eid2;

	v3_t ev = sub(phd->verts[e->v[1]], phd->verts[e->v[0]]);
	v3_t ev2 = sub(phd2->verts[e2->v[1]], phd2->verts[e2->v[0]]);
	v3_t axis = cross(ev, ev2);
	float a=dot(axis,axis);
	if(a<0.0001) return 0;

	axis = scale(axis, 1.0f/sqrt(a));
	float d = dot(axis, sub(phd->verts[e->v[0]], phd2->verts[e2->v[0]]));
	if(fabs(d)>CONTACT_EPSILON) return 0;

	plane_t pln[2];
	pln[0].n = cross(ev, axis);
	pln[0].d = dot(phd->verts[e->v[0]], pln[0].n);

	
}


/* TODO: use GJK/SAT to find contact points */
static int contact(phd_t *phd)
{
	int cont=0;
	float intxbnd[2][3];
	if(!bbintx2(intxbnd, meh_xp.pbnd, phd->bnd)) return 0;
	
	for(int i=0; i<meh_xp.phd->nverts; i++){
		v3_t *v = meh_xp.phd->verts + i; 
		if(!pntinbnd(v, intxbnd)) continue;
		for(int j=0; j<phd->npolys; j++){
			poly_t *poly = phd->polys + j;
			float d = dot(*v, poly->plane.n) - poly->plane.d;
			if(d<-CONTACT_EPSILON || d>CONTACT_EPSILON) continue;
			if(pntinpoly(v, phd, j)){
				meh_contacts[meh_ncontacts++] = *v;
				cont=1;
			}
		}
		
	}

	for(int i=0; i<phd->nverts; i++){
		v3_t *v = phd->verts + i;
		if(!pntinbnd(v, intxbnd)) continue;
		for(int j=0; j<meh_xp.phd->npolys; j++){
			poly_t *poly = meh_xp.phd->polys + j;
			float d = dot(*v, poly->plane.n) - poly->plane.d;
			if(d<-CONTACT_EPSILON || d>CONTACT_EPSILON) continue;
			if(pntinpoly(v, meh_xp.phd, j)){
				meh_contacts[meh_ncontacts++] = *v;
				cont=1;
			}
		}
	}

#if 0
	/* meh_xp.phd vs. intxbnd */
	for(int i=0; i<meh_xp.phd->nedges; i++){
	}

	/* phd vs. intxbnd */
	for(int i=0; i<phd->nedges; i++){
	}

	for(int i=0; i<nme; i++){
		for(int j=0; j<ne; j++){
			
		}
	}
#endif


	for(int i=0; i<meh_xp.phd->nedges; i++){
		edge_t *me = meh_xp.phd->edges + i;
		v3_t *mv[2] = {meh_xp.phd->verts+me->v[0], meh_xp.phd->verts+me->v[1]};

		for(int j=0; j<phd->nedges; j++){
			edge_t *e = phd->edges + j;

			v3_t *v[2] = {phd->verts+e->v[0], phd->verts+e->v[1]};
			float d = dot(me->pln.n, *v[0]);
			
			if((d<me->pln.d)^(dot(me->pln.n, *v[1])>me->pln.d)) continue;

			v3_t ev = sub(*v[1], *v[0]);
			float denom = dot(me->pln.n, ev);
			if(denom < 0.0001) continue;

			float t = (me->pln.d - dot(me->pln.n, *v[0])) / denom;
			if(t<0.0f || t>1.0f) continue;

			v3_t mev = sub(*mv[1], *mv[0]);
			v3_t intx = add(*v[0], scale(ev, t));
			float d0 = dot(sub(intx, *mv[0]), mev);
			float d1;
			if(d0<0 || (d1=dot(sub(intx, *mv[1]), mev))>0) continue;
			
			float k = -d0/d1;
			v3_t dist = sub(intx, add(*mv[0], scale(mev, k/(k+1))));
			if(dot(dist,dist) < CONTACT_EPSILON*.1){
				meh_contacts[meh_ncontacts++] = intx;
				cont=1;
			}
		}
	}
	meh_xp.stuck += cont;
	return cont;
}



#if 0
static int contact(phd_t *phd)
{
	/*
	  phd1->v vs. phd2->p
	  phd1->p vs. phd2->v
	  phd1->e vs. phd2->e
	*/
	float bnd[2][3];
	char v[32];
	int nv=0;
	float negx[3];
	float *x = &meh_xp.xlat;
	int cont=0;
	neg3(negx, ARRAY(meh_xp.xlat));
	xlatbnd4d(bnd, phd->bnd, negx);
	if(!bbintx(bnd, meh_xp.pbnd)) return 0;

	/* 1. */
	for(int i=0; i<meh_xp.phd->nverts; i++){
		if(pntinbnd(meh_xp.phd->verts+i, bnd)) v[nv++] = i;
	}
	for(int i=0; i<phd->npolys; i++){
		float polybnd[2][3];
		poly_t *poly = phd->polys + i;
#if EDGECACHE
#else
		plane_t edgeplns[16];
		
		for(int j=0; j<poly->nv; j++){ /* TODO: cache this! */
			v3_t *v = phd->verts + poly->v[j];
			v3_t *v2 = phd->verts + poly->v[(j+1)%poly->nv];
			edgeplns[j].n = cross(poly->plane.n, sub(*v2, *v));
			edgeplns[j].d = dot(*v, edgeplns[j].n);
		}
#endif

		poly2bnd(phd, i, polybnd);
		xlatbnd4d(polybnd, polybnd, negx);
		for(int j=0; j<nv; j++){
			v3_t *pnt = meh_xp.phd->verts + v[j];
			float tmin=1;
			v3_t intx;
			int stuck;
			if(!pntinbnd(pnt, polybnd)) continue;
			stuck = xlatvpln(x, pnt, &poly->plane, &tmin, &intx);

			if(stuck){
#if EDGECACHE
 				if(pntinpoly(&intx, phd, i)){
#else
 				if(pntinpoly(&intx, edgeplns, poly->nv)){ 
#endif
					if(meh_ncontacts < MAX_NCONTACTS){
						meh_contacts[meh_ncontacts++] = intx;
						cont=1;
					}
				}
			}
		}		
	}
	
	/* 2. */
	nv = 0;
	for(int i=0; i<phd->nverts; i++){
		if(pntinbnd(phd->verts+i, meh_xp.pbnd)) v[nv++] = i;
	}

	for(int i=0; i<meh_xp.phd->npolys; i++){
		float polybnd[2][3];
		poly_t *poly = meh_xp.phd->polys + i;
#if EDGECACHE
#else
		plane_t edgeplns[16];

		for(int j=0; j<poly->nv; j++){ /* TODO: cache this! */
			v3_t *v = meh_xp.phd->verts + poly->v[j];
			v3_t *v2 = meh_xp.phd->verts + poly->v[(j+1)%poly->nv];
			edgeplns[j].n = cross(poly->plane.n, sub(*v2, *v));
			edgeplns[j].d = dot(*v, edgeplns[j].n);
		}
#endif

		poly2bnd(meh_xp.phd, i, polybnd);
		xlatbnd4d(polybnd, polybnd, x);
		for(int j=0; j<nv; j++){
			v3_t *pnt = phd->verts + v[j];
			float tmin = 1;
			v3_t intx;
			int stuck;
			if(!pntinbnd(pnt, polybnd)) continue;
			stuck = xlatvpln(negx, pnt, &poly->plane, &tmin, &intx);
			if(stuck){
				intx = add(intx, meh_xp.xlat);
#if EDGECACHE
 				if(pntinpoly(&intx, meh_xp.phd, i)){ 
#else
 				if(pntinpoly(&intx, edgeplns, poly->nv)){ 
#endif
					if(meh_ncontacts < MAX_NCONTACTS){
						meh_contacts[meh_ncontacts++] = intx;
						cont=1;
					}
				}
			}
		}		
	}

	/* 3. */
	for(int i=0; i<meh_xp.phd->nedges; i++){
		float *d1,*d2,*e1,*e2;
		int *v;
		v = meh_xp.phd->edges[i].v;
		d1 = meh_xp.phd->verts + v[0];
		d2 = meh_xp.phd->verts + v[1];
		for(int j=0; j<phd->nedges; j++){
			v = phd->edges[j].v;
			e1 = phd->verts + v[0];
			e2 = phd->verts + v[1];
			float tmin = meh_xp.t;
			v3_t intx;
			int stuck = xlatee(x, d1, d2, e1, e2, &tmin, &intx);
			if(stuck){
				if(meh_ncontacts < MAX_NCONTACTS){
					meh_contacts[meh_ncontacts++] = intx;
					cont=1;
				}
			}
		}
	}
	meh_xp.stuck += cont;
	return cont;
}
#endif


static int endxlatphd()
{
	if(meh_xp.stuck){
	/* 	meh_xp.t *= 0.98; */ /* TODO: change this to relative threshold */
		v3_t xlat = scale(meh_xp.xlat, meh_xp.t);
		float d = length(xlat);
		meh_xp.t *= (d - CONTACT_EPSILON) / d;
		if(meh_xp.t < 0) meh_xp.t = 0;
		*meh_xp.x0 = scale(meh_xp.xlat, meh_xp.t);
		meh_xp.ccd->intx = meh_xp.intx;
		meh_xp.ccd->toi = meh_xp.t;
		return 1;
	}
	else{
		meh_xp.ccd->toi = meh_xp.t;
		return 0;
	}
}


void dbg_drawphdnormals(phd_t *phd)
{
#if 0
	glBegin(GL_LINES);
	for(int i=0; i<phd->npolys; i++){
		v3_t c = polycntr(phd, i);
		v3_t c2 = add(c, normalize(phd->polys[i].plane.n));
		glVertex3fv(&c);
		glVertex3fv(&c2);
	}
	glEnd();
#endif

	glColor3f(0,1,0);
	glBegin(GL_LINES);
	for(int i=0; i<phd->nedges; i++){
		edge_t *e = phd->edges + i;
		v3_t c = scale(add(phd->verts[e->v[0]], phd->verts[e->v[1]]),.5f);
		v3_t c2 = add(c, e->pln.n);
		glVertex3fv(&c);
		glVertex3fv(&c2);
	}
	glEnd();
}

void testxlatphd()
{

	extern phd_t *globphds[];
	static v3_t xlat = {.03, 0,0};
	ccdinfo_t ccd;

	beginxlatphd(&xlat, globphds[1], &ccd);
/* 	static int s=0; */
	static v3_t intx;
	xlatagstphd(globphds[0]);
	int stuck = endxlatphd();

/* 	if(rot.rad) rotphd(globphds[1], &rot); */

#if 1
/* 	if(!s) */{
		if(stuck){
/* 			s=1; */
/* 			if(rot.rad<0) rot.rad=0; */
			intx = ccd.intx;/*  add(meh_rp.rot.org, mattv3(meh_rp.zmat, ccd.intx)); */
/*			printf("shit deg=%f, ft=%f, intx=(%f %f %f), rot.org=(%f %f %f)\n", 
			       RAD2DEG(meh_rp.rot.rad), meh_rp.ft, intx.x, intx.y, intx.z, rot.org.x, rot.org.y, rot.org.z);
*/
		}
	}

#endif

	if(ccd.toi){
		xlatphd(globphds[1], xlat);
	}

/*    	if(stuck) xlat.x = .03;    */

/* 	printf("intx=(%f %f %f)\n", intx.x, intx.y, intx.z); */

#if 0
	if(stuck){
		if(rot.rad>0){
			rotphd(globphds[1], &rot);
			intx = ccd.intx;
		}
	}
#endif


#if 1
	glLineWidth(1);
	glColor3f(1,0,0);
	drawphd(globphds[0]);
	glColor3f(0,0,1);
	drawphd(globphds[1]);
	glLineWidth(1);

	glPointSize(2);
	glColor3f(0,1,1);
	glBegin(GL_POINTS);
	glVertex3fv(&intx);
	glEnd();
	glPointSize(1);
#endif

#if 0
	glColor3f(0,1,0);
	dbg_drawphdnormals(globphds[1]);
#endif

#if 0
	glPointSize(2);
	glColor3f(1,.5,.5);
	glBegin(GL_POINTS);
	glVertex3fv(&ii);
	glEnd();
#endif
	

#if 0
	float color[3] = {.2,.2,1};
	drawbox(meh_xp.pbnd, color);
#endif


#if 0
	rot_t rot2;
	v3_t pnt = {1,0,0};
	rot2.deg = 50;
	rot2.axis = vec(0,0,1);
	rot2.rad = DEG2RAD(rot2.deg);
	rot2.org = vec(0,0,0);
	float mat[3][3];
	rot2mat(&rot2, mat);
	v3_t pnt2 = matv3(mat, pnt);

	rot2.deg = 210;
	rot2.rad = DEG2RAD(rot2.deg);
	rot2mat(&rot2, mat);
	v3_t pnt3 = matv3(mat, pnt);

#endif

#if 0
	rot_t rot2;
	v3_t pnt = {1,0,0};
	rot2.deg = 100;
	rot2.axis = vec(0,0,1);
	rot2.rad = DEG2RAD(rot2.deg);
	rot2.org = vec(0,0,0);
	float mat[3][3];
	rot2mat(&rot2, mat);
	v3_t pnt2 = matv3(mat, pnt);
#endif
#if 0
	glLineWidth(5);
	glColor3f(.6, .3, .9);
	glBegin(GL_LINE_STRIP);
	glVertex3fv(&pnt);
	glVertex3fv(&pnt2);
	glVertex3fv(&pnt3);
	glEnd();
	glLineWidth(1);
#endif

}

void testcontact()
{
	extern phd_t *globphds[];
	static v3_t xlat = {.3, 0,0};
	v3_t *contacts;

	begincontact(globphds[1], &contacts);
/* 	static int s=0; */
	static v3_t intx;
	contact(globphds[0]);
	int n = endcontact();
	printf("KKKKKKKKKKKK n=%i\n", n);

	glPointSize(5);
	glColor3f(1,0,0);
	glBegin(GL_POINTS);
	for(int i=0; i<n; i++)
		glVertex3fv(contacts+i);
	glEnd();
	glPointSize(1);

}



typedef struct{
	v3_t o,u,r;
} rect3d_t;


static void rotpntbnd_rect3d(rot_t *rot)
{
}


static void rotpntbnd(rot_t *rot, v3_t p0, v3_t p1, float bnd[2][3])
{
	v3_t r,u,r1,proj;
	float radius;
	float c1,s1;
	float frad; /* f(rad) = -/+ cos(rad/2)^2    (f(rad) is increasing)*/
	r = sub(p0,rot->org);
	proj = scale(rot->axis, dot(r,rot->axis));
	r = sub(r, proj);
	radius = length(r);
/* 	u = cross(r,rot->axis); */
 	u = cross(rot->axis, r); 

 	r1 = sub(p1,rot->org); 
 	r1 = sub(r1, proj); 
#if 0
	c1 = dot(r1,r);
	s1 = dot(r1,u);
	frad = .5f * (1.0f + c1/radius);
	if(rot->deg<180) frad = -frad;
#endif

	for(int i=0; i<3; i++){
		float a = ((float*)&r)[i];
		float b = ((float*)&u)[i];
		int inc = b>0;
		int dec = !inc;
#if 0
		if(a){
			float t = b/a;
			float frads[3];
			float fradmax[2];
			fradmax[1] = sqrt(a*a + b*b);
			fradmax[0] = -fradmax[1];

			frads[1] = 1.0f / (1.0f + t*t);
			if(a*b>0){
				frads[0] = .5f * (1.0f + sqrt(frads[1]));
				frads[1] = -frads[1];
				frads[2] = frads[0] - 1;
			}
			else{
				frads[0] = .5f * (1.0f - sqrt(frads[1]));
				frads[2] = frads[0] + 1;
			}
			frads[0] = -frads[0];
			frads[2] = -frads[2];

			if(frad < frads[0]) bnd[inc][i] = ((float*)&r1)[i];
			else bnd[inc][i] = fradmax[inc];

			if(frad < frads[1]) bnd[dec][i] = a;
			else if(frad < frads[2]) bnd[dec][i] = ((float*)&r1)[i];
			else bnd[dec][i] = fradmax[dec];			
		}
		else{
			if(rot->rad<PI/2) bnd[inc][i] = ((float*)&r1)[i];
			else bnd[inc][i] = radius;

			if(rot->rad<PI) bnd[dec][i] = 0;
			else if(rot->rad<PI*3/2) bnd[dec][i] = ((float*)&r1)[i];
			else bnd[dec][i] = radius;
		}
#endif

#if 1
		float rad;
		float radmax[2];
		if(a){ 
			rad = atanf(b/a);
			if(rad<0) rad+=PI;
			radmax[1] = sqrt(a*a+b*b);
			radmax[0] = -radmax[1];
		}
		else { 
			rad = PI/2;
			radmax[0] = -radius;
			radmax[1] = radius;
		}

		if(rot->rad<rad){ bnd[inc][i] = ((float*)&r1)[i];}
		else bnd[inc][i] = radmax[inc];
		
		if(rot->rad<2*rad) bnd[dec][i]=a;
		else if(rot->rad<rad+PI) bnd[dec][i] = ((float*)&r1)[i];
		else bnd[dec][i] = radmax[dec];

/*		
		if(a){
			float radmax[2];
			radmax[1] = sqrt(a*a+b*b);
			radmax[0] = -radmax[1];
			if(rot->rad<rad){ bnd[inc][i] = ((float*)&r1)[i];}
			else bnd[inc][i] = radmax[inc];

			if(rot->rad<2*rad) bnd[dec][i]=a;
			else if(rot->rad<rad+PI) bnd[dec][i] = ((float*)&r1)[i];
			else bnd[dec][i] = radmax[dec];
		}
		else{
			if(rot->rad<PI/2) bnd[inc][i] = ((float*)&r1)[i];
			else bnd[inc][i] = radius;

			if(rot->rad<PI) bnd[dec][i] = 0;
			else if(rot->rad<PI*3/2) bnd[dec][i] = ((float*)&r1)[i];
			else bnd[dec][i] = radius;

		}
*/

#endif

		bnd[0][i] += ((float*)&rot->org)[i] + ((float*)&proj)[i];
		bnd[1][i] += ((float*)&rot->org)[i] + ((float*)&proj)[i];
	}
}

#define dot2(a,b) ((a)[0]*(b)[0]+(a)[1]*(b)[1])
#define cross2(a,b) ((a)[0]*(b)[1]-(a)[1]*(b)[0])
#define proj2(a,b) ((a)[0]*(b)[0]+(a)[1]*(b)[1])
#define sub2(c,a,b) {(c)[0]=(a)[0]-(b)[0];(c)[1]=(a)[1]-(b)[1];}
/* #define sub3(c,a,b) {(c)[0]=(a)[0]-(b)[0];(c)[1]=(a)[1]-(b)[1];(c)[2]=(a)[2]-(b)[2];} */
#define scale2(c,a,s) {(c)[0]=(s)*((a)[0]); (c)[1]=(s)*((a)[1]);}
#define mad2(c,a,k,b) {(c)[0]=(k)*(a)[0]+(b)[0]; (c)[1]=(k)*(a)[1]+(b)[1];} /* c=ka+b (a,b,c are v2, k is scalar) */

static int rotvpln(float v[3],  plane_t *p, float *ft, float itx[3])
{
	float a,b,c,r2;
	float *n = &p->n;
	float d;
	float i0[2],i1[2];

	r2 = dot2(v,v);
	a = dot2(n,n);
	if(!a) return 0;

	b = n[2]*v[2] - p->d;
	c = b*b - n[0]*n[0]*r2;
	b = n[1] * b;
	
	d = b*b - a*c;
	if(d<0){
		if(d<-0.00001) return 0;
		d = 0;
	}

	if(d<0.00001){
		if(n[0]){
/* 			i0[0] = (p->d - n[1]*i0[1] - n[2]*v.z) / n[0]; */
			return 0;
		}
 		else{
			float d2;
			i0[1] = -b/a;
			d2 = r2 - i0[1]*i0[1];
			if(d2<0) return 0;

			i1[0] = sqrt(d2);
			i0[0] = -i1[0];
			i1[1] = i0[1];
		}
	}
	else{
		float sqrtd = sqrt(d);
		i0[1] = (-b + sqrtd) / a;
		i1[1] = (-b - sqrtd) / a;
		/*
		  if n[0]==0, then a=n[1]*n[1], c=b*b/(n[1]*n[1])   => d = b*b - a*c = b*b - n[1]*n[1]*b*b/(n[1]*n[1]) == 0
		  so n[0]!=0
		 */
		i0[0] = (p->d - n[1]*i0[1] - n[2]*v[2]) / n[0];
		i1[0] = (p->d - n[1]*i1[1] - n[2]*v[2]) / n[0];
	}

	/* find the first intersection */
	if(1/* cross2(n,i0)>0 */){
 		int gt180 = cross2(v, i0)<0; 
 		float len = .5f*(dot2(i0,i0)+dot2(v,v)); 
/* 		float len = sqrt(dot2(i0,i0)*dot2(v,v)); */
		float cosine = dot2(i0,v)/len;
		float ft2 = gt180 ? -1.0f-cosine: 1.0f+cosine;
		if(ft2 > *ft){
			*ft = ft2;
			itx[0] = i0[0];
			itx[1] = i0[1];
			itx[2] = v[2];
			return 1;
		}
	}
	/* else  */if(1/* cross2(n,i1)>0 */){ /* collide with back face: ignored */
 		int gt180 = cross2(v, i1)<0; 
 		float len = .5f*(dot2(i1,i1)+dot2(v,v)); 
/* 		float len = sqrt(dot2(i1,i1)*dot2(v,v)); */
		float cosine = dot2(i1,v)/len;
		float ft2 = gt180 ? -1.0f-cosine : 1.0f+cosine;
		if(ft2 > *ft){
			*ft = ft2;
			itx[0] = i1[0];
			itx[1] = i1[1];
			itx[2] = v[2];
			return 1;
		}
	}

	return 0;
}

static int rotee(float d1[3], float d2[3], float e1[3], float e2[3], float *ft, float itx[3])
{
	float d[3], e[3], v[2], i1[3], i2[3];
	float a,b,c,t,t1,t2;
	float dz2, dd, ee;
	float b2ac;

	sub3(d,d2,d1);
	sub3(e,e2,e1);
	dz2 = d[2]*d[2];

	{
		float k1,k2,b1,b2;	/* f1(t) = k1 t + b1;  f2(t) = k2 t + b2 */
		k1 = d[0]*e[2];
		b1 = d[0]*(e1[2]-d1[2]) + d[2]*d1[0];
		k2 = d[1]*e[2];
		b2 = d[1]*(e1[2]-d1[2]) + d[2]*d1[1];
		a = dot2(e,e)*dz2 - k1*k1 - k2*k2;
		b = dot2(e,e1)*dz2 - k1*b1 - k2*b2;
		c = dot2(e1,e1)*dz2 - b1*b1 - b2*b2;
	}

	if(fabs(a)<0.00001) return 0; /* intersection is the edge */

	b2ac = b*b-a*c;

	/*
	  calculated by using xmaxima, if d[2]==0, then b*b-a*c=0;
	*/

	if(b2ac<0){
		if(b2ac<-0.00001) return 0;
		else b2ac=0;
	}

	int intx=0;
	if(b2ac<0.00001){		/* d[2]==0 */
		t = -b/a;
		if(t>=0 && t<=1){ /* between e1 and e2 */
			float r2;
			float ft2;
			itx[0] = e[0]*t + e1[0];
			itx[1] = e[1]*t + e1[1];
			itx[2] = e[2]*t + e1[2];

			r2 = dot2(itx,itx);
			a = dot2(d,d);
			b = dot2(d,d1);
			c = dot2(d1,d1)-r2;
			b2ac = b*b-a*c;
			if(!b2ac){
				float td = -b/a;
				if(td>=0 && td<=1){
					mad2(v,d,td,d1);
					if(cross2(v,itx)>=0){ /* <PI */
						ft2 = 1.0f + dot2(itx,v)/r2;
					}
					else{
						ft2 = -1.0f - dot2(itx,v)/r2;
					}
					if(ft2 > *ft){
						*ft = ft2;
						intx = 1;
					}
				}
			}
			else{
				float sqrtb2ac = sqrt(b2ac);
				float td[2] = {(-b+sqrtb2ac)/a, (-b-sqrtb2ac)/a};
				for(int i=0; i<2; i++){
					if(td[i]>=0 && td[i]<=1){
						mad2(v,d,td[i],d1);
						if(cross2(v,itx)>=0){
							ft2 = 1.0f + dot2(itx,v)/r2;
						}
						else{
							ft2 = -1.0f - dot2(itx,v)/r2;
						}
						if(ft2 > *ft){
							*ft = ft2;
							intx = 1;
						}
					}
				}
			}
		}
	}
	else{
		float sqrtb2ac = sqrt(b2ac);
		t1 = (-b+sqrtb2ac)/a;
/*  		assert(t1);  */
		if(t1>=0 && t1<=1){
			float td;
			i1[2] = e[2]*t1 + e1[2];
			td = (i1[2] - d1[2]) / d[2];

			if(td>=0 && td<=1){ /* i1 is in between d1 and d2 */
				printf("a=%f, b=%f, c=%f, t1=%f, td=%f\n",a,b,c,t1, td);
				float v[2];
				float ft2;

				mad2(i1,e,t1,e1); /* i1 = t1*e + e1 */
				mad2(v,d,td,d1);  /* v = td*d + d1 */
				float r2 = 0.5f*(dot2(i1,i1)+dot2(v,v));
				if(cross2(v,i1)>=0){ /* <2*PI */
					ft2 = 1.0f + dot2(i1,v)/r2;
				}
				else{
					ft2 = -1.0f - dot2(i1,v)/r2;
				}

				printf("dot2(v,v)=%f, dot2(i1,i1)=%f ft2=%f, ft=%f JJJ\n", dot2(v,v),r2, ft2, *ft);

				if(ft2 > *ft){ /* found ealier hit point */
					*ft = ft2;
					intx = 1;
					itx[0] = i1[0];
					itx[1] = i1[1];
					itx[2] = i1[2];
				}
			}
		}

		t2 = (-b-sqrtb2ac)/a;
		if(t2>=0 && t2<=1){
			float td;
			i2[2] = e[2]*t2 + e1[2];
			td = (i2[2] - d1[2]) / d[2];

			if(td>=0 && td<=1){ /* i1 is in between d1 and d2 */
				float v[2];
				float ft2;

				mad2(i2,e,t2,e1);
				mad2(v,d,td,d1);
				float r2 = .5f*(dot2(i2,i2)+dot2(v,v)); /* dot2(i2,v)/dot2(i2,i2) may > 1 */
				if(cross2(v,i2)>=0){ /* <2*PI */
					ft2 = 1.0f + dot2(i2,v)/r2;
				}
				else{
					ft2 = -1.0f - dot2(i2,v)/r2;
				}
				if(ft2 > *ft){ /* found ealier hit point */
					*ft = ft2;
					intx = 1;
					itx[0] = i2[0];
					itx[1] = i2[1];
					itx[2] = i2[2];
					printf("KKK\n");
				}
			}
		}
	}

 	if(intx){
		if(*ft>2.0f){
/* 			*ft=2.0f; */
			printf("ft=%f, d[2]=%f\n", *ft, d[2]);
		}
		assert(*ft<=2.0); 
	}
	return intx;
}


void testee()
{
#if 0
	v3_t v0={-2.3, .5, 0.1},v1={-0.5, .5, 0.1};
	v3_t obs0={-1.5, .7, -0.5},obs1={3, .9, 0.3};
	rot_t rot;
	rot.axis=vec(0,0,1);
	rot.org=vec(0,0,0);
	rot.deg=250;
	rot.rad=DEG2RAD(rot.deg);
	float ft= rot.rad<PI ? 1.0f + cosf(rot.rad) : -1.0f - cosf(rot.rad);
	v3_t itx;
	int intx;

/*
	float ft0 = ft;
	for(int i=0; i<10000000; i++){
		ft = ft0;
		intx = rotee(&v0,&v1,&obs0,&obs1,&ft,&itx);
	}
	extern int done;
	done=1;
*/

	intx = rotee(&v0,&v1,&obs0,&obs1,&ft,&itx);
	float c,rad;
	if(ft>=0) rad=acosf(ft-1.0f);
	else rad=2*PI-acosf(-1.0f-ft);
	rot.rad = rad;
	float rotmat[3][3];
	rot2mat(&rot, rotmat);
	v3_t w0,w1;
	w0 = matv3(rotmat, v0);
	w1 = matv3(rotmat, v1);

	glLineWidth(5);
	glBegin(GL_LINES);
	glColor3f(0,1,0);
	glVertex3fv(&v0);
	glVertex3fv(&v1);
	glColor3f(1,0,0);
	glVertex3fv(&obs0);
	glVertex3fv(&obs1);
	glColor3f(0,0,1);
	glVertex3fv(&w0);
	glVertex3fv(&w1);
	glEnd();
	glLineWidth(1);

	glColor3f(1,1,1);
	glPointSize(5);
	glBegin(GL_POINTS);
	glVertex3fv(&itx);
	glEnd();
	glPointSize(1);

	printf("intx=%i, ft=%f\n",intx,ft);
#endif
}


static inline void rotpln2z(float zmat[3][3], v3_t *org, plane_t *pln)
{
	v3_t v = matv3(zmat, sub(scale(pln->n, pln->d), *org));
	pln->n = matv3(zmat, pln->n);
	pln->d = dot(pln->n,v);
}


extern int globrebuildpntbnd;
void debug_drawrotpnt_n_bnd()
{
	static rot_t rot;
	v3_t p0, p1;
	float rotmat[3][3];
	float rotbnd[2][3];
	float zmat[3][3];
	v3_t v;
	v3_t intx;
	testee();

	if(globrebuildpntbnd){
		rot.axis.x = rand();
		rot.axis.y = rand();
		rot.axis.z = rand();
		rot.axis = normalize(rot.axis);
		rot.org = vec(0,0,0);
		rot.deg = rand()%360;
		rot.rad = DEG2RAD(rot.deg);
		globrebuildpntbnd=0;
	}

#if 0
	if(globrebuildpntbnd){
		rot.axis = normalize(vec(1,0,1));
		rot.org = vec(0,0,0);
		rot.deg = rand()%360;
		rot.rad = DEG2RAD(rot.deg);
		globrebuildpntbnd=0;
	}		
#endif

	p0.x = p0.y = p0.z = 1;
	
	v3_t proj = scale(rot.axis, dot(p0,rot.axis));
	v3_t rotorg = add(rot.org, proj);
	v3_t r,u;
	r = sub(p0, rotorg);
	u = cross(rot.axis,r);

 	rot2mat(&rot, rotmat); 
  	p1 = add(rot.org, matv3(rotmat, sub(p0,rot.org)));  

/*  	p1 = add(rotorg, add(scale(r,cosf(rot.rad)),scale(u,sinf(rot.rad))));  */
	rotpntbnd(&rot, p0, p1, rotbnd);
	v = add(rot.org, rot.axis);
	
	v3_t z0,z1;
	rot2zmat(&rot,zmat);
	z0=matv3(zmat, sub(p0,rot.org));
	z1=matv3(zmat, sub(p1,rot.org));
	v3_t oldz0=z0,oldz1=z1;
	plane_t pln={0,1,0,.5};
	rotpln2z(zmat, &rot.org, &pln);

	float ft = rot.rad>PI ? -1.0f-cosf(rot.rad) : 1.0f+cosf(rot.rad);

	int hit;
#if 0
	float ft0=ft;
	for(int i=0; i<10000000;i++){
		ft=ft0;
		hit=rotvpln(&z0, &pln, &ft, &z1);
	}
	extern int done;
	done=1;
#endif

/* 	printf("ft=%f\n",ft); */

	hit=rotvpln(&z0, &pln, &ft, &z1);
	if(hit){
		z1=add(rot.org,mattv3(zmat, z1));
/* 		printf("hit=(%f %f %f)\n",z1.x,z1.y,z1.z); */
	}
	z0=matv3(zmat,rot.axis);

	/* draw axis */
	glBegin(GL_LINES);
	glColor3f(1,0,0);
	glVertex3fv(&rot.org);
	glColor3f(0,0,1);
	glVertex3fv(&v);
	glEnd();

	/* draw trace */
	glColor3f(0,0,.2);
	glBegin(GL_LINE_STRIP);
	for(float rad=0; rad<rot.rad; rad+=0.02){
		v3_t p;
		p = add(rotorg, add(scale(r,cosf(rad)),scale(u,sinf(rad))));
		glVertex3fv(&p);
	}
	glEnd();

	/* draw points */
	glPointSize(5);
	glBegin(GL_POINTS);
	glColor3f(0,1,1);
	glVertex3fv(&p0);
	glColor3f(0,1,0);
	glVertex3fv(&p1);
	if(hit){
		glColor3f(1,0,0);
		glVertex3fv(&z1);
	}

	glColor3f(0,0,0);
	glVertex3fv(&oldz0);
	glColor3f(1,0,1);
	glVertex3fv(&oldz1);

	glEnd();
	glPointSize(1);

	/* draw bounds */
	float bndcolor[3] = {0,0,1};
	drawbox(rotbnd, bndcolor);

	glColor3f(1,1,1);
	glBegin(GL_LINES);
	glVertex3f(0,0,0);
	glVertex3fv(&z0);
	glEnd();

}


int rotobj(rot_t *rot, obj_t *obj, sctr_t *sctr, ccdinfo_t *ccdinfo)
{
	/*
	  1. find all obstacles (brush/obj)
	  2. transform the object and obstacles into a collision coordinate system
	     (1) translate center of rotation to (0,0,0)
	     (2) rotate rotation axis to (0,0,1)

	 */


#if 0
	v3_t objverts[64];
	v3_t obsverts[64];
	float zmat[3][3];
	float rotbnd[2][3];
	int n;

	/* 1. */
	rotationalbnd(rot, rotbnd, obj->phd);
	n = rngqry(sctr->sp, rotbnd, meh_touchers);	
	if(!n){
		return 0;
	}

	/* 2. */
	rot2zmat(rot->axis, zmat);
	for(int i=0; i<obj->phd->nverts; i++){
		objverts[i] = matv3(zmat, sub(obj->phd->verts[i], rot->org));
	}
	for(int i=0; i<n; i++){
		stuff_t *stuff = SP_GETOWNER(sctr->sp, meh_touchers[i]);
		phd_t *phd = stuff->common.phd;
		for(int j=0; j<phd->nverts; j++){
			obsverts[j] = matv3(zmat, sub(phd->verts[i], rot->org));
		}
	}
#endif

}

static struct{
	rot_t *rot0;
	phd_t *phd;
	rot_t rot;
	float zmat[3][3];
	float rotmat[3][3];
	v3_t zverts[32];
	float vbnds[6*32];	/* bnd of verts */
	float pbnd[2][3];	/* bnd of phd */
	float ft;
	v3_t intx;
	ccdinfo_t *ccd;
	int stuck;
} meh_rp;

#define RADFUNC(rad) (((rad)>PI ? -1.0f-cosf(rad) : 1.0f+cosf(rad)))

static void beginrotphd(rot_t *rot, phd_t *phd, ccdinfo_t *ccd)
{
	meh_rp.phd = phd;
	meh_rp.rot0 = rot;
	meh_rp.ccd = ccd;
	meh_rp.stuck=0;

	if(rot->rad<0){
		meh_rp.rot.rad = -meh_rp.rot.rad;
		meh_rp.rot.deg = -meh_rp.rot.deg;
		meh_rp.rot.axis = scale(rot->axis, -1);
		meh_rp.rot.org = rot->org;
	}
	else{
		meh_rp.rot = *rot;
	}

	rot2zmat(&meh_rp.rot, meh_rp.zmat);
	rot2mat(&meh_rp.rot, meh_rp.rotmat);
	for(int i=0; i<phd->nverts; i++){
		v3_t r, p;
		r = sub(phd->verts[i], rot->org);
		p = add(rot->org, matv3(meh_rp.rotmat, r));
		meh_rp.zverts[i] = matv3(meh_rp.zmat, r);
		rotpntbnd(&meh_rp.rot, phd->verts[i], p, meh_rp.vbnds+6*i);
	}

	meh_rp.pbnd[0][0] = meh_rp.pbnd[0][1] = meh_rp.pbnd[0][2] = 1e10f;
	meh_rp.pbnd[1][0] = meh_rp.pbnd[1][1] = meh_rp.pbnd[1][2] = -1e10f;
	for(int i=0; i<6*phd->nverts; i+=6){
		float *bnd=meh_rp.vbnds+i;
		for(int j=0; j<3; j++){
			if(bnd[j]<meh_rp.pbnd[0][j]){ meh_rp.pbnd[0][j] = bnd[j];}
			if(bnd[j+3]>meh_rp.pbnd[1][j]){ meh_rp.pbnd[1][j] = bnd[j+3];}
		}
	}

	if(meh_rp.rot.rad>2*PI){
		meh_rp.ft = -2.0f;
	}
	else{
		meh_rp.ft = RADFUNC(meh_rp.rot.rad);
	}
}


int rotagstphd(phd_t *phd)
{
	/*
	  Three types of tests
	  -- phd1->vert vs. phd2->poly
	  -- phd1->poly vs. phd2->vert
	  -- phd1->edge vs. phd2->edge

	  Mid-phase:
	*/

	/* brute force */

	/* 1 */
	char bbitx[32];
	char bbitx2[32];
	plane_t edgeplns[16];
	v3_t zverts[32];
	int coll=0;

	for(int i=0; i<meh_rp.phd->nverts; i++){
		bbitx[i] = bbintx(meh_rp.vbnds+6*i, phd->bnd);
	}

	for(int i=0; i<phd->nverts; i++){
		zverts[i] = matv3(meh_rp.zmat, sub(phd->verts[i], meh_rp.rot.org));
	}
#if 1
	for(int i=0; i<phd->npolys; i++){
		poly_t *p = phd->polys + i;
		float pbnd[2][3];
		plane_t pln = p->plane;
		rotpln2z(meh_rp.zmat, &meh_rp.rot.org, &pln);
		poly2bnd(phd, i, pbnd);

#if EDGECACHE
#else
		for(int j=0; j<p->nv; j++){
			v3_t *v,*v2;
			v = zverts + p->v[j];
			v2 = zverts + p->v[(j+1)%p->nv];
			edgeplns[j].n = cross(pln.n, sub(*v2, *v));
			edgeplns[j].d = dot(edgeplns[j].n, *v);
		}
#endif

		for(int j=0; j<meh_rp.phd->nverts; j++){
			if(!bbitx[j]) continue;
			if(!bbintx(meh_rp.vbnds+6*j, pbnd)) continue;
			v3_t intx,intx2;
			float ft = meh_rp.ft;
/*
			printf("damn dot(pln.n,v)=%f, pln.d=%f, pln.n=(%f %f %f), v0=(%f, %f, %f)\n",
			       dot(pln.n, meh_rp.zverts[j]), pln.d, pln.n.x, pln.n.y, pln.n.z,
			       meh_rp.phd->verts[j].x, meh_rp.phd->verts[j].y, meh_rp.phd->verts[j].z);
*/
			int stuck = rotvpln(meh_rp.zverts+j, &pln, &ft, &intx);
/* 			printf("j=%i stuck=%i KKKKKKKKKKKKKKKKKKKKKKKK\n",j,stuck); */


			if(stuck){


#if EDGECACHE
  				intx2 = add(meh_rp.rot.org, mattv3(meh_rp.zmat, intx));    
  				if(pntinpoly(&intx2, phd, i)){  
#else
  				if(pntinpoly(&intx, edgeplns, p->nv)){  
#endif
					meh_rp.ft = ft;
					meh_rp.intx = intx; /* in z space */
					coll=1;
				}
			}
		}
	}
#endif

#if 1
	/* 2 */
	for(int i=0; i<phd->nverts; i++){
		bbitx[i] = pntinbnd(phd->verts+i, meh_rp.pbnd);
	}
	for(int i=0; i<meh_rp.phd->npolys; i++){
		poly_t *p = meh_rp.phd->polys + i;
		plane_t pln = p->plane;
		rotpln2z(meh_rp.zmat, &meh_rp.rot.org, &pln);

#if EDGECACHE
#else
		for(int j=0; j<p->nv; j++){ 
			v3_t *v,*v2;
			v = meh_rp.zverts + p->v[j];
			v2 = meh_rp.zverts + p->v[(j+1)%p->nv];
			edgeplns[j].n = cross(pln.n, sub(*v2, *v));
			edgeplns[j].d = dot(edgeplns[j].n, *v);
		}
#endif

		pln.n.x = -pln.n.x;
		
		for(int j=0; j<phd->nverts; j++){
			if(!bbitx[j]) continue;
			v3_t intx, intx2;
			float ft = meh_rp.ft;
			v3_t v = zverts[j];
			v.x = -v.x;
			int stuck = rotvpln(&v, &pln, &ft, &intx);
			if(stuck){
				intx.x = -intx.x;
#if EDGECACHE
 				intx2 = add(meh_rp.rot.org, mattv3(meh_rp.zmat, intx));   
 				if(pntinpoly(&intx2, meh_xp.phd, i)){ 
#else
  				if(pntinpoly(&intx, edgeplns, p->nv)){  
#endif
					meh_rp.ft = ft;
					meh_rp.intx = intx;
					coll=1;
				}
			}
		}
	}
#endif

	/* 3 */
	for(int i=0; i<meh_rp.phd->nedges; i++){
		float *d1, *d2, *e1, *e2;
		int *v;
		v = meh_rp.phd->edges[i].v;
		d1 = meh_rp.zverts + v[0];
		d2 = meh_rp.zverts + v[1];
		for(int j=0; j<phd->nedges; j++){
			v3_t intx;
			float ft = meh_rp.ft;
			v = phd->edges[j].v;
			e1 = zverts + v[0];
			e2 = zverts + v[1];
			int stuck = rotee(d1,d2,e1,e2,&ft,&intx);
			if(stuck){
				printf("yeah KKKKKKKKK\n");


/*  				if(ft>meh_rp.ft){  */
					meh_rp.ft = ft;
					meh_rp.intx = intx;
/*  				}  */
				coll=1;
			}
		}
	}
	meh_rp.stuck += coll;
	return coll;
}


int endrotphd()
{
#if 0
	if(meh_rp.ft>=2){
/* 		meh_rp.ccd->intx = vec(0,0,0); */
		meh_rp.rot0->rad = 0;
		return 1;
	}
	if(meh_rp.ft<=-2){
		meh_rp.rot0->rad = meh_rp.rot.rad;
		return 0;
	}
#endif
#if 0
	if(meh_rp.ft>2.0f) meh_rp.ft=2.0f;
	else if(meh_rp.ft<-2.0f) meh_rp.ft=-2.0f;
#endif

	if(meh_rp.stuck){	
		meh_rp.rot0->rad = meh_rp.ft>0 ? acosf(meh_rp.ft - 1.0f) : 2*PI - acosf(-meh_rp.ft - 1.0f);
		v3_t r = mattv3(meh_rp.zmat, meh_rp.intx);
 		meh_rp.ccd->intx = add(meh_rp.rot.org, r);   

		float l = meh_rp.rot0->rad * length(r);
		meh_rp.rot0->rad *= (l - CONTACT_EPSILON)/ l;

/* 		meh_rp.ccd->intx = meh_rp.intx; */

/*  		meh_rp.rot0->rad *= 0.998;  */




/*  		meh_rp.rot0->rad -= 0.0003; */
		if(meh_rp.rot0->rad<0) meh_rp.rot0->rad = 0;
	}
	else{
		return 0;
	}

	return meh_rp.stuck;
}


void testrotphd()
{
	static float deg=1;
	extern phd_t *globphds[];
	rot_t rot;
	ccdinfo_t ccd;
	rot.org = bndcntr(globphds[1]->bnd);
	rot.axis = normalize(vec(.5,1.3, -.3));
	rot.deg = deg;
	rot.rad = DEG2RAD(rot.deg);
	printf("rad=%f\n",rot.rad);

	beginrotphd(&rot, globphds[1], &ccd);
/* 	static int s=0; */
	static v3_t intx;
	rotagstphd(globphds[0]);
	int stuck = endrotphd();

/* 	if(rot.rad) rotphd(globphds[1], &rot); */

#if 1
/* 	if(!s) */{
		if(stuck){
			deg=0;
/* 			s=1; */
/* 			rot.rad -= 0.0003; */
/* 			if(rot.rad<0) rot.rad=0; */

			intx = ccd.intx;/*  add(meh_rp.rot.org, mattv3(meh_rp.zmat, ccd.intx)); */
			printf("shit deg=%f, ft=%f, intx=(%f %f %f), rot.org=(%f %f %f)\n", 
			       RAD2DEG(meh_rp.rot.rad), meh_rp.ft, intx.x, intx.y, intx.z, rot.org.x, rot.org.y, rot.org.z);

		}
	}

	if(rot.rad){
		rotphd(globphds[1], &rot);
	}

	printf("stuck=%i, rot.rad=%f\n", stuck, rot.rad);
#endif

/* 	printf("intx=(%f %f %f)\n", intx.x, intx.y, intx.z); */

#if 0
	if(stuck){
		if(rot.rad>0){
			rotphd(globphds[1], &rot);
			intx = ccd.intx;
		}
	}
#endif

#if 1

	glLineWidth(1);
	glColor3f(1,0,0);
	drawphd(globphds[0]);
	glColor3f(0,0,1);
	drawphd(globphds[1]);
	glLineWidth(1);

	glPointSize(5);
	glColor3f(0,1,1);
	glBegin(GL_POINTS);
	glVertex3fv(&intx);
	glEnd();
	glPointSize(1);
#endif

#if 0
	glPointSize(2);
	glColor3f(1,.5,.5);
	glBegin(GL_POINTS);
	glVertex3fv(&ii);
	glEnd();
#endif
	

#if 0
	float color[3] = {.2,.2,1};
	drawbox(meh_rp.pbnd, color);
#endif


#if 0
	rot_t rot2;
	v3_t pnt = {1,0,0};
	rot2.deg = 50;
	rot2.axis = vec(0,0,1);
	rot2.rad = DEG2RAD(rot2.deg);
	rot2.org = vec(0,0,0);
	float mat[3][3];
	rot2mat(&rot2, mat);
	v3_t pnt2 = matv3(mat, pnt);

	rot2.deg = 210;
	rot2.rad = DEG2RAD(rot2.deg);
	rot2mat(&rot2, mat);
	v3_t pnt3 = matv3(mat, pnt);

#endif

#if 0
	rot_t rot2;
	v3_t pnt = {1,0,0};
	rot2.deg = 100;
	rot2.axis = vec(0,0,1);
	rot2.rad = DEG2RAD(rot2.deg);
	rot2.org = vec(0,0,0);
	float mat[3][3];
	rot2mat(&rot2, mat);
	v3_t pnt2 = matv3(mat, pnt);
#endif

#if 0
	glLineWidth(5);
	glColor3f(.6, .3, .9);
	glBegin(GL_LINE_STRIP);
	glVertex3fv(&pnt);
	glVertex3fv(&pnt2);
	glVertex3fv(&pnt3);
	glEnd();
	glLineWidth(1);
#endif

}


/*
  find contact points between two phds
  
  just another translational ccd?
  or some static checking? 

  which one is faster?

  this test shouldn't be O(n^2)

  for a phd, three contact points are enough.


  or, use persistent contact manifold?

 */



int addlight()
{
}

int rmlight()
{
}

int xlatlight()
{
}


/* which sector the point is in? */
sctr_t *inwhichsctr(v3_t pnt)
{
	int n=pntqry(map.sp, pnt, meh_touchers);
	for(int i=0;i<n;i++){
		sctr_t *s = SP_GETOWNER(map.sp, meh_touchers[i]);
		if(pntinphd(pnt, s->phd)) return s;
	}

	return NULL;
}

static void drawsctr(sctr_t *s)
{
	for(int i=0;i<s->nwalls;i++){
		wall_t *w = s->walls + i;
		for(brush_t *b=w->brushlist; b; b=b->next){
			drawphd(b->phd);
		}
	}
	for(brush_t *b=s->brushlist; b; b=b->next){
		drawphd(b->phd);
	}
	glLineWidth(3);
	for(int i=0;i<s->nwalls;i++){
		wall_t *w = s->walls + i;
		for(metaptl_t *mp=w->mptllist; mp; mp=mp->next){
			glColor3f(0,1,0);
			glBegin(GL_POLYGON);
			for(int j=0;j<mp->nverts;j++){
				glVertex3fv(mp->verts+j);
			}
			glEnd();
		}
	}
	glLineWidth(1);
}


static void drawmap()
{
	for(sctr_t *s=map.sctrlist; s; s=s->next){
		drawphd(s->phd);
	}

	glLineWidth(3);
	for(sctr_t *s=map.sctrlist; s; s=s->next){
		for(int i=0;i<s->nwalls;i++){
			wall_t *w = s->walls + i;
			for(metaptl_t *mp=w->mptllist; mp; mp=mp->next){
				glColor3f(0,1,0);
				glBegin(GL_POLYGON);
				for(int j=0;j<mp->nverts;j++){
					glVertex3fv(mp->verts+j);
				}
				glEnd();
			}
		}
	}
	glLineWidth(1);

}


static void integrate()
{
}

/* simulate the object */
int sim(obj_t *obj)
{
}

#if 0
/* take a snapshot of the map, regardless of whether or not other things like physics are finished. */
int snapshot(pnt_t *p, v3_t *dir)
{
}
#endif

void render(cam_t *cam)
{
	sctr_t *s = inwhichsctr(cam->pos);
	if(s){
		glColor3f(0,0,0);
		drawsctr(s);
	}
	else{
		glColor3f(1,1,1);
		drawmap();
	}

}


int loadsctr()
{
}

int stosctr()
{
}


/* create a new empty map (the previous map is ) */
void initmap()
{
	if(map.sp){
		freesp(map.sp);
		/* remove sectors */
	}

	map.sp = allocsp(32);
	map.sctrlist = NULL;
	map.nsctrs = 0;
}


void shutdownmap()
{
	if(map.sp){
		freesp(map.sp);
		rmallsctrs();
	}
}


/* 
   the sectors are stored in one file, separated from the content in each sector.
   sectors, portals and light volumes are all loaded to RAM; brushes, objects and textures are streamed in later when needed.
 */
int loadmap()
{
}



sctr_t *debug_getsctr(int id)
{
	if(id>=map.nsctrs) return NULL;
	sctr_t *s=map.sctrlist;
	while(id){
		s=s->next;
		id--;
	}
	return s;
}
