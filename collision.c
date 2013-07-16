#include "collision.h"


// Sweep and Prune
#define MULTI_SP 0

#ifdef MULTI_SP
/*
  Spatial Query

  To simplify thing a bit, only one SP operation can be executed at a time.
  Num of object operations are not limited.
 */
#include "types.h"

// An query object is something you use to communicate with the spatial query system.
// You should add a query_t typed variable to your struct if you want to add it to the spatial database

/*
  The design of query object:
  
  Each object has:
  () a list of sub-SPs it's in
  () a list of overlapping objects

  The sub-SP doesn't store a list of objects in it -- it stores their AABBs which store pointers to them.
  A pair is two linked list nodes
 */

// handle to sub-SP -- one sub-SP may have multiple instances (used by query object)
typedef struct sp_h_s{
	struct sp_h_s *prev, *next;
	sp_t *sp;
	ushort bid; // index to the bound in SP
	ushort id; // SP's id
} sp_h_t;

typedef struct{
	struct queryobj_h_s *ol;	// list of objects that touch this object
	sp_h_t *sl;			// list of SPs that touch this object
	int flag;
//	ushort ntouchobj;
	uchar type;
	float bnd[6];
} queryobj_t;


typedef struct queryobj_h_s{
	struct queryobj_h_s *next, *twin;
	queryobj_t *obj;
} queryobj_h_t;


/*
  Each multi-SP has a number of sub-SPs inside. You can add/remove/move/reshape them freely.
  Each sector or something in the engine stores a pointer to the corresponding sub-SP.

  Six types of operation are allowed:
  Add/remove/reshape-move sub-SPs;
  Add/remove/reshape-move query objects;


 */


/*
  Each SP is allocated/freed by malloc/free.
 */

#define ITER(n) for(int i=0; i<n; i++)
#define ITER1(n) for(int j=0; i<n; i++)

#define MAX_NUM_SP 4096
#define INVALID_SP_ID 0xffff

enum {SPOP_NONE, SPOP_ADD, SPOP_DEL, SPOP_MOVE};
/*
  0 --------------------------------- 4096-1
  first(top) <------------------------last(bottum)
   /\
   |
  pop 
 */
typedef struct{
	sp_t			*sp_ptrs[64*64+MAX_NUM_SP];
	ushort			idstack[MAX_NUM_SP]; // stores unused ids
	ushort			top; // first non-empty slot
	ushort			ns;
} sppool_t;


static void initsppool(sppool_t *pool)
{
	memset(pool->sp_ptrs, NULL, sizeof(pool->sp_ptrs));

	// pre-fill sppool_t::sp_ptrs[0 ~ 64*64-1] with SP-grid
	int n=0;
	ITER(64){
		ITER1(64){
			world.spgrid[i][j].sp.id = n;
			pool->sp_ptrs[n++] = &world.spgrid[i][j].sp;
		}
	}

	ITER(MAX_NUM_SP) pool->idstack[i]=i+64*64;

	pool->top = 0;
	pool->ns = 0;
}

static void shutdownsppool(sppool_t *pool)
{
//	for(int i=0; i<MAX_NUM_SP; i++)
	ITER(MAX_NUM_SP){
		if(pool->sp_ptrs[i])
			free(pool->sp_ptrs[i]);
	}
}



// internal use only -- called by flushsp()
static sp_t *allocsp(sppool_t *pool, ushort *id)
{
	if(pool->top >= MAX_NUM_SP) return INVALID_SP_ID;

	int i=pool->idstack[pool->top];
	assert(i<MAX_NUM_SP && !pool->sp_ptrs[i]);

	sp_t *sp;
	sp = pool->sp_ptrs[i] = malloc(sizeof(sp_t));
	pool->top++;
	pool->ns++;
	sp->id = *id = i;

	return sp;
}


// internal use only -- called by flushsp()
static void freesp(sppool_t *pool, ushort id)
{
	if(id >= MAX_NUM_SP) return;
	if(!pool->sp_ptrs[id]) return;

	free(pool->sp_ptrs[id]);
	pool->sp_ptrs[id] = NULL;
	pool->idstack[--pool->top] = id;
	poll->ns--;
}


/*
  Query objects

  There're many of them (32k)
 */

#define MAX_NUM_OBJ 32768
#define INVALID_OBJ_ID 0xffff

typedef struct{
	queryobj_t		obj[MAX_NUM_OBJ]; // 
	ushort			idstack[MAX_NUM_OBJ]; // 64k
	uint			bitmap[MAX_NUM_OBJ/32]; // 4k   bit==0 => unused
	int			top;
} objpool_t;

static void initobjpool(objpool_t *pool)
{
	memset(pool->bitmap, 0, sizeof(pool->bitmap));

//	for(int i=0; i<MAX_NUM_OBJ; i++) pool->idstack[i] = i;	

	ITER(MAX_NUM_OBJ) pool->idstack[i]=i;
	pool->top=0;
}


static void setbit(uint bitmap[], ushort bit)
{
	ushort i=bit/32;
	bitmap[i] |= 1<<(bit%32);
}

static void clrbit(uint bitmap[], ushort bit)
{
	ushort i=bit/32;
	bitmap[i] &= ~(1<<(bit%32));
}

static int bitis0(uint bitmap[], ushort bit)
{
	return !(bitmap[bit/32] & (1<<(bit%32)))
}

// internal use only -- called by flushsp()
// both pointer and index are returned => pointer can be used for direct access, index can be used for sorting
static queryobj_t *allocobj(objpool_t *pool, ushort *id)
{
	if(pool->top >= MAX_NUM_OBJ) return INVALID_OBJ_ID;
	ushort i=pool->idstack[pool->top];
	setbit(pool->bitmap, i);
	pool->top++;
	*id = i;
	return pool->obj+i;
}

// internal use only -- called by flushsp()
static void freeobj(objpool_t *pool, ushort id)
{
	if(id >= MAX_NUM_OBJ) return;
	if(bitis0(pool->bitmap, id)) return;
	pool->idstack[--pool->top] = id;
	clrbit(pool->bitmap, id);
}


//#define DEFSTACK(type,name) struct{type *buf; int n; int cap; int elementsize;} name;



/**
   Dynamic Array
    _____________
   [_|_________|_]
   -1 0         n-1
    |           |
   \|/         \|/
   num         sentinel (e.g. 0xffff)
   of
   items


 */

#if 0
typedef struct{
//	sppool_t sppool;
//	objpool_t objpool;

	sp_t *splist;
	int principleaxis;

//	pm_t *pm;

//	pool_t sp_h_pool;
//	pool_t obj_h_pool;

	struct {
		struct {
			int type;
			ushort id;
			float bnd[6];
		} sp;
		
		struct {
			// for ushort stack, uint &stack[-2] stores num of elements


			stacku16_t *d, *m, *a; // storing object ids is enough -- other info is stored in query object.
			// new/moved bounds are stored in corresponding objects (endpnt[] and bnd[])
			/*
			  using __SENTINEL__ value to tag the end of the array (in this case 0xffff).
			  if the array is full, realloc twice the memory for the array
			 */

		} obj;
	} op;
		
	// delete object from SP
	// delspbid_buf[] should only contain (bound_id + SP_id)
	// delobj_buf[] is used to store objcects being deleted from some SP

	struct {
		stacku32_t *bid_sp;
		stacku32_t *obj_sp;
		stacku32_t *obj_obj; // same SP
		stacku32_t *obj_jbo; // obj and jbo come from different SPs, hence "jbo"
	
		stacku16_t *mend;
		stacku16_t *nearsp;

		uint sprank[MAX_NUM_SP];
		ushort spid[MAX_NUM_SP];
		float spmin[MAX_NUM_SP];
		queryobj_t *q[MAX_NUM_OBJ]; // stores result of proximity query

		uchar prtab[8][8];
	} lut; // look up tables


	// each MSP has its own coordinate system
	
} msp_t;
#endif


typedef struct msp_s{
	struct msp_s		*prev, *next;
	sp_t			*splist;
	int			pax; //principleaxis;
	aabbnode_t		*aabbtree;	
} msp_t;


struct {
	/* Memory */
	sppool_t sppool;
	objpool_t objpool;
	
	pm_t *pm;

	pool_t sp_h_pool;
	pool_t obj_h_pool;


	/* Look up tables */
	struct {
		stacku32_t *bid_sp;
		stacku32_t *obj_sp;
		stacku32_t *obj_obj; // same SP
		stacku32_t *obj_jbo; // obj and jbo come from different SPs, hence "jbo"
	
		stacku16_t *mend;
		stacku16_t *nearsp;

		uint sprank[MAX_NUM_SP];
		ushort spid[MAX_NUM_SP];
		float spmin[MAX_NUM_SP];
		queryobj_t *q[MAX_NUM_OBJ]; // stores result of proximity query

		uchar prtab[8][8];
	} lut; // look up tables


	/* Operations */
	// add/del/move an MSP
	/*
	  When an MSP is moved, all objects in it is deactivated (and won't be automatically activated after the operation --
	  objects will float in the air until hit by something. )

	 */


	/*
	  You can move multiple MSPs in a frame (in single player mode you can move only one MSP at a time).
	  
	  or

	  For each MSP:
	    you can move multiple SPs

	  and

	  You can move any object at any time


	 */



	/*
	  


	 */


	// add/del/move an SP in MSP

	// add/del/move objects


	struct {
		struct {
			ushort id;
		} msp;

		struct {
			int type;
			ushort id;
			float bnd[6];
		} sp;
		
		struct {
			// for ushort stack, uint &stack[-2] stores num of elements


			stacku16_t *d, *m, *a; // storing object ids is enough -- other info is stored in query object.
			// new/moved bounds are stored in corresponding objects (endpnt[] and bnd[])
			/*
			  using __SENTINEL__ value to tag the end of the array (in this case 0xffff).
			  if the array is full, realloc twice the memory for the array
			 */

		} obj;
	} op;


	
	/*
	  A (MGSP) multi-grid SP for objects that are not fully contained in some MSPs.

	  The game world is 32k*32k;
	  Each grid is 256*256 => 128^2 cells
	  Each grid is 512*512 => 64^2 cells --- yeah!


	 */


	struct {
		sp_t			sp;
		stacku32_t		*mspnodes; // nodes of MSPs' AABB trees
	} spgrid[64][64];


	/*
	  MSPs
	 */
	msp_t			*msplist;
	

	/*
	  32k*32k terrain
	 */

	/*
	  Virtual texturing 
	 */


} world;




#define SENTINEL_U16 0xffff
#define SENTINEL_U32 0xffffffff

static inline int isprvalid(ushort a, ushort b)
{
	queryobj_t *A=OBJPTR(a);
	queryobj_t *B=OBJPTR(b);
	return cmsp->lut.prtab[A->type][B->type];
}

static inline void pushu16(stacku16_t **stackptr, ushort dat)
{
	ushort *stack=*stackptr;
	uint *pn=stack-2;
	uint n = *pn;
	if(stack[n]==SENTINEL_U16){ //stack is full. enlarge it
		stack = *stackptr = realloc(stack-2, (n+1)<<2)+2; //(2+2*n)*sizeof(ushort));
		memset(stack+n, 0);
		stack[2*n-1]=SENTINEL_U16;
		pn=stack-2;
	}

	stack[n]=dat;
	(*pn)++;
}

static inline void pushu32(stacku32_t **stackptr, uint dat)
{
	uint *stack=*stackptr;
	uint *pn=stack-1;
	uint n = *pn;
	if(stack[n]==SENTINEL_U32){ //stack is full. enlarge it
		stack = *stackptr = realloc(stack-1, (n+1)<<2)+1; //(2+2*n)*sizeof(ushort));
		memset(stack+n, 0);
		stack[2*n-1]=SENTINEL_U32;
		pn=stack-1;
	}

	stack[n]=dat;
	(*pn)++;
}

static inline uint getnumu16(stacku16_t *stack)
{
	return *((uint*)(stack-1));
}

static inline uint getnumu32(stacku32_t *stack)
{
	return stack[-1];
}


static void resetstacku16(stacku16t *stack) // **stackptr to make the interface more consistent
{
	*((uint*)(stack-2))=0;
}

static void resetstacku32(stacku32t *stack) // **stackptr to make the interface more consistent
{
	*((uint*)(stack-1))=0;
}

#define OBJPTR(id) (cmsp->objpool.obj+(id))
#define SPPTR(id) (cmsp->sppool.sp_ptrs[id])

static msp_t *cmsp=NULL;

static inline void delprb_objjbo(ushort id1, ushort id2)
{
	uint pr= ( (id1|0x8000) <<16) | id2; // set sign bit to 1 => neg
	pushu32(&cmsp->lut.obj_jbo, pr);
}

static inline void addprb_objjbo(ushort id1, ushort id2)
{
	uint pr = PACK(id1,id2);
	pushu32(&cmsp->lut.obj_jbo, pr);
}

static inline void delprb_objobj(ushort id1, ushort id2)
{
	uint pr=( (id1|0x8000) <<16) | id2; // set sign bit to 1 => neg
	pushu32(&cmsp->lut.obj_obj, pr);
}

static inline void addprb_objobj(ushort id1, ushort id2)
{
	uint pr = PACK(id1,id2);
	pushu32(&cmsp->lut.obj_obj, pr);
}


/* whenever sp doesn't have enough mem for new bnds, just allocate! */
static void prepspace(sp_t *s, int nbnds)
{
	if(nbnds > s->nfree){
		int n = nbnds + s->nbnds;
		if(n > s->cap){ /* mem is too small, need enlargement */
			s->cap = nextp2(n);
			int size = 2*s->cap*sizeof(endpnt_t);
			s->endpnts[0] = realloc(s->endpnts[0], size);
			s->endpnts[1] = realloc(s->endpnts[1], size);
			s->endpnts[2] = realloc(s->endpnts[2], size);
			s->bndeps = realloc(s->bndeps, 6*s->cap*sizeof(ushort));
			s->bnds = realloc(s->bnds, s->cap*sizeof(bnd_t));
		}
		
		/*  | s->nbnds | s->nfree | s->cap - s->nbnds - s->nfree |  <--s->cap */
		
		int ntot = s->nbnds + s->nfree; /* total num of bnds allocated */
		int id=ntot*6;
		for(int i=2*(s->nbnds+s->nfree); i<2*(s->nbnds+nbnds); i+=2){ /* set values for new endpoints */
			s->endpnts[0][i].id = id; 
			id+=6;
		}
		s->nfree = 0;
	}
	else{
		s->nfree-=nbnds;
	}

}


//static void addbnds(sp_t *s, endpnt_t *endpnts[3], ushort *bndids[], void **owners, int nbnds)
static void addbnds(sp_t *s, ushort obj[], int no)
{
	prepspace(s, no);

	endpnt_t *endpnts[3];
	ITER(3){
		endpnts[i]=alloca(sizeof(2*endpnt_t*no));
		ITER1(no){
			int k=j*2;
			queryobj_t *o=OBJPTR(obj[j]);
			endpnts[k].value=o->bnd[i];
			endpnts[k+1].value=o->bnd[i+3];
		}
	}

	/* copy values of the available endpoints to endpoints of the newly added bounds */
	endpnt_t *p=s->endpnts[0]+2*s->nbnds;
	ITER(no){	
		ushort id=p->id;
		*bndids[i] = id/6;
		ITER1(6){
			endpnts[j/2][2*i+(j&1)].id = id++;
		}

		bnd_t *b=s->bnds+p->id/6;
		b->owner = owners[i];
		p+=2;
	}

	int np=2*no;
	/* sort the new endpoints, and update the links pointing to endpoints after sorting (to be used in the following step) */

	for(int a=0;a<3;a++){
		if (np>2) qsort(endpnts[a], np, sizeof(endpnt_t), cmpep);
		ITER(np){
			s->bndeps[endpnts[a][i].id] = i;
		}
	}

	// add pairs between the newly added bnds (self pruning through x axis)
	if (np>2) {
		ITER(np) {
			endpnt_t *min = endpnts[0] + i; /* endpnts[] are already sorted*/
			if (ISMIN(min->flag)) {
				endpnt_t *max = endpnts[0] + s->bndeps[min->id+1];
				bnd_t *b = s->bnds + min->id/6;
				for (endpnt_t *p=min+1; p<max; p++) { /* min<p<max */
					if (ISMIN(p->flag)/*  && sp->validpr(min->flag, p->flag) */) {
						/* endpnts[i] are sorted => comparing indices is enough */
						ushort *pid = s->bndeps+p->id; /* PLACE 1 */
						ushort *minid = s->bndeps+min->id;
						if (pid[2]<minid[3] && minid[2]<pid[3] &&
						    pid[4]<minid[5] && minid[4]<pid[5]) {
							bnd_t *b1=s->bnds + p->id/6;
							addprb_objjbo(b->oid, b1->oid);
						}
					}
				}
			}
		}
	}


	/* move the inserted bnds to the right positions */
	/* y,z axis */
	endpnt_t *newp,*oldp,*tail;
	for (int a=1;a<3;a++) {
		newp = endpnts[a]+np-1;
		oldp = s->endpnts[a]+2*s->nbnds-1;
		tail = oldp+2*nbnds;
		while (newp>=endpnts[a]) {
			if (oldp>=s->endpnts[a] && oldp->value>newp->value) {
				*tail = *oldp--;
			}
			else {
				*tail = *newp--;
			}
			s->bndeps[tail->id] = tail-s->endpnts[a];
			--tail;
		}
	}

	/* x axis, plus box pruning */
	newp = endpnts[0]+np-1;
	oldp = s->endpnts[0]+2*s->nbnds-1;
	tail = oldp+2*nbnds;
 	int stab=0; 

	endpnt_t *old0, *new0;
	old0=s->endpnts[0];
	new0=endpnts[0];
	// == insert a sorted array into another
	int newismin=0; //ISMIN(newp->flag);
	float newmin;
	float newmax=newp->value;
	while (newp>=new0) {
		if (oldp>=old0 && newp->value<oldp->value) { /* newp<oldp */
			if (ISMIN(oldp->flag)) {
				if (newismin /*&& sp->validpr(oldp->flag, newp->flag) */) {
					float oldmax=s->endpnts[0][s->bndeps[oldp->id+1]].value;
					if (oldmax<newmax) { // oldp2<newp2 => newp<oldp<oldp2<newp2
						ushort *oldid, *newid;
						oldid=s->bndeps+oldp->id;
						newid=s->bndeps+newp->id;
						if (oldid[2]<newid[3] && newid[2]<oldid[3] && 
						    oldid[4]<newid[5] && newid[4]<oldid[5]) {
							bnd_t *ob=s->bnds+oldp->id/6;
							bnd_t *nb=s->bnds+newp->id/6;
							addprb_objjbo(ob->oid, nb->oid);
						}
					}
				}
				stab--; 
			}
			else{
 				stab++; 
			}
			*tail=*oldp--;
		}
		else{ 		/* oldp<newp */
			if (stab) { /* the endpoint is contained in some intervals */
				int st=stab;
				endpnt_t *p=oldp;
				while (st) {
					if (ISMIN(p->flag)) {
						endpnt_t *p2=s->endpnts[0]+s->bndeps[p->id+1];
						if (newp->value <= p2->value) { /* p<newp<p2 */
							if (newismin || // p<newp<p2<newp2 or p<newp<newp2<p2
							    (newmin<p->value)) { // newp<p<newp2<p2
								ushort *pid, *newid;
								pid = s->bndeps + p->id/6*6;
								newid = s->bndeps + newp->id/6*6;
								if (pid[2]<newid[3] && newid[2]<pid[3] && 
								    pid[4]<newid[5] && newid[4]<pid[5]) {
									bnd_t *ob=s->bnds+p->id/6;
									bnd_t *nb=s->bnds+newp->id/6;
									addprb_objjbo(ob->oid, nb->oid);
								}
							}
							st--;
						}
					}
					p--;
				}
			}
			*tail=*newp--;
			if(ISMIN(newp->flag)){
				newismin=1;
				newmin=newp->value;
			}
			else{
				newismin=0;
				newmax=newp->value;
			}
		}
		s->bndeps[tail->id] = tail-s->endpnts[0];
		--tail;
	}
	
	s->nbnds+=nbnds;
}


/*
  No pair ops
*/
static void rmbnds(sp_t *s, ushort d[], ushort nd)
{
	// tag bounds to be deleted with 0xffff
	ushort first[3]={0xffff,0xffff,0xffff};
	ITER(nd) {	
		s->bnds[d[i]].flag = 0xffff;
		ushort *b=s->bndeps + 6*d[i];

		ITER(3) {
			if(b[j*2]<first[j]) first[j]=b[j*2];
		}
	}
	
	// shift bounds
	ITER(3){
		ushort offset=1;
		for(ushort j=first[i]+1; j<2*s->nbnds; j++){
			int p = s->endpnts[i][j].id;
			int b = p/6;
			if(s->bnds[b].flag != INVALID_ID){
				s->endpnts[i][j-offset] = s->endpnts[i][j];
				s->bndeps[p]=j-offset;
			}
			else{
				offset++;
			}
		}
	}


	/* link endpoints to unused bounds for later use. 
	   === now endpoints[i].id is used to store __bnd__ id, NOT index to the bndeps[] array ===
	   only endpoints[0][s->bndeps[6*d[i]]].id is set. the other five endpoints' .id field is not set.
	*/
	s->nbnds-=nd;
	s->nfree+=nd;
	endpnt_t *p=s->endpnts[0]+s->nbnds*2;
	for(ushort i=0;i<nd; i++,p+=2){
		p->id=d[i];
	}
}


static void updendpnt(sp_t *s, endpnt_t *pnt, float oldvalue)
{
	int ismax = !ISMIN(pnt->flag);
	int d = ((pnt->value > oldvalue)<<1)-1; /* delta = -1 or 1*/
	ushort a[3];
	a[0] = (pnt->id % 6) / 2;
	a[1] = 2*(a[0]+1) % 6;
	a[2] = 2*(a[0]+2) % 6;
	int bid=pnt->id/6;
	ushort *oldpid = s->bndeps + bid*6;/* (pnt->id & 0xfffe); */
	endpnt_t oldp=*pnt;
	float value2 = s->endpnts[a[0]][s->bndeps[pnt->id^1]].value; /* updated value */
	int inc=d>0;
	endpnt_t *first = s->endpnts[a[0]]-1;
	endpnt_t *last = s->endpnts[a[0]]+2*s->nbnds;
	int ds=((ismax==inc)<<1)-1;
	ushort oid=s->bnds[bid].oid;

	for(endpnt_t *p=pnt+d; p>=first && p<=last; p+=d){
		if(p==first){	/* TODO: make this less ugly */
			s->bndeps[oldp.id]=0;
			p[1]=oldp;
			break;
		}
		else if(p==last){
			s->bndeps[oldp.id]=2*s->nbnds-1;
			p[-1]=oldp;
			break;
		}

		if((oldp.value>p->value) == inc){
			if(!ismax == !ISMIN(p->flag)){ 		/* max vs. min / min vs. max */
				ushort *pid=s->bndeps + p->id/6*6;
				int overlap = oldpid[a[1]] < pid[a[1]+1] && 
					      oldpid[a[1]+1] > pid[a[1]] &&
					      oldpid[a[2]] < pid[a[2]+1] &&
					      oldpid[a[2]+1] > pid[a[2]];
				int lt=value2 < s->endpnts[a[0]][s->bndeps[p->id^1]].value;
				if(overlap && (lt == ismax)){ 	/* max,< / min,> */
					/* if(sp->validpr(oldp.flag, p->flag) */{
						ushort o=s->bnds[p->id/6].oid;
						if(ismax == inc) /* max,<,inc / min,>,dec */{
							addprb_objobj(oid,o);
						}
						else  /* max,<,dec / min,>,inc */{
							delprb_objobj(oid,o);
						}
					}
				}
			}

			s->bndeps[p->id]-=d;
			p[-d]=p[0];
		}
		else{	/* done */
			s->bndeps[oldp.id]=p-d-s->endpnts[a[0]];
			p[-d]=oldp;
			break;
		}
	}
}



void updbnd(sp_t *s, int bid, float bnd[2][3])
{
	ushort *pid = s->bndeps+bid*6;
	float old[6];
	int p[6];
	for(int i=0;i<6;i+=2){
		int a=i/2;
		int dec = bnd[0][a] < s->endpnts[0][pid[0]].value; /* bnd[0][a] is correct: dec means min bounds decrease */
		float *min = &(s->endpnts[a][pid[i]].value);
		float *max = &(s->endpnts[a][pid[i+1]].value);
		old[i]=*min;
		old[i+1]=*max;
		*min=bnd[0][a];
		*max=bnd[1][a];
		p[i]=(i+1)^dec;
		p[i+1]=i^dec;
	}

	for(int i=0;i<6;i++){
		updendpnt(s, s->endpnts[i/2]+pid[p[i]], old[p[i]]);
	}
}


static void delobj()
{
	int n=getnumu16(cmsp->op.obj.d);
	ITER(n){
		ushort d = cmsp->op.obj.d[i];
		queryobj_t*o = OBJPTR(d);
		sp_h_t *next;
		for (sp_h_t *s=o->sl; s; s=next) {
			next=s->next;
			pushu32(&cmsp->lut.bid_sp, *((uint*)&s->bid));

			// delete the splist node
			poolfree(&cmsp->sp_h_pool, s);
		}

		// delete all pairs involving this object
		unlinkmates(cmsp->pm, d);
	}
}


static void mendlink()
{
	int n = getnumu16(cmsp->lut.mend);
	ushort *mend = cmsp->lut.mend;

	// remove duplicated entries in mendlink[]

	rsort16(mend, n);
	int realn = 0;
	uint  m = INVALID_LINK;
	ITER(n){
		if (mend[i] != m) {
			mend[realn++] = m;
			m = mend[i];
		}
	}

	static uint t=0; // timestamp
	ITER(realn){
		t++;
		uint m = mend[i];
		queryobj_t *o = OBJPTR(m);
		for (sp_h_t *s=o->sl; s; s=s->next) {
			s->flag = t;
		}
		
		int common=0;
		for (uint p=cmsp->pm->first[m]; p!=INVALID_LINK; p=cmsp->pm->next[p]) { // p == partner
			queryobj_t *po = OBJPTR(cmsp->pm->link[p].obj);
			for (sp_h_t *s=po->sl; s; s=s->next) {
				if (s->flag == t) {
					common=1;
					goto commonlabel;
				}
			}

		commonlabel:
			if (!common) {
				delprb_objjbo(m, cmsp->pm->link[p].obj);
			}
		}
	}

	resetstacku16(cmsp->lut.mend);
}


#define LO16(bitfield) ((bitfield)&0xffff)
#define HI16(bitfield) ((bitfield)>>16)


static void batdel()
{
	// sort delspbid_buf[] by SP
	int n = getnumu32(cmsp->lut.bid_sp);
	uint *bid_sp = cmsp->lut.bid_sp;

	rsortlo16(bid_sp,n);

	int i = 0;
	ushort *b = alloca(32*1024);
	int done=0;
	do{
		ushort nb=0;
		ushort sp=LO16(bid_sp[i]);
		ushort nextsp;
		do{
			b[nb++] = HI16(bid_sp[i++]);
			if(i==n){
				done=1;
				break;
			}
			
			nextsp = LO16(bid_sp[i]);
		}while(nextsp == sp);

		rmbnds(SPPTR(sp), b, nb); // no delpr_buf_() calls
	}while(!done);

	resetstacku32(cmsp->lut.bid_sp);
}


static void moveobj()
{
	int n=getnumu16(cmsp->op.obj.m);
	ITER(n){
		queryobj_t *o = OBJPTR(cmsp->op.obj.m[i]);
		for(sp_h_t *s=o->sl; s; s=s->next){
			updbnd(s->sp, o->bid, o->bnd); // contains addprb_objobj()/delprb_objobj() calls
		}
	}
}


// no need for this -- it's done in boxpruning()
static void addobj()
{
	/*
	  for each object, find MSPs that it touches
	 */

	// 1. rasterize the object into the grid
	


	// 2. for each grid it touches{ for each MSP in that grid {test collision between the obj & MSP}}


}


static void delsp()
{
	if(cmsp->op.sp.type != SPOP_DEL) return;

	
	// tag the SP to delete
	sp_t *sp = SPPTR(cmsp->op.sp.id);
	sp->flag = 0xffffffff;

	ITER(sp->nbnds){
		// remove SP from the object's splist
		sp_h_t *h=sp->bnds[i].h;
		ushort oid=sp->bnds[i].obj;
		queryobj_t *o=OBJPTR(oid);
		
		if (h->prev) h->prev->next = h->next;
		else o->sl = o->next;
		
		if (h->next) h->next->prev = h->prev;
		
		// add object to mendlink[] if it's in some SP
		int insomesp=0;
		for (sp_h_t *s=o->sl; s; s=s->next) {
			if (s->sp->flag != 0xffffffff) {
				insomesp=1;
				break;
			}
		}
		
		if (insomesp) {
			pushu16(&cmsp->lut.mend, oid);
		}

		// remove the SP from the multiSP
		freesp(&cmsp->sp, sp);
	}
}

#define MAX_OBJ_SIZE 64

// for each moved/newly added SP, find nearby (old) SPs
static void findnearsps()
{
	// if no one moves, skip this step
	if(cmsp->op.sp.type != SPOP_MOVE || cmsp->op.sp.type != SPOP_ADD) return;

	// just test this SP against all other SPs -- brute force ftw!
	// cache the result?

	float bnd[6];
	ITER(3){
		bnd[i]=cmsp->op.sp.bnd[i]-MAX_OBJ_SIZE;
		bnd[i+3]=cmsp->op.sp.bnd[i+3]+MAX_OBJ_SIZE;
	}

	for(sp_t *sp=cmsp->splist; sp; sp=sp->next){
		if(bndintx(bnd, sp->bnd)){
			if(sp->id != cmsp->op.sp.id){
				pushu16(&cmsp->lut.nearsp, sp->id);
			}
		}
	}
}



static void updspaabb()
{
	int n=cmsp->op.sp.nm;
	ushort *m=cmsp->op.sp.m;
	ITER(n){
		sp_t *sp=SPPTR(m[i]);
		sp->flag = 0xffffffff;
		memcpy(sp->bnd, cmsp->op.sp.bnd, sizeof(sp->bnd));
	}
}


static void updprincipleaxis(msp_t *msp)
{
	float max_var=-1;
	int max_axis=-1;
	int n=msp->sppool.ns*2;
	float *values = alloca(n*sizeof(float));
	for(int a=0; a<3; a++){
		float *v=values;
		float avg=0;
		int i=0;
		for(sp_t *sp=msp->splist; sp; sp=sp->next){
			v[i]=sp->bnd[a];
			v[i+1]=sp->bnd[a+3];
			avg += v[i]+v[i+1];
			i+=2;
		}

		avg /= n;
		float var=0;
		for(i=0; i<n; i++){
			float k=v[i]-avg;
			var += k*k;
		}

		var /= n-1;

		if(var>max_var){
			max_var=var;
			max_axis=a;
		}
	}

	msp->principleaxis = max_axis;
}

/*
  To update adjancent SPs for moved objects, several approaches can be used:
  1. Box pruning
  2. Dynamic AABB tree

  To tell which one is better, just implement both of them, then compare!
 */
static void boxpruning()
{
	// choose the principle axis by finding the axis of the largest variance
	int type=cmsp->op.sp.type;
	if(type==SPOP_MOVE || type==SPOP_DEL){ // delay the updating to addsp() if type==SPOP_ADD
		updprincipleaxis(cmsp);
		int i=0;
		for(sp_t *sp=cmsp->splist; sp; sp=sp->next){
			cmsp->lut.spid[i]=sp->id;
			cmsp->lut.spmin[i++]=sp->bnd[cmsp->principleaxis];
		}

		// sort SP's min values
		rank32f(cmsp->lut.spmin, cmsp->lut.sprank, cmsp->sppool.ns);
	}

	// sort obj's min values -- do this each frame since they change very often
	int nm=cmsp->op.obj.nm;
	int n=nm + cmsp->op.obj.na;


	int n=cmsp->op.obj.nm + cmsp->op.obj.na;
	float *objmin=alloca(n*sizeof(float));
	uint *objrank=alloca(n*sizeof(uint));
	ushort *objid=alloca(n*sizeof(ushort));

	for(int i=0; i<n; i++){
		ushort m = cmsp->op.obj.m[i];
		objmin[i]=(OBJPTR(m))->bnd + cmsp->principleaxis;
		objid[i]=m;
	}
	for(int i=nm; i<nm+na; i++){
		ushort a = cmsp->op.obj.a[i-nm];
		objmin[i]=(OBJPTR(a))->bnd + cmsp->principleaxis;
		objid[i]=a;
	}

	ushort movedsp = cmsp->op.sp.id;

	rank32f(objmin, objrank, n);

	// box pruning to find obj-SP intersections
	uint *lastobj=objrank+n;
	uint *obj=objrank;
	uint *cursp=cmsp->lut.sprank;
	uint *lastsp=cmsp->lut.sprank+cmsp->sppool.ns;
	float *spmin=cmsp->lut.spmin;
	int a1=(cmsp->principleaxis+1)%3;
	int a2=(cmsp->principleaxis+2)%3;
	while(obj<lastobj && cursp<lastsp){
		uint sid=*cursp++;
		while(obj<lastobj && objmin[*obj]<spmin[sid]) obj++;
		uint *objinsp=obj;
		ushort spid=cmsp->lut.spid[sid];
		if(spid==movedsp) continue; // skip moved SP

		sp_t *s=SPPTR(spid);
		float spmax=s->bnd[cmsp->principleaxis+3];
		while(objinsp<lastobj && objmin[*objinsp]<=spmax){
			uint oid=objid[*objinsp++];
			queryobj_t *o=OBJPTR(oid);
			if(o->bnd[a1]<s->bnd[a1+3] && o->bnd[a2]<s->bnd[a2+3] &&
			   o->bnd[a1+3]>s->bnd[a1] && o->bnd[a2+3]>s->bnd[a2]){
				pushu32(&cmsp->lut.obj_sp, PACK(oid,spid));
			}
		}
	}

	uint *sp=cmsp->lut.sprank;
	uint *curobj=objrank;
	while(sp<lastsp && curobj<lastobj){
		uint oid=*curobj++;
		while(sp<lastsp && spmin[*sp]<objmin[oid]) sp++;
		uint *spinobj=sp;
		queryobj_t *o=OBJPTR(objid[oid]);
		float objmax=o->bnd[cmsp->principleaxis+3];
		while(spinobj<lastsp && spmin[*spinobj]<=objmax){
			uint sid=cmsp->lut.spid[*spinobj++];
			if(sid==movedsp) continue; // skip moved SP
			sp_t *s=SPPTR(sid);

			if(o->bnd[a1]<s->bnd[a1+3] && o->bnd[a2]<s->bnd[a2+3] &&
			   o->bnd[a1+3]>s->bnd[a1] && o->bnd[a2+3]>s->bnd[a2]){
				pushu32(&cmsp->lut.obj_sp, PACK(oid,sid));
			}
		}
	}

	
	SPPTR(cmsp->op.sp.id)->flag = 0;
}


static inline int bndinbnd(float in[6], float bnd[6])
{
	return in[0]>=bnd[0] && in[1]>=bnd[1] && in[2]>=bnd[2] &&
		in[3]<=bnd[3] && in[4]<=bnd[4] && in[5])<=bnd[5];
}


#define PACK(hi,lo) (((hi)<<16)|(lo))

// updates moved objects' splist
static void updsplist()
{
	int n=cmsp->op.obj.nm;
	int nm=0;

	ITER(n){
		ushort m=cmsp->op.obj.m[i];
		queryobj_t *o = OBJPTR(m);

		if(o->sl && !o->sl->next){ // the object is only in one SP; This optimizatin can only be used
			                                         // if SPs don't intersect with each other
			sp_t *s = o->sl->s;
			if(bndinbnd(o->bnd, s->bnd)) continue;
		}

		// only process those problematic moved objects
		cmsp->op.obj.m[nm++] = m;
	}

	boxpruning();	

	// for each moved object
	int done=0;
	uint *obj_sp=cmsp->lut.obj_sp;
	int n=getnumu32(cmsp->lut.obj_sp);
	rsorthi16(obj_sp,n);
	int i=0;
	int pr;
	ushort movedsp=cmsp->op.sp.id;
	do {
		ushort curobj=HI16(obj_sp[i]);
		int ntouches=0;
		pr=i++;
		do {
			sp_t *sp=SPPTR(curobj);
			sp->flag = 1;
			ntouches++;

			if (i==n) {
				done=1;
				break;
			}
			next=HI16(obj_sp[i]);
		} while (next==curobj);

		queryobj_t *obj=OBJPTR(curobj);
		sp_h_t *next;
		for(sp_h_t *sh=obj->sl; sh; sh=next){
			next=sh->next;
			if(sh->id==movedsp) continue; // skip moved SP

			if(sh->sp->flag != 1){ // only in old SP
				// add pair to bid_sp
				uint pr= PACK(sh->bid,sh->id);
				pushu32(&cmsp->lut.bid_sp, pr);

				// add obj to mendlink[]
				pushu16(&cmsp->lut.mend, curobj);

				// remove the SP from object's splist
				if(sh->prev){
					sh->prev->next = sh->next;
				}
				else{
					obj->sl = sh->next;
				}
				if(sh->next){
					sh->next->prev = sh->prev;
				}

				poolfree(&cmsp->sppool, sh);
			}
			else{
				sh->sp.flag = 0;
			}
		}

		for (int i=pr; i<pr+ntouches; i++) {
			ushort sid= LO16(obj_sp[i]);
			sp_t *s= SPPTR(sid);
			if (s->flag==1) {
				// add SP to object's SPlist
				sp_h_t *sh=poolalloc(&cmsp->sppool);
				sh->prev=NULL;
				sh->next=obj->sl;
				if (obj->sl) {
					obj->sl->prev = sh;
				}
				obj->sl=sh;

				// add to obj_sp[]
				pushu32(&cmsp->lut.obj_sp, PACK(curobj,sid));

				s->flag=0; // So an SP will be added to the object's list only once -- s->flag==1 only at the first time
					   // it's encountered.
			}
		}
	} while (!done);
}


/*
  Add dynamic/static object to (remove dynamic/static object from) the moved SP

  This step may:
  1. Add dynamic obj to moved SP
  2. Remove dynamic obj from moved SP

  But updsplist() may also do 1&2.

  So only one of them should do the above things.

 */
static void movesp()
{
	if (cmsp->op.sp.type != SPOP_MOVE) return;

	// find _static_ objects that the new (moved) SP contains
	int size=0;
	int ns=getnumu16(cmsp->lut.nearsp);
	ITER(ns){
		sp_t *sp=SPPTR(cmsp->lut.nearsp[i]);
		size += sp->nbnds;
	}
	ushort *objinnew=alloca(size*sizeof(ushort));

	int n=0;
	float *sb=cmsp->op.sp.bnd;
	ITER(ns){
		sp_t *s=SPPTR(cmsp->lut.nearsp[i]);
		ITER1(s->nbnds){
			ushort oid = s->bnds[j].oid;
			float *ob=OBJPTR(oid)->bnd;
			if (bndintx(sb, ob)) {
				objinnew[n++]=oid;
			}
		}
	}

	enum {NOT_ENCOUNTERED, IN_OLD, IN_NEW, IN_OLD_AND_NEW};
	
	// for each object in the old SP
	sp_t *s=SPPTR(cmsp->op.sp.id);
	ITER(s->nbnds){	// each object will be encountered once
		queryobj_t *obj=OBJPTR(sp->bnds[i].oid);
		obj->flag=IN_OLD;
	}

	// for each object in the new SP
	ITER(n){// object may be encountered multiple times
		queryobj_t *obj=OBJPTR(objinnew[i]);
		if(obj->flag==NOT_ENCOUNTERED){ // only in new SP
			obj->flag=IN_NEW;

			// add to obj_sp[]
			uint pr= PACK(objinnew[i], cmsp->op.sp.id);
			pushu32(&cmsp->lut.obj_sp, pr);
		}
		else if(obj->flag==IN_OLD){ // already in old SP => in both old and new SP
			obj->flag=IN_OLD_AND_NEW;
		}
	}


	ITER(s->nbnds){ // each object will be encountered once
		queryobj_t *obj=OBJPTR(sp->bnds[i].oid);
		if(obj->flag==IN_OLD){ // only in old SP
			uint pr= PACK(i, cmsp->op.sp.id);
			pushu32(&cmsp->lut.bid_sp. pr); // add to bid_sp[]
			pushu16(&cmsp->lut.mend, sp->bnds[i].oid); // add to mend[]
		}
		obj->flag=NOT_ENCOUNTERED; // clear flags
	}

	// clear flags
	ITER(n){
		OBJPTR(objinnew[i])->flag=NOT_ENCOUNTERED;
	}

	resetstacku16(cmsp->lut.nearsp);

}

static void addsp()
{
	if(cmsp->op.sp.type!=SPOP_ADD) return;

	// find objects that touch this SP

	// reserve enough space for objinnew[] -- count!
	int size=0;
	int ns=getnumu16(cmsp->lut.nearsp);
	ITER(ns){
		sp_t *s=SPPTR(cmsp->lut.nearsp[i]);
		size += s->nbnds;
	}
	usort *objinnew=alloca(size*sizeof(ushort));


	// perform AABB test to find objects the SP actually touches
	int n=0;
	float *bnd = cmsp->op.sp.bnd;
	ITER(ns){
		sp_t *s=SPPTR(cmsp->lut.nearsp[i]);
		ITER1(s->nbnds){
			ushort oid=s->bnds[j].oid;
			float *objbnd=OBJPTR(oid)->bnd;
			if (bndintx(bnd, objbnd)) {
				objinnew[n++]=oid;
			}
		}
	}


	// new SP
	ushort spid;
	sp_t *sp=allocsp(&cmsp->sppool, &spid);
	sp->prev=NULL;
	sp->next=cmsp->splist;
	if(sp->next) sp->next->prev=sp;
	cmsp->splist=sp;
	
	ITER(n){
		pushu32(&cmsp->lut.obj_sp, PACK(objinnew[i], spid));
	}

	updprincipleaxis(cmsp);

	resetstacku16(cmsp->lut.nearsp);
}


static void batadd()
{
	int n=cmsp->lut.obj_sp.n;
	uint *obj_sp=cmsp->lut.obj_sp;
	rsortlo16(obj_sp, n);

	int i=0;
	ushort *obj=alloca(MAX_NUM_OBJ);
	int done=0;
	do {
		ushort no=0;
		ushort cursp=obj_sp[i] & 0xffff;
		ushort nextsp;
		do {
			obj[no++] = obj_sp[i++]>>16;
			if (i==n) {
				done=1;
				break;
			}
			
			nextsp = obj_sp[i] & 0xffff;
		} while(nextsp == cursp);

		addbnds(SPPTR(cursp), obj, no); // TODO -- rewrite addbnds() ?
	} while (!done);

	resetstacku32(cmsp->lut.obj_sp);
}


#define CLRSB(b) ((b)&0x7fffffff)
#define GETSB(b) ((b)&0x80000000)

static void flushobjjbo()
{
	// sort by pairs
	int n=cmsp->lut.obj_jbo.n;
	uint *obj_jbo=cmsp->lut.obj_jbo;
	rsort32_ignore_signbit(obj_jbo, n); // TODO: implement this
	
	// count addpr() ops
	// count delpr() ops

	int i=0;
	int done=0;
	do {
		int na,nd;
		na=nd=0;
		uint cur=CLRSB(obj_jbo[i]);
		uint next;
		do {
			int sign=GETSB(obj_jbo[i++]);
			if(sign) nd++;
			else na++;

			if(i==n){
				done=1;
				break;
			}

			next=CLRSB(obj_jbo[i]);
		} while (next==cur);

		ushort p[2];
		*((uint*)p) = cur;

		if (isprvalid(p[0],p[1])){
			if (na && !nd) {
				addpr(cmsp->pm, p[1], p[0]);
			}
			else if (!na && nd) {
				delpr(cmsp->pm, p[1], p[0]);
			}
		}
	}while(!done);

	resetstacku32(cmsp->lut.obj_jbo);
}

static void flushobjobj()
{
	// sort by pairs
	int n = getnumu32(cmsp->lut.obj_jbo);
	uint *obj_jbo = cmsp->lut.obj_jbo;
	rsort32_ignore_signbit(obj_jbo, n); // TODO: implement this
	
	// count addpr() ops
	// count delpr() ops

	int i=0;
	int done=0;
	do {
		int na,nd;
		na=nd=0;
		uint cur=CLRSB(obj_obj[i]);
		uint next;
		do {
			int sign=GETSB(obj_obj[i++]);
			if(sign) nd++;
			else na++;

			if (i==n) {
				done=1;
				break;
			}

			next=CLRSB(obj_obj[i]);
		} while (next==cur);

		ushort p[2];
		*((uint*)p) = cur;

		if (isprvalid(p[0],p[1])){
			if (na > nd){
				addpr(cmsp->pm, p[0], p[1]);
			}
			else if (na < nd){
				delpr(cmsp->pm, p[0], p[1]);
			}
		}
	} while (!done);

	resetstacku32(cmsp->lut.obj_obj);
}

/*
  Each object (cube/wall/lightvolume) stores a pointer (or handle) to the query object

 */

void mspaddobj(float bnd[6], int type, ushort *h)
{
	queryobj_t *obj=allocobj(&cmsp->objpool, h);
	memcpy(obj->bnd, bnd, sizeof(bnd));
	obj->type = type;
	pushu16(&cmsp->op.obj.a, *h);
}

void mspdelobj(ushort h)
{
	pushu16(&cmsp->op.obj.d, h);
}

void mspmoveobj(float bnd[6], ushort h)
{
	queryobj_t *o = OBJPTR(h);
	memcpy(o->bnd, bnd, sizeof(bnd));
	pushu16(&cmsp->op.obj.m, h);
}

// this OP should be delayed -- later in findnearsps() we should test whether it intersect with other SPs
void mspaddsp(float bnd[6], ushort *h)
{
	cmsp->op.sp.type = SPOP_ADD;
	memcpy(cmsp->op.sp.bnd, bnd, sizeof(bnd));
}

void mspdelsp(ushort h)
{
	cmsp->op.sp.type = SPOP_DEL;
	cmsp->op.sp.id = h;
}

void mspmovesp(float bnd[6], ushort h)
{
	cmsp->op.sp.type = SPOP_MOVE;
	memcpy(cmsp->op.sp.bnd, bnd, sizeof(bnd));
	cmsp->op.sp.id=h;
}

void mspflush()
{
	delobj();
	mendlink();		/*  */
	batdel();		/*  */
	updaabb();		/*  */
	moveobj();		/*  */
//	addobj();

	delsp();
	findnearsps();
	updspaabb();
	updsplist();		/*  */
	movesp();
	addsp();

	batadd();		/*  */
	flushobjjbo();		/*  */
	flushobjobj();		/*  */
}


queryobj_t *mspquery(ushort h, int *num)
{
	queryobj_t *obj = OBJPTR(h);
	if(!obj) return 0;
	
	int n=0;
	for (queryobj_h_t *o=obj->ol; o; o=o->next) {
		cmsp->lut.q[n]=o->obj;
		n++;
	}

	*num=n;
	return cmsp->lut.q;
}


void loadmsp(char *name)
{
}

void storemsp()
{
}

void newmsp(char *name)
{
}


/*
  World contains MSP (building) & other stuff;
  MSP contains SP (room);
  SP contains AABB (object);


  The engine should be as simple as possible.





 */







/**
  multiple MSPs is a must

  The current MSP should be automatically chosen according to the view point.
  
  Each MSP should have an AABB tree;
  All MSPs should be stored in a SP;
  
  When the player is in editmode (a.k.a. shipyard), no AABB tree is used 
  => don't update until the ship leaves the shipyard.
  => the shape of the ship won't change
  => compute the AABB tree once


  Multiplayer game:
  Multiplayer physics is hard

  Interactions between players


      ______
    \|/     |
  players <-+
     /|\
      |
     \|/
  objects <-+
    /|\     |
     |______|


  To each player, other players are equal to talking NPCs.
  Other players can move and talk freely.
  NPCs can move but can not talk.

  Movements of other players may generate complex events which are very difficult to synchronize--each local event may have global effects
  => use NPC to simulate other players' movements.

  Just let NPCs to emulate other players' movements:
  A player moves --> network lag --> the NPC moves

  In the remote player P's own simulation, world state is A. => the event (P moves in A) is simulated by the remote player.
  In the local player's simulation, world state is B. => the event (NPC_P moves in B) is simulated by the local player.

  And the outcome of the two simulations may be different.
  
  Broardcasting the events.


  Direct interactions between players:
  Player A kicked player B.
  In player B's simulation, NPC_A approached to kicked him. The path that NPC_A took may differ from the one player A took.


  1. A kick/punch/stab B =>
  In A's simulation, these actions are triggered by key press/mouse clicks;
  In B's simulation, NPC_A did all these actions.
  The two sets of actions may differ from each other -- player's orientation / hit point may be different.

  2. A throw an object at B


  When in play mode, the ship is less destructible than in edit mode:
  You cannot dig holes on the floor or roof -- you can only dig holes on the wall => avoiding lots of objects falling downstairs.
  
 */


#else


#define BUFSIZE 32
sp_t *allocsp(int nbnds)
{
    nbnds = nextp2(nbnds);
    sp_t *s = malloc(sizeof(sp_t));
    int size = 2*nbnds*sizeof(endpnt_t);
    s->endpnts[0] = malloc(size);
    s->endpnts[1] = malloc(size);
    s->endpnts[2] = malloc(size);
    s->bndeps = malloc(6*sizeof(unsigned short)*nbnds);
    s->bnds = calloc(nbnds,sizeof(bnd_t));
    s->nbnds = s->navailbnds = 0;
    s->cap = nbnds;
    s->dbnds = malloc(nbnds*sizeof(unsigned short));
    s->ndbnds = 0;
    s->dcap = nbnds;
    s->mm = allocmm();

    size = 2*BUFSIZE*sizeof(endpnt_t);
    s->addbuf.endpnts[0] = malloc(size);
    s->addbuf.endpnts[1] = malloc(size);
    s->addbuf.endpnts[2] = malloc(size);
    s->addbuf.nbnds = 0;
    s->addbuf.cap = BUFSIZE;
    s->addbuf.bids = malloc(BUFSIZE*sizeof(unsigned short*));
    s->addbuf.owners = malloc(BUFSIZE*sizeof(void*));
    s->rmbuf.bids = malloc(BUFSIZE*sizeof(unsigned short));
    s->rmbuf.nbnds = 0;
    s->rmbuf.cap = BUFSIZE;
    return s;
}

void freesp(sp_t *s)
{
    for(int a=0;a<3;a++){
	if(s->endpnts[a]) free(s->endpnts[a]);
	if(s->addbuf.endpnts[a]) free(s->addbuf.endpnts[a]);
    }
    if(s->bndeps) free(s->bndeps);
    if(s->bnds) free(s->bnds);
    if(s->addbuf.bids) free(s->addbuf.bids);
    if(s->addbuf.owners) free(s->addbuf.owners);
    if(s->rmbuf.bids) free(s->rmbuf.bids);
    freemm(s->mm);
    free(s);
}

inline void *getowner(const sp_t *s, unsigned short bid)
{
    return s->bnds[bid].owner;
}


static int cmpep(endpnt_t *p1, endpnt_t *p2)
{
    if(p1->value < p2->value) return -1;
    else if(p1->value > p2->value) return 1;
    else{
	int ismin[2];
	ismin[0] = p1->flag & ISMIN;
	ismin[1] = p2->flag & ISMIN;
	if(ismin[0]==ismin[1]) return 0;
	else if(ismin[0]) return -1;
	else return 1;
    }
}


static inline void addtouchers(sp_t *s, unsigned short id, unsigned short id2, unsigned short flag, unsigned short flag2)
{
    assert(id!=id2);
    bnd_t *b=s->bnds+id;
    bnd_t *b2=s->bnds+id2;

    assert(b->ns+b->nd<=32);
    assert(b2->ns+b2->nd<=32);

    if(flag & ISDYNAMIC){
	assert(b->touchers);

	if(flag2 & ISDYNAMIC){
	    b->touchers[b->nd++]=id2;
	    b2->touchers[b2->nd++]=id;
	}
	else{
	    b->touchers[31-b->ns++]=id2;
	}
    }
    else if(flag2 & ISDYNAMIC){
	assert(b2->touchers);

	b2->touchers[31-b2->ns++]=id;
    }
}


static inline void rmtouchers(sp_t *s, unsigned short id, unsigned short id2, unsigned short flag, unsigned short flag2)
{
    assert(id!=id2);

    bnd_t *b=s->bnds+id;
    bnd_t *b2=s->bnds+id2;
    if(flag & ISDYNAMIC){
	if(flag2 & ISDYNAMIC){
	    for(int i=0;i<b->nd;i++){
		if(b->touchers[i]==id2){
		    b->touchers[i]=b->touchers[--b->nd];
		    break;
		}
	    }
	    for(int i=0;i<b2->nd;i++){
		if(b2->touchers[i]==id){
		    b2->touchers[i]=b2->touchers[--b2->nd];
		    break;
		}
	    }
	}
	else{
	    for(int i=31;i>31-b->ns;i--){
		if(b->touchers[i]==id2){
		    b->touchers[i]=b->touchers[32-b->ns--];
		    break;
		}
	    }
	}

    }
    else if(flag2 & ISDYNAMIC){
	for(int i=31;i>31-b2->ns;i--){
	    if(b2->touchers[i]==id){
		b2->touchers[i]=b2->touchers[32-b2->ns--];
		break;
	    }
	}
    }
}


static void rm1bnd(sp_t *s, unsigned short bid)
{
    bnd_t *b=s->bnds+bid;
    int id=bid*6;;
    int d=s->endpnts[0][s->bndeps[id]].flag & ISDYNAMIC;


    {
	endpnt_t *p0,*p1;
	p0=s->endpnts[0]+s->bndeps[id];
	p1=s->endpnts[0]+s->bndeps[id+1];
	for(endpnt_t *p=p0+1;p<p1;p++){
	    p->stab--;
	}
    }


    for(int a=0;a<3;a++){
	endpnt_t *p0,*p1;
	p0=s->endpnts[a]+s->bndeps[id++];
	p1=s->endpnts[a]+s->bndeps[id++];
	for(endpnt_t *p=p0+1; p<p1; p++){
	    p[-1]=p[0];
	    s->bndeps[p->id]--;
	}
	for(endpnt_t *p=p1+1; p<s->endpnts[a]+2*s->nbnds; p++){
	    p[-2]=p[0];
	    s->bndeps[p->id]-=2;
	}
    }
    s->nbnds--;
    s->navailbnds++;
    id=bid*6;
    s->endpnts[0][s->nbnds*2].id=id;

    /* remove pairs */
    if(d){
	for(int i=0;i<b->nd;i++){
	    bnd_t *t=s->bnds+b->touchers[i];
	    if(s->endpnts[0][s->bndeps[b->touchers[i]*6]].flag & ISDYNAMIC){
		for(int j=0;j<t->nd;j++){
		    if(t->touchers[j]==bid){
			t->touchers[j]=t->touchers[--t->nd];
			break;
		    }
		}
	    }
	}
	for(int i=0;i<s->ndbnds;i++){ /* TODO: store index to s->dbnds[]? */
	    if(s->dbnds[i]==bid){
		s->dbnds[i]=s->dbnds[--s->ndbnds];
		break;
	    }
	}
#if 1
	if(b->touchers){
	    freeblk(s->mm,b->touchers);
	    b->touchers = NULL;
	}
#endif
    }
    else{	/* removing static objects is more expensive than removing dynamic ones */
	for(int i=0;i<s->ndbnds;i++){
	    bnd_t *b=s->bnds+s->dbnds[i];
	    for(int j=31;j>31-b->ns;j--){
		if(b->touchers[j]==bid){
		    b->touchers[j]=b->touchers[32-b->ns--];
		    break;
		}
	    }
	}
    }

/* 	b->ns = 0xffff; */

}



static void addbnds(sp_t *s, endpnt_t *endpnts[3], unsigned short *bndids[], void **owners, int nbnds)
{
    if(nbnds > s->navailbnds){
	int n = nbnds + s->nbnds;
	if(n > s->cap){ /* mem is too small, need enlargement */
	    s->cap = nextp2(n);
	    int size = 2*s->cap*sizeof(endpnt_t);
	    s->endpnts[0] = realloc(s->endpnts[0], size);
	    s->endpnts[1] = realloc(s->endpnts[1], size);
	    s->endpnts[2] = realloc(s->endpnts[2], size);
	    s->bndeps = realloc(s->bndeps, 6*s->cap*sizeof(unsigned short));
	    s->bnds = realloc(s->bnds, s->cap*sizeof(bnd_t));
	}

	/* || s->nbnds | s->navailbnds | s->cap - s->nbnds - s->navailbnds ||<--s->cap */

	int ntot = s->nbnds + s->navailbnds; /* total num of bnds allocated */
	int id=ntot*6;
	for(int i=2*(s->nbnds+s->navailbnds); i<2*(s->nbnds+nbnds); i+=2){ /* set values for new endpoints */
	    s->endpnts[0][i].id = id; 
	    id+=6;
	}
	s->navailbnds = 0;
    }
    else{
	s->navailbnds-=nbnds;
    }

    /* copy values of the available endpoints to endpoints of the newly added bounds */
    endpnt_t *p=s->endpnts[0]+2*s->nbnds;
    for(int i=0; i<nbnds; i++){
	unsigned short id=p->id;
	*bndids[i] = id/6;
	for(int j=0;j<6;j++){
	    endpnts[j/2][2*i+(j&1)].id = id++;
	}

	bnd_t *b=s->bnds+p->id/6;
	b->owner = owners[i];
	if(endpnts[0][2*i].flag & ISDYNAMIC){
	    if(s->ndbnds==s->dcap){
		s->dcap*=2;
		s->dbnds = realloc(s->dbnds, s->dcap*sizeof(unsigned short));
	    }
	    s->dbnds[s->ndbnds++] = p->id/6;
	    b->touchers = allocblk(s->mm);
	    b->ns = b->nd = 0;
	}
	else{
	    b->touchers = NULL;
	}
	p+=2;
    }

    int np=2*nbnds;
    /* sort the new endpoints, and update the links pointing to endpoints after sorting (to be used in the following step) */
    for(int a=0;a<3;a++){
	if(np>2) qsort(endpnts[a], np, sizeof(endpnt_t), cmpep);
	for(int i=0;i<np;i++){
	    s->bndeps[endpnts[a][i].id] = i;
	}
    }

    /* add pairs between the newly added bnds (self pruning through x axis) */
    if(np>2){
	for(int i=0; i<np; i++){
	    if(endpnts[0][i].flag & ISMIN){
		endpnt_t *min = endpnts[0] + i; /* endpnts[] are already sorted*/
		endpnt_t *max = endpnts[0] + s->bndeps[min->id+1];
		int d = min->flag & ISDYNAMIC;
		bnd_t *b=s->bnds + min->id/6;
		for(endpnt_t *p=min+1; p<max; p++){ /* min<p<max */
		    if((d || (p->flag & ISDYNAMIC)) && (p->flag & ISMIN)){
			/* endpnts[i] are sorted => comparing indices is enough */
			unsigned short *pid = s->bndeps+p->id; /* PLACE 1 */
			unsigned short *minid = s->bndeps+min->id;
			if(pid[2]<minid[3] && minid[2]<pid[3] && pid[4]<minid[5] && minid[4]<pid[5]){
			    addtouchers(s,p->id/6,min->id/6,p->flag,min->flag);
			}
		    }
		}
	    }
	}
    }


    /* move the inserted bnds to the right positions */
    /* y,z axis */
    endpnt_t *newp,*oldp,*tail;
    for(int a=1;a<3;a++){
	newp = endpnts[a]+np-1;
	oldp = s->endpnts[a]+2*s->nbnds-1;
	tail = oldp+2*nbnds;
	while(newp>=endpnts[a]){
	    if(oldp>=s->endpnts[a] && oldp->value>newp->value){
		*tail = *oldp--;
	    }
	    else{
		*tail = *newp--;
	    }
	    s->bndeps[tail->id] = tail-s->endpnts[a];
	    --tail;
	}
    }

    /* x axis, plus box pruning */
    newp = endpnts[0]+np-1;
    oldp = s->endpnts[0]+2*s->nbnds-1;
    tail = oldp+2*nbnds;
    int stab=0; 
    while(newp>=endpnts[0]){
	if(oldp>=s->endpnts[0] && newp->value<oldp->value){ /* newp<oldp */
	    if(oldp->flag & ISMIN){
		if(newp->flag & ISMIN){
		    endpnt_t *oldp2=s->endpnts[0]+s->bndeps[oldp->id+1];
		    float n2 = s->endpnts[0][s->bndeps[newp->id+1]].value;
		    if(oldp2->value<n2){ /* oldp2<newp2 => newp<oldp<oldp2<newp2 */
			unsigned short *oldid, *newid;
			oldid=s->bndeps+oldp->id;
			newid=s->bndeps+newp->id;
			if(oldid[2]<newid[3] && newid[2]<oldid[3] && 
			   oldid[4]<newid[5] && newid[4]<oldid[5]){
			    addtouchers(s, oldp->id/6, newp->id/6, oldp->flag, newp->flag);
			}
		    }
		}
		stab--; 
	    }
	    else{
		stab++; 
	    }
	    *tail=*oldp--;
	}
	else{ 		/* newp>oldp */
	    if(stab){ /* the endpoint is contained in some intervals */
		int st=stab;
		endpnt_t *p=oldp;
		while(st){
		    if(p->flag & ISMIN){
			endpnt_t *p2=s->endpnts[0]+s->bndeps[p->id+1];
			if(p->value <= newp->value && p2->value >= newp->value){ /* p<newp<p2 */
			    if((newp->flag & ISMIN)|| /* 1,2 */
			       (endpnts[0][s->bndeps[newp->id-1]].value<p->value)){ /* 3 */
				unsigned short *pid, *newid;
				pid = s->bndeps + p->id/6*6;
				newid = s->bndeps + newp->id/6*6;
				if(pid[2]<newid[3] && newid[2]<pid[3] && 
				   pid[4]<newid[5] && newid[4]<pid[5]){
				    addtouchers(s, p->id/6, newp->id/6, p->flag, newp->flag);
				}
			    }
			    st--;
			}
		    }
		    p--;
		}
	    }
	    *tail=*newp--;
	}
	s->bndeps[tail->id] = tail-s->endpnts[0];
	--tail;
    }
	
    s->nbnds+=nbnds;


    /*  */
    {
	int st=0;
	endpnt_t *p;
	for(p=endpnts[0]+2*nbnds-1; p>endpnts[0]; p--){
	    endpnt_t *real = s->endpnts[0] + s->bndeps[p->id];
	    endpnt_t *next = s->endpnts[0] + s->bndeps[p[-1].id];
	    if(real->flag & ISMIN){
		real->stab = real[1].stab-1;
		--st;
	    }
	    else{
		if(real == s->endpnts[0]+2*s->nbnds-1){
		    real->stab=1;
		}
		else{
		    real->stab = real[1].stab+1;
		}
		++st;
	    }
	    if(!st) continue;
	    for(endpnt_t *stab = real-1; stab>next; stab--){
		stab->stab += st;
	    }
	}
	endpnt_t *first = s->endpnts[0] + s->bndeps[p->id];
	first->stab = first[1].stab-1;
    }

}


static int cmpushort(unsigned short *a, unsigned short *b)
{
    return ((*a>*b)<<1)-1;
}

static void rmbnds(sp_t *s, unsigned short bids[], unsigned short nbnds)
{
    if(nbnds==1){
	rm1bnd(s,bids[0]);
	return;
    }


    /* update stabbing numbers */
    /* three versions:
     */
#define BNDEPSBUF_SIZE 2048


    int deferred=0;
    if(nbnds>s->nbnds/2 || nbnds>BNDEPSBUF_SIZE/2){	/* brute force 1: TODO: find a better threshold other than 1/2 */
	deferred=1;
    }
    else if(nbnds<30){	/* brute force 2: TODO: replace the magic number "30"*/
	for(int i=0;i<nbnds;i++){
	    endpnt_t *p0=s->endpnts[0]+s->bndeps[bids[i]*6];
	    endpnt_t *p1=s->endpnts[0]+s->bndeps[bids[i]*6+1];
	    for(endpnt_t *p=p0+1;p<p1;p++){
		p->stab--;
	    }
	}
    }
    else{			/* optimaaaaaarz */
	static unsigned short buf[BNDEPSBUF_SIZE];
	unsigned short *b=buf;
	for(int i=0;i<nbnds;i++){
	    unsigned short *p=s->bndeps+bids[i]*6;
	    *b++ = *p++;
	    *b++ = *p;
	}
	qsort(buf,2*nbnds,sizeof(unsigned short),cmpushort); /* TODO: write some faster shit other than qsort? e.g. radix*/
	int stab=0;
	for(int i=2*nbnds-1;i>=0;i--){
	    unsigned char ismin = s->endpnts[0][buf[i]].flag & ISMIN;
	    if(ismin) --stab; else ++stab;
	    if(!stab) continue;
	    for(int j=buf[i]-1; j>buf[i-1]; j--){
		s->endpnts[0][j].stab -= stab;
	    }
	}
    } 

    /* 1. set flags of the endpoints which are going to be deleted to "invalid"(0xffff) */
    /* 2. shift endpoints */
    unsigned short minp[3]={0xffff,0xffff,0xffff};
    for(int a=0;a<3;a++){
	for(int i=0;i<nbnds;i++){
	    unsigned short *p = s->bndeps+6*bids[i];
	    s->endpnts[a][p[2*a]].flag = INVALID_ID;
	    s->endpnts[a][p[2*a+1]].flag = INVALID_ID;
	    /* find min endpoints */
	    if(p[a]<minp[a]) minp[a]=p[a<<1];
	    /* for a gonna-be-deleted bound, .ns is useless, so use it to tag this bound as "garbage" */
	    s->bnds[bids[i]].ns=0xffff;
	}
	unsigned short offset=1;
	for(unsigned short i=minp[a]+1; i<2*s->nbnds; i++){
	    if(s->endpnts[a][i].flag == INVALID_ID){
		offset++;
	    }
	    else{
		s->endpnts[a][i-offset] = s->endpnts[a][i];
		s->bndeps[s->endpnts[a][i].id]=i-offset;
	    }
	}
    }



    /* 3. link endpoints to unused bounds for later use. 
       now endpoints[i].id is used to store bnd id, NOT index to the bndeps[] array.
       only endpoints[0][s->bndeps[6*bids[i]]].id is set. the other five endpoints' .id field is not set.
    */
    s->nbnds-=nbnds;
    s->navailbnds+=nbnds;
    endpnt_t *p=s->endpnts[0]+s->nbnds*2;
    for(unsigned short i=0;i<nbnds; i++,p+=2){
	p->id=bids[i];
	/* delete colliding pairs */
	if(s->bnds[bids[i]].touchers){ /* for each deleted dynamic bound */
	    bnd_t *b=s->bnds+bids[i];
	    for(unsigned short j=0;j<b->nd;j++){ /* for each toucher */
		bnd_t *t=s->bnds+b->touchers[j];
		if(t->ns != 0xffff) t->nd |= 0x8000; /* set the highest bit to 1 => t->touchers changes */
	    }
	}
    }
    /* for each dynamic bound, check whether its touchers are valid. */
    for(unsigned short i=0;i<s->ndbnds;i++){
	bnd_t *b=s->bnds+s->dbnds[i];
	if(b->ns == 0xffff){ /* the bound is going to be removed */
	    freeblk(s->mm,b->touchers);
	    s->dbnds[i--]=s->dbnds[s->ndbnds--];
	    continue;
	}
	{
	    unsigned short j=31;
	    while(j>31-b->ns){
		bnd_t *t=s->bnds + b->touchers[j];
		if(t->ns == 0xffff){
		    b->touchers[j]=b->touchers[32-b->ns--];
		}
		else j--;
	    }
	}
	if(b->nd & 0x8000){
	    b->nd ^= 0x8000;
	    short j=0;
	    while(j<b->nd){
		bnd_t *t=s->bnds + b->touchers[j];
		if(t->ns == 0xffff){
		    b->touchers[j]=b->touchers[--b->nd];
		}
		else j++;
	    }
	}
    }

    if(deferred){		/* update the stabbing numbers from scratch if necessary */
	int stab=0;
	for(endpnt_t *p=s->endpnts[0]+s->nbnds*2-1; p>=s->endpnts[0]; p--){
	    if(p->flag & ISMIN){
		p->stab = --stab;
	    }
	    else{
		p->stab = ++stab;
	    }
	}
    }

}


static void updendpnt(sp_t *s, endpnt_t *pnt, float oldvalue)
{
    int ismax = !(pnt->flag & ISMIN);
    int d = ((pnt->value > oldvalue)<<1)-1; /* delta = -1 or 1*/
    unsigned short a[3];
    a[0] = (pnt->id % 6) / 2;
    a[1] = 2*(a[0]+1) % 6;
    a[2] = 2*(a[0]+2) % 6;
    int bid=pnt->id/6;
    unsigned short *oldpid = s->bndeps + bid*6;/* (pnt->id & 0xfffe); */
    endpnt_t oldp=*pnt;
    float value2 = s->endpnts[a[0]][s->bndeps[pnt->id^1]].value;
    int inc=d>0;
    endpnt_t *first = s->endpnts[a[0]]-1;
    endpnt_t *last = s->endpnts[a[0]]+2*s->nbnds;
    int ds=((ismax==inc)<<1)-1;
    for(endpnt_t *p=pnt+d; p>=first && p<=last; p+=d){
	if(p==first){	/* TODO: make this less ugly */
	    s->bndeps[oldp.id]=0;
	    p[1]=oldp;
	    break;
	}
	else if(p==last){
	    s->bndeps[oldp.id]=2*s->nbnds-1;
	    p[-1]=oldp;
	    break;
	}

	if((oldp.value>p->value) == inc){
	    if(!ismax == !(p->flag & ISMIN)){ 		/* max vs. min / min vs. max */
		unsigned short *pid=s->bndeps + p->id/6*6;
		int overlap = oldpid[a[1]] < pid[a[1]+1] && 
		    oldpid[a[1]+1] > pid[a[1]] &&
		    oldpid[a[2]] < pid[a[2]+1] &&
		    oldpid[a[2]+1] > pid[a[2]];
		int lt=value2 < s->endpnts[a[0]][s->bndeps[p->id^1]].value;
		if(overlap && (lt == ismax)){ 	/* max,< / min,> */
		    int d = oldp.flag & ISDYNAMIC;
		    int d2 = p->flag & ISDYNAMIC;
		    if(ismax == inc) /* max,<,inc / min,>,dec */{
			addtouchers(s,bid,p->id/6,d,d2);
		    }
		    else  /* max,<,dec / min,>,inc */{
			rmtouchers(s,bid,p->id/6,d,d2);		
		    }
		}
	    }

	    /* TEST--> */
	    p->stab += ds; /* TODO: skip this for y,z axis? */
	    /* <--TEST */


	    s->bndeps[p->id]-=d;
	    p[-d]=p[0];
	}
	else{	/* done */
	    s->bndeps[oldp.id]=p-d-s->endpnts[a[0]];
	    p[-d]=oldp;
	    break;
	}
    }
}



void updbnd(sp_t *s, int bid, float bnd[2][3])
{
    unsigned short *pid = s->bndeps+bid*6;
    float old[6];
    int p[6];
    for(int i=0;i<6;i+=2){
	int a=i/2;
	int dec = bnd[0][a] < s->endpnts[0][pid[0]].value; /* bnd[0][a] is correct: dec means min bounds decrease */
	float *min = &(s->endpnts[a][pid[i]].value);
	float *max = &(s->endpnts[a][pid[i+1]].value);
	old[i]=*min;
	old[i+1]=*max;
	*min=bnd[0][a];
	*max=bnd[1][a];
	p[i]=(i+1)^dec;
	p[i+1]=i^dec;
    }

    for(int i=0;i<6;i++){
	updendpnt(s, s->endpnts[i/2]+pid[p[i]], old[p[i]]);
    }

    /* update stabbing numbers for the updated bound */
    endpnt_t *p0,*p1;
    p0=s->endpnts[0]+s->bndeps[bid*6];
    p1=s->endpnts[0]+s->bndeps[bid*6+1];
    if(p1==s->endpnts[0]+s->nbnds*2-1){
	p1->stab=1;
    }
    else{
	p1->stab=p1[1].stab+1;
    }
    p0->stab=p0[1].stab-1;
}


/* DON'T EVER call updbnd() during two updbnd_def() */
void updbnd_def(sp_t *s, int bid, float bnd[2][3])
{
    unsigned short *pid = s->bndeps+bid*6;
    s->endpnts[0][pid[0]].value = bnd[0][0];
    s->endpnts[0][pid[1]].value = bnd[1][0];
    s->endpnts[1][pid[2]].value = bnd[0][1];
    s->endpnts[1][pid[3]].value = bnd[1][1];
    s->endpnts[2][pid[4]].value = bnd[0][2];
    s->endpnts[2][pid[5]].value = bnd[1][2];
}

/* usage: loop{updbnd_def()}; updallbnds() */
void updallbnds(sp_t *s)
{
    int np = 2*s->nbnds;

    /* sort endpoints & update bndeps*/
    for(int a=0; a<3; a++){
	qsort(s->endpnts[a], np, sizeof(endpnt_t), cmpep);
	for(int i=0; i<np; i++){
	    s->bndeps[s->endpnts[a][i].id] = i;
	}
    }
	
    /* clear old neighborhoods */
    for(int i=0; i<s->ndbnds; i++){
	bnd_t *b = s->bnds + s->dbnds[i];
	b->ns = b->nd = 0;
    }

    /* prune along x axis */
    for(int i=0; i<np; i++){
	if(s->endpnts[0][i].flag & ISMIN){
	    endpnt_t *min = s->endpnts[0] + i; /* endpnts[] are already sorted*/
	    endpnt_t *max = s->endpnts[0] + s->bndeps[min->id+1];
	    int d = min->flag & ISDYNAMIC;
	    bnd_t *b=s->bnds + min->id/6;
	    for(endpnt_t *p=min+1; p<max; p++){ /* min<p<max */
		if((d || (p->flag & ISDYNAMIC)) && (p->flag & ISMIN)){
		    /* endpnts[i] are sorted => comparing indices is enough */
		    unsigned short *pid = s->bndeps+p->id; /* PLACE 1 */
		    unsigned short *minid = s->bndeps+min->id;
		    if(pid[2]<minid[3] && minid[2]<pid[3] && pid[4]<minid[5] && minid[4]<pid[5]){
			addtouchers(s,p->id/6,min->id/6,p->flag,min->flag);
		    }
		}
	    }
	}
    }

    /* update stabbing numbers */
    int st=0;
    endpnt_t *p;
    for(p=s->endpnts[0]+np-1; p>s->endpnts[0]; p--){
	endpnt_t *next = s->endpnts[0] + s->bndeps[p[-1].id];
	if(p->flag & ISMIN){
	    p->stab = p[1].stab-1;
	    --st;
	}
	else{
	    if(p == s->endpnts[0]+2*s->nbnds-1){
		p->stab=1;
	    }
	    else{
		p->stab = p[1].stab+1;
	    }
	    ++st;
	}
	if(!st) continue;
	for(endpnt_t *stab = p-1; stab>next; stab--){
	    stab->stab += st;
	}
    }
    endpnt_t *first = s->endpnts[0] + s->bndeps[p->id];
    first->stab = first[1].stab-1;
}


static inline int qrydbnd(const sp_t *s, int bid, unsigned short touchers[])
{
    int n=0;
    bnd_t *b=s->bnds+bid;
    for(int i=0;i<b->nd;i++){
	touchers[n++] = b->touchers[i];
    }
    for(int i=31;i>31-b->ns;i--){
	touchers[n++] = b->touchers[i];
    }
    return n;
}



static inline int qrysbnd(const sp_t *s, int bid, unsigned short touchers[])
{
    int n=0;
    endpnt_t *p0,*p1;
    p0=s->endpnts[0]+s->bndeps[6*bid];
    p1=s->endpnts[0]+s->bndeps[6*bid+1];
    unsigned short *qid=s->bndeps + p0->id;
    for(endpnt_t *p=p0+1; p<p1; p++){
	if(p->flag & ISMIN){
	    unsigned short *pid=s->bndeps+p->id;
	    if(qid[2]<pid[3] && pid[2]<qid[3] && 
	       qid[4]<pid[5] && pid[4]<qid[5]){
		touchers[n++]=p->id/6;
	    }
	}
    }
    int st=p0->stab;
    endpnt_t *p=p0-1;	/* no need to check whether p<s->endpnts[0], since if p0==s->endpnts[0], then p0->stab==0 */
    while(st){
	if(p->flag & ISMIN){
	    if(s->bndeps[p->id+1]>qid[0]){
		unsigned short *pid=s->bndeps+p->id;
		if(qid[2]<pid[3] && pid[2]<qid[3] && 
		   qid[4]<pid[5] && pid[4]<qid[5]){
		    touchers[n++]=p->id/6;
		}
		st--;
	    }
	}
	p--;
    }
    return n;
}


int bndqry(const sp_t *s, int bid, unsigned short touchers[])
{
    bnd_t *b=s->bnds+bid;
    if(b->touchers) 
	qrydbnd(s,bid,touchers);
    else
	qrysbnd(s,bid,touchers);
}


static int testupd(const sp_t *s, endpnt_t *pnt, unsigned short pid[6], float value, float value2, unsigned short touchers[], int n)
{
    int ismax = !(pnt->flag & ISMIN);
    int d = ismax ? 1 : -1;

    unsigned short a[3];
    a[0] = pnt->id % 6 / 2;
    a[1] = 2*(a[0]+1) % 6;
    a[2] = 2*(a[0]+2) % 6;
    int bid=pnt->id/6*6;
    endpnt_t *first=s->endpnts[a[0]]-1;
    endpnt_t *last = s->endpnts[a[0]]+2*s->nbnds;

    for(endpnt_t *p=pnt+d; p>first && p<last; p+=d){
	if((value>p->value)==ismax){ /* max increase / min decrease */
	    unsigned short *id=s->bndeps + p->id/6*6;
	    int overlap = pid[a[1]] < id[a[1]+1] && 
		pid[a[1]+1] > id[a[1]] &&
		pid[a[2]] < id[a[2]+1] &&
		pid[a[2]+1] > id[a[2]];
	    if(!ismax == !(p->flag & ISMIN)){ /* max vs. min / min vs. max */
		int lt = value2 < s->endpnts[a[0]][s->bndeps[p->id^1]].value;
		if(overlap && (lt == ismax)){
		    touchers[n++]=p->id/6;
		}
	    }
	}
	else{	/* done */
	    pid[2*a[0]+ismax]=p-d-s->endpnts[a[0]];  /* for later use */
	    break;
	}
    }
    return n;
}



int bndqry4(const sp_t *s, int bid, const float bnd[2][3], unsigned short touchers[])
{
    /* the given bound should contain the old bound */
    int n=bndqry(s,bid,touchers);
    unsigned short pid[6];
    for(int i=0;i<6;i++){
	pid[i]=s->bndeps[6*bid+i];
    }
    for(int i=0;i<6;i++){
	int a=i/2;
	int ismax=i%2;
	endpnt_t *p=s->endpnts[a]+pid[i];
	if(fabs(bnd[ismax][a]-p->value)>1.0e-6f){
	    n = testupd(s, p, pid, bnd[ismax][a], bnd[!ismax][a], touchers, n);
	}
    }
    return n;
}


void flushsp(sp_t *s)
{
    if(s->rmbuf.nbnds){
	rmbnds(s, s->rmbuf.bids, s->rmbuf.nbnds);
	if(s->rmbuf.cap > BUFSIZE){
	    free(s->rmbuf.bids);
	    s->rmbuf.bids = malloc(BUFSIZE*sizeof(unsigned short));
	    s->rmbuf.cap = BUFSIZE;
	}
	s->rmbuf.nbnds = 0;
    }
    if(s->addbuf.nbnds){
	addbnds(s, s->addbuf.endpnts, s->addbuf.bids, s->addbuf.owners, s->addbuf.nbnds);
	if(s->addbuf.cap > BUFSIZE){
	    free(s->addbuf.endpnts[0]);
	    free(s->addbuf.endpnts[1]);
	    free(s->addbuf.endpnts[2]);
	    free(s->addbuf.bids);
	    free(s->addbuf.owners);
	    int size = 2*BUFSIZE*sizeof(endpnt_t);
	    s->addbuf.endpnts[0] = malloc(size);
	    s->addbuf.endpnts[1] = malloc(size);
	    s->addbuf.endpnts[2] = malloc(size);
	    s->addbuf.bids = malloc(BUFSIZE*sizeof(unsigned short*));
	    s->addbuf.owners = malloc(BUFSIZE*sizeof(void*));
	    s->addbuf.nbnds = 0;
	    s->addbuf.cap = BUFSIZE;
	}
	s->addbuf.nbnds = 0;
    }
}


void addbnd(sp_t *s, float bnd[2][3], int isdynamic, void *owner, unsigned short *bid)
{
    if(s->addbuf.cap == s->addbuf.nbnds){
	s->addbuf.cap*=2;
	for(int a=0;a<3;a++){
	    s->addbuf.endpnts[a] = realloc(s->addbuf.endpnts[a], 2*s->addbuf.cap*sizeof(endpnt_t));
	}
	s->addbuf.bids = realloc(s->addbuf.bids, s->addbuf.cap*sizeof(unsigned short*));
	s->addbuf.owners = realloc(s->addbuf.owners, s->addbuf.cap*sizeof(void*));
    }
    isdynamic *= ISDYNAMIC;
    for(int a=0;a<3;a++){
	endpnt_t *p = s->addbuf.endpnts[a] + 2*s->addbuf.nbnds;
	p[0].value = bnd[0][a];
	p[0].flag = ISMIN | isdynamic;
	p[1].value = bnd[1][a];
	p[1].flag = isdynamic;
    }
    s->addbuf.bids[s->addbuf.nbnds]=bid;
    s->addbuf.owners[s->addbuf.nbnds++]=owner;
}

void rmbnd(sp_t *s, unsigned short bid)
{
    if(s->rmbuf.cap == s->rmbuf.nbnds){
	s->rmbuf.cap *=2;
	s->rmbuf.bids = realloc(s->rmbuf.bids, s->rmbuf.cap*sizeof(unsigned short));
    }
    s->rmbuf.bids[s->rmbuf.nbnds++]=bid;
}


void getbnd(sp_t *s, unsigned short bid, float bnd[2][3])
{
    unsigned short *pid = s->bndeps+6*bid;
    bnd[0][0]=s->endpnts[0][pid[0]].value;
    bnd[1][0]=s->endpnts[0][pid[1]].value;
    bnd[0][1]=s->endpnts[1][pid[2]].value;
    bnd[1][1]=s->endpnts[1][pid[3]].value;
    bnd[0][2]=s->endpnts[2][pid[4]].value;
    bnd[1][2]=s->endpnts[2][pid[5]].value;
}


static inline endpnt_t *searchpnt(const endpnt_t *base, int np, float x)
{
    int low=0;
    int high=np-1;
    int mid=high>>1;
    while(low+1<high){
	float midx=base[mid].value;
	if(x==midx) return base+mid;
	if(x<midx){
	    high=mid;
	}
	else{
	    low=mid;
	}
	mid=(low+high)>>1;
    }
    return base+low;
}


int rngqry(const sp_t *s, float bnd[2][3], unsigned short touchers[])
{
    if(!s->nbnds) return 0;
    int n=0;
    endpnt_t *p0,*p1;
    endpnt_t *last = s->endpnts[0]+s->nbnds*2-1;
    if(bnd[0][0]>=last->value || bnd[1][0]<=s->endpnts[0]->value) return 0;
    if(bnd[0][0]<=s->endpnts[0]->value){
	p0=s->endpnts[0];
    }
    else{
	p0 = searchpnt(s->endpnts[0], s->nbnds*2, bnd[0][0])+1;
    }
    if(bnd[1][0]>=s->endpnts[0][s->nbnds*2-1].value){
	p1=s->endpnts[0]+s->nbnds*2-1;
    }
    else{
#if 1
	p1 = searchpnt(p0, s->endpnts[0]+s->nbnds*2-p0, bnd[1][0]);
#endif
#if 0
	p1 = searchpnt(s->endpnts[0], s->nbnds*2, bnd[1][0]);
#endif
    }
    /* TODO: search bnd's index for y and z axis? */
    for(endpnt_t *p=p0; p<=p1; p++){
	if(p->flag & ISMIN){
	    unsigned short *pid=s->bndeps + p->id/6*6;
	    if(s->endpnts[1][pid[2]].value<bnd[1][1] &&
	       s->endpnts[1][pid[3]].value>bnd[0][1] &&
	       s->endpnts[2][pid[4]].value<bnd[1][2] &&
	       s->endpnts[2][pid[5]].value>bnd[0][2]){
		touchers[n++]=p->id/6;
	    }
	}
    }
    endpnt_t *p=p0-1;
    int st=p0->stab;
    while(st){
	if(p->flag & ISMIN){
	    if(s->endpnts[0][s->bndeps[p->id+1]].value>bnd[0][0]){
		unsigned short *pid=s->bndeps + p->id/6*6;
		if(s->endpnts[1][pid[2]].value<bnd[1][1] &&
		   s->endpnts[1][pid[3]].value>bnd[0][1] &&
		   s->endpnts[2][pid[4]].value<bnd[1][2] &&
		   s->endpnts[2][pid[5]].value>bnd[0][2]){
		    touchers[n++]=p->id/6;
		}
		st--;
	    }
	}
	p--;
    }
    return n;
}


int pntqry(const sp_t *s, v3_t pnt, unsigned short touchers[])
{
    if(pnt.x<=s->endpnts[0]->value) return 0;
    if(pnt.x>=s->endpnts[0][s->nbnds*2-1].value) return 0;
    int n=0;
    endpnt_t *p=searchpnt(s->endpnts[0], s->nbnds*2, pnt.x);
    int st=p[1].stab;
    while(st){
	if(p->flag & ISMIN){
	    if(s->endpnts[0][s->bndeps[p->id+1]].value > pnt.x){
		unsigned short *pid=s->bndeps + p->id/6*6;
		if(s->endpnts[1][pid[2]].value < pnt.y &&
		   s->endpnts[1][pid[3]].value > pnt.y &&
		   s->endpnts[2][pid[4]].value < pnt.z &&
		   s->endpnts[2][pid[5]].value > pnt.z){
		    touchers[n++]=p->id/6;
		}
		st--;
	    }
	}
	p--;
    }
    return n;
}



static inline int pntinbnd(v3_t p, float bnd[2][3])
{
    return p.x<bnd[1][0] && p.x>bnd[0][0] &&
	p.y<bnd[1][1] && p.y>bnd[0][1] &&
	p.z<bnd[1][2] && p.z>bnd[0][2];
		
}



/*

  min                    max
  org --- | -----------end        |
  0    t[0]            1        t[1]
  |--------  ---------|
  \/
  dir



  min             max
  |       org --- | -----------end
  t[0]      0    t[1]            1
  |--------  ---------|
  \/
  dir



  0  1  min     max
  0    +min 1   max
  0    +min    -max   1
  min 0 1 max
  min 0  -max   1
  min     max 0 1


  1  0  min     max
  1     min- 0  max
  1     min-    max+  0
  min 1 0 max
  min 1   max+  0
  min     max 1 0


*/
static inline int clip1d(float org, float end, float min, float max, float *t0, float *t1)
{
    float dir=end-org;
    float t[2];
    int ret;
    if(dir==0){
	if(org<min || org>max) return -1;
	return 2;
    }
    else if(dir>0){
	ret = 0;	/* hit min plane */
	t[0] = (min-org)/dir;
	t[1] = (max-org)/dir;
    }
    else{
	ret = 1;	/* hit max plane */
	t[0] = (max-org)/dir;
	t[1] = (min-org)/dir;
    }
    if(t[0]>*t1 || t[1]<*t0) return -1;
    if(t[1]<*t1){
	*t1=t[1];
    }
    if(t[0]>*t0){
	*t0=t[0];
	return ret;	/* new first hit plane */
    }
    else
	return 2;	/* not the first hit */
}

int bndclipray(const float bnd[2][3], float begin[3], float end[3])
{
    float t[2]={0,1};
    int hitplane = -1;
    for(int a=0;a<3;a++){
	int hit = clip1d(begin[a],end[a],bnd[0][a],bnd[1][a],t,t+1);
	if(hit<0) return -1;
	if(hit<2) hitplane=a*2+hit;
    }
    for(int a=0;a<3;a++){
	float d=end[a]-begin[a];
	end[a] = begin[a]+t[1]*d;
	begin[a] += t[0]*d;
    }
    return hitplane;
}



/* bounds in which the ray starts are not considered as hit */
void initray(const sp_t *s, v3_t begin, v3_t end, ray_t *r)
{
    int max=s->nbnds*2-1;
    float bnd[2][3] = {s->endpnts[0][0].value,s->endpnts[1][0].value,s->endpnts[2][0].value,
		       s->endpnts[0][max].value,s->endpnts[1][max].value,s->endpnts[2][max].value};

    int in = pntinbnd(begin,bnd);
    if(!in){
	int hit = bndclipray(bnd,&begin,&end);
	if(hit<0){
	    r->done = 1;
	    return;
	}
	else{
	    r->hit = hit;
	}
    }
    else{
	r->hit = -2;
    }

    glPointSize(10);
    glBegin(GL_POINTS);
    glColor3f(1,0,0);
    glVertex3fv(&begin);
    glColor3f(0,0,1);
    glVertex3fv(&end);
    glEnd();
    glPointSize(1);

    r->t = 0;
    r->o = begin;
    r->dir = sub(end,begin);
    for(int i=0;i<3;i++){
	if(((float *)&r->dir)[i]==0){
	    ((float *)&r->dir)[i]=1.0e-20f;
	}

	if(r->hit/2 == i){
	    int ismax = r->hit%2;
	    r->id[i] = ismax*(s->nbnds*2-2);
	}
	else{
	    r->id[i] = searchpnt(s->endpnts[i], s->nbnds*2, ((float *)&begin)[i]) - s->endpnts[i];
	}

	r->delta[i] = ((float *)&r->dir)[i]>0 ? 1 : -1;
#if 1
	if(r->delta[i]<0){
	    r->id[i]++;
	}
#endif
    }
    r->done = 0;
}

int rayqry(const sp_t *s, ray_t *r, unsigned short touchers[])
{
    if(r->done) return 0;
    if(r->t>=1.0f) return 0;
    if(r->hit>=0){
	int a[3];
	a[0] = r->hit/2;
	a[1] = (a[0]+1)%3;
	a[2] = (a[0]+2)%3;
	int ismax = r->hit%2;
	endpnt_t *p = s->endpnts[a[0]] + ismax*(s->nbnds*2-1);
	unsigned short *pid = s->bndeps + p->id/6*6;
	r->hit=-2;
	float *o = (float *)&r->o;
	if(o[a[1]] > s->endpnts[a[1]][pid[2*a[1]]].value &&
	   o[a[1]] < s->endpnts[a[1]][pid[2*a[1]+1]].value &&
	   o[a[2]] > s->endpnts[a[2]][pid[2*a[2]]].value &&
	   o[a[2]] < s->endpnts[a[2]][pid[2*a[2]+1]].value){
	    touchers[0] = p->id/6;
	    return 1;
	}
    }

    int facing[3] = {ISMIN*(r->delta[0]>0), ISMIN*(r->delta[1]>0), ISMIN*(r->delta[2]>0)};
    float co[3] = {r->o.x*r->t, r->o.y*r->t, r->o.z*r->t};
    int id[3] = {r->id[0], r->id[1], r->id[2]};

    /* for each axis, find the next front-facing plane. */
    int iii=0;
    while(1){
	float d[3];
	float t[3];
	float tmin;
	int amin;
	iii++;
	assert(iii<10000);

	/* next x */
	if(id[0]+r->delta[0]<0 || id[0]+r->delta[0]>=s->nbnds*2) return 0;
	d[0] = s->endpnts[0][id[0]+r->delta[0]].value;
	t[0] = (d[0] - r->o.x)/r->dir.x;

	/* next y */
	if(id[1]+r->delta[1]<0 || id[1]+r->delta[1]>=s->nbnds*2) return 0;
	d[1] = s->endpnts[1][id[1]+r->delta[1]].value;
	t[1] = (d[1] - r->o.y)/r->dir.y;

	/* next z */
	if(id[2]+r->delta[2]<0 || id[2]+r->delta[2]>=s->nbnds*2) return 0;
	d[2] = s->endpnts[2][id[2]+r->delta[2]].value;
	t[2] = (d[2] - r->o.z)/r->dir.z;

	if(t[0]<t[1]){
	    if(t[0]<t[2]){
		tmin = t[0];
		amin = 0;
	    }
	    else{
		tmin = t[2];
		amin = 2;
	    }
	}
	else{
	    if(t[2]>t[1]){
		tmin = t[1];
		amin = 1;
	    }
	    else{
		tmin = t[2];
		amin = 2;
	    }
	}

	if(tmin>=1.0f) {
	    r->done = 1;
	    return 0;
	}
	if(id[amin]>=2*s->nbnds || id[amin]<0){
	    r->done = 1;
	    return 0;
	}
	id[amin] += r->delta[amin];
	endpnt_t *p = s->endpnts[amin] + id[amin];
	if((p->flag & ISMIN) == facing[amin]){
	    int bid=p->id/6;
	    unsigned short *pid = s->bndeps + bid*6;
	    int a[2] = {(amin+1)%3, (amin+2)%3};
	    float min[2]={s->endpnts[a[0]][pid[2*a[0]]].value, s->endpnts[a[1]][pid[2*a[1]]].value};
	    float max[2]={s->endpnts[a[0]][pid[2*a[0]+1]].value,s->endpnts[a[1]][pid[2*a[1]+1]].value};
	    float *o = (float *)&r->o;
	    float *d = (float *)&r->dir;
	    float value[2] = {o[a[0]]+d[a[0]]*tmin, o[a[1]]+d[a[1]]*tmin};
	    if(value[0]>min[0] && value[0]<max[0] && value[1]>min[1] && value[1]<max[1]){
		touchers[0] = bid;
		r->t = tmin;
		r->id[0]=id[0];
		r->id[1]=id[1];
		r->id[2]=id[2];
		return 1;
	    }
	}

    }
    return 0;
}


void printbnd(sp_t *s, unsigned short bid)
{
    float bnd[2][3];
    getbnd(s,bid,bnd);
    printf("(%f,%f,%f)-(%f,%f,%f)\n", bnd[0][0],bnd[0][1],bnd[0][2],bnd[1][0],bnd[1][1],bnd[1][2]);
}

void printbndinfo(sp_t *s, unsigned short bid)
{
    unsigned short *pid=s->bndeps+6*bid;
    printf("bnd%u: (%u %u)(%u %u)(%u %u)\n",bid,pid[0],pid[1],pid[2],pid[3],pid[4],pid[5]);
}

void printspinfo(sp_t *s)
{
    printf("nbnds=%i\nnavailbnds=%i\ncap=%u\nndbnds=%i\ndcap=%u\n",s->nbnds,s->navailbnds,s->cap,s->ndbnds,s->dcap);
}

static void drawbndwire(sp_t *s, unsigned short bid, float color[3])
{
    static float oldcolor[3];

    float bnd[2][3];
    getbnd(s,bid,bnd);

    if(oldcolor[0]!=color[0] || oldcolor[1]!=color[1] || oldcolor[2]!=color[2]){
	glColor3fv(color);
	oldcolor[0]!=color[0];
	oldcolor[1]!=color[1];
	oldcolor[2]!=color[2];
    }

    glBegin(GL_LINES);
	
    /* x */
    glVertex3f(bnd[0][0], bnd[0][1], bnd[0][2]);
    glVertex3f(bnd[1][0], bnd[0][1], bnd[0][2]);
	
    glVertex3f(bnd[0][0], bnd[0][1], bnd[1][2]);
    glVertex3f(bnd[1][0], bnd[0][1], bnd[1][2]);
	
    glVertex3f(bnd[0][0], bnd[1][1], bnd[0][2]);
    glVertex3f(bnd[1][0], bnd[1][1], bnd[0][2]);
	
    glVertex3f(bnd[0][0], bnd[1][1], bnd[1][2]);
    glVertex3f(bnd[1][0], bnd[1][1], bnd[1][2]);

    /* y */
    glVertex3f( bnd[0][0], bnd[0][1], bnd[0][2]);
    glVertex3f(bnd[0][0], bnd[1][1], bnd[0][2] );
	
    glVertex3f(bnd[1][0],bnd[0][1], bnd[0][2] );
    glVertex3f(bnd[1][0],bnd[1][1], bnd[0][2] );
	
    glVertex3f(bnd[0][0],bnd[0][1], bnd[1][2] );
    glVertex3f(bnd[0][0],bnd[1][1], bnd[1][2] );
	
    glVertex3f( bnd[1][0],bnd[0][1], bnd[1][2]);
    glVertex3f(bnd[1][0],bnd[1][1], bnd[1][2] );

    /* z */
    glVertex3f(bnd[0][0], bnd[0][1],bnd[0][2] );
    glVertex3f(bnd[0][0], bnd[0][1],bnd[1][2] );
	
    glVertex3f(bnd[0][0], bnd[1][1],bnd[0][2]);
    glVertex3f(bnd[0][0], bnd[1][1],bnd[1][2]);
	
    glVertex3f( bnd[1][0], bnd[0][1],bnd[0][2]);
    glVertex3f( bnd[1][0], bnd[0][1],bnd[1][2]);
	
    glVertex3f( bnd[1][0], bnd[1][1],bnd[0][2]);
    glVertex3f( bnd[1][0], bnd[1][1],bnd[1][2]);


    glEnd();
}


void drawbox(float bnd[2][3], float color[3])
{
    static float oldcolor[3];
    if(oldcolor[0]!=color[0] || oldcolor[1]!=color[1] || oldcolor[2]!=color[2]){
	glColor3fv(color);
	oldcolor[0]!=color[0];
	oldcolor[1]!=color[1];
	oldcolor[2]!=color[2];
    }

    glBegin(GL_QUAD_STRIP);
    glVertex3f(bnd[0][0],bnd[0][1],bnd[1][2]);
    glVertex3f(bnd[1][0],bnd[0][1],bnd[1][2]);
    glVertex3f(bnd[0][0],bnd[1][1],bnd[1][2]);
    glVertex3f(bnd[1][0],bnd[1][1],bnd[1][2]);
    glVertex3f(bnd[0][0],bnd[1][1],bnd[0][2]);
    glVertex3f(bnd[1][0],bnd[1][1],bnd[0][2]);
    glVertex3f(bnd[0][0],bnd[0][1],bnd[0][2]);
    glVertex3f(bnd[1][0],bnd[0][1],bnd[0][2]);
    glEnd();

    glBegin(GL_QUAD_STRIP);
    glVertex3f(bnd[1][0],bnd[1][1],bnd[0][2]);
    glVertex3f(bnd[1][0],bnd[1][1],bnd[1][2]);
    glVertex3f(bnd[1][0],bnd[0][1],bnd[0][2]);
    glVertex3f(bnd[1][0],bnd[0][1],bnd[1][2]);
    glVertex3f(bnd[0][0],bnd[0][1],bnd[0][2]);
    glVertex3f(bnd[0][0],bnd[0][1],bnd[1][2]);
    glVertex3f(bnd[0][0],bnd[1][1],bnd[0][2]);
    glVertex3f(bnd[0][0],bnd[1][1],bnd[1][2]);
    glEnd();
}

void drawbnd(sp_t *s, unsigned short bid, float color[3])
{
    float bnd[2][3];
    getbnd(s,bid,bnd);
    drawbox(bnd,color);
}

static unsigned short bndbuf[1024*10];

#if 0

void drawsp(sp_t *s, short activeid)
{
    float dynamic[3]={0,0,1};
    float stat[3]={0,0,.5};
    float dactive[3]={1,0,0};
    float sactive[3]={.5,0,0};
    float dtouch[3]={0,1,0};
    float stouch[3]={0,.5,0};
    float touch4d[3]={1,0,1};
    float active4d[3]={.8,.8,0};
    float range[3]={.3,.5,.7};
    float pntcolor[3]={0,.5,.5};


    extern int glob4d;
    extern int globntouchers;
    extern unsigned short globtouchers[];
    extern int globlight,globpntenable;
    extern v3_t globpnt;
    extern int globnormal;

    for(unsigned short i=0;i<2*s->nbnds;i++){
/* 		if(s->bnds[i].ns == 0xffff) continue; */

	endpnt_t *p=s->endpnts[0]+i;
	int bid=p->id/6;
	if(p->flag & ISMIN){
	    if(p->flag & ISDYNAMIC){
		drawbnd(s,bid,dynamic);
	    }
	    else{
		drawbnd(s,bid,stat);
	    }
	}
    }

    if(activeid>=0){
	int n=bndqry(s,activeid,bndbuf);
	for(int i=0;i<n;i++){
	    if(s->endpnts[0][s->bndeps[bndbuf[i]*6]].flag & ISDYNAMIC){
		drawbnd(s,bndbuf[i],dtouch);
	    }
	    else{
		drawbnd(s,bndbuf[i],stouch);
	    }
	}
		
	if(globnormal){
	    if(s->endpnts[0][s->bndeps[activeid*6]].flag & ISDYNAMIC){
		drawbnd(s,activeid,dactive);
	    }
	    else{
		drawbnd(s,activeid,sactive);
	    }
	}
	else{
	    if(glob4d){
		extern float globbnd4d[2][3];
		drawbox(globbnd4d, active4d);
	    }
	    else if(globlight){
		extern float globlightbnd[2][3];
		float light[3]={.8,.8,0};
		drawbox(globlightbnd,light);
	    }
	    else{
		glPointSize(10);
		glBegin(GL_POINTS);
		glColor3f(.8,.8,0);
		glVertex3fv((float *)&globpnt);
		glEnd();
		glPointSize(1);
	    }
	}
    }

    float *color;
    if(glob4d){
	color=touch4d;
    }
    else if(globlight){
	color=range;
    }
    else if(globpntenable){
	color=pntcolor;
    }
    else{
#include "camera.h"
	extern cam_t *cam;
	extern int globmvr,globmvl,globmvu,globmvd,globmvf,globmvb;
	ray_t ray;
#if 0
	static v3_t begin={0,0,-1};
	static v3_t end={0,0,1};

	if(globmvr) {begin.x+=.1;end.x+=.1;}
	if(globmvl) {begin.x-=.1;end.x-=.1;}
	if(globmvu) {begin.y+=.1;end.y+=.1;}
	if(globmvd) {begin.y-=.1;end.y-=.1;}
	if(globmvf) {begin.z-=.1;end.z-=.1;}
	if(globmvb) {begin.z+=.1;end.z+=.1;}
	initray(s, begin, end, &ray);
#endif


#if 1
	v3_t begin=cam->pos;
	v3_t end=add(cam->pos,scale(cam->axis,-200));
	initray(s, cam->pos, end, &ray);
#endif
	int n=0;
	float raycolor[3] = {1,1,1};
	do{
	    n = rayqry(s,&ray,globtouchers);
	    for(int i=0;i<n;i++){
		drawbnd(s,globtouchers[i],raycolor);
	    }
	    if(n) break;
	}while(n);
	glBegin(GL_LINES);
	glColor3f(1,1,1);
	glVertex3fv(&begin);
	glColor3f(0,0,0);
	glVertex3fv(&end);
	glEnd();
	return;
    }

    for(int i=0;i<globntouchers;i++){
	drawbnd(s,globtouchers[i],color);
    }
	



}

#endif

void xlatbnd(float bnd[2][3], float dx,float dy,float dz)
{
    bnd[0][0]+=dx;
    bnd[1][0]+=dx;
    bnd[0][1]+=dy;
    bnd[1][1]+=dy;
    bnd[0][2]+=dz;
    bnd[1][2]+=dz;
}

void xlatbnd4d(float dest[2][3], float bnd[2][3], float d[3])
{
#if 0
    if(d[0]>0){
	dest[1][0] = bnd[1][0]+d[0];
	dest[0][0] = bnd[0][0];
    }
    else{
	dest[0][0] = bnd[0][0]+d[0];
	dest[1][0] = bnd[1][0];
    }
    if(d[1]>0){
	dest[1][1] = bnd[1][1]+d[1];
	dest[0][1] = bnd[0][1];
    }
    else{
	dest[0][1] = bnd[0][1]+d[1];
	dest[1][1] = bnd[1][1];
    }
    if(d[2]>0){
	dest[1][2] = bnd[1][2]+d[2];
	dest[0][2] = bnd[0][2];
    }
    else{
	dest[0][2] = bnd[0][2]+d[2];
	dest[1][2] = bnd[1][2];
    }
#endif

#if 1
    for(int i=0; i<3; i++){
	int gtz = d[i]>0;
	dest[gtz][i] = bnd[gtz][i] + d[i];
	dest[!gtz][i] = bnd[!gtz][i];
    }
#endif
}


void scalebnd(float bnd[2][3], float sx, float sy, float sz)
{
    float center[3];
    float width[3];
    float s[3]={sx,sy,sz};
    for(int i=0;i<3;i++){
	float c=(bnd[0][i]+bnd[1][i])/2.0f;
	float w=(bnd[1][i]-bnd[0][i])/2.0f;
	bnd[0][i]=c-w*s[i];
	bnd[1][i]=c+w*s[i];
    }
}

void printstab(sp_t *s)
{
    for(int i=0;i<s->nbnds*2;i++){
	endpnt_t *p=s->endpnts[0]+i;
	if(p->flag & ISMIN){
	    printf("[  ");
	}
	else{
	    printf("]  ");
	}
    }
    printf("\n");
    int prevstab=1;
    for(int i=0;i<s->nbnds*2;i++){
	endpnt_t *p=s->endpnts[0]+i;
	printf("%u  ",p->stab);
#if 1
	assert(abs(prevstab-p->stab)<2 || (printf("p->id=%u,s->bndeps[p->id]=%u\n",p->id,s->bndeps[p->id]),freesp(s),0));
	assert((prevstab!=p->stab) || (printf("p->id=%u,s->bndeps[p->id]=%u\n",p->id,s->bndeps[p->id]),freesp(s),0));
#endif
	prevstab=p->stab;
    }
    printf("\n");
}

void xlatsp(sp_t *s, v3_t x)
{
    for(int i=0; i<s->nbnds*2; i++){
	s->endpnts[0][i].value += x.x;
	s->endpnts[1][i].value += x.y;
	s->endpnts[2][i].value += x.z;
    }
}


void setdbnd(sp_t *s, int bid)
{
    unsigned short *pid;
    unsigned short touchers[32];
    int n;
    endpnt_t *p;

    pid = s->bndeps + bid*6;
    p = s->endpnts[0]+pid[0];
    if(p->flag & ISDYNAMIC) return;

    n=bndqry(s, bid, touchers);
    p->flag |= ISDYNAMIC;
    s->bnds[bid].touchers = allocblk(s->mm);

    for(int i=0; i<n; i++){
	endpnt_t *t = s->endpnts[0] + s->bndeps[touchers[i]*6];
	addtouchers(s, bid, touchers[i], p->flag, t->flag);
    }

    if(s->ndbnds >= s->dcap){
	s->dcap *= 2;
	s->dbnds = realloc(s->dbnds, s->dcap*sizeof(unsigned short));
    }
    s->dbnds[s->ndbnds++] = bid;
}

#endif


// OBB-OBB collision detection
#define proj dot

static inline volatile long long rdtsc() {
    register long long tsc __asm__("eax");
//	__asm__ __volatile__ (".byte 15, 49" : : : "eax", "edx");
    __asm__ volatile ("rdtsc" : "=A" (tsc));

    return tsc;
} 

#define BEGIN_TIMING long long begin__ = rdtsc();
#define END_TIMING long long end__ = rdtsc(); printf("%qd cycles\n", end__-begin__);

typedef struct{
    float			proj[3];
    float			face;
} projinfo_t;



static float signtab[] = {-1,1};

typedef struct{
    int type;
    int edge[2];
    int face[2];
    v3_t axis;
    float depth;
    float minaxis;
} bbinfo_t;



enum {INTERSECT, SEPARATE};
float epsilon = 0.001f;
static inline int axisproj_edge(obb_t *b0, obb_t *b1, int a0, int a1, v3_t *t, float *depth,
				float dp[3][3], float dpf[3][3], float adotc[3][2]/* v3_t *_axis  */) {
    float			r0, r1, r=0;
    v3_t			axis;
    int			o0, o1;

    o0 = (a0+1)%3;
    o1 = (a0+2)%3;

    r0 = b0->extent[o0] * dpf[o1][a1]
	+ b0->extent[o1] * dpf[o0][a1];

    r = dp[o1][a1]*adotc[o0][1] - dp[o0][a1]*adotc[o1][1];

    o0 = (a1+1)%3;
    o1 = (a1+2)%3;

    r1 = b1->extent[o0] * dpf[a0][o1]
	+ b1->extent[o1] * dpf[a0][o0];

    r += dp[a0][o1]*adotc[o0][0] - dp[a0][o0]*adotc[o1][0];

    r=fabsf(r);

    *depth = r0+r1-r;
    if (*depth<0) {
	return SEPARATE;
    }
    else {
	float len=sqrtf(1.0f - dp[a0][a1]*dp[a0][a1]);
	*depth /= len;
	return INTERSECT;
    }
}

/*
  axis is b0's axis[a]

*/


static inline int axisproj_face(obb_t *b0, obb_t *b1, int a, v3_t *t, float *depth, projinfo_t *p)
{
    float			r0, r1, r;
    v3_t			axis;

    r0 = b0->extent[a];
    axis = b0->axis[a];

#if 0
    p->proj[0] = dot(axis, b1->axis[0]);
    p->proj[1] = dot(axis, b1->axis[1]);
    p->proj[2] = dot(axis, b1->axis[2]);

#else
    p->proj[0] = axis.x*b1->axis[0].x  +  axis.y*b1->axis[0].y  + axis.z*b1->axis[0].z;
    p->proj[1] = axis.x*b1->axis[1].x  +  axis.y*b1->axis[1].y  + axis.z*b1->axis[1].z;
    p->proj[2] = axis.x*b1->axis[2].x  +  axis.y*b1->axis[2].y  + axis.z*b1->axis[2].z;
#endif

    p->face = proj(*t, axis);

    r1 = b1->extent[0]*fabsf(p->proj[0]) + b1->extent[1]*fabsf(p->proj[1]) + b1->extent[2]*fabsf(p->proj[2]);

    r = fabsf(p->face);


    /* DEBUG! */
    if(r*.95f>r0+r1)
	return SEPARATE;
    else{
	*depth = r0+r1-r;
	return INTERSECT;
    }


}


static void closestfaces(int minaxis, projinfo_t *p, bbinfo_t *b)
{
    int			maxaxis;
    float			maxproj, proj;

    maxproj = fabsf(p->proj[0]);
    maxaxis = 0;

    for(int i=1; i<3; i++){
	proj = fabsf(p->proj[i]);
	if(proj>maxproj){
	    maxproj = proj;
	    maxaxis = i;
	}
    }

    int neg=minaxis>2;
    minaxis %= 3;
    b->face[neg] = minaxis + 3*((p->face<0)^neg);
    b->face[!neg] = maxaxis + 3*((p->proj[maxaxis]*p->face>=0)^neg);
}


/*
  return: boolean

  type0: face-vert
  type1: edge-edge
*/
enum{TYPE_FV, TYPE_EE};

static int boxboxtest(obb_t *b0, obb_t *b1,float dp[3][3], bbinfo_t *bbinfo)
{
    float			depth[15];
    projinfo_t		projinfo[6];
    int			minaxis;

    v3_t t = sub(b1->center, b0->center);
    float *d = depth;

    for(int i=0; i<3; i++) {
	if(axisproj_face(b0, b1, i, &t, d++, projinfo+i) == SEPARATE)
	    return SEPARATE;
    }

    for(int i=0; i<3; i++) {
	if(axisproj_face(b1, b0, i, &t, d++, projinfo+i+3) == SEPARATE)
	    return SEPARATE;
    }

    int naxes=15;
    float dpf[3][3];
    float adotc[3][2];
    for(int i=0; i<3; i++) {
	for(int j=0; j<3; j++) {
	    float d = dot(b0->axis[i],b1->axis[j]);
	    float fd = fabsf(d);
	    if(0.9995f<fd){
		naxes=6;
	    }
	    dp[i][j]=d;
	    dpf[i][j]=fd;
	}
	adotc[i][0] = dot(b1->axis[i], b0->center);
	adotc[i][1] = dot(b0->axis[i], b1->center);
    }

    if(naxes==6)
	goto omit_edges;


    for(int i=0;i<3;i++){
	for(int j=0;j<3;j++){
	    if(axisproj_edge(b0, b1, i, j, &t, d++, dp, dpf, adotc) == SEPARATE)
		return SEPARATE;
	}
    }

  omit_edges:

    minaxis=0;
    for(int i=1; i<naxes; i++){
	if(depth[i]<depth[minaxis]){
	    minaxis=i;
	}
    }

    if(minaxis<6){
	closestfaces(minaxis, projinfo + minaxis, bbinfo);
	bbinfo->type = TYPE_FV;
    }
    else {
	int i,j;
	i = (minaxis-6)/3;
	j = (minaxis-6)%3;
	bbinfo->axis = normalize(cross(b0->axis[i], b1->axis[j]));
	bbinfo->edge[0]=i;
	bbinfo->edge[1]=j;
	bbinfo->type = TYPE_EE;
    }
    bbinfo->depth = depth[minaxis];
    bbinfo->minaxis = minaxis;

    return INTERSECT;
}



static inline float dot3(float a[3], float b[3])
{
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}


static void rotate2box(v3_t *w, obb_t *box)
{
    float x,y,z;
    x = dot3(box->axis, w);
    y = dot3(box->axis+1, w);
    z = dot3(box->axis+2, w);

    w->x=x; w->y=y; w->z=z;
}

static inline void interpolate(float vint[3], float vtx[3], float vtx2[3], float t)
{
    for(int i=0; i<3; i++)
	vint[i] = vtx[i] + t*(vtx2[i]-vtx[i]);
}


static void transpoly(obb_t *host, obb_t *box, int face, float R[3][3], int transpose, v3_t poly[])
{
    int o0 = (face+1)%3;
    int o1 = (face+2)%3;

    /* the other box's space */
    float axis0[3], axis1[3];

/*
  axis0[o0]=box->extent[o0];
  axis1[o1]=box->extent[o1];

*/

    /* T is in world space; translate it into the other box's space */
	
    v3_t T = sub(box->center, host->center);
    float s=signtab[face<3];
    face%=3;
    T = add(T, scale(box->axis[face], s*box->extent[face]));
    rotate2box(&T, host);


    /* box -> world -> the other box */
    if(transpose){
	for(int i=0; i<3; i++){
	    axis0[i] = R[o0][i]*box->extent[o0];
	    axis1[i] = R[o1][i]*box->extent[o1];
	}
    }
    else{
	for(int i=0; i<3; i++){
	    axis0[i] = R[i][o0]*box->extent[o0];
	    axis1[i] = R[i][o1]*box->extent[o1];
	}
    }

    /*
      T - (a0 + a1)  -- (0,0)
      T + a0 - a1    -- (1,0)
      T + a0 + a1    -- (1,1)
      T - (a0 - a1)  -- (0,1)
    */

    v3_t a0adda1 = {axis0[0]+axis1[0], axis0[1]+axis1[1], axis0[2]+axis1[2]};
    v3_t a0suba1 = {axis0[0]-axis1[0], axis0[1]-axis1[1], axis0[2]-axis1[2]};

    poly[0] = sub(T, a0adda1);
    poly[1] = add(T, a0suba1);
    poly[2] = add(T, a0adda1);
    poly[3] = sub(T, a0suba1);

}



static inline void copyv3(float dst[3], float src[3])
{
    dst[0]=src[0];
    dst[1]=src[1];
    dst[2]=src[2];
}


/*
  for host face:


  (axis<<2) | (signbit<<1) | signbit   (-1 ==> signbit ==1)

  face0:
  (-extenty, extenty, -extentz, extentz) => (1001, 1000, 0110, 0100) == (9,8,6,4)
  face3:
  (-extenty, extenty, -extentz, extentz) => (1011, 1010, 0111, 0101) == (11,10,7,5)

  face1:
  (-extentz, extentz, -extentx, extentx) => (0001, 0000, 1010, 1000) == (1,0,10,8)
  face4:
  (-extentz, extentz, -extentx, extentx) => (0011, 0010, 1011, 1001) == (3,2,11,9)

  face2:
  (-extentx, extentx, -extenty, extenty) => (0101, 0100, 0010, 0000) == (5,4,2,0)
  face5:
  (-extentx, extentx, -extenty, extenty) => (0111, 0110, 0011, 0001) == (7,6,3,1);

*/


static unsigned hosttab[][4] = {
    {9<<12, 8<<12,  6<<12,  4<<12},
    {1<<12,  0<<12,   10<<12, 8<<12},
    {5<<12,  4<<12,   2<<12,  0<<12},
    {11<<12,  10<<12,   7<<12,  5<<12},
    {3<<12,  2<<12,   11<<12,  9<<12},
    {7<<12,  6<<12,   3<<12,  1<<12}
};


/*
  for guest face: (3 2 4 1)
  face0:
  (-extentz, extenty, extentz, -extenty)
  face3:
  (-extentz, extenty, extentz, -extenty)

  face1:
  (-extentx, extentz, extentx, -extentz)
  face4:
  (-extentx, extentz, extentx, -extentz)

  face2:
  (-extenty, extentx, extenty, -extentx)
  face5:
  (-extenty, extentx, extenty, -extentx)

*/

static unsigned guesttab[][4] = {
    {6<<12, 8<<12,  4<<12,  9<<12},
    {10<<12,  0<<12,   8<<12, 1<<12},
    {2<<12,  4<<12,   0<<12,  5<<12},
    {7<<12,  10<<12,   5<<12,  11<<12},
    {11<<12,  2<<12,   9<<12,  3<<12},
    {3<<12,  6<<12,   1<<12,  7<<12}
};


static void guestfeats(int id, int face, unsigned feats[4])
{
    for(int i=0; i<4; i++){
	unsigned help=guesttab[face][i]>>12;
	feats[i] = guesttab[face][i] | id;
    }
}

static void hostfeats(int id, int face, unsigned feats[4])
{
    for(int i=0; i<4; i++){
	feats[i] = hosttab[face][i] | id;
    }
}



void getfeats(unsigned feature, unsigned feats[])
{
    feats[0] = (feature & 0x0fff0000)>>16;
    feats[1] = (feature & 0xf0000000)>>28;
    feats[2] = (feature & 0xfff);
    feats[3] = (feature & 0xf000)>>12;
}


#define OBJID_MASK 0xfffffff
static int clippoly(obb_t *host, int face, int hostid, float vtx[][3], unsigned ft[])
{
//BEGIN_TIMING
    /*
      host edge id => higher 16 bits
      guest edge id => lower l6 bits
    */

    int nv=4;
    int nv2=0;
    float vtx2[8][3];
    unsigned hf[2][2];
    hostfeats(hostid,face,hf);

    unsigned ft2[8];

    for(int k=1; k<3; k++){
	int a=(face+k)%3;

	/* clipped by -host->extent[a] */
	nv2=0;
	for(int i=0; i<nv; i++){
	    int in=vtx[i][a]>-host->extent[a];
	    if(in){
		copyv3(vtx2[nv2], vtx[i]);
		ft2[nv2++]=ft[i];
	    }
		
	    int nexti=(i+1)%nv;
	    int nextin=vtx[nexti][a]>-host->extent[a];
	    if(in!=nextin){
		float t= (-host->extent[a]-vtx[i][a])/(vtx[nexti][a]-vtx[i][a]);
		interpolate(vtx2[nv2], vtx[i], vtx[nexti], t);
		if(nextin){
		    ft2[nv2++]=ft[i];
		}
		else{
		    ft2[nv2++]=hf[k-1][0];
		}
	    }
	}

	/* clipped by host->extent[a] */
	nv=0;
	for(int i=0; i<nv2; i++){
	    int in=vtx2[i][a]<host->extent[a];
	    if(in){
		copyv3(vtx[nv], vtx2[i]);
		ft[nv++]=ft2[i];
	    }

	    int nexti=(i+1)%nv2;
	    int nextin=vtx2[nexti][a]<host->extent[a];
	    if(in!=nextin){
		float t= (host->extent[a]-vtx2[i][a])/(vtx2[nexti][a]-vtx2[i][a]);
		interpolate(vtx[nv], vtx2[i], vtx2[nexti], t);
		if(nextin){
		    ft[nv++]=ft2[i];
		}
		else{
		    ft[nv++]=hf[k-1][1];
		}
	    }		
	}
    }

//END_TIMING

    /* contact features */
    if(nv){
	for(int i=0; i<nv; i++){
	    int previ = (i-1+nv)%nv;
	    if(ft[i]<ft[previ]){
		ft2[i] = (ft[previ]<<16) | ft[i];
	    }
	    else{
		ft2[i] = (ft[i]<<16) | ft[previ];
	    }
	}
	memcpy(ft, ft2, nv*sizeof(unsigned));
    }

    return nv;
}



v3_t obb2world(obb_t *obb, float vtx[3])
{
    v3_t w;

    w.x = obb->axis[0].x*vtx[0] + obb->axis[1].x*vtx[1] + obb->axis[2].x*vtx[2] + obb->center.x;
    w.y = obb->axis[0].y*vtx[0] + obb->axis[1].y*vtx[1] + obb->axis[2].y*vtx[2] + obb->center.y;
    w.z = obb->axis[0].z*vtx[0] + obb->axis[1].z*vtx[1] + obb->axis[2].z*vtx[2] + obb->center.z;

    return w;
}


static int facecontacts(obb_t *box, int face, float vtx[][3], unsigned features[], int nv, collinfo_t *col, 
			int id0, int id1, float invm0, float invm1)
{
    int n=0;
    int o[2]={(face+1)%3, (face+2)%3};
    int back=face>=3;
    float s=signtab[back];
    int f=face%3;
    for(int i=0; i<nv; i++){
	float depth = box->extent[f]/*  *.95f  */+ s*vtx[i][f];
	if(depth>=0/*  -.01*box->extent[f] */){
	    collinfo_t *c=col+n;
	    /* transform to world space */
	    c->tangent[0] = box->axis[o[back]];
	    c->tangent[1] = box->axis[o[!back]];
	    c->normal = scale(box->axis[f],-s);
	    c->depth = -depth;
	    c->poi[1] = obb2world(box, vtx[i]);
	    c->poi[0] = c->poi[1];//add(c->poi[1], scale(c->normal, c->depth));
	    c->id[0] = id0;
	    c->id[1] = id1;
	    c->feature = features[i];
	    n++;
	}
    }
	
    if(n){
	float fakemass;
	if (invm0 == 0 && invm1 == 0) {
	    fakemass = 1000000.0f;
	}
	else {
	    fakemass =  invm0 > invm1 ? 1.0 / invm0 / (float)n : 1.0 / invm1 / (float)n;
	}
	for(int i=0; i<n; i++){
	    col[i].fakemass = fakemass;
	}
    }

#if 0
    for(int i=0; i<n; i++){
	unsigned f[4];
	getfeats(col[i].feature, f);
	printf("feat b%i e%i b%i e%i\n", f[0], f[1], f[2], f[3]);
    }
#endif

    return n;
}


static int edgecontact(obb_t *b0, obb_t *b1, bbinfo_t *bbinfo, float R[3][3], collinfo_t *col,
		       int id0, int id1, float invm0, float invm1)
{

//	BEGIN_TIMING;

    /* transform box1's edge into box0's space */
    int a0=bbinfo->edge[0], a1=bbinfo->edge[1];
    v3_t axis=bbinfo->axis;
    float depth=bbinfo->depth; 

    float e1dir[3] = {R[0][a1], R[1][a1], R[2][a1]};	/* edge direction */

    /* find mid point of the edge */
    int o0 = (a1+1)%3;
    int o1 = (a1+2)%3;

    float T[3] = {b1->center.x-b0->center.x, 
		  b1->center.y-b0->center.y, 
		  b1->center.z-b0->center.z};

    /* transform T from world space to box0's space */
    rotate2box(T, b0);

    /* axis0 X axis1 in box0's space:

       axis0 = (0..1..0)
       axis1 = R (0..1..0) = (R0i, R1i, R2i)
       axis0 X axis1 = (a0.y*a1.z - a0.z*a1.y,  =  a0[1]*R2i - a0[2]*R1i
       a0.z*a1.x - a0.x*a1.z,  =  a0[2]*R0i - a0[0]*R2i
       a0.x*a1.y - a0.y*a1.x)  =  a0[0]*R1i - a0[1]*R0i

       axis0 = (1,0,0) => axis0 X axis1 = (0,-R2i,R1i)
       axis0 = (0,1,0) => axis0 X axis1 = (R2i,0,-R0i)
       axis0 = (0,0,1) => axis0 X axis1 = (-R1i,R0i,0)

    */

    int o00 = (a0+1)%3;
    int o01 = (a0+2)%3;

/*
  float axis[3];
  axis[a0] = 0;
  axis[o00] = -R[o01][a1];
  axis[o01] = R[o00][a1];

  samedir = dot(a0xa1, T)>=0

  dot(a0xa1, T) == -R[o01][a1]*T[o00] + R[o00][a1]*T[o01]

*/
    int samedir = (-R[o01][a1]*T[o00] + R[o00][a1]*T[o01])>=0.0f;
    float s00 = signtab[R[o01][a1]<0.0f];
    float s01 = signtab[R[o00][a1]>=0.0f];
    float s10 = signtab[R[a0][o1]<0.0f];
    float s11 = signtab[R[a0][o0]>=0.0f];

    if(!samedir){s00=-s00; s01=-s01; s10=-s10; s11=-s11;}

    float mid[3];
    for(int i=0; i<3; i++){
	mid[i] = s10*b1->extent[o0]*R[i][o0] + s11*b1->extent[o1]*R[i][o1] + T[i];
    }

    /* now edge1 is in box0's space */

    float d=-s00*b0->extent[o00]*R[o00][a1]-s01*b0->extent[o01]*R[o01][a1];
    float ndotmid = mid[a0]*R[a0][a1]- dot3(e1dir,mid);
    float ndote = R[a0][a1]*R[a0][a1]-1.0f;
    if(fabsf(ndote)<epsilon) return 0;

    float t=(d-ndotmid)/ndote;

    float p[3] = {mid[0]+e1dir[0]*t, mid[1]+e1dir[1]*t, mid[2]+e1dir[2]*t};
    col->poi[1] = obb2world(b0, p);
    col->depth = -depth;//depth;
    col->normal = axis;
    if(!samedir){
	col->normal = scale(col->normal, -1);
    }
//	col->poi[0] = add(col->poi[1], scale(col->normal, depth));
    col->poi[0] = col->poi[1];
    col->tangent[0] = b0->axis[a0];
    col->tangent[1] = b1->axis[a1];
    col->id[0] = id0;
    col->id[1] = id1;
    float fakemass;
    if (invm0 == 0 && invm1 == 0) {
	fakemass = 1000000.0f;
    }
    else {
	fakemass =  invm0 > invm1 ? 1.0 / invm0 : 1.0 / invm1;
    }
    col->fakemass = fakemass;

    /* contact features */
    unsigned ft0 = (a0<<14) | ((s00<0)<<13) | ((s01<0)<<12) | id0;
    unsigned ft1 = (a1<<14) | ((s10<0)<<13) | ((s11<0)<<12) | id1;
    col->feature = (ft0>ft1) ? (ft0<<16)|ft1 : (ft1<<16)|ft0;


//	END_TIMING;

    return 1;
}

void printfeats(unsigned ft[], int n)
{
    printf("feature\n");
    for(int i=0; i<n; i++){
	printf("f%i: ", i);
	unsigned f[4];
	getfeats(ft[i], f);
	printf("b%i e%i b%i e%i\n", f[0], f[1], f[2], f[3]);
    }
}


int boxboxintersect(solid_t *s0, solid_t *s1, int id0, int id1, collinfo_t *col) {
    float R[3][3];
    bbinfo_t bbinfo;
    obb_t *b[2] = {&s0->obb, &s1->obb};
    int id[2]={id0,id1};
    float invm[2]={s0->rigid.invmass, s1->rigid.invmass};

    if(boxboxtest(b[0], b[1], R, &bbinfo) == SEPARATE) {
	return 0;
    }

    /* face vs. face, in box0's space */

    /* b0 is box, b1 is floor */
    if (bbinfo.type==TYPE_FV) {
	/* 
	   transform box1 to box0's space
		   
	   b1 -> world
	   world -> b0

	   R = 
	   / x0.x1  x0.y1  x0.z1 \
	   | y0.x1  y0.y1  y0.z1 |
	   \ z0.x1  z0.y1  z0.z1 /

	   "." means dot product  (computed by boxboxtest())
		   
	*/
	int *face = bbinfo.face;		
	v3_t poly[8];
	int A=bbinfo.minaxis>2;
	int B=!A;

	transpoly(b[A], b[B], face[B], R, A, poly);

	unsigned ft[8];
	guestfeats(id[B], face[B], ft);
	int n=clippoly(b[A], face[A], id[A], poly, ft);

//		printfeats(ft, n);

	return facecontacts(b[A], face[A], poly, ft, n, col, id[A], id[B], invm[A], invm[B]);
    }

    /* edge0 vs. edge1, edge0's space */
    else {
	return edgecontact(b[0], b[1], &bbinfo, R, col,  id[0], id[1], invm[0], invm[1]);
    }

}

