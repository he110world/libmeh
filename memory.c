#include "memory.h"

// internal data structures

typedef struct {
    int allocated;
    uint slot32;
    uint slot[32];
} idmngr1k_t;


typedef struct poolelem_s{
    struct poolelem_s	*next;
    int			data; // not necessarily int
} poolelem_t;

typedef struct poolblk_s{
    struct poolblk_s	*next;
    poolelem_t		elems[0]; // variable sized
} poolblk_t;

typedef struct pool_s{
    int			blk_ne; // num of elements each block has
    int			blksz; // size of the block in bytes
    int			elemsz;
    poolblk_t		*blklist;
    poolelem_t		*freelist;
    int			total, used;
    int			allocated;
} pool_t;

typedef struct bdblk_s{
    ushort			prev, next;
    int			flag;
} bdblk_t;

typedef struct{
    int blksz;
    int nblks;
    int nlevels;
    bdblk_t *tree;
    ushort *freelists;
    int allocated;
} bdsys_t;


// can manage (alloc/free) 1k indices
// bit == 1 => available

/*
  
 */
#define CHECK_IDMNGR(m)	{if ((m)<0 || (m)>=arr_len(G.id1ks) || !G.id1ks[(m)].allocated) return -1;}

void id1k_gen(int *m)
{
    if (!m) return;
    int sz=arr_len(G.id1ks);
    idmngr1k_t *idmngr=NULL;
    for (int i=0; i<sz; ++i) {
	if (!G.id1ks[i].allocated) {
	    idmngr=G.id1ks+i;
	    break;
	}
    }
    if (!idmngr) {
	arr_pushn(G.id1ks,NULL,1);
	idmngr=G.id1ks+sz;
    }

    memset(idmngr, 0xff, sizeof(idmngr1k_t));
    *m=sz;
}

int id1k_del(int *m)
{
    if (!m) return 0;
    CHECK_IDMNGR(*m);
    return G.id1ks[*m].allocated=0;
}

int id1k_alloc(int m)
{
    CHECK_IDMNGR(m);

    idmngr1k_t *idmngr=G.id1ks+m;
    int s32 = nlz(idmngr->slot32);
    if(s32==32) return -1;

    uint *bf = idmngr->slot + s32;
    int s = nlz(*bf);
    *bf ^= 0x80000000>>s;
	
    if(!*bf) idmngr->slot32 ^= 0x80000000>>s32;

    return (s32<<5) + s;
}

int id1k_free(int m, int id)
{
    CHECK_IDMNGR(m);

    idmngr1k_t *idmngr=G.id1ks+m;	

    id &= 1023;
    int s32 = id>>5;
    int s = id & 31;

    if(idmngr->slot[s32] |= 0x80000000>>s) idmngr->slot32 |= 0x80000000>>s32;
}

int id1k_allocated(int m, int id)
{
    CHECK_IDMNGR(m);

    idmngr1k_t *idmngr=G.id1ks+m;		

    id &= 1023;
    int s32 = id>>5;
    int s = id & 31;

    return !(idmngr->slot[s32] & (0x80000000>>s));
}

void id1k_killall()
{
    arr_kill(G.id1ks);
}

#undef CHECK_IDMNGR



///////////////////////////////////


/*
  Memory Pool:
  Block allocator for linked list

  Allocate large blocks;
  Allocate list nodes from large blocks => fewer malloc()/free() + less fragmental

  small constant sized blocks
*/

#define CHECK_POOL(p) {	if ((p)<0 || (p)>=arr_len(G.pools) || !G.pools[(p)].allocated) return NULL;}

void pool_gen(int *p, int blk_ne, int elemsz)
{
    if (!p) return;
    pool_t *pool=NULL;
    int sz=arr_len(G.pools);
    for (int i=0; i<sz; ++i) {
	if (!G.pools[i].allocated) {
	    pool=G.pools+i;
	    break;
	}
    }

    if (!pool) {
	arr_pushn(G.pools,NULL,1);
	pool=G.pools+sz;
    }
	
    pool->blk_ne = blk_ne;
    pool->blksz = blk_ne*(elemsz+sizeof(poolelem_t*)) + sizeof(poolblk_t*);
    pool->elemsz = elemsz;
    pool->blklist = NULL;
    pool->freelist = NULL;
    pool->total = pool->used = 0;
    pool->allocated=1;
}

int pool_del(int *p)
{
    if (!p) return NULL;
    CHECK_POOL(*p);
	
    pool_t *pool=G.pools+*p;
    poolblk_t *next;
    for(poolblk_t *b=pool->blklist; b; b=next){
	next=b->next;
	free(b);
    }

    pool->allocated=0;
    return NULL;
}

void *pool_alloc(int p)
{
    CHECK_POOL(p);

    pool_t *pool=G.pools+p;
    if(!pool->freelist){
	poolblk_t *b = malloc(pool->blksz);
	b->next = pool->blklist;
	pool->blklist = b;
	pool->total += pool->blk_ne;

	// add newly allocated elements to free list
	char *mem = b->elems;
	int size = pool->elemsz + sizeof(poolelem_t*);
	for(int i=0; i<pool->blk_ne; i++){
	    poolelem_t *e=mem;
	    e->next = pool->freelist;
	    pool->freelist = e;
	    mem+=size;
	}
    }
    pool->used ++;
    poolelem_t *e = pool->freelist;
    pool->freelist = e->next;
    return &e->data;

}

int pool_free(int p, char *dat)
{
    CHECK_POOL(p);

    pool_t *pool=G.pools+p;
    poolelem_t *e = dat-sizeof(poolelem_t*);
    e->next = pool->freelist;
    pool->freelist = e;
    pool->used--;

    return 0;
}

void pool_killall()
{
    int sz=arr_len(G.pools);
    for (int i=0; i<sz; ++i) {
	pool_del(&i);
    }
}

#undef CHECK_POOL



/*--------------------------------------
  Buddy System 
  -------------------------------------*/

enum{ NULL_BLK=0xffff, MAX_MEMSIZE=16*1024*1024};


/*
  blk's (usage) flag: 
  0 => not on free list && not allocated
  1 => on list
  2 => allocated
*/

enum{UNINIT, ALLOCATED, ONLIST};

#define CHECK_BD(b) {if ((b)<0 || (b)>=arr_len(G.bds) || !G.bds[(b)].allocated) return -1;}

void bd_gen(int *id, int blksz, int nblks)
{
    int sz=arr_len(G.bds);
    bdsys_t *b=NULL;
    for (int i=0; i<sz; ++i) {
	if (!G.bds[i].allocated) {
	    b=G.bds+i;
	    break;
	}
    }
    if (!b) {
	arr_pushn(G.bds,NULL,1);
	b=G.bds+sz;
    }
	
    blksz = nextp2(blksz);
    nblks = nextp2(nblks);
    if(blksz*nblks > MAX_MEMSIZE) return NULL;

    b->blksz = blksz;
    b->nblks = nblks;
    b->nlevels = lg2(nblks)+1;
    b->tree = malloc((2*nblks-1) * sizeof(bdblk_t));

    b->freelists = malloc(b->nlevels*sizeof(ushort));
    memset(b->freelists, NULL_BLK, b->nlevels*sizeof(ushort));
    b->freelists[b->nlevels-1] = 0;
    b->tree->prev = b->tree->next = NULL_BLK;
    b->tree->flag = ONLIST;
    b->allocated=1;
	
}


int bd_del(int *id)
{
    if (!id) return -1;

    CHECK_BD(*id);

    bdsys_t *b=G.bds+*id;
    free(b->tree);
    free(b->freelists);
    b->allocated=0;

    return 0;
}

static inline void bd_freelist_pop(bdsys_t *b, int i, int flag)
{
    ushort node = b->freelists[i]; // pop a unused node from freelist
    ushort next = b->tree[node].next;
    b->tree[node].flag = flag;
    if(next != NULL_BLK){ b->tree[next].prev = NULL_BLK; }
    b->freelists[i] = next;
}

static inline void bd_freelist_push(bdsys_t *b, ushort node, int i) // push a node to freelist
{
    b->tree[node].prev = NULL_BLK;
    b->tree[node].next = b->freelists[i];
    b->tree[node].flag = ONLIST;
    if(b->freelists[i] != NULL_BLK){ b->tree[b->freelists[i]].prev = node; }
    b->freelists[i] = node;
}

static inline void bd_freelist_del(bdsys_t *b, ushort node, int i)
{
    ushort prev = b->tree[node].prev;
    ushort next = b->tree[node].next;
    if(prev != NULL_BLK) {b->tree[prev].next = next;}
    else {b->freelists[i] = next;}
    if(next != NULL_BLK) {b->tree[next].prev = prev;}
    b->tree[node].flag = UNINIT;
}


static inline int bd_getaddr(bdsys_t *b, ushort node, ushort level)
{
    /*
      nblks = 2^(b->nlevels-1 - level) = b->nblks / (2^level)

      addr = b->nblks * b->blksz * (node+1-nblks) / nblks
      = b->nblks * b->blksz * ( (node+1)/nblks - 1 )
      = b->nblks * b->blksz * (node+1) * (2^level) / b->nblks - b->nblks * b->blksz
      = b->blksz * (node+1) * (2^level) - b->blksz * b->nblks
      = b->blksz * ( (node+1)*(2^level) - b->nblks )

    */

    return b->blksz * (((node+1)<<level) - b->nblks);
}

static inline int bd_getnode(bdsys_t *b, int blk, int level)
{
    /*
      node = 2^(b->nlevels-1-level) - 1 + offset/(2^level)
      = 2^(b->nlevels-1) / (2^level) + offset / (2^level) - 1
      = (2^(b->nlevels-1) + offset) / (2^level) -1
      = (b->nblks + offset) / (2^level) - 1
    */

    return (b->nblks + blk) / (1<<level) - 1;
}

/*
  1. find the best fit level
  2.1 remove from free list
  or 
  2.2 split larger block, add blocks to free list



  15:  0
  14:  1 2
  13:  3 4 5 6
  12:  7 8 9 10 11 12 13 14
  11:  15...           30
  10:  31 32...        62
  9:   63 64...        126
  8:   127 128...      254
  7:   255 256...      510
  6:   511 512...      1022
  5:   1023 1024...    2046
  4:   2047 2048...    4094
  3:   4095 4096...    8190
  2:   8191 8192 ...   16382
  1:   16383 16384 ... 32766
  0:   32767 32768 ... 65534

*/

int bd_alloc(int id, int size, int *finalsize)
{
    CHECK_BD(id);

    bdsys_t *b=G.bds+id;
    /* 1. */
    *finalsize = size = size<=b->blksz ? b->blksz : nextp2(size);
    int level = lg2(size/b->blksz);

    /* 2.1
       about the for loop below:
       break out from (i<b->nlevels) when: i==b->nlevels
       break out from (b->freelists[i] == NULL_BLK) when: (i<b->nlevels) && (b->freelists[i] != NULL_BLK)
    */
    int i;
    for(i=level; (i<b->nlevels) && (b->freelists[i] == NULL_BLK); ++i);
    if(i==b->nlevels) return -1;

    /* 2.2 */
    ushort node = b->freelists[i];
    bd_freelist_pop(b,i, i==level);
    if(i>level){ // available block > wanted block, so split
	for(i=i-1, node=(node<<1)+1; i>=level; --i, node=(node<<1)+1){
	    int flag = i==level;
	    bd_freelist_push(b,node+1,i);
	    b->tree[node].flag = flag; /* flag == UNINIT or ALLOCATED ==> 0/1 */
	    if(flag) break;
	}
    }
    return bd_getaddr(b,node,i);
}


#define BUDDY(node) (((node)<<1)-((node)^1))  /* buddy(2k+1)=2k+2; buddy(2k)=2k-1 */
/*
  1. find the highest level L that the block can be at
  2. find the block in level <= L (block.flag == 1)
  3. put the block on the free list / merge the block with its buddy
  {
  if buddy is on the list:
  3.1 remove the buddy from the list
  3.2 goto the higher level
  else
  3.3 put the block on the list
  }
*/
int bd_free(int id, int addr)
{
    CHECK_BD(id);

    bdsys_t *b=G.bds+id;
    int blk = addr / b->blksz;

    /* 1. */
    int level = nfact2(blk);
    if(level<0){ level = b->nlevels-1; }

    /* 2. */
    int i,node;
    for(i=level, node=bd_getnode(b,blk,i);  
	(i>=0) && (b->tree[node].flag != ALLOCATED); 
	--i, node=(node<<1)+1);

    if(i<0) return -1;

    /* 3. */
    b->tree[node].flag = UNINIT;
    for(int buddy=BUDDY(node); (buddy>0) && (b->tree[buddy].flag == ONLIST); ++i){
	bd_freelist_del(b,buddy,i);
	node=(node-1)>>1;
	buddy=BUDDY(node);
    }
    bd_freelist_push(b,node,i);

    return 0;
}

#undef BUDDY


void bd_killall()
{
    int sz=arr_len(G.bds);
    for (int i=0; i<sz; ++i) {
	if (G.bds[i].allocated) {
	    bd_del(&i);
	}
    }
}


/*-------------------------------------------------
  Test whether the addr is already reused by others. If not, recover it.

  TODO: implement this some days later!
  ------------------------------------------------*/
int bd_reusable(const bdsys_t *b, int addr)
{
    return 0;
}



// Pair manager

#define PM_HASH 0

#if PM_HASH
#define INVALID_ID 0xffffffff

static unsigned hash(unsigned a, unsigned b)
{
    unsigned key = (b<<16) | a;
    key = ~key + (key << 15);
    key = key ^ (key >> 12);
    key = key + (key << 2);
    key = key ^ (key >> 4);
    key = key * 2057;
    key = key ^ (key >> 16);
    return key;
}

/* new pairmngr_t */
pairmngr_t *allocpm(int cap)
{
    pairmngr_t *pm = malloc(sizeof(pairmngr_t));
    cap = nextp2(cap);
    pm->cap = cap;
    pm->mask = cap-1;
    pm->pair = malloc(cap*sizeof(pair_t));
    pm->next = malloc(cap*sizeof(unsigned));
    pm->hashtab = malloc(cap*sizeof(unsigned));
    pm->npairs = 0;
    memset(pm->hashtab, INVALID_ID, cap*sizeof(unsigned));
    memset(pm->next, INVALID_ID, cap*sizeof(unsigned));
    return pm;
}

void freepm(pairmngr_t *pm)
{
    free(pm->pair);
    free(pm->next);
    free(pm->hashtab);
    free(pm);
}

pair_t *findpr(pairmngr_t *pm, unsigned short a, unsigned short b)
{
    if(a>b){
	unsigned short tmp=a;
	a=b;
	b=tmp;
    }
    unsigned h = hash(a,b) & pm->mask;
    unsigned i=pm->hashtab[h];
    while(i!=INVALID_ID && !(a==pm->pair[i].id[0] && b==pm->pair[i].id[1])){
	i=pm->next[i];
    }
    if(i==INVALID_ID) return NULL;
    return pm->pair+i;
}


/* suppose a<=b */
static pair_t *findpr1(pairmngr_t *pm, 
		       unsigned short a, 
		       unsigned short b, 
		       unsigned h)
{
    unsigned i=pm->hashtab[h];
    while(i!=INVALID_ID && !(a==pm->pair[i].id[0] && b==pm->pair[i].id[1])){
	i=pm->next[i];
    }
    if(i==INVALID_ID) return NULL;
    return pm->pair+i;
}

/* add pair */
pair_t *addpr(pairmngr_t *pm, unsigned short a, unsigned short b)
{
    if(a>b){
	unsigned short tmp=a;
	a=b;
	b=tmp;
    }
    unsigned h = hash(a,b) & pm->mask;
    pair_t *p=findpr1(pm,a,b,h);
    if(p) return p;

    if(pm->npairs >= pm->cap){
	int c=pm->cap<<1;
	pm->mask = c-1;
	pm->pair = realloc(pm->pair, c*sizeof(pair_t));
	free(pm->next);
	free(pm->hashtab);
	pm->next = malloc(c*sizeof(unsigned));
	pm->hashtab = malloc(c*sizeof(unsigned));
	for(int i=0; i<pm->npairs; i++){
	    pair_t *p = pm->pair+i;
	    unsigned h = hash(p->id[0], p->id[1]) & pm->mask; // mask changes, so we have to rehash
	    pm->next[i] = pm->hashtab[h];
	    pm->hashtab[h] = i;
	}
	memset(pm->next+pm->cap, INVALID_ID, pm->cap*sizeof(unsigned));
	memset(pm->hashtab+pm->cap, INVALID_ID, pm->cap*sizeof(unsigned));
	pm->cap<<=1;
    }

    p = pm->pair + pm->npairs;
    p->id[0]=a;
    p->id[1]=b;
    pm->next[pm->npairs] = pm->hashtab[h];
    pm->hashtab[h] = pm->npairs++;
    return p;
}

/* delete pair */
void delpr(pairmngr_t *pm,unsigned short a, unsigned short b)
{
    if(a>b){
	unsigned short tmp=a;
	a=b;
	b=tmp;
    }

    /* find the pair */
    unsigned h = hash(a,b) & pm->mask;
    pair_t *p=findpr1(pm,a,b,h);
    if(!p) return;

    unsigned d = p-pm->pair; /* idx to the found pair */
    unsigned i = pm->hashtab[h]; /* first pair that hashes to h */

    /* find the previous pair */
    unsigned prev=INVALID_ID;
    while(i!=d){
	prev=i;
	i=pm->next[i];
    }
    if(prev!=INVALID_ID)
	pm->next[prev]=pm->next[i];
    else 
	pm->hashtab[h]=pm->next[i];

    /* move tail to the removed pair */
    pair_t *tail = pm->pair + pm->npairs - 1;
    if(tail == p){
	pm->npairs--;
	return;
    }

    h = hash(tail->id[0],tail->id[1]) & pm->mask;
    i = pm->hashtab[h];
    unsigned t = tail-pm->pair;
    prev = INVALID_ID;
    while(i!=t){
	prev=i;
	i=pm->next[i];
    }
    if(prev!=INVALID_ID)
	pm->next[prev]=pm->next[t];
    else
	pm->hashtab[h]=pm->next[t];
    pm->pair[d]=pm->pair[t];
    pm->next[d]=pm->hashtab[h];
    pm->hashtab[h]=d;
    pm->npairs--;
}

#else
/* (Hopefully) improved pair manager -- random access, no hashing, and you can get each object's "partners" immediately. */

/*
  num partners   [32k]
  idx to next[]  [32k]
  next           [???k]
  free list      [??k]    

  next[i]:   next|twin
 */

typedef struct{
	uint			prev, next;
	ushort			obj;
} link_t;

typedef struct{
	ushort			nmates[MAX_NUM_OBJ];
	uint			first[MAX_NUM_OBJ];
	link_t			*link;
	int			link_cap, nl;

	uint			*free; // 2i being in the freelist means that 2i+1 is in the list too.
	uint			free_cap, nf;
	uint			new; // links whose new<= id <cap are newly allocated (unused)
} pm_t;


pm_t *newpm()
{
	pm_t *pm=calloc(1, sizeof(pm_t));
	return pm;
}

void delpm(pm_t *pm)
{
	free(pm->link);
	free(pm->free);
}

void addpr(pm_t *pm, ushort A, ushort B)
{
	// do nothing if pair(A,B) exists
	if(pm->nmates[A] > pm->nmates[B]){
		ushort tmp=A;
		A=B;
		B=tmp;
	}
	for(uint i=pm->first[A]; i!=INVALID_LINK; i=pm->link[i].next){
		if(pm->link[i].obj == B) return; // the pair already exists
	}

	
	// now this is a new pair.
	// but first we need to make sure there's enough space to store it.
	uint a,b;
	if(pm->nf){
		a=pm->free[--pm->nf];
		b=a ^ 1;
	}
	else{	 // free list is empty => links are used up => new links
		if(pm->new >= pm->link_cap){
			pm->new = pm->link_cap;
			pm->link_cap = nextp2(pm->link_cap);
			pm->link = realloc(pm->link, pm->link_cap);
		}
		a=pm->link + pm->new;
		b=a+1;
		pm->new++;
	}

	link_t *lnka = pm->link + a;
	link_t *lnkb = pm->link + b;
	lnka->prev=INVALID_LINK; 
	lnka->next=pm->first[A];  
	if(pm->first[A] != INVALID_LINK) pm->link[pm->first[A]].prev = a;
	lnka->obj=A;  pm->first[A]=a;  pm->npartners[A]++;

	lnkb->prev=INVALID_LINK; 
	lnkb->next=pm->first[B];  
	if(pm->first[B] != INVALID_LINK) pm->link[pm->first[B]].prev = b;
	lnkb->obj=B;  pm->first[B]=b;  pm->npartners[B]++;
	
}


void delpr(pm_t *pm, ushort A, ushort B)
{
	if(pm->nmates[A] > pm->nmates[B]){
		ushort tmp=A;
		A=B;
		B=tmp;
	}
	
	for(uint i=pm->first[A]; i!=INVALID_LINK; i=pm->link[i].next){
		if(pm->link[i].obj == B) break;
	}

	if(i==INVALID_LINK) return; // didn't find the pair

	uint prevA=pm->link[i].prev;
	uint nextA=pm->link[i].next;
	uint prevB=pm->link[i^1].prev;
	uint nextB=pm->link[i^1].next;

	if(prevA != INVALID_LINK) link[prevA].next = nextA;
	else pm->first[A] = nextA;

	if(prevB != INVALID_LINK) link[prevB].next = nextB;
	else pm->first[B] = nextB;

	if(nextA != INVALID_LINK) link[nextA].prev = prevA;
	if(nextB != INVALID_LINK) link[nextB].prev = prevB;
	
	--pm->nmates[A];
	--pm->nmates[B];

	if(pm->nf >= pm->free_cap){
		pm->free_cap=nextp2(pm->free_cap);
		pm->free = realloc(pm->free, pm->free_cap);
	}
	pm->free[pm->nf++] = i;
}

// unlink all mates of the object from it
void unlinkmates(pm_t *pm, ushort obj)
{
	
}

int getmates(pm_t *pm, ushort obj, ushort mates[])
{
	if(obj>=32*1024) return 0;
	if(pm->nmates[obj]==0) return 0;

	ushort *m=mates;
	for(uint i=pm->first[obj]; i!=INVALID_LINK; i=pm->link[i].next){
		*m++ = pm->link[i].obj;
	}

	return pm->nmates[obj];
}


#endif // PM_HASH

#if 0

int main(int argc, char *argv[])
{
    pairmngr_t *pm = allocpm(128);

    printf("0x%x\n", addpr(pm, 10, 100));
    printf("0x%x\n", 	addpr(pm, 3, 9));
    printf("0x%x\n", 	addpr(pm, 666, 50));
    printf("0x%x\n", 	addpr(pm, 666, 50));
    printf("0x%x\n", 	addpr(pm, 666, 50));
    pair_t *p = findpr(pm, 10, 100);
    if(p){
	printf("pair = (%u %u)\n", p->id[0], p->id[1]);
    }

    p = findpr(pm, 666, 50);
    if(p){
	printf("pair = (%u %u)\n", p->id[0], p->id[1]);
    }

//	delpr(pm, 666, 50);
    delpr(pm, 10, 100);
    delpr(pm, 3, 9);

    p = findpr(pm, 666, 50);
    if(p){
	printf("found pair = (%u %u)\n", p->id[0], p->id[1]);
    }


    p = findpr(pm, 10, 100);
    if(p){
	printf("found pair = (%u %u)\n", p->id[0], p->id[1]);
    }


    p = findpr(pm, 3, 9);
    if(p){
	printf("found pair = (%u %u)\n", p->id[0], p->id[1]);
    }

    return 0;
}
#endif
