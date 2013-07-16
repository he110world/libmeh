#ifndef MEMORY_H
#define MEMORY_H

// Memory management routines

extern void id1k_gen(int *m);
extern int id1k_del(int *m);
extern int id1k_alloc(int m);
extern int id1k_free(int m, int id);
extern int id1k_allocated(int m, int id);
extern void id1k_killall();

extern void pool_gen(int *p, int blk_ne, int elemsz);
extern int pool_del(int *p);
extern void *pool_alloc(int p);
extern int pool_free(int p, char *dat);
extern void pool_killall();

extern void bd_gen(int *id, int blksz, int nblks);
extern int bd_del(int *id);
extern int bd_alloc(int id, int sz, int *finalsz);
extern int bd_free(int id, int addr);
extern void bd_killall();

// pair manager
typedef struct{
	ushort			id[2];
} pair_t;

typedef struct{
	pair_t			*pair;
	unsigned			*next;
	unsigned			*hashtab;
	unsigned			cap;
	unsigned			npairs;
	unsigned			mask;
} pairmngr_t;


extern pairmngr_t *allocpm(int cap);
extern void freepm(pairmngr_t *pm);
extern pair_t *findpr(pairmngr_t *pm, ushort a, ushort b);
extern pair_t *addpr(pairmngr_t *pm, ushort a, ushort b);
extern void delpr(pairmngr_t *pm,ushort a, ushort b);


#endif
