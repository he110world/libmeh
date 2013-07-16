#ifndef UTIL_H
#define UTIL_H

extern unsigned nextp2(unsigned);
extern int	nlz(unsigned);
#define lg2(x) (31-nlz(x))
#define nfact2(x) (31-nlz((x) & -(x)))
// nfact2(x) -- 800=25*32=25*2^5 => nfact2(800)=5
// 600=75*8=75*2^3 => nfact2(600)=3

#define ARRAY_INIT_CAP 16
#define ARR_CKSUM 0x4da3a98b
typedef struct {
	int	nbytes, cap;
//	unsigned	cksum;
} arr__t;
extern void *arr__new(int);
extern int arr__addnb(arr__t**,void*,int);
extern int arr__popnb(arr__t*,int);

#define arr_push(head,val) \
	((head) \
	 ? arr__addnb(&(head), NULL, sizeof(*(head))), ((head)[((arr__t*)(head)-1)->nbytes/sizeof(*(head))-1]=(val)) \
	 : ((head)=arr__new(sizeof(*(head))), *(head)=(val)))

#define arr_pushn(head,ptr,n)			\
	((head) \
	 ? arr__addnb(&(head), (ptr), (n)*sizeof(*(head)))		\
	 : ((head)=arr__new((n)*sizeof(*(head))), ((ptr) ? mymemcpy((head), (ptr), (n)*sizeof(*(head))) : 0)  ))
			       
#define arr_len(head) \
	((head) \
	 ? ((arr__t*)(head)-1)->nbytes/sizeof(*head) \
	 : 0)
/*
#define arr_pop(head) \
	((head) \
	 ? (head)[arr__popnb((head), sizeof(*(head))) / sizeof(*(head))] \
	 : NULL)
*/

#define arr_popall(head) \
	((head) \
	 ? ((arr__t*)(head)-1)->nbytes=0 \
	 : 0)

#define arr_popn(head,n) \
	((head) \
	 ? ((arr__t*)(head)-1)->nbytes-= (n)*sizeof(*head) \
	 : 0)

#define arr_pop(head) arr_popn(head,1)

#define arr_kill(head) \
	((head) \
	 ? free((arr__t*)(head)-1), head=NULL \
	 : 0)

extern int flen(FILE*);
extern int flenbyname(const char*nam);
extern char *readtextfile(const char*nam, char*mem);
#define loadtext_alloca(nam)  readtextfile((nam), alloca(flenbyname((nam))+1))
#define loadtext_malloc(nam)  readtextfile((nam), malloc(flenbyname((nam))+1))

#endif
