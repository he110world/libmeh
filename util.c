#include "util.h"

/*------------------------------------------------------
  bit utility
  -----------------------------------------------------*/
uint nextp2(uint x){
    x = x-1; 
    x = x|(x>>1); 
    x = x|(x>>2); 
    x = x|(x>>4); 
    x = x|(x>>8); 
    x = x|(x>>16); 
    return x+1; 
}

/* number of leading zeros */
#define LE 1
int nlz(uint k) {
    union {
	uint i[2]; 
	double d; 
    }u; 
    u.d = (double)k + 0.5; 
    return 1054 - (u.i[LE] >> 20); 
}


void *arr__new(int nbytes)
{
    int cap = nextp2(nbytes);
    cap = cap > 16 ? cap : 16;
    arr__t *arr = malloc(sizeof(arr__t) + cap);
    arr->nbytes = nbytes;
    arr->cap = cap;
    return arr+1;	
}


int arr__addnb(arr__t **head, void *ptr, int nbytes)
{
    assert(*head);

    arr__t *base=*head-1;
    int avail = base->cap - base->nbytes;
    int size = base->nbytes + nbytes;
    if(nbytes > avail){
	int cap = nextp2(size);
	void *p=realloc(base, cap+sizeof(arr__t)); // DEBUGGED: cap => cap+sizeof(arr__t) -- don't forget the header!
	base=p;

	//		base = realloc(base, cap+sizeof(arr__t)); // DEBUGGED: cap => cap+sizeof(arr__t) -- don't forget the header!
	*head = base+1;
	base->cap = cap;
    }
	
    if (ptr) {
	memcpy((char*)*head+base->nbytes, ptr, nbytes);
    }

    base->nbytes=size;
    return 0;
}


int arr__popnb(arr__t *head, int nbytes)
{
    arr__t *base=head-1;
    if (base->nbytes>=nbytes) {
	return base->nbytes-=nbytes;
    }

    return 0;
}





/*--------------------------------------
  Length of file
  -------------------------------------*/
int flen(FILE *f)
{
    if(!f) return 0;
    int cur=ftell(f);
    rewind(f);
    int beg=ftell(f);
    fseek(f,0,SEEK_END);
    int end=ftell(f);
    fseek(f,cur-beg,SEEK_SET);
    return end-beg;
}


int flenbyname(const char *nam)
{
    FILE *f=fopen(nam, "r");
    if(!f) return 0;
    int len=flen(f);
    fclose(f);
    return len;
}


char *readtextfile(const char *nam, char *mem)
{
    FILE *f=fopen(nam, "r");
    if(!f) return NULL;
    int size=flen(f);
    fread(mem, 1, size, f);
    fclose(f);
    mem[size]='\0';
    return mem;
}




