#ifndef MEH_H__
#define MEH_H__

#include "math2.h"
#include "sp.h"
#include "cd.h"
#include "rigid.h"
#include "memutil.h"
#include "render.h"

#define BRUSH_TYPE 1
#define OBJ_TYPE 2
#define NULL_BRUSH_TYPE 0xff

typedef struct brush_s{
	unsigned short type;
	struct brush_s *prev,*next;
	phd_t *phd;
	phdr_t *phdr;
	unsigned short spid;
	short wallid;
} brush_t;


#define OBJ_CONVEX 1
#define OBJ_PHANTOM 2		/* phantom objects are those objects that intersect with other brushes and objects */

typedef struct obj_s{
	unsigned short type;
	struct obj_s *prev, *next;
	phd_t *phd;
	unsigned short spid;
	rigid_t rb;
	unsigned short convex;
	unsigned int flag;
	phd_t *localphd;
} obj_t;


/* portals are rectanglular */
typedef struct{
	
} ptl_t;


typedef struct{
	struct brush_s *brushlist;
	struct metaptl_s *mptllist;
	int phdr_totnverts;
	int phdr_totnvids;	/* totol num of indices to verts in phdr on this wall */
} wall_t;


#define MAX_NWALLS 32
typedef struct sctr_s{
	struct sctr_s		*prev,*next;
	sp_t			*sp;
	phd_t			*phd;
	wall_t			walls[MAX_NWALLS];
	short			nwalls;
	unsigned short		spid;
	brush_t			*brushlist;
	avl_t			*planetree;  /* store plane equations in an AVL tree */
	v3_t			rightvec[2]; /* used only if the plane normal = {0, +/-1, 0}  (at most two of them in a convex phd)*/
	
	vbo_t			*vbo;
	ibo_t			*ibo;
	int			phdr_totnverts;
	int			phdr_totnvids; /* total num of indices to verts in phdr */
} sctr_t;

typedef struct metaptl_s{
	struct metaptl_s *prev,*next,*twin;
	v3_t *verts;		/* verts[32] */
	sctr_t *sctr;
	int nverts;
	int ccw;
	int wallid;
} metaptl_t;



typedef struct{
} pnt_t;



typedef union{
	struct {
		unsigned short type;
		void *prev, *next;
		phd_t *phd;
		unsigned short spid;
	} common;
	brush_t brush;
	obj_t obj;
	unsigned short type;
} stuff_t;


extern int phd2sctr(phd_t *phd);
extern void rmsctr(sctr_t *s);
extern int xlatsctr_def(sctr_t *sctr, v3_t x, v3_t *xdone); /* deffered updating */
extern void xlatsctr_upd(sctr_t *s, v3_t x);
extern int xlatsctr(sctr_t *sctr, v3_t dx);
extern void initmap();
extern void shutdownmap();
extern sctr_t *inwhichsctr(v3_t pnt);
extern void render();
extern sctr_t *debug_getsctr(int id);

/* engine variables. you can access them directly or use some variable managers like quake's cvar. */
extern float fov;


#endif
