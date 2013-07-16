#ifndef GEOMETRY_H
#define GEOMETRY_H

/*
  The engine is data driven

 */


// Arculated bodies
uint ab_load();
void ab_render();

uint mesh_load();

// Brush
typedef struct{
	int			v[2];
	int			p[2];
	plane_t			pln;
} edge_t;

typedef struct{
	plane_t			plane;
	int			nv;
	short			*v;
	short			*e;
} poly_t;


typedef struct{
	unsigned int		flag;
	int			npolys, nverts, nedges;
	v3_t			*verts; /* no malloc for verts/edges/polys; just mem=malloc(), and assign it to V/E/P*/
	edge_t			*edges;
	poly_t			*polys;

	void			*mem; /* |<--verts-->|<--edges-->|<--polys-->| */
	int			cap;
	float			bnd[2][3];
} phd_t;			/* polyhedron, not Ph.D. :) */

typedef struct{
	v3_t org;
	v3_t axis;
	float deg;
	float rad;
} rot_t;


extern phd_t *allocphd();
extern void freephd(phd_t *phd);
extern void bnd2planes(float bnd[2][3], plane_t planes[]);
extern int splitphd(plane_t *split, phd_t *phd);
extern int splitphd2(phd_t *phd, plane_t *split, phd_t **front, phd_t **back);
extern void addplns(phd_t *phd, plane_t planes[], int np);
extern int phdsub(phd_t *p, phd_t *p2, phd_t *result[]);
extern void xlatphd(phd_t *p, v3_t x);
extern int phdqry(const phd_t *p, const phd_t *p2);
extern void drawphd(phd_t *phd);
extern int splitpoly(v3_t vert[], int nverts, plane_t *split);
extern void splitpoly2(plane_t *split, v3_t vert[], v3_t front[], v3_t back[], int nverts, int *nfronts, int *nbacks);
extern int pntinphd(v3_t pnt, phd_t *phd);
extern int phdabovepln(phd_t *phd, plane_t *pln);
extern phd_t *poly2phd(phd_t *phd, int pid, float thichness);

#define PHD_INTX_PLN 0
#define PHD_ABOVE_PLN 1
#define PHD_BELOW_PLN 2

#define SAME_PLN(p,p2) (dot((p)->n, (p2)->n)>0.99 && fabs((p)->d - (p2)->d)<0.01)



#endif
