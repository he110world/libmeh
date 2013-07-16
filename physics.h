#ifndef PHYSICS_H
#define PHYSICS_H


// Finite Element Method
typedef struct{
	v3_t			x0;
	v3_t			x;
	v3_t			v;
	v3_t			f;
	mat3x3_t			*K;
	int			*entry; /* non-zero entries */
	int			ne;	/* num of non-zero entries */
	int			maxne;	/* similar to stl vector */
	float			mass;	/* lumped mass */
	int			fixed;

	/* temp vars for conjugate gradient */
	v3_t			Ap;
	v3_t			p;
	v3_t			z;
	v3_t			residual;
	float			Cinv[3];

} node_t;


/* tetrahedron element */
typedef struct{
	int			nodeidx[4];
	node_t			*node[4];
	float			B[6][12];
	float			P[12][6];
	float			K[12][12];
	float			R[3][3];
	float			V;		/* volume */
	int			entryidx[4][4];
	float			plastic_strain[6];

	/* material */
	float			yield;
	float			maxstrain;
	float			creep;
	float			density;
} element_t;



typedef struct{
	int			ne;			/* num of elements */
	int			nn;			/* num of nodes */
	element_t		*element;
	node_t			*node;

	float			youngs_modulus;
	float			poissons_ratio;

	float			E[6][6];

	float			damping; 		/* C = damping M --> mass damping */
} femobj_t;

extern femobj_t *fem_newobj(float vert[], int nv, int tetra[], int nt);
extern void fem_delobj(femobj_t *obj);
extern void fem_setmaterial(femobj_t *obj, 
			    float young, float poisson, float damping,
			    float yield, float maxstrain, float creep, float density) ;
extern void fem_simulate(femobj_t *obj, float dt);
extern void fem_drawobj(femobj_t *obj);



// Mass-Spring
typedef struct{
	float m;
	float invm;
	v3_t v;
	v3_t f;
	v3_t x;
	int lock;
} mass_t;


typedef struct{
	float stiff;
	float damp;
	float creep;
	float maxdeform;
	float maxforce;
} elastic_t;


typedef struct{
	mass_t *point0,*point1;
	float l;		/* length */
	float l0;		/* origional (rest) length */
	elastic_t *elastic;		/* material */
	int collide;			/* 0: not collidable -- diagonals(volume,facet), concave edges
					 * 1: planar edges
					 * 2: convex edges
					 */
} spring_t;


typedef struct{
	mass_t **point;
	int np;
	int *convexedge;
	int nce;
} facet_t;


typedef struct{
	mass_t *point;
	spring_t *spring;
	facet_t *facet;
	int np;
	int ns;
	int nf;
} obj_t;

extern void initobj();
extern void drawtetra();
extern void simtetra();
extern obj_t *makecube(v3_t c0, v3_t c1, float m, elastic_t *e);
extern void delobj(obj_t *obj);
extern void drawobj_wire(obj_t *obj);
extern void simobj(obj_t *obj);
extern void setelastic(elastic_t *e, float stiff, float damp, float creep, float maxdeform, float maxforce);
extern obj_t *voxmesh(int vox[],int nx,int ny,int nz,v3_t size,float m, elastic_t *e, int hollow);
extern void moveobj(obj_t *obj,v3_t offset);

#endif
