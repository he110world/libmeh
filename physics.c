#include "physics.h"


// Finite Element Method
/*

  -- The stiffness matrix K is a (sparse) matrix of size N_n*N_n, where N_n is the number of nodes in the mesh.

  -- Each node stores one row of K, so adding/removing nodes (adding/removing rows&columns of K) would be easy.

  -- Each row is a linked list (because each row has different number of non-zero entries, and during the initialization

  of K we have to add these entries to each row; and ease the updating of the mesh connectivity )

  -- Everytime a node fractures, a new node is generated (hence adds 3 dofs to the system). So the dof of the system
  (or number of nodes) should be able to change. => linked list
  But using linked list avoids random access to the nodes => each element (tetrahedron) has to stores four pointers to
  the nodes, not four indices.
*/

#include <stdlib.h>
#include "math2.h"
#include "fem2.h"
#include <GL/gl.h>

void print_v3(v3_t v);

static void compute_K_entries(femobj_t *obj)
{
    for(element_t *e=obj->element; e<&(obj->element[obj->ne]); e++){

	for(int i=0; i<4; i++){
	    node_t *n = e->node[i];

	    for(int j=0; j<4; j++){
		int found=0;
		for(int k=0; k<n->ne; k++){
		    if(n->entry[k] == e->nodeidx[j]){
			e->entryidx[i][j] = k;
			found = 1;
			break;
		    }
		}

		if(!found){
		    if(n->ne >= n->maxne){
			if(n->ne == 0){
			    n->entry = malloc(sizeof(int));
			    n->maxne = 1;
			}
			else{
			    n->entry = realloc(n->entry, 2*n->maxne*sizeof(int));
			    n->maxne *= 2;
			}
		    }

		    n->entry[n->ne] = e->nodeidx[j];
		    e->entryidx[i][j] = n->ne++;

		}
	    }
	}
    }

    for(int i=0; i<obj->nn; i++){
	obj->node[i].K = calloc(obj->node[i].ne, sizeof(mat3x3_t));
    }
}


static void compute_E(femobj_t *obj)
{
    memset(obj->E, 0, 36*sizeof(float));

    float e = obj->youngs_modulus;
    float v = obj->poissons_ratio;
    float factor = e/((1+v)*(1-2*v));

    obj->E[0][0] = obj->E[1][1] = obj->E[2][2] = factor*(1-v);
    obj->E[3][3] = obj->E[4][4] = obj->E[5][5] = factor*(.5-v);
    obj->E[0][1] = 	obj->E[1][0] = 	obj->E[0][2] = 	obj->E[2][0] =	obj->E[1][2] = 	obj->E[2][1] = 	factor*v;

    printf("E = \n");
    for(int i=0; i<6; i++){
	for(int j=0; j<6; j++){
	    printf("%f\t ",obj->E[i][j]);
	}

	printf("\n");
    }


}



static void compute_Be(element_t *e)
{
    float a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4;
    float x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
    float z21,z31,z41,z32,z42,z43;
    float y21,y31,y41,y32,y42,y43;
    float x21,x31,x41,x32,x42,x43;
    float six_times_V;

    x1=e->node[0]->x.x;	x2=e->node[1]->x.x;	x3=e->node[2]->x.x;	x4=e->node[3]->x.x;
    y1=e->node[0]->x.y;	y2=e->node[1]->x.y;	y3=e->node[2]->x.y;	y4=e->node[3]->x.y;
    z1=e->node[0]->x.z;	z2=e->node[1]->x.z;	z3=e->node[2]->x.z;	z4=e->node[3]->x.z;

    x21=x2-x1;	x31=x3-x1;	x41=x4-x1;	x32=x3-x2;	x42=x4-x2;	x43=x4-x3;
    y21=y2-y1;	y31=y3-y1;	y41=y4-y1;	y32=y3-y2;	y42=y4-y2;	y43=y4-y3;
    z21=z2-z1;	z31=z3-z1;	z41=z4-z1;	z32=z3-z2;	z42=z4-z2;	z43=z4-z3;

    a1 = y2*z43 - y3*z42 + y4*z32;
    a2 = -y1*z43 + y3*z41 - y4*z31;
    a3 = y1*z42 - y2*z41 + y4*z21;
    a4 = -y1*z32 + y2*z31 - y3*z21;

    b1 = -x2*z43 + x3*z42 - x4*z32;
    b2 = x1*z43 - x3*z41 + x4*z31;
    b3 = -x1*z42 + x2*z41 - x4*z21;
    b4 = x1*z32 - x2*z31 + x3*z21;

    c1 = x2*y43 - x3*y42 + x4*y32;
    c2 = -x1*y43 + x3*y41 - x4*y31;
    c3 = x1*y42 - x2*y41 + x4*y21;
    c4 = -x1*y32 + x2*y31 - x3*y21;

    six_times_V = x21*(y31*z41-y41*z31) + y21*(z31*x41-z41*x31) + z21*(x31*y41-x41*y31);

    /* ********************************************* */
    memset(e->B, 0, 6*12*sizeof(float));

    e->B[0][0] = a1;	e->B[0][3] = a2;	e->B[0][6] = a3;	e->B[0][9] = a4;
    e->B[1][1] = b1;	e->B[1][4] = b2;	e->B[1][7] = b3;	e->B[1][10] = b4;
    e->B[2][2] = c1;	e->B[2][5] = c2;	e->B[2][8] = c3;	e->B[2][11] = c4;

    e->B[3][0] = b1;	e->B[3][1] = a1;
    e->B[3][3] = b2;	e->B[3][4] = a2;
    e->B[3][6] = b3;	e->B[3][7] = a3;
    e->B[3][9] = b4;	e->B[3][10] = a4;

    e->B[4][1] = c1;	e->B[4][2] = b1;
    e->B[4][4] = c2;	e->B[4][5] = b2;
    e->B[4][7] = c3;	e->B[4][8] = b3;
    e->B[4][10] = c4;	e->B[4][11] = b4;

    e->B[5][2] = a1;	e->B[5][3] = c2;
    e->B[5][5] = a2;	e->B[5][6] = c3;
    e->B[5][8] = a3;	e->B[5][9] = c4;
    e->B[5][11] = a4;	e->B[5][0] = c1;


    /* *********************************************** */
    for(int i=0; i<6; i++)
	for(int j=0; j<12; j++)
	    e->B[i][j] /= six_times_V;

    e->V = six_times_V/6.0;
}


/* Pe = Ve Be^T E */
void compute_Pe(element_t *e, float E[6][6])
{
    memset(e->P, 0, sizeof(float)*12*6);

    for(int i=0; i<12; i++){
	for(int j=0; j<6; j++){
	    for(int k=0; k<6; k++){
		e->P[i][j] += e->B[k][i]*E[k][j];
	    }
	}
    }

    for(int i=0; i<12; i++){
	for(int j=0; j<6; j++){
	    e->P[i][j] *= e->V;
	}
    }
}

/* Ke = Pe Be */
static void compute_Ke(element_t *e)
{
    memset(e->K, 0, sizeof(float)*12*12);

    for(int i=0; i<12; i++){
	for(int j=0; j<12; j++){
	    for(int k=0; k<6; k++){
		e->K[i][j] += e->P[i][k]*e->B[k][j];
	    }
	}
    }
}


static void compute_M(femobj_t *obj)
{
    int i=0;
    for(element_t *e=obj->element; e<&(obj->element[obj->ne]); e++){
	float node_mass = e->density * e->V / 4; /* element_mass / 4 */

	for(int i=0; i<4; i++){
	    e->node[i]->mass += node_mass;
	}
    }
}




static void compute_Re(element_t *e)
{
    v3_t e01, e02, e03;
    v3_t e1, e2, e3;
    v3_t n01, n02, n03;

    e01 = sub(e->node[1]->x0, e->node[0]->x0);
    e02 = sub(e->node[2]->x0, e->node[0]->x0);
    e03 = sub(e->node[3]->x0, e->node[0]->x0);

    e1 = sub(e->node[1]->x, e->node[0]->x);
    e2 = sub(e->node[2]->x, e->node[0]->x);
    e3 = sub(e->node[3]->x, e->node[0]->x);

    n01 = scale(cross(e02, e03), 1.0/6.0/e->V);
    n02 = scale(cross(e03, e01), 1.0/6.0/e->V);
    n03 = scale(cross(e01, e02), 1.0/6.0/e->V);

    e->R[0][0] = e1.x*n01.x + e2.x*n02.x + e3.x*n03.x;
    e->R[0][1] = e1.x*n01.y + e2.x*n02.y + e3.x*n03.y;
    e->R[0][2] = e1.x*n01.z + e2.x*n02.z + e3.x*n03.z;

    e->R[1][0] = e1.y*n01.x + e2.y*n02.x + e3.y*n03.x;
    e->R[1][1] = e1.y*n01.y + e2.y*n02.y + e3.y*n03.y;
    e->R[1][2] = e1.y*n01.z + e2.y*n02.z + e3.y*n03.z;

    e->R[2][0] = e1.z*n01.x + e2.z*n02.x + e3.z*n03.x;
    e->R[2][1] = e1.z*n01.y + e2.z*n02.y + e3.z*n03.y;
    e->R[2][2] = e1.z*n01.z + e2.z*n02.z + e3.z*n03.z;

    orthonormalization3x3(e->R);
}


static void update_plastic_strain(element_t *e, float dt)
{
    float dof[12];
    float total_strain[6];
    float elastic_strain[6];
    float elastic_norm;
    float plastic_norm;

    /* Re^-1 x - x0 */
    for(int i=0; i<4; i++){
	for(int j=0; j<3; j++){
	    dof[3*i + j] = e->R[0][j] * e->node[i]->x.x 
		+ e->R[1][j] * e->node[i]->x.y
		+ e->R[2][j] * e->node[i]->x.z

		+ ((float *)&(e->node[i]->x0))[j];
	}
    }

    /* Be (Re^-1 x - x0) */
    for(int i=0; i<6; i++){
	total_strain[i] = 0;

	for(int j=0; j<12; j++){
	    total_strain[i] += e->B[i][j]*dof[j];
	}
    }

    /* elastic_strain = total_strain - plastic_strain */
    elastic_norm = plastic_norm = 0;
    for(int i=0; i<6; i++){
	elastic_strain[i] = total_strain[i] - e->plastic_strain[i];

	elastic_norm += elastic_strain[i] * elastic_strain[i];
	plastic_norm += e->plastic_strain[i] * e->plastic_strain[i];
    }

    elastic_norm = sqrt(elastic_norm);
    plastic_norm = sqrt(plastic_norm);


    if(elastic_norm > e->yield){
	for(int i=0; i<6; i++){
	    e->plastic_strain[i] += dt * e->creep * elastic_strain[i];
	}
    }

    if(plastic_norm > e->maxstrain){
	for(int i=0; i<6; i++){
	    e->plastic_strain[i] *= e->maxstrain / plastic_norm;
	}
    }

}


static void clear_K(femobj_t *obj)
{
    for(int i=0; i<obj->nn; i++){
	memset(obj->node[i].K, 0, sizeof(mat3x3_t)*obj->node[i].ne);
    }
}

/* call clear_K() first */
static void assemble_K(femobj_t *obj)
{
    for(element_t *e=obj->element; e<&(obj->element[obj->ne]); e++){
	for(int i=0; i<4; i++){
	    node_t *n = e->node[i];

	    for(int j=0; j<4; j++){
		float KR[3][3];

#if 1
		/* K = sum{ Re Ke Re^T }*/

		/* Ke Re^T */

		/* ************************************************ */
		/* begin VERY SLOW */
		for(int m=0; m<3; m++){
		    for(int n=0; n<3; n++){
			KR[m][n] = e->K[3*i+m][3*j] * e->R[n][0]
			    + e->K[3*i+m][3*j+1] * e->R[n][1]
			    + e->K[3*i+m][3*j+2] * e->R[n][2];
		    }
		}
		/* end VERY SLOW */
		/* **************************************************** */
#endif
		/* Re (Ke Re^T) */
		int entry = e->entryidx[i][j];
		mat3x3_t *Knode = &e->node[i]->K[entry];

#if 1
		for(int m=0; m<3; m++){
		    for(int n=0; n<3; n++){
			Knode->data[m][n] += e->R[m][0]*KR[0][n]
			    + e->R[m][1]*KR[1][n]
			    + e->R[m][2]*KR[2][n];
		    }
		}
#endif


#if 0
		for(int m=0; m<3; m++){
		    for(int n=0; n<3; n++){
			Knode->data[m][n] += e->K[3*i+m][3*j+n];
		    }
		}
#endif

	    }
	}
    }
}


static void clear_forces(femobj_t *obj)
{
    v3_t zero={0,0,0};
	
    for(int i=0; i<obj->nn; i++){
	obj->node[i].f = zero;
    }
}


/* call clear_forces() first */
static void assemble_forces(femobj_t *obj)
{
    for(element_t *e=obj->element; e<&(obj->element[obj->ne]); e++){
	for(int i=0; i<4; i++){
	    v3_t f={0,0,0};
	    int row=3*i;

	    /* factor Re out => save 3 mat vec productions */
	    /* f0 = -Re Ke x0 */


	    /* ***************************************************** */
	    /* begin VERY SLOW */
	    for(int j=0; j<4; j++){
		int col=3*j;
		v3_t *x0 = &(e->node[j]->x0);


		/* Ke X0 */
		f.x -= e->K[row][col]*x0->x + e->K[row][col+1]*x0->y + e->K[row][col+2]*x0->z;
		f.y -= e->K[row+1][col]*x0->x + e->K[row+1][col+1]*x0->y + e->K[row+1][col+2]*x0->z;
		f.z -= e->K[row+2][col]*x0->x + e->K[row+2][col+1]*x0->y + e->K[row+2][col+2]*x0->z;
	    }
	    /* end VERY SLOW */
	    /* ******************************************************* */

#if 1
	    /* fplastic = Re Pe elastic_strain */
	    for(int j=0; j<6; j++){
		f.x += e->P[row][j] * e->plastic_strain[j];
		f.y += e->P[row+1][j] * e->plastic_strain[j];
		f.z += e->P[row+2][j] * e->plastic_strain[j];
	    }
#endif
	    /* f = f0 + fplastic */
	    e->node[i]->f = add(e->node[i]->f, mat3x3_v3(e->R, f)); 

#if 0
	    e->node[i]->f = add(e->node[i]->f, f);  
#endif
	}
    }
}


static void assemble_motion(femobj_t *obj, float dt)
{
    for(int i=0; i<obj->nn; i++){
	node_t *node = &obj->node[i];

	/* K x */
	v3_t Kx={0,0,0};
	for(int j=0; j<node->ne; j++){
	    Kx = add( mat3x3_v3( node->K[j].data, obj->node[node->entry[j]].x ), Kx );
	}

	/* K' = M + dt C + dt^2 K */
	float dt2 = dt*dt;
	for(int j=0; j<node->ne; j++){
	    int col = node->entry[j];
	    mat3x3_t *K = &node->K[j];

	    for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
		    K->data[i][j] *= dt2;
		}
	    }

	    if(i == col){
		float mass_damp = node->mass 
		    + obj->damping * node->mass * dt;

		K->data[0][0] += mass_damp;
		K->data[1][1] += mass_damp;
		K->data[2][2] += mass_damp;
	    }
	}


	/* ***************************** */
/* 		print_v3(node->f); */


	/* f' = Mv -dt(Kx + f0 + fplastic -fext) */
	node->f = sub( scale(node->v, node->mass), scale(add(Kx, node->f), dt) ); 
	/* 		node->f = scale(node->v, node->mass); *//* correct if no K involved */
    }

}

/* rmin : the smallest ||r||^2 */
/* preconditioned conjugate gradient solver */
static void pcg(femobj_t *obj, float dt, int num_iterates, float rmin)
{
    /* solve K' v = f' */

    /*
      v = 0

      residual = f'

      C^-1 = inv( diag( K ) )

      z = C^-1 r

      p = z

      z^T r
    */

    float zr = 0;
    for(int i=0; i<obj->nn; i++){
	node_t *node = &obj->node[i];

#if 1
	/* ********************TESTING*********************** */
	if(node->fixed)
	    continue;
	/* ***************END TESTING************************ */
#endif

	node->v = vec(0,0,0);
	node->residual = node->f;

	for(int j=0; j<node->ne; j++){
	    if(i == node->entry[j]){
		node->Cinv[0] = 1.0 / node->K[j].data[0][0];
		node->Cinv[1] = 1.0 / node->K[j].data[1][1];
		node->Cinv[2] = 1.0 / node->K[j].data[2][2];
	    }
	}


	node->z.x = node->Cinv[0]*node->residual.x;
	node->z.y = node->Cinv[1]*node->residual.y;
	node->z.z = node->Cinv[2]*node->residual.z;

	node->p = node->z;

	zr += dot( node->z, node->residual ); 
    }

    int hello;

    for(int iter=0; iter<num_iterates; iter++){
	float pAp = 0;

	/*
	  z^T r
	  alpha = --------
	  p^t A p
	*/

	for(int row=0; row<obj->nn; row++){
	    node_t *node = &obj->node[row];

#if 1
	    /* *****************TESTING******************* */
	    if(node->fixed){
		continue;
	    }
	    /* ****************END TESTING**************** */
#endif

	    /* A p */
	    node->Ap = vec(0,0,0);
	    for(int i=0; i<node->ne; i++){
		int col = node->entry[i];
		v3_t *p = &obj->node[col].p;
		mat3x3_t *K = &node->K[i];

		node->Ap.x += K->data[0][0] * p->x
		    + K->data[0][1] * p->y
		    + K->data[0][2] * p->z;

		node->Ap.y += K->data[1][0] * p->x
		    + K->data[1][1] * p->y
		    + K->data[1][2] * p->z;

		node->Ap.z += K->data[2][0] * p->x
		    + K->data[2][1] * p->y
		    + K->data[2][2] * p->z;

	    }

	    pAp += dot(node->p, node->Ap);

	}
/* 		printf("pAp=%f\n", pAp); */

	float alpha = zr/pAp;
	float zr_old = zr;

	zr = 0;
	for(node_t *node=obj->node; node<&obj->node[obj->nn]; node++){

#if 1
	    /* ********************TESTING************************* */
	    if(node->fixed)
		continue;
	    /* ******************END TESTING*********************** */
#endif

	    /*  v = v + alpha p	*/
	    node->v = add( node->v, scale( node->p, alpha ));

	    /* residual = residual - alpha A p */
	    node->residual = sub( node->residual, scale( node->Ap, alpha ) );

	    /* z = C^-1 residual */
	    node->z.x = node->Cinv[0]*node->residual.x;
	    node->z.y = node->Cinv[1]*node->residual.y;
	    node->z.z = node->Cinv[2]*node->residual.z;

	    zr += dot( node->z, node->residual );
	}

	/* if r is small enough, quit */
	if(zr < rmin){
	    /* 			printf("zr%i = %f\n", iter, zr); */
	    return;
	}

	/*
	  z^T r
	  beta = ----------------
	  {z^T r}_old
	*/
	float beta = zr / zr_old;

	/* p = z + beta p */
	for(node_t *node=obj->node; node<&obj->node[obj->nn]; node++){


#if 1
	    /* ************TESTING************** */
	    /* fixed node */
	    if(node->fixed)
		continue;
	    /* **********END TESTING************ */
#endif

	    node->p = add( node->z, scale( node->p, beta ) );

	}

	hello = iter;
    }

/*  	printf("zr%i = %f\n", hello, zr); */

}



/* rmin : the smallest ||r||^2 */
/* (NON-preconditioned) conjugate gradient solver */
static void cg(femobj_t *obj, float dt, int num_iterates, float rmin)
{
    /* solve K' v = f' */

    /*
      v = 0

      residual = f'

      p = r

      r^T r
    */

    float rr = 0;
    for(int i=0; i<obj->nn; i++){
	node_t *node = &obj->node[i];

	node->v = vec(0,0,0);
	node->residual = node->f;

	node->p = node->residual;

	rr += dot( node->residual, node->residual ); 
    }

    int hello;

    for(int iter=0; iter<num_iterates; iter++){
	float pAp = 0;

	/*
	  r^T r
	  alpha = --------
	  p^t A p
	*/

	for(int row=0; row<obj->nn; row++){
	    node_t *node = &obj->node[row];

	    /* A p */
	    node->Ap = vec(0,0,0);
	    for(int i=0; i<node->ne; i++){
		int col = node->entry[i];
		v3_t *p = &obj->node[col].p;
		mat3x3_t *K = &node->K[i];

		node->Ap.x += K->data[0][0] * p->x
		    + K->data[0][1] * p->y
		    + K->data[0][2] * p->z;

		node->Ap.y += K->data[1][0] * p->x
		    + K->data[1][1] * p->y
		    + K->data[1][2] * p->z;

		node->Ap.z += K->data[2][0] * p->x
		    + K->data[2][1] * p->y
		    + K->data[2][2] * p->z;

	    }

	    pAp += dot(node->p, node->Ap);

	}
/* 		printf("pAp=%f\n", pAp); */

	float alpha = rr/pAp;
	float rr_old = rr;

	rr = 0;
	for(node_t *node=obj->node; node<&obj->node[obj->nn]; node++){
	    /*  v = v + alpha p	*/
	    node->v = add( node->v, scale( node->p, alpha ));

	    /* residual = residual - alpha A p */
	    node->residual = sub( node->residual, scale( node->Ap, alpha ) );

	    rr += dot( node->residual, node->residual );
	}

	/* if r is small enough, quit */
	if(rr < rmin) {
/* 			printf("rr%i = %f\n", iter, rr); */
	    break;
	}
	/* 		printf("r = %f\n", sqrt(r_norm));  */
	/*
	  r^T r
	  beta = ----------------
	  {r^T r}_old
	*/
	float beta = rr / rr_old;

	/* p = z + beta p */
	for(node_t *node=obj->node; node<&obj->node[obj->nn]; node++){
	    node->p = add( node->residual, scale( node->p, beta ) );
	}

	hello = iter;
    }

/* 	printf("zr%i = %f\n", hello, rr); */

}


static void update_pos(femobj_t *obj, float dt)
{
    for(node_t *node=obj->node; node<&obj->node[obj->nn]; node++){

#if 1
	/* *******************TESTING******************* */
	if(node->fixed)
	    continue;

	/* ****************END TESTING****************** */
#endif

	node->x = add( node->x, scale( node->v, dt ) );
    }
}


static void compute_stress(element_t *e)
{

}


void fracture()
{
}

void post_process()
{
}

void apply_gravity(femobj_t *obj)
{
#if 1
    for(int i=0; i<obj->nn; i++){
	obj->node[i].f = add(obj->node[i].f,  vec(0,obj->node[i].mass*98,0));
    }
#endif

}


void fem_setvelocity(femobj_t *obj, v3_t v)
{
    for(node_t *n=obj->node; n<&obj->node[obj->nn]; n++){
	n->v = v;
    }
}



femobj_t *fem_newobj(float vert[], int nv, int tetra[], int nt)
{
    femobj_t *obj;

    obj = malloc(sizeof(femobj_t));

    obj->node = calloc(nv, sizeof(node_t));
    obj->element = calloc(nt, sizeof(element_t));

    /* nodes */
    {			
	v3_t *v = (v3_t *)vert;
	for(int i=0; i<nv; i++, v++){
	    obj->node[i].x0 = obj->node[i].x = *v;
	}
    }

    /* elements */
    {
	int *t=tetra;
	for(int i=0; i<nt; i++, t+=4){
	    element_t *e = obj->element + i;

/* 			memcpy(e->nodeidx, t, 4*sizeof(int)); */
	    e->nodeidx[0] = t[1];
	    e->nodeidx[1] = t[0];
	    e->nodeidx[2] = t[2];
	    e->nodeidx[3] = t[3];


	    for(int j=0; j<4; j++){
		e->node[j] = obj->node + e->nodeidx[j];
	    }
	}
    }
		
    obj->nn=nv;
    obj->ne=nt;

    return obj;
}


void fem_init(femobj_t *obj)
{
    /* properties */
    compute_K_entries(obj);

    compute_E(obj);

    for(int i=0; i<obj->ne; i++){
	element_t *e = &obj->element[i];

	compute_Be(e);
	compute_Pe(e, obj->E);
	compute_Ke(e);
    }



    compute_M(obj);


}


void fem_delobj(femobj_t *obj)
{
    free(obj->element);

    for(int i=0; i<obj->nn; i++){
	free(obj->node[i].K);
	free(obj->node[i].entry);
    }
    free(obj->node);
    free(obj);

}


void fem_setmaterial(femobj_t *obj, 
		     float young, float poisson, float damping,
		     float yield, float maxstrain, float creep, float density) 
{
    obj->youngs_modulus = young;
    obj->poissons_ratio = poisson;
    obj->damping = damping;

    for(int i=0; i<obj->ne; i++){
	element_t *e = &obj->element[i];

	e->yield = yield;
	e->maxstrain = maxstrain;
	e->creep = creep;
	e->density = density;
    }
}


void fem_lock(femobj_t *obj, int idx)
{
    obj->node[idx].fixed = 1;
    obj->node[idx].mass = 1e100;
}

void fem_unlock(femobj_t *obj, int idx)
{
    obj->node[idx].fixed = 0;
}


void print_K(femobj_t *obj);

void fem_simulate(femobj_t *obj, float dt)
{
    for(element_t *e=obj->element; e<&obj->element[obj->ne]; e++){
	compute_Re(e);

#if 0
	/* no stiffness warping */
	memset(e->R, 0, sizeof(float)*9);
	e->R[0][0] = e->R[1][1] = e->R[2][2] = 1;
#endif

	update_plastic_strain(e, dt);   
    }


    clear_K(obj);
    assemble_K(obj);

    clear_forces(obj);
    assemble_forces(obj);

    apply_gravity(obj);  

    assemble_motion(obj, dt);

    pcg(obj, dt, 30, 1e-10*obj->youngs_modulus);
    update_pos(obj, dt);  


}

void fem_drawobj(femobj_t *obj)
{

    glColor3f(0,0,1);

    glBegin(GL_LINES);
    for(int i=0; i<obj->ne; i++){
	element_t *e = &obj->element[i];

	for(int i=0; i<3; i++){
	    for(int j=i+1; j<4; j++){
		glVertex3fv(&e->node[i]->x);
		glVertex3fv(&e->node[j]->x);
	    }
	}
    }
    glEnd();

    /* locked nodes */

    glColor3f(0,1,0);
    glPointSize(8);

    glBegin(GL_POINTS);
    for(int i=0; i<obj->nn; i++){
	node_t *n = &obj->node[i];

	if(n->fixed){
	    glVertex3fv(&obj->node[i].x);
	}
    }
    glEnd();
    glPointSize(1);

}


femobj_t *test, *test2;


void print_K_entries(femobj_t *obj)
{
    for(int i=0; i<obj->nn; i++){
	node_t *n = &obj->node[i];

	printf("row %i    (", i);

	for(int j=0; j<n->ne; j++){
	    printf("%i, ", n->entry[j]);
	}
	printf(")\n\n");

    }

}


void print_mat3x3(float m[3][3])
{
    for(int i=0; i<3; i++){
	printf("| ");
	for(int j=0; j<3; j++){
	    printf("%f ", m[i][j]);
	}
	printf(" |\n");
    }
    printf("\n");
}

void print_K(femobj_t *obj)
{
    for(int i=0; i<obj->nn; i++){
	node_t *n = &obj->node[i];
		
	printf("row %i\n", i);
	for(int j=0; j<n->ne; j++){
	    print_mat3x3( n->K[j].data );
	}
	printf("--------------------------------\n");
    }
}

void print_v3(v3_t v)
{
    printf("(%f, %f, %f)\n", v.x, v.y, v.z);
}


#define SOFT 30000.0
#define HARD 60000.0
#define VERY_HARD 120000.0

/*
  Stiffness - zr_threshold relationship

  30000 -- 0.000003
  60000 -- 0.000006
  120000 -- 0.000012
  240000 -- 0.000024
  ...

  stiffness -- stiffness * 10^-10

  (stiffness = young's modulus)
*/
 

void fem_translate(femobj_t *obj, v3_t v)
{
    for(int i=0; i<obj->nn; i++){
	obj->node[i].x0 = add(obj->node[i].x0, v);
	obj->node[i].x = add(obj->node[i].x, v);
    }
}	 

void testfem()
{
    extern float dog_vert[];
    extern int dog_tetra[];
    extern int num_vert;
    extern int num_tetra;

    test = fem_newobj(dog_vert, num_vert,
		      dog_tetra, num_tetra);
    fem_setmaterial(test, 10000, 0.33, 2, 
		    1, 0, 0, 2);
    fem_init(test);
    /* 	fem_setvelocity(test, vec(0,-1,0)); */
#if 1
    fem_lock(test, 10); 
    fem_lock(test, 20);
    fem_lock(test, 30);
#endif

    test2 = fem_newobj(dog_vert, num_vert,
		       dog_tetra, num_tetra);
    fem_setmaterial(test2, 10000, 0.33, 2,
		    1, 0, 0, 2);
    fem_init(test2);
    fem_translate(test2, vec(0,-1.5,0));

#if 1
    fem_lock(test2, 120);
    fem_lock(test2, 130); 
    fem_lock(test2, 110);
/* 	fem_lock(test2, 50); */
#endif

    fem_simulate(test, 0.01);

    fem_simulate(test2, 0.01);
}





// Mass-Spring


mass_t tetra[4];
spring_t spring[6];


void setmasspoint(mass_t *point, float m, v3_t x)
{
    point->m=m;
    point->x=x;
    point->v=point->f=vec(0,0,0);
}



void setelastic(elastic_t *e, float stiff, float damp, float creep, float maxdeform, float maxforce)
{
    e->stiff=stiff;
    e->damp=damp;
    e->creep=creep;
    e->maxdeform=maxdeform;
    e->maxforce=maxforce;
}

void setspring(spring_t *s, mass_t *point0, mass_t *point1, elastic_t *e)
{
    s->point0=point0;
    s->point1=point1;

    s->l0=length(sub(point0->x,point1->x));

    s->elastic=e;
}


/* semi-implicit euler */
void integrate(mass_t point[], int n, float dt)
{
    int i;
    mass_t *p;

    for(i=0;i<n;i++){
	p=point+i;

	p->v=add(p->v,scale(p->f,dt/p->m)); /* v(t+dt) = v(t) + a(t)*dt */
	p->x=add(p->x,scale(p->v,dt));	   /* x(t+dt) = x(t) + v(t+dt)*dt */
    }

}


/* spring force */
void springf(spring_t *s)
{
    v3_t xrel, vrel, f;
    mass_t *p0, *p1;
    float deform, v;		/* dl: change in length; dv: change in velocity */
    elastic_t *e;
    float len;

    p0=s->point0;
    p1=s->point1;

    xrel=sub(p0->x,p1->x);
    vrel=sub(p0->v,p1->v);

    len=length(xrel);
    s->l=len;

    v=dot(vrel,xrel)/len; /* project vrel to the direction of xrel ( xrel/|xrel| ) */
    deform=(len - s->l0)/len; /* 1 - L0/L */

    e=s->elastic;
    f=scale(xrel, e->stiff*deform + e->damp*v);


#if 1
    if(fabs(deform) > e->maxdeform) s->l0=(1.0-e->creep)*s->l0 + e->creep*len;
#endif

#if 0
    if(length(f)>e->maxforce) s->l0=(1.0-e->creep)*s->l0 + e->creep*len;
#endif


    p0->f=sub(p0->f,f);
    p1->f=add(p1->f,f);
}


/* external force */
void extf(mass_t point[], int n)
{
    int i;

    for(i=0;i<n;i++){
	point[i].f=vec(0,-9.8*point[i].m,0);
    }
}

/* simple collision detection and response (with ground) */
int simplecol(mass_t point[], int n)
{
    int i;
    mass_t *p;
    v3_t v;
    float lv;
    int nc;
    float pd;

    p=point;
    nc=0;
    for(i=0;i<n;i++){
	if(p->x.y<0.2)
	    p->lock=1;
	else
	    p->lock=0;

	/* testing! */
	if(p->x.y<-0.0001){
	    v=normalize(p->v);
	    v=scale(v,-p->x.y/v.y);
	    p->x=add(p->x,v); 
	    p->v.y=-p->v.y;
	    p->v=scale(p->v,.3);
	    nc++;
	}
	p++;
    }

    return nc;
}

void rigidify(spring_t spring[], int ns)
{
    int i;
    float diff;
    v3_t d;
    float len;
    mass_t *p0,*p1;

    for(i=0;i<ns;i++){
	p0=spring[i].point0;
	p1=spring[i].point1;
		
	if(p0->lock && p1->lock)
	    continue;

	d=sub(p1->x, p0->x);
	len=length(d);
	diff=(spring[i].l0-len)/len;

	if(!p0->lock && !p1->lock){
	    p0->x = sub(p0->x,scale(d,.5*diff));
	    p1->x = add(p1->x,scale(d,.5*diff));
	}
	else if(p1->lock){
	    p0->x = sub(p0->x,scale(d,diff));
	}
	else {
	    p1->x = add(p1->x,scale(d,diff));
	}
    }
}


void simulate(mass_t point[], int np, spring_t spring[], int ns, float dt)
{
    int i;
    static int nc=1;

    extf(point,np);

#if 1
    for(i=0;i<ns;i++)
	springf(spring+i);
#endif

    integrate(point,np,dt);

#if 0
	
    for(i=0;i<5;i++)
	rigidify(spring,ns);

#endif

    nc=simplecol(point,np);
		
}


void printv(v3_t v)
{
    printf("(%f %f %f)\n",v.x,v.y,v.z);
}


extern elastic_t elastic;
void initobj()
{
    int i;

    setmasspoint(tetra,.5,vec(5,5+20,2));
    setmasspoint(tetra+1,.5,vec(3,5+20,0));
    setmasspoint(tetra+2,.5,vec(7,3+20,0));
    setmasspoint(tetra+3,.5,vec(5,8+20,1));

    for(i=0;i<4;i++)
	printv(tetra[i].x);

    setelastic(&elastic,500,1,0.01,0.02,10);

    setspring(spring,tetra,tetra+1,&elastic);
    setspring(spring+1,tetra,tetra+2,&elastic);
    setspring(spring+2,tetra,tetra+3,&elastic);
    setspring(spring+3,tetra+1,tetra+2,&elastic);
    setspring(spring+4,tetra+1,tetra+3,&elastic);
    setspring(spring+5,tetra+2,tetra+3,&elastic);



}



obj_t *makecube(v3_t c0, v3_t c1, float m, elastic_t *e)
{
    obj_t *obj;
    float v[2][3];
    int i,j,k;
    mass_t *p;
    spring_t *s;

    /* if the cube is too small, then enlarge it */
    if(fabs(c0.x-c1.x)<0.1) c0.x=c1.x+1.0;
    if(fabs(c0.y-c1.y)<0.1) c0.y=c1.y+1.0;
    if(fabs(c0.z-c1.z)<0.1) c0.z=c1.z+1.0;


    v[0][0]=c0.x<c1.x ? c0.x : c1.x;
    v[1][0]=c0.x<c1.x ? c1.x : c0.x;
    v[0][1]=c0.y<c1.y ? c0.y : c1.y;
    v[1][1]=c0.y<c1.y ? c1.y : c0.y;
    v[0][2]=c0.z<c1.z ? c0.z : c1.z;
    v[1][2]=c0.z<c1.z ? c1.z : c0.z;


    obj=malloc(sizeof(obj_t));

    obj->point=malloc(sizeof(mass_t)*8);
    obj->np=8;

    obj->spring=malloc(sizeof(spring_t)*28);
    obj->ns=28;

	
    /* mass points */
    for(i=0;i<2;i++)
	for(j=0;j<2;j++)
	    for(k=0;k<2;k++){
		p=&(obj->point[4*i+2*j+k]);

		p->x.x=v[i>0][0];
		p->x.y=v[j>0][1];
		p->x.z=v[k>0][2];

		p->m=m/8.0;

		p->v=p->f=vec(0,0,0);
	    }


    /* springs */
    k=0;
    for(i=0;i<7;i++)
	for(j=i+1;j<8;j++){
	    s=&(obj->spring[k]);

	    s->point0=obj->point+i;
	    s->point1=obj->point+j;

	    s->l0=length(sub(s->point0->x,s->point1->x));

	    s->elastic=e;

	    k++;
	}


    return obj;
}


void delobj(obj_t *obj)
{
    int i;

    free(obj->point);
    free(obj->spring);

    for(i=0;i<obj->nf;i++){
	free(obj->facet[i].point);
	free(obj->facet[i].convexedge);
    }
    free(obj->facet);

    free(obj);

}


int 
veceq(v3_t a, v3_t b)
{
    return a.x==b.x && a.y==b.y && a.z==b.z;
}


void 
drawobj_wire(obj_t *obj)
{
    int i,j;
    v3_t zero;

#if 0
    glBegin(GL_LINES);
    glColor3f(1,0,0);

    for(i=0;i<obj->ns;i++){
	glVertex3fv(&(obj->spring[i].point0->x));
	glVertex3fv(&(obj->spring[i].point1->x));

    }
    glEnd();
#endif

#if 0
    glBegin(GL_LINES);
    glColor3f(0,1,0);

    for(i=0;i<obj->ns;i++){
	if(obj->spring[i].collide==2){
	    glVertex3fv(&(obj->spring[i].point0->x));
	    glVertex3fv(&(obj->spring[i].point1->x));

	}
    }
    glEnd();

#endif	/* 0 */

#if 1
    glBegin(GL_LINES);
    glColor3f(0,0,1);
    for(i=0;i<obj->nf;i++){
	for(j=0;j<obj->facet[i].np;j++){
	    glVertex3fv(&(obj->facet[i].point[j]->x));
	    glVertex3fv(&(obj->facet[i].point[(j+1)%obj->facet[i].np]->x));
	}
    }
    glEnd();
#endif

#if 0
    glBegin(GL_LINES);
    glColor3f(1,1,0);
    for(i=0;i<obj->nf;i++){
	for(j=0;j<obj->facet[i].nce;j++){
	    glVertex3fv(&(obj->facet[i].point[obj->facet[i].convexedge[j]]->x));
	    glVertex3fv(&(obj->facet[i].point[(obj->facet[i].convexedge[j]+1)%obj->facet[i].np]->x));
	}
    }

    glEnd();
#endif


}


void drawtetra()
{
    int i;
    v3_t *x;

    glBegin(GL_LINES);
    for(i=0;i<6;i++){
#if 0
	printv(spring[i].point0->x);
	printv(spring[i].point1->x);
#endif

	glVertex3fv((float *)&(spring[i].point0->x));
	glVertex3fv((float *)&(spring[i].point1->x));
    }
    glEnd();

}


void simtetra()
{
    int i;

    usleep(1000*16);
    for(i=0;i<5;i++)
	simulate(tetra,4,spring,6,0.02);
}


void simobj(obj_t *obj)
{
    int i;

    usleep(1000*16);
    for(i=0;i<5;i++)
	simulate(obj->point,obj->np,obj->spring,obj->ns,0.01/3);

}


void
moveobj(obj_t *obj, v3_t offset)
{
    int i;

    for(i=0;i<obj->np;i++){
	obj->point[i].x = add(obj->point[i].x,offset);
    }
}


/* 26 neighbours */
int *hollowvox(int vox[],int nx,int ny,int nz, int *n)
{
    int x,y,z;
    int osx,osy,osz;
    int i,j,k;
    int id;
    int ni;			/* n internal voxels */
    int *newvox;

    newvox=malloc(nx*ny*nz*sizeof(int));
    memcpy(newvox,vox,nx*ny*nz*sizeof(int));

    /* the outer voxels (x==0 || y==0 || z==0) can't be internal voxels. so x,y,z starts from 1 and ends n-1 */

    ni=0;
    osx=ny*nz;
    osy=nz;
    for(x=1;x<nx-1;x++)
	for(y=1;y<ny-1;y++)
	    for(z=1;z<nz-1;z++){

		id=x*osx + y*osy + z;

		if(vox[id]){
		    for(i=-1;i<=1;i++)
			for(j=-1;j<=1;j++)
			    for(k=-1;k<=1;k++){
				if(!vox[id+ i*osx + j*osy + k])
				    goto non_internal;
			    }
					
		    newvox[id]=0;
		    ni++;
		}

	      non_internal:
		continue;
	    }

    if(!ni){
	free(newvox);
	newvox=NULL;
    }

    *n=ni;

    return newvox;
}


int collidetab[]={0,2,2,1,2,1,0,0,2,0,1,0,1,0,0,0};

/* size == size of a voxel */
obj_t *voxmesh(int vox0[],int nx,int ny,int nz,v3_t size,float m, elastic_t *elastic, int hollow)
{
    int nvox,nv,ne;		/* nvox: num of non-empty vox */
    int *v;			/* vertex */
    int *e[3];	        /* edge */
    int *f[3];		/* facet */
    int i,j;
    int x,y,z;
    obj_t *obj;
    int a,b,c,d;
    mass_t *p;
    int v0,v1;
    spring_t *s;
    int ns;
    int fv0,fv1,fv2,fv3;		/* facet vertex index */
    int ev0,ev1;			/* edge vertex index */
    int osx,osy,osz;		/* index offsets */
    int posx,posy,posz;
    int eosx,eosy,eosz;
    int vos[3];
    int cv[8];		/* cube vertex index */
    int *vox;
    int ni;
    int bv;			/* bit flags indicating whether a bounding voxel (there're four of them) of an edge is solid */
    int bvi[4];
    int nf;
    facet_t *facet;
    int ce[4];
    int nc;
    int ix,iy,iz;

    if(hollow){
	vox=hollowvox(vox0,nx,ny,nz,&ni);
	if(!vox){
	    vox=vox0;
	    hollow=0;
	}
    }
    else{
	vox=vox0;
    }


    obj=malloc(sizeof(obj_t));

    v=malloc((nx+1)*(ny+1)*(nz+1)*sizeof(int));
    memset(v,-1,(nx+1)*(ny+1)*(nz+1)*sizeof(int)); /* set all vertex index to -1 (set to 0 => big bang! -- 1st vox's id==0!)*/

    e[0]=calloc(nx*(ny+1)*(nz+1), sizeof(int)); /* edgex */
    e[1]=calloc((nx+1)*ny*(nz+1), sizeof(int));
    e[2]=calloc((nx+1)*(ny+1)*nz, sizeof(int));

    f[0]=calloc((nx+1)*ny*nz, sizeof(int));
    f[1]=calloc(nx*(ny+1)*nz, sizeof(int));
    f[2]=calloc(nx*ny*(nz+1), sizeof(int));

    /* lock x,y => lock x */
    nvox=0;
    for(x=0;x<nx;x++)
	for(y=0;y<ny;y++)
	    for(z=0;z<nz;z++){

		i = x*ny*nz + y*nz + z;

		if(vox[i]){

		    /* mass points */
		    for(a=0;a<2;a++)
			for(b=0;b<2;b++)
			    for(c=0;c<2;c++){
				j=(x+a)*(ny+1)*(nz+1) + (y+b)*(nz+1) + (z+c);
				v[j]=1;
			    }


		    for(a=0;a<3;a++)
			for(b=0;b<2;b++){
			    /* edges */
			    for(c=0;c<2;c++){
				if(a==0){
				    j=x*(ny+1)*(nz+1) + (y+b)*(nz+1) + (z+c);
				}
				else if(a==1){
				    j=(x+b)*ny*(nz+1) + y*(nz+1) + (z+c);
				}
				else{
				    j=(x+b)*(ny+1)*nz + (y+c)*nz + z;
				}
				e[a][j]=1;
			    }


			    /* facets */
			    if(a==0)
				j=(x+b)*ny*nz + y*nz + z;
			    else if(a==1)
				j=x*(ny+1)*nz + (y+b)*nz + z;
			    else 
				j=x*ny*(nz+1) + y*(nz+1) + z+b;

			    f[a][j]=1;
			}

		    nvox++;
		}
	    }


    /* mass points */
    nv=0;
    for(i=0; i<(nx+1)*(ny+1)*(nz+1); i++){
	if(v[i]>0){
	    v[i]=nv++;
	}
    }

    obj->point=malloc(sizeof(mass_t)*nv);
    obj->np=nv;

    /* positions of mass points */
    for(x=0;x<nx+1;x++)
	for(y=0;y<ny+1;y++)
	    for(z=0;z<nz+1;z++){
		i=x*(ny+1)*(nz+1) + y*(nz+1) + z;

		if(v[i]>=0){
		    p=&(obj->point[v[i]]);
		    p->x.x = x*size.x;
		    p->x.y = y*size.y;
		    p->x.z = z*size.z;
					
		    p->f=p->v=vec(0,0,0);
		    p->m=m;

		}


	    }

    /* edge */
    ne=0;
    for(i=0; i<nx*(ny+1)*(nz+1); i++){if(e[0][i]) ne++;}
    for(i=0; i<(nx+1)*ny*(nz+1); i++){if(e[1][i]) ne++;}
    for(i=0; i<(nx+1)*(ny+1)*nz; i++){if(e[2][i]) ne++;}

    /* facet diagonals */
    for(i=0; i<(nx+1)*ny*nz; i++){if(f[0][i]) ne+=2;}
    for(i=0; i<nx*(ny+1)*nz; i++){if(f[1][i]) ne+=2;}
    for(i=0; i<nx*ny*(nz+1); i++){if(f[2][i]) ne+=2;}

    /* ===springs=== */
    obj->ns=ne+4*nvox;
    obj->spring=malloc(sizeof(mass_t)*obj->ns); /* 4*nvox == num of diagonal springs */

    /* edge springs */
    ns=0;
    nc=0;
    osx=ny*nz;
    osy=nz;
    for(x=0;x<nx;x++)
	for(y=0;y<ny+1;y++)
	    for(z=0;z<nz+1;z++){
		i=x*(ny+1)*(nz+1) + y*(nz+1) + z;

		if(e[0][i]){
		    ev0=i;
		    ev1=ev0+(ny+1)*(nz+1);

		    s=&(obj->spring[ns++]);
		    s->point0=&(obj->point[v[ev0]]);
		    s->point1=&(obj->point[v[ev1]]);
		    s->l0=length(sub(obj->point[v[ev0]].x,obj->point[v[ev1]].x));


		    /* find collidable edge */
		    bv=15;
		    if(y<=0) bv&=12;
		    else if(y>=ny) bv&=3;
		    if(z<=0) bv&=10;
		    else if(z>=nz) bv&=5;

		    if((bv&1) && (!vox0[x*osx + (y-1)*osy + (z-1)])) bv&=~1;
		    if((bv&2) && (!vox0[x*osx + (y-1)*osy + z])) bv&=~2;
		    if((bv&4) && (!vox0[x*osx + y*osy + (z-1)])) bv&=~4;
		    if((bv&8) && (!vox0[x*osx + y*osy + z])) bv&=~8;

		    s->collide=collidetab[bv];

		    if(s->collide==2){
			e[0][i]=2;
			nc++;
		    }

		}
	    }

    for(x=0;x<nx+1;x++)
	for(y=0;y<ny;y++)
	    for(z=0;z<nz+1;z++){
		i=x*ny*(nz+1) + y*(nz+1) + z; 

		if(e[1][i]){
		    ev0=x*(ny+1)*(nz+1) + y*(nz+1) + z;
		    ev1=ev0+nz+1;

		    s=&(obj->spring[ns++]);
		    s->point0=&(obj->point[v[ev0]]);
		    s->point1=&(obj->point[v[ev1]]);
		    s->l0=length(sub(obj->point[v[ev0]].x,obj->point[v[ev1]].x));

		    /* find collidable edge */
		    bv=15;
		    if(z<=0) bv&=12;
		    else if(z>=nz) bv&=3;
		    if(x<=0) bv&=10;
		    else if(x>=nx) bv&=5;

		    if((bv&1) && (!vox0[(x-1)*osx + y*osy + (z-1)])) bv&=~1;
		    if((bv&2) && (!vox0[x*osx + y*osy + (z-1)])) bv&=~2;
		    if((bv&4) && (!vox0[(x-1)*osx + y*osy + z])) bv&=~4;
		    if((bv&8) && (!vox0[x*osx + y*osy + z])) bv&=~8;

		    s->collide=collidetab[bv];					

		    if(s->collide==2){
			e[1][i]=2;
			nc++;
		    }

		}
	    }

    for(x=0;x<nx+1;x++)
	for(y=0;y<ny+1;y++)
	    for(z=0;z<nz;z++){
		i=x*(ny+1)*nz + y*nz + z;

		if(e[2][i]){
		    ev0=x*(ny+1)*(nz+1) + y*(nz+1) + z;
		    ev1=ev0+1;

		    s=&(obj->spring[ns++]);
		    s->point0=&(obj->point[v[ev0]]);
		    s->point1=&(obj->point[v[ev1]]);
		    s->l0=length(sub(obj->point[v[ev0]].x,obj->point[v[ev1]].x));

		    /* find collidable edge */
		    bv=15;
		    if(x<=0) bv&=12;
		    else if(x>=nx) bv&=3;
		    if(y<=0) bv&=10;
		    else if(y>=ny) bv&=5;

		    if((bv&1) && (!vox0[(x-1)*osx + (y-1)*osy + z])) bv&=~1;
		    if((bv&2) && (!vox0[(x-1)*osx + y*osy + z])) bv&=~2;
		    if((bv&4) && (!vox0[x*osx + (y-1)*osy + z])) bv&=~4;
		    if((bv&8) && (!vox0[x*osx + y*osy + z])) bv&=~8;

		    s->collide=collidetab[bv];

		    if(s->collide==2){
			e[2][i]=2;
			nc++;
		    }
		}
	    }
    printf("ns=%i\n", ns);
    printf("nc=%i\n", nc);


    /* facet springs */
    /* x */
    nf=0;
    for(x=0;x<nx+1;x++)
	for(y=0;y<ny;y++)
	    for(z=0;z<nz;z++){
		i=x*ny*nz + y*nz + z;

		if(f[0][i]){
		    fv0=x*(ny+1)*(nz+1) + y*(nz+1) + z;
		    fv1=fv0+nz+1;
		    fv2=fv1+1;
		    fv3=fv0+1;
					
		    s=&(obj->spring[ns++]);
		    s->point0=&(obj->point[v[fv0]]);
		    s->point1=&(obj->point[v[fv2]]);
		    s->l0=length(sub(s->point0->x,s->point1->x));
					
		    s=&(obj->spring[ns++]);
		    s->point0=&(obj->point[v[fv1]]);
		    s->point1=&(obj->point[v[fv3]]);
		    s->l0=length(sub(s->point0->x,s->point1->x));


		    /* find outer facets */
		    /*       
			     ----axis direction---->
			     _______
			     |
			     _____|    2

			     _______
			     /
			     3  |_______
					  

		    */
		    if(x==0){
			f[0][i]=2;
			nf++;
		    }
		    else if(x==nx){
			f[0][i]=3;
			nf++;
		    }
		    else {
			v0=vox0[i-ny*nz];
			v1=vox0[i];
			if(!v0 && v1){
			    f[0][i]=2;
			    nf++;
			}
			else if(v0 && !v1){
			    f[0][i]=3;
			    nf++;
			}
		    }


		}
	    }

    /* y */
    for(x=0;x<nx;x++)
	for(y=0;y<ny+1;y++)
	    for(z=0;z<nz;z++){
		i=x*(ny+1)*nz + y*nz + z;

		if(f[1][i]){
		    fv0=x*(ny+1)*(nz+1) + y*(nz+1) + z;
		    fv1=fv0+(ny+1)*(nz+1);
		    fv2=fv1+1;
		    fv3=fv0+1;
					
		    s=&(obj->spring[ns++]);
		    s->point0=&(obj->point[v[fv0]]);
		    s->point1=&(obj->point[v[fv2]]);
		    s->l0=length(sub(s->point0->x,s->point1->x));

		    s=&(obj->spring[ns++]);
		    s->point0=&(obj->point[v[fv1]]);
		    s->point1=&(obj->point[v[fv3]]);
		    s->l0=length(sub(s->point0->x,s->point1->x));


		    /* find outer facets */
		    if(y==0){
			f[1][i]=2;
			nf++;
		    }
		    else if(y==ny){
			f[1][i]=3;
			nf++;
		    }
		    else {
			j=x*ny*nz + y*nz + z;
			v0=vox0[j-nz];
			v1=vox0[j];
			if(!v0 && v1){
			    f[1][i]=2;
			    nf++;
			}
			else if(v0 && !v1){
			    f[1][i]=3;
			    nf++;
			}
		    }


		}
	    }

    /* z */
    for(x=0;x<nx;x++)
	for(y=0;y<ny;y++)
	    for(z=0;z<nz+1;z++){
		i=x*ny*(nz+1) + y*(nz+1) + z;

		if(f[2][i]){
		    fv0=x*(ny+1)*(nz+1) + y*(nz+1) + z;
		    fv1=fv0+(ny+1)*(nz+1);
		    fv2=fv1+nz+1;
		    fv3=fv0+nz+1;
					
		    s=&(obj->spring[ns++]);
		    s->point0=&(obj->point[v[fv0]]);
		    s->point1=&(obj->point[v[fv2]]);
		    s->l0=length(sub(s->point0->x,s->point1->x));

		    s=&(obj->spring[ns++]);
		    s->point0=&(obj->point[v[fv1]]);
		    s->point1=&(obj->point[v[fv3]]);
		    s->l0=length(sub(s->point0->x,s->point1->x));


		    /* find outer facets */
		    if(z==0){
			f[2][i]=2;
			nf++;
		    }
		    else if(z==nz){
			f[2][i]=3;
			nf++;
		    }
		    else {
			j=x*ny*nz + y*nz + z;
			v0=vox0[j-1];
			v1=vox0[j];
			if(!v0 && v1){
			    f[2][i]=2;
			    nf++;
			}
			else if(v0 && !v1){
			    f[2][i]=3;
			    nf++;
			}
		    }
		}
	    }

    printf("nf=%i\n",nf);

    /* facets */
    obj->facet = malloc(nf*sizeof(facet_t));
    obj->nf=nf;
    facet=obj->facet;

    /* facet x */

    osx=ny*nz;
    osy=nz;
    vos[0]=(ny+1)*(nz+1);
    vos[1]=nz+1;
    vos[2]=1;


    for(i=0;i<(nx+1)*ny*nz;i++){
	x=i/osx;
	y=(i-x*osx)/osy;
	z=i-x*osx-y*osy;
		
	ix=x*(ny+1)*(nz+1) + y*(nz+1) + z;
	iy=x*ny*(nz+1) + y*(nz+1) + z;
	iz=x*(ny+1)*nz + y*nz + z;

	if(f[0][i]==2 || f[0][i]==3){
	    facet->point=malloc(4*sizeof(mass_t*));
	    facet->np=4;
			
	    /* point */
	    j=x*(ny+1)*(nz+1) + y*(nz+1) + z;
	    if(f[0][i]==2) {a=0;b=1;c=2;d=3;}
	    else {a=3;b=2;c=1;d=0;}
			
	    facet->point[a]=&(obj->point[v[j]]);
	    facet->point[b]=&(obj->point[v[j+vos[1]]]);
	    facet->point[c]=&(obj->point[v[j+vos[1]+vos[2]]]);
	    facet->point[d]=&(obj->point[v[j+vos[2]]]);
			
	    /* edge */
	    nc=0;
	    if(f[0][i]==2) {a=0;b=2;c=3;d=1;}
	    else {a=2;b=0;c=3;d=1;}
			
	    /* y */
	    if(e[1][iy]==2) ce[nc++]=a;
	    if(e[1][iy+1]==2) ce[nc++]=b;
			
	    /* z */
	    if(e[2][iz]==2) ce[nc++]=c;
	    if(e[2][iz+nz]==2) ce[nc++]=d;
			
	    facet->nce=nc;
	    if(nc>0){
		facet->convexedge=malloc(nc*sizeof(int));
		memcpy(facet->convexedge,ce,nc*sizeof(int));
	    }
			
	    facet++;
	}
    }



    /* facet y */
    for(i=0;i<nx*(ny+1)*nz;i++){
	x=i/((ny+1)*nz);
	y=(i-x*(ny+1)*nz)/nz;
	z=i-x*(ny+1)*nz-y*nz;

	if(f[1][i]==2 || f[1][i]==3){
	    facet->point=malloc(4*sizeof(mass_t*));
	    facet->np=4;
			
	    /* point */
	    j=x*(ny+1)*(nz+1) + y*(nz+1) + z;
	    if(f[1][i]==2) {a=0;b=1;c=2;d=3;}
	    else {a=3;b=2;c=1;d=0;}
			
	    facet->point[a]=&(obj->point[v[j]]);
	    facet->point[b]=&(obj->point[v[j+vos[2]]]);
	    facet->point[c]=&(obj->point[v[j+vos[2]+vos[0]]]);
	    facet->point[d]=&(obj->point[v[j+vos[0]]]);
			
	    /* edge */
	    nc=0;
	    if(f[1][i]==2) {a=0;b=2;c=3;d=1;}
	    else {a=2;b=0;c=3;d=1;}
			
	    /* z */
	    j=x*(ny+1)*nz + y*nz + z;
	    if(e[2][j]==2) ce[nc++]=a;
	    if(e[2][j+(ny+1)*nz]==2) ce[nc++]=b;
			
	    /* x */
	    j=x*(ny+1)*(nz+1)+y*(nz+1)+z;
	    if(e[0][j]==2) ce[nc++]=c;
	    if(e[0][j+1]==2) ce[nc++]=d;


	    facet->nce=nc;
	    if(nc>0){
		facet->convexedge=malloc(nc*sizeof(int));
		memcpy(facet->convexedge,ce,nc*sizeof(int));
	    }
			
	    facet++;
	}
    }


    /* facet z */
    for(i=0;i<nx*ny*(nz+1);i++){
	x=i/(ny*(nz+1));
	y=i-x*ny*(nz+1);
	z=i-x*ny*(nz+1)-y*(nz+1);

	if(f[2][i]==2 || f[2][i]==3){
	    facet->point=malloc(4*sizeof(mass_t*));
	    facet->np=4;

	    /* point */
	    j=x*(ny+1)*(nz+1) + y*(nz+1) + z;
	    if(f[2][i]==2) {a=0;b=1;c=2;d=3;}
	    else {a=3;b=2;c=1;d=0;}

	    facet->point[a]=&(obj->point[v[j]]);
	    facet->point[b]=&(obj->point[v[j+vos[0]]]);
	    facet->point[c]=&(obj->point[v[j+vos[0]+vos[1]]]);
	    facet->point[d]=&(obj->point[v[j+vos[1]]]);

	    /* edge */
	    nc=0;
	    if(f[2][i]==2) {a=0;b=2;c=3;d=1;}
	    else {a=2;b=0;c=3;d=1;}

	    /* x */
	    j=x*(ny+1)*(nz+1)+y*(nz+1)+z;
	    if(e[0][j]==2) ce[nc++]=a;
	    if(e[0][j+nz+1]==2) ce[nc++]=b;

	    /* y */
	    j=x*ny*(nz+1) + y*(nz+1) + z;
	    if(e[1][j]==2) ce[nc++]=c;
	    if(e[1][j+ny*(nz+1)]==2) ce[nc++]=d;

	    facet->nce=nc;
	    if(nc>0){
		facet->convexedge=malloc(nc*sizeof(int));
		memcpy(facet->convexedge,ce,nc*sizeof(int));
	    }

	    facet++;
	}
    }


    printf("obj->np=%i\n",obj->np);

    /* volume diagonals */
    osx=(ny+1)*(nz+1);
    osy=nz+1;
    osz=1;

	
    for(i=0;i<nx*ny*nz;i++){
	x=i/(ny*nz);
	y=(i-x*ny*nz)/nz;
	z=i-x*ny*nz-y*nz;
		
	if(vox[i]){
	    cv[0]=x*osx + y*osy +z;
			
	    for(j=1;j<8;j++){
		cv[j]=cv[0] + ((j&4)>>2)*osx + ((j&2)>>1)*osy + (j&1); 
		/* stupid mistake: j&1 -- & is of extremely low priority. without (),
		   it will be equal to (... +j) & 1 */
	    }
			
	    for(j=0;j<4;j++){
		s=&(obj->spring[ns++]);
		s->point0=&(obj->point[v[cv[j]]]);
		s->point1=&(obj->point[v[cv[7-j]]]);
		s->l0=length(sub(s->point0->x,s->point1->x));
				
	    }
	}
    }
	
	



    /* clean up */
    free(v);
    for(i=0;i<3;i++){
	free(e[i]);
	free(f[i]);
    }
    if(hollow)
	free(vox);			

    /* elastic */
    for(i=0;i<obj->ns;i++){
	obj->spring[i].elastic=elastic;
	if(obj->spring[i].point0==NULL)
	    printf("spring%i point0==NULL\n", i);
	if(obj->spring[i].point1==NULL)
	    printf("spring%i point1==NULL\n", i);

    }


    return obj;

}






/* 
   Damping that produces a damping force proportional to the mass's velocity is commonly 
   referred to as "viscous damping" 
   
   
   references:
   http://en.wikipedia.org/wiki/Viscoelasticity
   http://en.wikipedia.org/wiki/Creep_(deformation)
   http://thayer.dartmouth.edu/defmech/
   
*/


/* 

   "semi-implicit euler"

   Explicit Euler is not stable:

   x = x + dt*v
   v = v + dt*a

   But, if you swap the order, you get Semi-Implicit Euler which is just as stable as Verlet and simpler:

   v = v + dt*a
   x = x + dt*v

   Erin


   see also:
   http://www.cs.unc.edu/~coombe/comp259/hw1/
   http://www.gamedev.net/community/forums/ViewReply.asp?id=1826867



   I have analyzed the stability of these integrators before. I have pasted my analysis below. Beware that you need to understand eigenvalues. (Symplectic Euler == Semi-Implicit Euler)

   -------------------------------------------------------------------------------
   Verlet has a stability limit like any explicit integrator, but it tends to
   be higher than many integrators (especially for the cost). And when Verlet
   is stable, it conserves energy on average.

   Consider a = -2*x:

   With Verlet I get: 1 -1 -1 1 1 -1 -1 1 1 -1 (stable)

   With the ballistics formula I get: 1 0 -2 -2 2 6 2 -10 -14 6
   (unstable)

   It is not hard to prove that Verlet and symplectic Euler are stable up to
   a = -4*x with dt = 1. The ballistic solver is _never_ stable for a positive stiffness. Here's the stability analysis to prove it:

   Verlet:
   x2 = 2*x1 - y1 - k*h^2*x1
   y2 = x1

   Growth matrix:
   A = [2-k*h^2 -1]
   [1 0]

   |lambda| <= 1 if k*h^2 <= 4


   Ballistic:
   x2 = x1 + h*v1 - 0.5*k*h^2*x1
   v2 = v1 - k*h*x1

   Growth matrix:
   A = [1-0.5*k*h^2 h]
   [-k*h 1]

   |lambda| <= 1 if k = 0


   Symplectic Euler:
   v2 = v1 - k*h*x1
   x2 = x1 + h*v2
   = x1 + h*v1 - k*h^2*x1

   Growth matrix:
   A = [1-k*h^2 h]
   [-k*h 1]

   |lambda| <= 1 if k*h^2 <= 4


   So for oscillatory motion, the exact ballistic solver is no better than
   explicit Euler (i.e. don't use it).

   Erin Catto


*/


// example code
#if 0
extern void drawedges();
extern void drawtri();
extern void draw_dog();

extern femobj_t *test, *test2;

cam_t *cam;

obj_t *cube;

int bx[4] = {400, -100, 200, 200};
int by[4] = {-200, -200, -100, -300};
void draw()
{
    int i;
    static int k = 0;
    extern int stop;

    glClearColor(1,1,1,1);
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );


    usecam(cam);
    glTranslatef(-1.5,-1.5,-10);


    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glLineWidth(1);

    glColor3f(1,0,0);


    if(!stop)
	simobj(cube);

    drawobj_wire(cube);

    glDisable(GL_DEPTH_TEST);


    SDL_GL_SwapBuffers();
}


elastic_t elastic={5000,1,0.01,0.01,10};

int smallcube[]={
    0,1,
    1,1,

};

#if 1
int cubevox[]={
    1,1,1,
    1,1,1,
    1,1,1,

    1,1,1,
    1,1,1,
    1,1,1,

    1,1,1,
    1,1,1,
    1,1,1,

};
#endif


#if 1
int smallchair[]={
    1,0,1,
    1,1,1,
    1,0,0,

    0,0,0,
    1,1,1,
    1,0,0,

    1,0,1,
    1,1,1,
    1,0,0,

};

int chair[]={
    1,0,0,1,
    1,0,0,1,
    1,1,1,1,
    1,0,0,0,
    1,0,0,0,
    1,0,0,0,

    0,0,0,0,
    0,0,0,0,
    1,1,1,1,
    1,0,0,0,
    1,0,0,0,
    1,0,0,0,

    0,0,0,0,
    0,0,0,0,
    1,1,1,1,
    1,0,0,0,
    1,0,0,0,
    1,0,0,0,

    1,0,0,1,
    1,0,0,1,
    1,1,1,1,
    1,0,0,0,
    1,0,0,0,
    1,0,0,0,


};


int bigcube2[]={
    0,0,0,1,0,
    1,0,1,0,1,
    0,1,0,1,0,
    1,0,1,0,1,
    0,1,0,1,0,

    1,0,1,0,1,
    0,1,1,1,0,
    1,1,1,1,1,
    0,1,1,1,0,
    1,0,1,0,1,

    0,1,0,1,0,
    1,1,1,1,1,
    0,1,1,1,0,
    1,1,1,1,1,
    0,1,0,1,0,

    1,0,1,0,1,
    0,1,1,1,0,
    1,1,1,1,1,
    0,1,1,1,0,
    1,0,1,0,1,

    0,1,0,1,0,
    1,0,1,0,1,
    0,1,0,1,0,
    1,0,1,0,1,
    0,1,0,1,0,

};

int bigcube[]={
    0,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,

    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,

    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,

    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,

    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,

};


int circle[] = {
    0,0,1,1,0,0,
    0,1,1,1,1,0,
    1,1,1,1,1,1,
    1,1,1,1,1,1,
    0,1,1,1,1,0,
    0,0,1,1,0,0	
};


#endif

extern float dx,dy,dz,du,dv,da;
void motion()
{
    extern int stop;
#if 1
    if(stop)
	usleep(1000*16);
#endif

    movecam(cam,dx,dy,dz);
    rotcam(cam,du,dv,da);
}

int hello[]={1,1,1,1};
void Init()
{
/* 	setelastic(&elastic,500,1,0.03,0.01,10); */

    int i;

    extern void testfem();
    testfem();

    cam=monitor2cam(17,4.0/3.0,20);

/* 	 cube=voxmesh(&hello,1,1,2,vec(1,1,1),1,&elastic,0);  */
/* 	cube=voxmesh(cubevox,3,3,3,vec(1,1,1),1,&elastic,0);  */
    /*      	cube=voxmesh(chair,4,6,4,vec(.5,.5,.5),.2,&elastic,0);      */     
//     	cube= voxmesh(circle,6,6,1,vec(.3,.3,.3),.5,&elastic,0);  
    cube=voxmesh(bigcube,5,5,5,vec(.3,.3,.3),.5,&elastic,0);   

    moveobj(cube,vec(0,10,0));

    printf("%i mass points\n",cube->np);
    printf("%i springs %i\n",cube->ns);

/* 	cube=voxmesh(smallcube,2,2,1,vec(1,1,1),1,&elastic); */

/* 	cube=makecube(vec(0,10+30,0),vec(2,12+30,2),8,&elastic); */

    InLoopInput();
    InLoop(motion);
    InLoop( draw );
}


#endif



// Rigid Body Dynamics
/* Rigid body dynamics simulation */
/**********************************/


static inline volatile long long rdtsc() {
    register long long tsc __asm__("eax");
//	__asm__ __volatile__ (".byte 15, 49" : : : "eax", "edx");
    __asm__ volatile ("rdtsc" : "=A" (tsc));

    return tsc;
} 

#define BEGIN_TIMING long long begin__ = rdtsc();
#define END_TIMING long long end__ = rdtsc(); printf("%qd cycles\n", end__-begin__);

#define BEGIN_TIMING2 begin__ = rdtsc();
#define END_TIMING2 end__ = rdtsc(); printf("%qd cycles\n", end__-begin__);

static float clamp(float a, float low, float high)
{
    return max(low, min(a,high));
}



typedef struct {
//	int edgeid[2];
//	float lambda[3];	/* each contact => three lambdas */
    unsigned feature;
    int idx;

} contactcache_t;

typedef struct{
    float Pn;
    float Pt[2];
    float Pnb;
    float nmass, tmass[2];
    float bias;
    v3_t tangent[2];
} contactdata_t;

/***********************************************************/
/* fixme: don't ... don't use boxhell. use something else. */
/***********************************************************/
#define MAX_NUM_BOX 100
#define MAX_NUM_JACK 10000
#define PRETTY_BIG 100000000000.0f
solid_t boxhell[MAX_NUM_BOX];
struct{
    float lambda0[MAX_NUM_JACK];
    float lambda1[MAX_NUM_JACK];
    float *lambda, *oldlambda;
    int nl, oldnl;
    contactcache_t cache0[MAX_NUM_JACK];
    contactcache_t cache1[MAX_NUM_JACK];
    contactcache_t *cache, *oldcache;
    int nc, oldnc;
    collinfo_t col[MAX_NUM_JACK];
    contactdata_t contact0[MAX_NUM_JACK];
    contactdata_t contact1[MAX_NUM_JACK];
    contactdata_t *contact, *oldcontact;
} meh_physics;


float d[MAX_NUM_JACK];
int numjacks = 0;
int numboxes = 0;
float beta = 5;
float g = 9.8;
float mu = 0.3;
int numiters = 10;
//contactcache_t ccache0[MAX_NUM_BOX * MAX_NUM_BOX * CONTACT_CACHE_SIZE];
//contactcache_t ccache1[MAX_NUM_BOX * MAX_NUM_BOX * CONTACT_CACHE_SIZE];
//contactcache_t *cachenow = ccache0;
//contactcache_t *cachethen = ccache1;
//int cacheusage0[MAX_NUM_BOX * MAX_NUM_BOX];
//int cacheusage1[MAX_NUM_BOX * MAX_NUM_BOX];
//int *cacheusagenow = cacheusage0;
//int *cacheusagethen = cacheusage1;

//#define ACCESS_CACHE(cache_, i, j, k)cache_[i*MAX_NUM_BOX*CONTACT_CACHE_SIZE + j*CONTACT_CACHE_SIZE + k]
//#define ACCESS_CACHE_USAGE(cu, i, j)cu[i*MAX_NUM_BOX + j]

void swapcache(){
//	swapptr(&cachenow, &cachethen);
//	swapptr(&cacheusagenow, &cacheusagethen);

    void *tmp = meh_physics.lambda;
    meh_physics.lambda = meh_physics.oldlambda;
    meh_physics.oldlambda = tmp;
    meh_physics.oldnl = meh_physics.nl;
    meh_physics.nl = 0;

    tmp = meh_physics.cache;
    meh_physics.cache = meh_physics.oldcache;
    meh_physics.oldcache = tmp;
    meh_physics.oldnc = meh_physics.nc;
    meh_physics.nc = 0;

    tmp = meh_physics.contact;
    meh_physics.contact = meh_physics.oldcontact;
    meh_physics.oldcontact = tmp;
}


/******************/
/* Update X and R */
/******************/
/*********************************************************/
/* Works ok because nothing wrong when in free motion,   */
/* and it cannot distinguish whether it's in free motion */
/* or not                                                */
/*********************************************************/
static void updateX(float dt){
    for(int k=0; k<numboxes; k++){
	solid_t *box = boxhell+k;
	rigid_t *rb = &(box->rigid);
	obb_t *obb = &(box->obb);
	rb->x = add(rb->x, scale(add(rb->bv, rb->v), dt)); /* x_new = x_old + v_new*dt */
	rb->q = addq(rb->q, scaleq(crossq(addq(rb->bomega,rb->omega), rb->q), dt/2.0f)); /* q_new = q_old +  */
	rb->q = normalizeq(rb->q);
	q2mat(rb->q, rb->R);
	obb->center = rb->x;
	matt(obb->axis, rb->R); /* each obb->axis[i] is a column of R, so u need to transpose R*/
    }
}

void draw_force(rigid_t *rb, v3_t f){
    v3_t v;

    v = add(rb->x, f);
    glBegin(GL_LINES);
    glColor3f(0,1,0);
    glVertex3fv(&(rb->x));
    glVertex3fv(&v);
    glEnd();

}


void printv(v3_t v, char *str){
    printf("%s(%f, %f, %f)\n", str, v.x, v.y, v.z);
}


v3_t clampv(v3_t v0, v3_t vmin, v3_t vmax){
    v3_t v;

    v.x = v0.x<vmax.x ? v0.x : vmax.x;
    v.y = v0.y<vmax.y ? v0.y : vmax.y;
    v.z = v0.z<vmax.z ? v0.z : vmax.z;
    v.x = v0.x>vmin.x ? v0.x : vmin.x;
    v.y = v0.y>vmin.y ? v0.y : vmin.y;
    v.z = v0.z>vmin.z ? v0.z : vmin.z;

    return v;
}

static void updateV(float dt){
    for(int i=0; i<numboxes; i++){
	rigid_t *rigid = &(boxhell[i].rigid);
	rigid = &boxhell[i].rigid;
	rigid->v = add(scale(rigid->F, dt*rigid->invmass), rigid->v);
	rigid->omega.v = add(scale(matv(rigid->invI, rigid->torque), dt), rigid->omega.v);
	rigid->bv = rigid->bomega.v = vec(0,0,0);
    }
}



void compute_invI(){
    int i;
    static m3_t Rt, IbodyRt;
    rigid_t *rb;

    for(i=0; i<numboxes; i++){
	rb = &(boxhell[i].rigid);
	matt(Rt, rb->R);
	matmat(IbodyRt, rb->invIbody, Rt);
	matmat(rb->invI, rb->R, IbodyRt);
    }
}


void addobb(v3_t x, float R[3][3], v3_t F, float invmass, v3_t v, v3_t omega, v3_t extent, int visible){
    rigid_t *rb;
    obb_t *obb;
    int i;

    if(numboxes > MAX_NUM_BOX){
	return;
    }

    boxhell[numboxes].visible = visible;

    rb = &(boxhell[numboxes].rigid);
    obb = &(boxhell[numboxes].obb);
    rb->x = x;
    rb->v = v;
    rb->torque = vec(0, 0, 0);
    rb->omega.v = omega;
    rb->omega.w = 0;
    memcpy(rb->R, R, sizeof(rb->R));


    rb->invmass = invmass;
    memset(rb->invIbody, 0, sizeof(rb->invIbody));
    if(invmass){
	rb->invIbody[0][0] = 3.0 /(extent.y*extent.y + extent.z*extent.z)* invmass;
	rb->invIbody[1][1] = 3.0 /(extent.x*extent.x + extent.z*extent.z)* invmass;
	rb->invIbody[2][2] = 3.0 /(extent.x*extent.x + extent.y*extent.y)* invmass;
    }

    rb->q = mat2q(R);
    rb->F = F;
    obb->center = rb->x;
    matt(obb->axis, rb->R);

    obb->extent[0] = extent.x;
    obb->extent[1] = extent.y;
    obb->extent[2] = extent.z;
    numboxes++;

}


/*
  currently three types of joints:

  ball-socket (3 linear constraints)
  universal (3 linear constraints + 1 angular constraint)
  hinge (3 linear + 2 anglar)

  joint limit

*/

/*
  compute jacobian for the ball-and-socket joint;
*/
typedef struct pjoint_s{
    struct pjoint_s *prev, *next;
    solid_t *s1, *s2;
    v3_t localpos1, localpos2;
    v3_t r1, r2;
    v3_t P;
    v3_t bias;
    float M[3][3];
} pjoint_t;

static pjoint_t *meh_jointlist=NULL;


/*
  J is one row of the Jacobian matrix
  (a complete constraint may consist of many rows)
*/
void addjoint(solid_t *s1, solid_t *s2, /* void (*compj)(pjoint_t*), int nrows,  */v3_t pos)
{
    pjoint_t *j = malloc(sizeof(pjoint_t));
    j->prev = NULL;
    j->next = meh_jointlist;
    if(j->next) j->next->prev = j;
    j->s1 = s1;
    j->s2 = s2;
    j->localpos1 = mattv3(s1->rigid.R, sub(pos, s1->rigid.x));
    j->localpos2 = mattv3(s2->rigid.R, sub(pos, s2->rigid.x));
    meh_jointlist = j;
    j->P = vec(0,0,0);
}


void testaddjoint()
{


#if 1
    solid_t *s1 = boxhell+9;
    solid_t *s2 = boxhell+10;
    addjoint(s1, s2, /* ballsocket_compj, 3,  */scale(add(s1->rigid.x,s2->rigid.x),.5));
#endif


#if 1
    s1 = boxhell+3;
    s2 = boxhell+10;
    addjoint(s1, s2, /* ballsocket_compj, 3,  */add(scale(add(s1->rigid.x,s2->rigid.x),.5), vec(2,0,0)));
#endif


#if 1
    s1 = boxhell+9;
    s2 = boxhell+8;
    addjoint(s1, s2, /* ballsocket_compj, 3,  */scale(add(s1->rigid.x,s2->rigid.x),.5));

    s1 = boxhell+8;
    s2 = boxhell+7;
    addjoint(s1, s2, /* ballsocket_compj, 3, */ scale(add(s1->rigid.x,s2->rigid.x),.5));
#endif

}


void deljoint(pjoint_t *j){
    if(j->prev){
	j->prev->next = j->next;
    }
    else{
	meh_jointlist = j->next;
    }

    if(j->next){
	j->next->prev = j->prev;
    }

//	free(j->J);
    free(j);
}

static inline float vmatv_inl(float v[3], float mat[3][3])
{
    return  (mat[0][0]*v[0] + mat[0][1]*v[1] + mat[0][2]*v[2])*v[0] +
	(mat[1][0]*v[0] + mat[1][1]*v[1] + mat[1][2]*v[2])*v[1] +
	(mat[2][0]*v[0] + mat[2][1]*v[1] + mat[2][2]*v[2])*v[2];
}





static float det3x3(float m[3][3]){
    return  m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2]) - 
	m[0][1]*(m[1][0]*m[2][2] - m[2][0]*m[1][2]) +
	m[0][2]*(m[1][0]*m[2][1] - m[2][0]*m[1][1]);
}



static void inv3x3(float ret[3][3], float m[3][3]) {
    float deti = 1.0f/det3x3(m);
	
    ret[0][0] = deti*(m[1][1]*m[2][2] - m[1][2]*m[2][1]);
    ret[0][1] = deti*(m[0][2]*m[2][1] - m[0][1]*m[2][2]);
    ret[0][2] = deti*(m[0][1]*m[1][2] - m[0][2]*m[1][1]);
	
    ret[1][0] = deti*(m[1][2]*m[2][0] - m[1][0]*m[2][2]);
    ret[1][1] = deti*(m[0][0]*m[2][2] - m[0][2]*m[2][0]);
    ret[1][2] = deti*(m[0][2]*m[1][0] - m[0][0]*m[1][2]);
	
    ret[2][0] = deti*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);
    ret[2][1] = deti*(m[0][1]*m[2][0] - m[0][0]*m[2][1]);
    ret[2][2] = deti*(m[0][0]*m[1][1] - m[0][1]*m[1][0]);
}    


static void screw(v3_t r, float m[3][3])
{
    m[0][0] = 0.0f;	m[0][1] = -r.z;	m[0][2] = r.y;
    m[1][0] = r.z;	m[1][1] = 0.0f;	m[1][2] = -r.x;
    m[2][0] = -r.y;	m[2][1] = r.x;	m[2][2] = 0.0f;
}

#define MAX_DIST 0.05

/*
  currently the mass matrix is diagnal
*/
void joint_prestep(pjoint_t *joint, float dt)
{
    rigid_t *b1 = &joint->s1->rigid;
    rigid_t *b2 = &joint->s2->rigid;

    joint->r1 = matv(b1->R, joint->localpos1);
    joint->r2 = matv(b2->R, joint->localpos2);


    float K[3][3];
    memset(K, 0, sizeof(K));
    K[0][0] = K[1][1] = K[2][2] = b1->invmass + b2->invmass;
	

    float screw1[3][3];
    screw(joint->r1, screw1);
    float tmp[3][3], tmp2[3][3];
    matmat(tmp, b1->invI, screw1);
    matmat(tmp2, screw1, tmp);

    for(int i=0; i<3; i++){
	for(int j=0; j<3; j++){
	    K[i][j] -= tmp2[i][j];
	}
    }

    float screw2[3][3];
    screw(joint->r2, screw2);
    matmat(tmp, b2->invI, screw2);
    matmat(tmp2, screw2, tmp);

    for(int i=0; i<3; i++){
	for(int j=0; j<3; j++){
	    K[i][j] -= tmp2[i][j];
	}
    }

    inv3x3(joint->M, K);

    v3_t w1 = add(b1->x, joint->r1);
    v3_t w2 = add(b2->x, joint->r2);

    v3_t C = sub(w2,w1);

    glPointSize(5);
    glBegin(GL_POINTS);
    glColor3f(1,0,0);
    glVertex3fv(&w1);
    glColor3f(0,1,0);
    glVertex3fv(&w2);
    glEnd();


    glBegin(GL_LINES);
    glVertex3fv(&w1);
    glVertex3fv(&b1->x);
    glVertex3fv(&w2);
    glVertex3fv(&b2->x);
    glEnd();

//	float len = length(d);

    /* TESTING */
//	joint->axis = normalize(C);
    /* END TESTING */

//	joint->axis = C;

    joint->bias = scale(C, -0.2f/dt);


#if 1

    b1->v = sub(b1->v, scale(joint->P, b1->invmass));
    b1->omega.v = sub(b1->omega.v, matv(b1->invI, cross(joint->r1, joint->P)));
    b2->v = add(b2->v, scale(joint->P, b2->invmass));
    b2->omega.v = add(b2->omega.v, matv(b2->invI, cross(joint->r2, joint->P)));
#endif


}



void test_updatex(rigid_t *body, v3_t *x, float R[3][3], float dt)
{
    *x = add(body->x, scale(add(body->bv, body->v), dt)); /* x_new = x_old + v_new*dt */
    q_t q = addq(body->q, scaleq(crossq(addq(body->bomega,body->omega), body->q), dt/2.0f)); /* q_new = q_old +  */
    q = normalizeq(q);
    q2mat(q, R);	
}

void joint_applyimpulse(pjoint_t *joint, float dt)
{
    rigid_t *body1 = &joint->s1->rigid;
    rigid_t *body2 = &joint->s2->rigid;

    v3_t dv = sub(add(body2->v, cross(body2->omega.v, joint->r2)), 
		  add(body1->v, cross(body1->omega.v, joint->r1)));

    v3_t zero={0,0,0};
    v3_t impulse = matv(joint->M, sub(joint->bias, dv));

#if 0
    printf("P =(%f %f %f)\n", impulse.x, impulse.y, impulse.z);
    printf("dv=(%f %f %f)\n", dv.x, dv.y, dv.z);
    printf("b =(%f %f %f)\n", joint->bias.x, joint->bias.y, joint->bias.z);
#endif


    body1->v = sub(body1->v, scale(impulse, body1->invmass));
    body1->omega.v = sub(body1->omega.v, scale(    matv(body1->invI, cross(joint->r1, impulse)), 1.f    ));
    body2->v = add(body2->v, scale(impulse, body2->invmass));
    body2->omega.v = add(body2->omega.v, scale(    matv(body2->invI, cross(joint->r2, impulse)), 1.f    ));



    /* TESTING */
    v3_t x1, x2;
    float R1[3][3], R2[3][3];
    test_updatex(body1, &x1, R1, dt);
    test_updatex(body2, &x2, R2, dt);

    v3_t w1 = add(x1, matv(R1, joint->localpos1));
    v3_t w2 = add(x2, matv(R2, joint->localpos2));
    v3_t C = sub(w2, w1);

#if 0
    printf("oldC =(%f %f %f)\n", joint->axis.x, joint->axis.y, joint->axis.z);
    printf("C =(%f %f %f)\n\n", C.x, C.y, C.z);
#endif


    joint->P = add(joint->P, impulse);

    /* END TESTING */

}



void draw_contact(collinfo_t *col){
    v3_t v;

    v = scale(col->normal, 10);
    v = add(col->poi[0], v);
#if 1
    glPointSize(5);
    glBegin(GL_POINTS);

    glColor3f(1,0,0);
    glVertex3fv(&(col->poi[0]));

    glColor3f(0,0,1);
    glVertex3fv(&(col->poi[1]));

    glEnd();
#endif

#if 1
    glBegin(GL_LINES);
    glColor3f(1,0,1);
    glVertex3fv((float *)&(col->poi[0]));
    glVertex3fv((float *)&v);
//	glColor3f(0,0,1);
//	glVertex3fv((float *)&(col->poi[1]));	
#endif
    glEnd();

}


/*
  friction


  In 3D, friction occurs in a plane. Since it is a plane, friction occurs in two directions (axes). 
  If you align the friction axes with the sliding motion, then you will get more realistic friction.

  The problem with using fixed directions for friction is that those directions will have weaker friction than diagonal directions. 
  Thus sliding boxes will steer towards the axial directions. This is bad.

  Twist friction is modeled by requiring the relative angular velocity about the contact normal to be zero. 
  The twist friction impulse torque needs to be clamped by:

  max_friction_torque = mu * lever_arm * accumulated_normal_impulse

  You should pre-compute the lever_arm as the average distance of all contact points on the manifold from the center of the manifold.


  Actually I am looking into improving the friction model in Erin GS solver. 
  At the moment I only solve the friction in the manifold center and apply an additional twist impulse to mimic the torsional friction. 
  Another trick (also suggested on the forum here) is to align one friction direction with the relative velocity in the contact plane.




*/

#define ALIGN_FRICTION_TO_VELOCITY 1

static void prestep(float dt)
{
    float penetration = 0.1f;
    float bias = 0.8f;

#if ALIGN_FRICTION_TO_VELOCITY
    for(int i=0; i<meh_physics.nc; i++){
	contactcache_t *c = meh_physics.cache + i;
	contactdata_t *cd = meh_physics.contact + i;
	collinfo_t *cinfo = meh_physics.col + c->idx;
	rigid_t *body1 = &boxhell[cinfo->id[0]].rigid;
	rigid_t *body2 = &boxhell[cinfo->id[1]].rigid;

	v3_t r1 = sub(cinfo->poi[0], body1->x);
	v3_t r2 = sub(cinfo->poi[0], body2->x);

	/* relative velocity between the two contact points, in world space. */
	v3_t vrel = sub(add(body2->v, cross(body2->omega.v, r2)),add(body1->v, cross(body1->omega.v, r1)));
	float vn = dot(vrel, cinfo->normal);
	v3_t tangent0 = normalize(sub(vrel, scale(cinfo->normal, vn)));
	v3_t tangent1 = normalize(cross(tangent0, cinfo->normal));
	cd->tangent[0] = tangent0;
	cd->tangent[1] = tangent1;
    }
#endif

    for(int i=0; i<meh_physics.nc; i++){
	contactcache_t *c = meh_physics.cache + i;
	contactdata_t *cd = meh_physics.contact + i;
	collinfo_t *cinfo = meh_physics.col + c->idx;
	rigid_t *body1 = &boxhell[cinfo->id[0]].rigid;
	rigid_t *body2 = &boxhell[cinfo->id[1]].rigid;

	v3_t r1 = sub(cinfo->poi[0], body1->x);
	v3_t r2 = sub(cinfo->poi[0], body2->x);

	/* normal impulse's denom */
//		float term1 = dot(cinfo->normal, cross(matv(body1->invI, cross(r1, cinfo->normal)), r1));
//		float term2 = dot(cinfo->normal, cross(matv(body2->invI, cross(r2, cinfo->normal)), r2));


	v3_t c1 = cross(r1, cinfo->normal);
	v3_t c2 = cross(r2, cinfo->normal);
	float term1 = dot(matv(body1->invI, c1),c1);
	float term2 = dot(matv(body2->invI, c2),c2);


	cd->nmass = 1.0f / (body1->invmass + body2->invmass + term1 + term2);


#if ALIGN_FRICTION_TO_VELOCITY
	v3_t tangent0 = cd->tangent[0];
	v3_t tangent1 = cd->tangent[1];
#else

	v3_t tangent0 = cd->tangent[0] = cinfo->tangent[0];
	v3_t tangent1 = cd->tangent[1] = cinfo->tangent[1];
#endif

	/* tangent[0] impulse's denom */
	term1 = dot(tangent0, cross(matv(body1->invI, cross(r1, tangent0)), r1));
	term2 = dot(tangent0, cross(matv(body2->invI, cross(r2, tangent0)), r2));
	cd->tmass[0] = 1.0f / (body1->invmass + body2->invmass + term1 + term2);

	/* tangent[1] impulse's denom */
	term1 = dot(tangent1, cross(matv(body1->invI, cross(r1, tangent1)), r1));
	term2 = dot(tangent1, cross(matv(body2->invI, cross(r2, tangent1)), r2));
	cd->tmass[1] = 1.0f / (body1->invmass + body2->invmass + term1 + term2);

	cd->bias = -bias * min(0.0f, cinfo->depth + penetration) / dt;

	/* accum impulse */
	v3_t P = add(scale(cinfo->normal, cd->Pn), 
		     add(scale(tangent0, cd->Pt[0]), 
			 scale(tangent1, cd->Pt[1])));

	body1->v = sub(body1->v, scale(P,body1->invmass));
	body1->omega.v = sub(body1->omega.v, matv(body1->invI, cross(r1, P)));

	body2->v = add(body2->v, scale(P,body2->invmass));
	body2->omega.v = add(body2->omega.v, matv(body2->invI, cross(r2, P)));
    }
}


static void applyimpulse(contactcache_t *c, contactdata_t *cd)
{
    collinfo_t *cinfo = meh_physics.col + c->idx;

    rigid_t *b1 = &boxhell[cinfo->id[0]].rigid;
    rigid_t *b2 = &boxhell[cinfo->id[1]].rigid;

    v3_t r1 = sub(cinfo->poi[0], b1->x);
    v3_t r2 = sub(cinfo->poi[0], b2->x);

    /* relative velocity */
    v3_t vrel = sub(add(b2->v, cross(b2->omega.v, r2)),add(b1->v, cross(b1->omega.v, r1)));

    /* normal impulse */
    float vn = dot(vrel, cinfo->normal);
    float dPn = cd->nmass * (-vn);
    float Pn0 = cd->Pn;
    cd->Pn = max(Pn0 + dPn, 0.0f);
    dPn = cd->Pn - Pn0;

    /* apply contact impulse */
    v3_t Pn = scale(cinfo->normal, dPn);
    b1->v = sub(b1->v, scale(Pn, b1->invmass));
    b1->omega.v = sub(b1->omega.v, matv(b1->invI, cross(r1, Pn)));
    b2->v = add(b2->v, scale(Pn, b2->invmass));
    b2->omega.v = add(b2->omega.v, matv(b2->invI, cross(r2, Pn)));

    /* split impulse */
    vrel = sub(add(b2->bv, cross(b2->bomega.v, r2)), 
	       add(b1->bv, cross(b1->bomega.v, r1)));
    float vnb = dot(vrel, cinfo->normal);
    float dPnb = cd->nmass * (-vnb + cd->bias);
    float Pnb0 = cd->Pnb;
    cd->Pnb = max(Pnb0 + dPnb, 0.0f);
    dPnb = cd->Pnb - Pnb0;

    v3_t Pb = scale(cinfo->normal, dPnb);
    b1->bv = sub(b1->bv, scale(Pb, b1->invmass));
    b1->bomega.v = sub(b1->bomega.v, matv(b1->invI, cross(r1, Pb)));
    b2->bv = add(b2->bv, scale(Pb, b2->invmass));
    b2->bomega.v = add(b2->bomega.v, matv(b2->invI, cross(r2, Pb)));


    /* friction1 */
    /* relative velocity */
    vrel = sub(add(b2->v, cross(b2->omega.v, r2)),add(b1->v, cross(b1->omega.v, r1)));
    v3_t tangent0 = cd->tangent[0];
    v3_t tangent1 = cd->tangent[1];

    float vt = dot(vrel, tangent0);
    float dPt = cd->tmass[0] * (-vt);
	
    float maxPt = 0.4f/*  friction */ * cd->Pn;
    float Pt0 = cd->Pt[0];
    cd->Pt[0] = clamp(Pt0+dPt, -maxPt, maxPt);
    dPt = cd->Pt[0] - Pt0;

    /* apply friction impulse */
    v3_t Pt = scale(tangent0, dPt);
    b1->v = sub(b1->v, scale(Pt, b1->invmass));
    b1->omega.v = sub(b1->omega.v, matv(b1->invI, cross(r1, Pt)));
    b2->v = add(b2->v, scale(Pt, b2->invmass));
    b2->omega.v = add(b2->omega.v, matv(b2->invI, cross(r2, Pt)));

    /* friction2 */
    /* relative velocity */
    vrel = sub(add(b2->v, cross(b2->omega.v, r2)),add(b1->v, cross(b1->omega.v, r1)));
    vt = dot(vrel, tangent1);
    dPt = cd->tmass[1] * (-vt);
    maxPt = /* friction */.4f * cd->Pn;
    float Pt1 = cd->Pt[1];
    cd->Pt[1] = clamp(Pt1+dPt, -maxPt, maxPt);
    dPt = cd->Pt[1] - Pt1;

    /* apply friction impulse */
    Pt = scale(tangent1, dPt);
    b1->v = sub(b1->v, scale(Pt, b1->invmass));
    b1->omega.v = sub(b1->omega.v, matv(b1->invI, cross(r1, Pt)));
    b2->v = add(b2->v, scale(Pt, b2->invmass));
    b2->omega.v = add(b2->omega.v, matv(b2->invI, cross(r2, Pt)));
}



/*
  store all contacts in an array (unsorted);
  sort the array according to contact feature (boxid1 edge1 boxid2 edge2)
  sweep the array over the old contact array, find same feature pairs and use the stored lambda;
*/


static void initdynamics()
{
    meh_physics.lambda = meh_physics.lambda0;
    meh_physics.oldlambda = meh_physics.lambda1;
    meh_physics.oldnl = meh_physics.nl = 0;

    meh_physics.cache = meh_physics.cache0;
    meh_physics.oldcache = meh_physics.cache1;
    meh_physics.oldnc = meh_physics.nc = 0;

    meh_physics.contact = meh_physics.contact0;
    meh_physics.oldcontact = meh_physics.contact1;
}

static int cmpcontact(void *p, void *p1)
{
    contactcache_t *c=p, *c2=p1;
    return c->feature>c2->feature ? 1 : -1;
}


void nextstep(float dt){
    int i, j, k;
    static collinfo_t col[10];
    int numcols;
    int l;
    static int firstrun = 1;


    if(firstrun){
	firstrun = 0;
	initdynamics();
    }


//	BEGIN_TIMING;

    /* TODO: sweep and prune */
    int n=0;
    for(int i=0; i<numboxes-1; i++){
	for(int j=i+1; j<numboxes; j++){
	    if(boxhell[i].rigid.invmass==0.0f && boxhell[j].rigid.invmass==0.0f) continue;
	    numcols = boxboxintersect(&(boxhell[i]), &(boxhell[j]), i, j, col);
	    for(int k=0; k<numcols; k++){
		meh_physics.col[n] = col[k];
		meh_physics.cache[n].feature = col[k].feature;
		meh_physics.cache[n].idx = n;

//				draw_contact(col+k);
		n++;
	    }
	}
    }

//	printf("collision cycles ");
//	END_TIMING;

    /* sort contact cache */
    qsort(meh_physics.cache, n, sizeof(contactcache_t), cmpcontact);
    meh_physics.nc = n;
    meh_physics.nl = 3*n;

    /* sweep against the old contact cache */
    contactcache_t *oc = meh_physics.oldcache;
    contactcache_t *c = meh_physics.cache;
    contactcache_t *oldend = meh_physics.oldcache + meh_physics.oldnc;
    contactcache_t *end = meh_physics.cache + meh_physics.nc;
    contactdata_t *contact = meh_physics.contact;
    contactdata_t *oldcontact = meh_physics.oldcontact;
    while(oc < oldend && c < end){
	if(c->feature == oc->feature){
	    contact->Pn = oldcontact->Pn;
	    contact->Pt[0] = oldcontact->Pt[0];
	    contact->Pt[1] = oldcontact->Pt[1];
	    contact->Pnb = oldcontact->Pnb;
	    c++;
	    contact++;
	}
	else if(c->feature < oc->feature){
	    contact->Pn = contact->Pt[0] = contact->Pt[1] = contact->Pnb = 0.0f;
	    c++;
	    contact++;
	}
	else{
	    oc++;
	    oldcontact++;
	}
    }
    if(c<=end){
	memset(contact, 0, meh_physics.contact+n-contact);
//		memset(c, 0, meh_physics.cache+n-c);
    }	

    compute_invI();

    updateV(dt);

//       	BEGIN_TIMING2;
/* TODO: add other constraints */


    prestep(dt);

    for(pjoint_t *j=meh_jointlist; j; j=j->next){
//		joint_presolve(j, dt);

	joint_prestep(j,dt);
//		joint_presolve3(j,dt);
    }

    for(int i=0; i<10; i++){
	for(int i=0; i<meh_physics.nc; i++){
	    applyimpulse(meh_physics.cache+i, meh_physics.contact+i);
	}

	for(pjoint_t *j=meh_jointlist; j; j=j->next){
	    joint_applyimpulse(j,dt);
	}
    }



//	printf("solver cycles ");
//	END_TIMING2;

    updateX(dt);
    swapcache();

    //	printf("%i boxes\n", numboxes);
}


void copy_vec_to_column2(v3_t v, m4_t m, int column){
    int i;

    i = 4 * column;
    m[i++] = v.x;
    m[i++] = v.y;
    m[i] = v.z;

}




int face2verttab[6][4] = {{4,5,7,6},
			  {6,7,3,2},
			  {3,7,5,1},
			  {0,2,3,1},
			  {5,4,0,1},
			  {4,6,2,0}
};


void drawobb(obb_t *b, int loc)
{
    static float s1[4]={-1,-1,1,1}, s2[4]={-1,1,1,-1};
    static float tc[8]={0,0,0,1,1,1,1,0};	
    float *extent = &b->extent;
    glBegin(GL_QUADS);
    for(int a=0; a<3; a++){
	int a1=(a+1)%3, a2=(a+2)%3;
	v3_t N = scale(b->axis[a], -1.0f);
	v3_t T = b->axis[a2];
	glNormal3fv(&N);
	glVertexAttrib3fv(loc, &T);

	/* min */
	float vert[3];
	vert[a] = -extent[a];
	for(int v=0; v<4; v++){
	    /* in object space */
	    vert[a1] = s1[v]*extent[a1];
	    vert[a2] = s2[v]*extent[a2];

	    /* to world space */
	    glTexCoord2fv(tc+v*2);
	    v3_t wldv = obb2world(b, vert);
	    glVertex3fv((float*)&wldv);
	}

	/* max */
	N = b->axis[a];
	T = scale(b->axis[a2], -1);
	glNormal3fv(&N);
	glVertexAttrib3fv(loc, &T);
	vert[a] = extent[a];
	for(int v=0; v<4; v++){
	    /* in object space */
	    vert[a1] = s1[3-v]*extent[a1];
	    vert[a2] = s2[3-v]*extent[a2];

	    /* to world space */
	    glTexCoord2fv(tc+v*2);
	    v3_t wldv = obb2world(b, vert);
	    glVertex3fv((float*)&wldv);
	}
		
    }
    glEnd();
}


void drawobbs(int loc)
{
    for(int i=0; i<numboxes; i++){
	if(boxhell[i].visible)
	    drawobb(&boxhell[i].obb, loc);
    }
}



void movebox(int id, float x, float y, float z)
{
    if(id<0 || id>numboxes-1)return;

    obb_t *box = &boxhell[id].obb;
    box->center.x += x;
    box->center.y += y;
    box->center.z += z;
}


