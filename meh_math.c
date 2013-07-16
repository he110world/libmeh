#include "meh_math.h"

v3_t add(v3_t a, v3_t b)
{
    v3_t v={a.x+b.x, a.y+b.y, a.z+b.z}; return v;
}


v3_t sub(v3_t a, v3_t b)
{
    v3_t v={a.x-b.x, a.y-b.y, a.z-b.z}; return v;
}


v3_t scale(v3_t a, float s)
{
    v3_t v={a.x*s, a.y*s, a.z*s}; return v;
}

float dot(v3_t a, v3_t b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

v3_t cross(v3_t a, v3_t b)
{
    v3_t v={a.y*b.z - b.y*a.z, a.z*b.x - b.z*a.x, a.x*b.y - b.x*a.y}; return v;
}

v3_t vec(float x, float y, float z)
{
    v3_t v={x,y,z};	return v;
}

v3_t normalize(v3_t v)
{
    float l2=v.x*v.x+v.y*v.y+v.z*v.z;
    if (!l2) return v;

    float l=sqrtf(l2);
    v.x/=l;
    v.y/=l;
    v.z/=l;
    return v;
}

/* using M[x][y] to index a element */
mat_t newmat(int r, int c)
{
    mat_t mat;
    float *m;
    float *mr;
    int i;

    mat = malloc(r*sizeof(float *));
    m = malloc(r*c*sizeof(float));
    for(i=0, mr=m; i<r; i++){
	mat[i] = mr;
	mr += c;
    }

    return mat;
}


v3_t matv3(float m[3][3], v3_t v){
    v3_t mv;
    mv.x = m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z;
    mv.y = m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z;
    mv.z = m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z;
    return mv;
}

/* mat^T v3 */
v3_t mattv3(float m[3][3], v3_t v){
    v3_t mv;
    mv.x = m[0][0]*v.x + m[1][0]*v.y + m[2][0]*v.z;
    mv.y = m[0][1]*v.x + m[1][1]*v.y + m[2][1]*v.z;
    mv.z = m[0][2]*v.x + m[1][2]*v.y + m[2][2]*v.z;
    return mv;
}

// rotate plane by matrix
void matpln(float m[3][3], plane_t *pln)
{
    pln->n = matv3(m, pln->n);
}

void makemat(int r, int c, mat_t mat, float m[])
{
    int i;
    float *mr;

    for(i=0, mr=m; i<r; i++){
	mat[i] = mr;
	mr += c;
    }
}


void delmat(mat_t mat)
{
    free(mat[0]);
    free(mat);
}


void eyemat(mat_t mat, int n)
{
    int i;

    memset(mat[0], 0, n*n*sizeof(float));
    for(i=0; i<n; i++){
	mat[i][i] = 1;
    }
}


v3_t matvec(m3_t m, v3_t v0){
    v3_t v;

    v.x = m[0]*v0.x + m[3]*v0.y + m[6]*v0.z;
    v.y = m[1]*v0.x + m[4]*v0.y + m[7]*v0.z;
    v.z = m[2]*v0.x + m[5]*v0.y + m[8]*v0.z;

    return v;
}


/* matrix transposition */
void matt(float dest[3][3], float src[3][3])
{
    float tmp;

    tmp = src[0][1];
    dest[0][1] = src[1][0];
    dest[1][0] = tmp;

    tmp = src[0][2];
    dest[0][2] = src[2][0];
    dest[2][0] = tmp;

    tmp = src[2][1];
    dest[2][1] = src[1][2];
    dest[1][2] = tmp;

    if(dest!=src){
	dest[0][0]=src[0][0];
	dest[1][1]=src[1][1];
	dest[2][2]=src[2][2];
    }

}


void matmat(float dest[3][3], float left[3][3], float right[3][3])
{
    assert(left!=right);
    for(int i=0; i<3; i++)
	for(int j=0; j<3; j++)
	    dest[i][j] = left[i][0]*right[0][j] + left[i][1]*right[1][j] + left[i][2]*right[2][j];
}

void cpmat(m3_t dest, m3_t src){
    int i;

    if(dest == src){
	return;
    }
    for(i=0; i<9; i++){
	dest[i] = src[i];
    }

}


float length(v3_t v){
    return sqrtf(v.x*v.x+v.y*v.y+v.z*v.z);
}


v4_t m4v4(float m[16], v4_t v0){
    v4_t v;

    v.x = m[0]*v0.x + m[4]*v0.y + m[8]*v0.z + m[12]*v0.w;
    v.y = m[1]*v0.x + m[5]*v0.y + m[9]*v0.z + m[13]*v0.w;
    v.z = m[2]*v0.x + m[6]*v0.y + m[10]*v0.z + m[14]*v0.w;
    v.w = m[3]*v0.x + m[6]*v0.y + m[11]*v0.z + m[15]*v0.w;

    return v;
}


/* Gram-Schimdt Orthonormalization */
void orthonormalization3x3(float m[3][3])
{
    v3_t row0, row1, row2;
    float l0, l1, l2;

    row0 = *((v3_t *)m[0]);
    row1 = *((v3_t *)m[1]);
    row2 = *((v3_t *)m[2]);

    l0 = length(row0);
    if(l0 > 0) row0 = scale(row0, 1.0/l0);

    row1 = sub(row1, scale(row1, dot(row0, row1)));
    l1 = length(row1);
    if(l1 > 0) row1 = scale(row1, 1.0/l1);

    row2 = cross(row0, row1);

    *((v3_t *)m[0]) = row0;
    *((v3_t *)m[1]) = row1;
    *((v3_t *)m[2]) = row2;
}


v3_t mat3x3_v3(float m[3][3], v3_t v)
{
    v3_t result;

    result.x = m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z;
    result.y = m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z;
    result.z = m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z;

    return result;
}


/*------------------------------------
  LU decomposition 
  -----------------------------------*/
int ludecomp(float A[], float LU[], int P[], int N){
    int i, j, k;
    float max, tmp;
    int who, tmpi;
    float sum;


    // init if necessary
    if(A != LU){
	for(i=0; i<N*N; i++){
	    LU[i] = A[i];
	}
    }

    // init permutation
    for(i=0; i<N; i++){
	P[i] = i;
    }

    // real stuff
    for(i=0; i<N; i++){

	// pivoting
	max = 0;
	who = -1;
	for(j=i; j<N; j++){
	    // find max element in current column
	    tmp = fabsf(MATRIX(LU, N, N, j, i));
	    if(tmp > max){
		max = tmp;
		who = j;
	    }
	}
	if(who == -1){
	    who = i;
	}

	// swap it with current row
	if(who != i){
	    for(j=0; j<N; j++){
		tmp = MATRIX(LU, N, N, i, j);
		MATRIX(LU, N, N, i, j)= MATRIX(LU, N, N, who, j);
		MATRIX(LU, N, N, who, j)= tmp;
	    }
	    tmpi = P[i];
	    P[i] = P[who];
	    P[who] = tmpi;
	}

	// row
	for(j=i; j<N; j++){
	    sum = 0;
	    for(k=0; k<i; k++){
		sum += MATRIX(LU, N, N, i, k)* MATRIX(LU, N, N, k, j);
	    }
	    MATRIX(LU, N, N, i, j)-= sum;
	}


	// column
	for(j=i+1; j<N; j++){
	    sum = 0;
	    for(k=0; k<i; k++){
		sum += MATRIX(LU, N, N, j, k)* MATRIX(LU, N, N, k, i);
	    }
	    MATRIX(LU, N, N, j, i)=(MATRIX(LU, N, N, j, i)- sum)/ MATRIX(LU, N, N, i, i);
	}
    }

}


/*------------------------------------------
  QR Decomposition
  -----------------------------------------*/
int qrdecomp(float A[], float Q[], float R[], int N){
    int row, col;
    float len;
    int i, j;
    float tmp;

    // init R
    for(j=0; j<N; j++){
	for(i=j+1; i<N; i++){
	    MATRIX(R, N, N, i, j)= 0;
	}
    }

    // real stuff
    for(col=0; col<N; col++){
	for(row=0; row<N; row++){
	    MATRIX(Q, N, N, row, col)= MATRIX(A, N, N, row, col);
	}
	for(j=0; j<col; j++){
	    len = 0;
	    for(i=0; i<N; i++){
		len += MATRIX(Q, N, N, i, j)* MATRIX(A, N, N, i, col);
	    }
	    MATRIX(R, N, N, j, col)= len;
	    for(i=0; i<N; i++){
		MATRIX(Q, N, N, i, col)-= len * MATRIX(Q, N, N, i, j);
	    }
	}
	len = 0;
	for(i=0; i<N; i++){
	    tmp = MATRIX(Q, N, N, i, col);
	    len += tmp*tmp;
	}
	len = sqrtf(len);
	if(len < EPSILON){
	    return 0;
	}
	for(i=0; i<N; i++){
	    MATRIX(Q, N, N, i, col)/= len;
	}
	MATRIX(R, N, N, col, col)= len;
    }

    return 1;
}


/*------------------------------------------------
  Forward and back substitution for LU
  -----------------------------------------------*/
static void lubacksub(float A[], float x[], float y[], float b0[], float b[], int P[], int N){
    float sum;
    int i, j;

    for(i=0; i<N; i++){
	b[i] = b0[P[i]];
    }
    for(i=0; i<N; i++){
	sum = 0;
	for(j=0; j<i; j++){
	    sum += MATRIX(A, N, N, i, j)* y[j];
	}
	y[i] = b[i] - sum;
    }
    for(i=N-1; i>-1; i--){
	sum = 0;
	for(j=N-1; j>i; j--){
	    sum += MATRIX(A, N, N, i, j)* x[j];
	}
	x[i] =(y[i] - sum)/ MATRIX(A, N, N, i, i);
    }
}


/*-------------------------------------------------------------------------
  Solve the linear algebra equation Ax = b using LU decomposition
  ------------------------------------------------------------------------*/
int lusolve(float A0[], float x[], float b0[], int N){
    int s, i, j;
    float tmp;
    float A[N*N], b[N], y[N];
    int P[N];
    memcpy(A, A0, sizeof(A));
    s = ludecomp(A, A, P, N);
    if(!s){
	return 0;
    }

    lubacksub(A, x, y, b0, b, P, N);
    return 1;
}


/*----------------------------------------------------
  Inversion of a matrix
  ---------------------------------------------------*/
int luinv(float A0[], float Ainv[], int N){
    int s;
    int i, j;
    int sz=N*sizeof(float);
    float A[N*sz], P[sz], x[sz], y[sz], b[sz], b0[sz];
    memcpy(A, A0, sizeof(A));

    // Decompose A to L and U
    s = ludecomp(A, A, P, N);
    if(!s){
	return 0;
    }

    // Solve each colume for Ainv
    for(j=0; j<N; j++){
	// Right hand side,(0,0,...,1,0,..,0)T
	for(i=0; i<N; i++){
	    b0[i] = 0;
	}
	b0[j] = 1;

	// Solve the jth colume of Ainv
	lubacksub(A, x, y, b0, b, P, N);

	// Copy colume to Ainv
	for(i=0; i<N; i++){
	    MATRIX(Ainv, N, N, i, j)= x[i];
	}
    }

    return 1;
}


void printmat(char *name, float A[], int N){
    int i, j;

    printf("%s =\n", name);
    for(i=0; i<N; i++){
	for(j=0; j<N; j++){
	    printf("%f\t", MATRIX(A, N, N, i, j));
	}
	printf("\n");
    }
    printf("\n");
}

void printvec(char *name, float v[], int N){
    int i;

    printf("%s =\n", name);
    for(i=0; i<N; i++){
	printf("%f\n", v[i]);
    }
    printf("\n");
}


/*-------------------------
  AB = C
  ------------------------*/
void matmul(float A[], float B[], float C0[], int N){
    float *C;
    int row, col, i;
    float sum;

    if(C0 == A || C0 == B){
	C = alloca(N*N*sizeof(float));
    }
    else {
	C = C0;
    }
    for(row=0; row<N; row++){
	for(col=0; col<N; col++){
	    sum = 0;
	    for(i=0; i<N; i++){
		sum += MATRIX(A, N, N, row, i)* MATRIX(B, N, N, i, col);
	    }
	    MATRIX(C, N, N, row, col)= sum;
	}
    }
    if(C != C0){
	for(i=0; i<N*N; i++){
	    C0[i] = C[i];
	}
    }
}


int eigenvalue(float A0[], float lambda[], int N){
    int iter, i;
    int row, col;
    float max;
    float tmp;
    int sz=N*N;
    float A[sz], Q[sz], QQ[sz], R[sz];
    memcpy(A,A0,sizeof(A));

    for(iter=0; iter<100; iter++){
	qrdecomp(A, Q, R, N);
	matmul(R, Q, A, N);
	max = 0;
	for(col=0; col<N; col++){
	    for(row=col+1; row<N; row++){
		tmp = fabsf(MATRIX(A, N, N, row, col));
		if(tmp > max){
		    max = tmp;
		}
	    }
	}
	if(max < 0.0001){
	    for(i=0; i<N; i++){
		lambda[i] = MATRIX(A, N, N, i, i);
	    }
	    return iter;
	}
    }
    for(i=0; i<N; i++){
	lambda[i] = MATRIX(A, N, N, i, i);
    }
    return 0;
}



/*-------------------------------------
  Copy A to B
  ------------------------------------*/
void cpmatrix(float A[], float B[], int N){
    int i;

    if(A==B){
	return;
    }
    for(i=0; i<N*N; i++){
	B[i] = A[i];
    }
}

/*----------------------------------------------------
  Eigenvector of symmetric matrix
  ---------------------------------------------------*/
int eigensym(float A0[], float lambda[], float U[], int N){
    int iter, i, j;
    int row, col;
    float max;
    float tmp;
    int sz=N*N;
    float A[sz], Q[sz], QQ[sz], R[sz];
    memcpy(A,A0,sizeof(A));
    for(i=0; i<N; i++){
	for(j=0; j<N; j++){
	    MATRIX(QQ, N, N, i, j)=(i==j)? 1 : 0;
	}
    }
    for(iter=0; iter<100; iter++){
	qrdecomp(A, Q, R, N);
	// The hack: eigenvectors = QQ = Qn*...*Q2*Q1
	matmul(QQ, Q, QQ, N);
	matmul(R, Q, A, N);
	max = 0;
	for(col=0; col<N; col++){
	    for(row=col+1; row<N; row++){
		tmp = fabsf(MATRIX(A, N, N, row, col));
		if(tmp > max){
		    max = tmp;
		}
	    }
	}
	if(max < 0.0001){
	    for(i=0; i<N; i++){
		lambda[i] = MATRIX(A, N, N, i, i);
	    }
	    cpmatrix(QQ, U, N);
	    return iter;
	}
    }
    for(i=0; i<N; i++){
	lambda[i] = MATRIX(A, N, N, i, i);
    }
    cpmatrix(QQ, U, N);

    return 0;

}



/**
   |Mij'|    |  c   -s  | | Mij|
   |    | =  |          | |    |
   |    |    |          | |    |
   |Mkl'|    |  s   c   | | Mkl|

*/
static void rotate(mat_t m, float c, float s, int i, int j, int k, int l)
{
    float mij;

    mij = c*m[i][j] - s*m[k][l];
    m[k][l] = s*m[i][j] + c*m[k][l];
    m[i][j] = mij;

}


/* only consider the upper triangle */
int maxinrow(mat_t M, int row, int c)
{
    int i;
    float max;
    int who;

    max = M[row][row+1];
    who = row+1;
    for(i=row+2; i<c; i++){
	if( M[row][i] > max ){
	    max = M[row][i];
	    who = i;
	}
    }

    return who;
}


void eyemat2(mat_t m, int n)
{
    for(int i=0; i<n; i++){
	for(int j=0; j<n; j++){
	    m[i][j] = (i==j);
	}
    }
}



/*----------------------------------------------------
  mir -- max (element) in a row
  ---------------------------------------------------*/
void jacobi_eigensys(int n, mat_t S, float e[], mat_t E, int mir[], int changed[])
{
    int nactives;
    int i;
    int m;
    int k,l;		// row, column to rotate
    float p,y,t,s,c;
    float tmp;
    float epsilon = 0.00001;

    eyemat2(E, n);
    nactives = n;

    for(i=0; i<n; i++){
	mir[i] = maxinrow(S,i,n);
	e[i] = S[i][i];
	changed[i]=1;
    }

    while(nactives>0){
	m=0;
	for(i=1; i<n; i++){
	    if(fabsf(S[i][mir[i]]) > fabsf(S[m][mir[m]])){
		m = i;
	    }
	}


	k = m;
	l = mir[k];

	p = S[k][l];
	y = (e[l]-e[k])*.5;
	t = fabsf(y)+sqrtf(p*p+y*y);
	s = sqrtf(p*p+t*t);
	c = t/s;
	s = fabsf(p/s);
	t = fabsf(p*p/t);

	S[k][l] = 0;

	/* the diagonal elements aren't stored in the matrix -- they're stored in the eigenvalue array directly
	 * --hence no copying needed
	 * */

	tmp = e[k];
	e[k] -= t;
	if(changed[k] && t<epsilon){
	    changed[k]=0;
	    nactives--;
	}
	else if(!changed[k] && t>epsilon){
	    changed[k]=1;
	    nactives++;
	}

	tmp = e[l];
	e[l] += t;
	if(changed[l] && l<epsilon){
	    changed[l]=0;
	    nactives--;
	}
	else if(!changed[l] && l>epsilon){
	    changed[l]=1;
	    nactives++;
	}


	/**
	   k       l
	   \  |  |    |  |
	   ____\|1_|____|1_|_____
	   |\ |    |  |
	   k ____|_\|_2__|__|__3__
	   |  |\   |  |
	   |  |  \ |2'|
	   __|__|___\|__|______
	   |  |  2 |\ |
	   l __|__|____|_\|__3___
	   |  |    |  |\
	   |  |    |  |  \

	   2' == 2

	*/

	for(i=0;i<k;i++)rotate(S,c,s,i,k,i,l);	 /* 1 */

	for(i=k+1;i<l;i++)rotate(S,c,s,k,i,i,l); /* 2 */

	for(i=l+1;i<n;i++)rotate(S,c,s,k,i,l,i); /* 3 */


	/* rotate the eigenvectors */
	for(i=0;i<n;i++)rotate(E,c,s,k,i,l,i);

	mir[k]=maxinrow(S,k,n);
	mir[l]=maxinrow(S,l,n);
    }
}


q_t addq(q_t q0, q_t q1)
{
    q_t q;

    q.w = q0.w + q1.w;
    q.v = add(q0.v, q1.v);

    return q;
}


q_t subq(q_t q0, q_t q1){
    q_t q;

    q.w = q0.w - q1.w;
    q.v = sub(q0.v, q1.v);

    return q;
}


q_t crossq(q_t q0, q_t q1){
    q_t q;

    q.w = q0.w * q1.w - dot(q0.v, q1.v);
    q.v = scale(q1.v, q0.w);
    q.v = add(q.v, scale(q0.v, q1.w));
    q.v = add(q.v, cross(q0.v, q1.v));
    return q;
}


float dotq(q_t q0, q_t q1){
    return q0.w*q1.w + dot(q0.v, q1.v);
}


q_t scaleq(q_t q0, float k){
    q_t q;

    q.w = q0.w * k;
    q.v = scale(q0.v, k);

    return q;
}

q_t normalizeq(q_t q0){
    float length;
    q_t q;

    length = sqrtf(q0.w*q0.w + dot(q0.v, q0.v));
    q = scaleq(q0, 1.0/length);

    return q;

}

static inline float dot_inl(const v3_t *a, const v3_t *b)
{
    return a->x*b->x + a->y*b->y + a->z*b->z;
}

static inline v3_t cross_inl(const v3_t *a, const v3_t *b)
{
    v3_t c;
    c.x = a->y*b->z - a->z*b->y;
    c.y = a->z*b->x - a->x*b->z;
    c.z = a->x*b->y - a->y*b->x;
    return c;
}

v3_t rotbyq(q_t q, v3_t v)
{
    q_t qv;
    qv.v = v;
    qv.w = 0;
    qv = crossq(q,qv);
    qv = crossq(qv, normalizeq(conjq(q)));
    return qv.v;
}

// Conjugate of a quaternion, that is: [v, w] => [-v, w]
q_t conjq(q_t q){
    q_t qc;

    qc.v = scale(q.v, -1);
    qc.w = q.w;

    return qc;
}


/*
  row major
*/
void q2mat(q_t q, float m[3][3]){
    float w,x,y,z;

    w = q.w;
    x=q.v.x;
    y=q.v.y;
    z=q.v.z;

    m[0][0] = 1 - 2*y*y - 2*z*z;	m[0][1] = 2*x*y - 2*w*z;	m[0][2] = 2*x*z + 2*w*y;
    m[1][0] = 2*x*y + 2*w*z;	m[1][1] = 1-2*x*x-2*z*z;	m[1][2] = 2*y*z - 2*w*x;
    m[2][0] = 2*x*z - 2*w*y;	m[2][1] = 2*y*z + 2*w*x;	m[2][2] = 1-2*x*x-2*y*y;
}


/* row major */
q_t mat2q(float m[3][3]){
    q_t q;
    float tr = m[0][0] + m[1][1] + m[2][2];
    float s;
    if(tr >= 0){
	s = sqrtf(tr + 1);
	q.w = 0.5 * s;
	s = 0.5 / s;
	q.v.x =(m[2][1] - m[1][2])* s;
	q.v.y =(m[0][2] - m[2][0])* s;
	q.v.z =(m[1][0] - m[0][1])* s;
    }
    else {
	int i = 0;

	if(m[1][1] > m[0][0]){
	    i = 1;
	}
	if(m[2][2] > m[i][i]){
	    i = 2;
	}

	switch(i){
	case 0:
	    s = sqrtf(m[0][0] - m[1][1] - m[2][2] + 1.0f);
	    q.v.x = 0.5f * s;
	    s = 0.5f / s;
	    q.v.y =(m[0][1] + m[1][0])* s;
	    q.v.z =(m[2][0] + m[0][2])* s;
	    q.w   =(m[2][1] - m[1][2])* s;
	    break;
	case 1:
	    s = sqrtf(m[1][1] - m[2][2] - m[0][0] + 1.0f);
	    q.v.y = 0.5f * s;
	    s = 0.5f / s;
	    q.v.z =(m[1][2] + m[2][1])* s;
	    q.v.x =(m[0][1] + m[1][0])* s;
	    q.w   =(m[0][2] - m[2][0])* s;
	    break;
	case 2:
	    s = sqrtf(m[2][2] - m[0][0] - m[1][1] + 1.0f);
	    q.v.z = 0.5f * s;
	    s = 0.5f / s;
	    q.v.x =(m[2][0] + m[0][2])* s;
	    q.v.y =(m[1][2] + m[2][1])* s;
	    q.w   =(m[1][0] - m[0][1])* s;
	    break;
	}
    }

    return q;
}


float angcos(float angle){
    return cosf(ANG_TO_RAD(angle));
}

float angsin(float angle){
    return sinf(ANG_TO_RAD(angle));
}


/* angle axis to quaternion */
q_t aa2q(float angle, v3_t axis){
    q_t q;
    float h;

    h = angle / 2;
    q.w = angcos(h);
    q.v = scale(axis, angsin(h));

    return q;
}


euler_t euler(float p, float y, float r){
    euler_t euler;

    euler.pitch = p;
    euler.yaw = y;
    euler.roll = r;

    return euler;
}

/*
 * world=>obj: yaw(y), pitch(x), roll(z)
 * obj=>world: roll, pitch, yaw.
 * go/strafe/climb offset should be in world space.
 *
 * Converts euler angles to obj=>world matrix(axes of the obj's frame)*/

void euler2axes(euler_t euler, v3_t axis[3]){
    /* 1. rot_around_z(roll degrees);
     * 2. rot_around_x(pitch degrees);
     * 3. rot_around_y(yaw degrees);
     */

    float sr, sp, sy;
    float cr, cp, cy;

    sr = angsin(euler.roll);
    sp = angsin(euler.pitch);
    sy = angsin(euler.yaw);

    cr = angcos(euler.roll);
    cp = angcos(euler.pitch);
    cy = angcos(euler.yaw);

    axis[0].x = sp*sr*sy + cr*cy;
    axis[0].y = sp*cr*sy - sr*cy;
    axis[0].z = cp * sy;

    axis[1].x = cp * sr;
    axis[1].y = cp * cr;
    axis[1].z = -sp;

    axis[2].x = sp*sr*cy - cr*sy;
    axis[2].y = sp*cr*cy + sr*sy;
    axis[2].z = cp * cy;

}

static float det4(const float m[16])
{
    return
	m[12]*m[9]*m[6]*m[3]-
	m[8]*m[13]*m[6]*m[3]-
	m[12]*m[5]*m[10]*m[3]+
	m[4]*m[13]*m[10]*m[3]+
	m[8]*m[5]*m[14]*m[3]-
	m[4]*m[9]*m[14]*m[3]-
	m[12]*m[9]*m[2]*m[7]+
	m[8]*m[13]*m[2]*m[7]+
	m[12]*m[1]*m[10]*m[7]-
	m[0]*m[13]*m[10]*m[7]-
	m[8]*m[1]*m[14]*m[7]+
	m[0]*m[9]*m[14]*m[7]+
	m[12]*m[5]*m[2]*m[11]-
	m[4]*m[13]*m[2]*m[11]-
	m[12]*m[1]*m[6]*m[11]+
	m[0]*m[13]*m[6]*m[11]+
	m[4]*m[1]*m[14]*m[11]-
	m[0]*m[5]*m[14]*m[11]-
	m[8]*m[5]*m[2]*m[15]+
	m[4]*m[9]*m[2]*m[15]+
	m[8]*m[1]*m[6]*m[15]-
	m[0]*m[9]*m[6]*m[15]-
	m[4]*m[1]*m[10]*m[15]+
	m[0]*m[5]*m[10]*m[15];
}

int invmat4(float i[16], const float m[16])
{
    float x=det4(m);
    if (x==0) return 0;

    i[0]= (-m[13]*m[10]*m[7] +m[9]*m[14]*m[7] +m[13]*m[6]*m[11]
	   -m[5]*m[14]*m[11] -m[9]*m[6]*m[15] +m[5]*m[10]*m[15])/x;
    i[4]= ( m[12]*m[10]*m[7] -m[8]*m[14]*m[7] -m[12]*m[6]*m[11]
	    +m[4]*m[14]*m[11] +m[8]*m[6]*m[15] -m[4]*m[10]*m[15])/x;
    i[8]= (-m[12]*m[9]* m[7] +m[8]*m[13]*m[7] +m[12]*m[5]*m[11]
	   -m[4]*m[13]*m[11] -m[8]*m[5]*m[15] +m[4]*m[9]* m[15])/x;
    i[12]=( m[12]*m[9]* m[6] -m[8]*m[13]*m[6] -m[12]*m[5]*m[10]
	    +m[4]*m[13]*m[10] +m[8]*m[5]*m[14] -m[4]*m[9]* m[14])/x;
    i[1]= ( m[13]*m[10]*m[3] -m[9]*m[14]*m[3] -m[13]*m[2]*m[11]
	    +m[1]*m[14]*m[11] +m[9]*m[2]*m[15] -m[1]*m[10]*m[15])/x;
    i[5]= (-m[12]*m[10]*m[3] +m[8]*m[14]*m[3] +m[12]*m[2]*m[11]
	   -m[0]*m[14]*m[11] -m[8]*m[2]*m[15] +m[0]*m[10]*m[15])/x;
    i[9]= ( m[12]*m[9]* m[3] -m[8]*m[13]*m[3] -m[12]*m[1]*m[11]
	    +m[0]*m[13]*m[11] +m[8]*m[1]*m[15] -m[0]*m[9]* m[15])/x;
    i[13]=(-m[12]*m[9]* m[2] +m[8]*m[13]*m[2] +m[12]*m[1]*m[10]
	   -m[0]*m[13]*m[10] -m[8]*m[1]*m[14] +m[0]*m[9]* m[14])/x;
    i[2]= (-m[13]*m[6]* m[3] +m[5]*m[14]*m[3] +m[13]*m[2]*m[7]
	   -m[1]*m[14]*m[7] -m[5]*m[2]*m[15] +m[1]*m[6]* m[15])/x;
    i[6]= ( m[12]*m[6]* m[3] -m[4]*m[14]*m[3] -m[12]*m[2]*m[7]
	    +m[0]*m[14]*m[7] +m[4]*m[2]*m[15] -m[0]*m[6]* m[15])/x;
    i[10]=(-m[12]*m[5]* m[3] +m[4]*m[13]*m[3] +m[12]*m[1]*m[7]
	   -m[0]*m[13]*m[7] -m[4]*m[1]*m[15] +m[0]*m[5]* m[15])/x;
    i[14]=( m[12]*m[5]* m[2] -m[4]*m[13]*m[2] -m[12]*m[1]*m[6]
	    +m[0]*m[13]*m[6] +m[4]*m[1]*m[14] -m[0]*m[5]* m[14])/x;
    i[3]= ( m[9]* m[6]* m[3] -m[5]*m[10]*m[3] -m[9]* m[2]*m[7]
	    +m[1]*m[10]*m[7] +m[5]*m[2]*m[11] -m[1]*m[6]* m[11])/x;
    i[7]= (-m[8]* m[6]* m[3] +m[4]*m[10]*m[3] +m[8]* m[2]*m[7]
	   -m[0]*m[10]*m[7] -m[4]*m[2]*m[11] +m[0]*m[6]* m[11])/x;
    i[11]=( m[8]* m[5]* m[3] -m[4]*m[9]* m[3] -m[8]* m[1]*m[7]
	    +m[0]*m[9]* m[7] +m[4]*m[1]*m[11] -m[0]*m[5]* m[11])/x;
    i[15]=(-m[8]* m[5]* m[2] +m[4]*m[9]* m[2] +m[8]* m[1]*m[6]
	   -m[0]*m[9]* m[6] -m[4]*m[1]*m[10] +m[0]*m[5]* m[10])/x;

    return 1;
}

