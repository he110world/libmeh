#ifndef MEH_MATH_H
#define MEH_MATH_H

#define PI 3.1415926535897932384626

typedef struct{
	float x,y;
} v2_t;

typedef struct {
	float x;
	float y;
	float z;
} v3_t;

typedef struct {
	float x;
	float y;
	float z;
	float w;
} v4_t;


typedef float **mat_t;

enum {X, Y, Z};
typedef struct {
	float _[6];
} v3__t;

typedef struct {
	float w;
	v3_t v;
} q_t;


typedef struct {
	float pitch;
	float yaw;
	float roll;
} euler_t;


typedef struct {
	v3_t n;
	float d;
} plane_t;

#if 0
typedef union{
	v3_t row[3];
	float _[3][3];
} mat3x3_t;
#endif

typedef struct{
	float			data[3][3];
} mat3x3_t;


typedef float m3_t[9];
typedef float m4_t[16];

extern v3_t add (v3_t v0, v3_t v1);
extern v3_t sub (v3_t v0, v3_t v1);
extern float dot (v3_t v0, v3_t v1);
extern v3_t cross (v3_t v0, v3_t v1);
extern v3_t normalize (v3_t v);
extern v3_t scale (v3_t v, float s);
extern v3_t vec (float x, float y, float z);
extern v3_t matvec (m3_t m, v3_t v);
extern void cpmat (m3_t dest, m3_t src);
extern void matt(float dest[3][3], float src[3][3]);
extern void matmat(float dest[3][3], float m[3][3], float r[3][3]);
extern float length (v3_t v);
extern int ludecomp (float A[], float LU[], int P[], int N);
extern int qrdecomp (float A[], float Q[], float R[], int N);
extern int lusolve (float A0[], float x[], float b0[], int N);
extern int luinv (float A0[], float Ainv[], int N);
extern int eigenvalue (float A0[], float lambda[], int N);
extern int eigensym (float A0[], float lambda[], float U[], int N);
extern void printmat (char *name, float A[], int N);
extern void printvec (char *name, float v[], int N);
extern void orthonormalization3x3(float m[3][3]);
extern v3_t mat3x3_v3(float m[3][3], v3_t v);
extern v3_t matv3(float m[3][3], v3_t v);
extern v3_t mattv3(float m[3][3], v3_t v);
extern q_t addq (q_t q0, q_t q1);
extern q_t subq (q_t q0, q_t q1);
extern q_t crossq (q_t q0, q_t q1);
extern float dotq (q_t q0, q_t q1);
extern q_t scaleq(q_t q0, float k);
extern q_t normalizeq (q_t q0);
extern q_t conjq (q_t q);
extern void q2mat(q_t q, float m[3][3]);
extern q_t mat2q(float m[3][3]);
extern v3_t rotbyq(q_t q, v3_t v);
extern void euler2axes (euler_t euler, v3_t axis[3]);
extern euler_t euler (float p, float y, float r);
extern int invmat4(float i[16], const float m[16]);

#define radsin(rad) sinf (rad)
#define radcos(rad) cosf (rad)

#define ANG_TO_RAD(ang) (ang) * 0.01745329252
#define RAD_TO_ANG(rad) (rad) * 57.29577951

#define DEG2RAD( d ) ( ( d ) *PI / 180.0f )
#define RAD2DEG( r ) ( ( r ) * 180.0f / PI )
#define SIND( d ) sinf( DEG2RAD( (d) ) )
#define COSD( d ) cosf( DEG2RAD( (d) ) )
#define TAND( d ) tanf( DEG2RAD( (d) ) )

#define MATRIX(m, w, h, i, j) m[(j)*(h) + (i)]
#define MATRIX_SAFE(m, w, h, i, j) m[(j)%(h)*(h) + (i)%(w)]

#define EPSILON 0.00001

#define max(a, b) ((a)>(b) ? (a) : (b))
#define min(a, b) ((a)<(b) ? (a) : (b))

#define WRAP(f, mini, maxi) ((f) >= (mini) ? ((f) < (maxi) ? (f) : (f) - (maxi) + (mini)) : (f) + (maxi) - (mini))
#define CLAMP(f, mini, maxi) ((f) < (mini) ? (mini) : ((f) > (maxi) ? (maxi) : (f))

#define v3_to_array (v) ((float *) & (v))
#define ARRAY(v) ((float *)&(v))

#define matv matv3
#define proj dot


#endif
