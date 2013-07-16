#ifndef BACKEND__H
#define BACKEND__H



// Memory management routines


extern void cvar_set(const char *nam, const char *str);
extern char *cvar_get(const char *nam, int *ok);
extern void cvar_seti(const char *nam, int i);
extern int cvar_geti(const char *nam, int *ok);
extern void cvar_setf(const char *nam, float f);
extern float cvar_getf(const char *nam, int *ok);
extern void cvar_del(char *nam);

extern float *cvec_get(char *nam);  // don't use sizeof() over cvec!


extern int fnt_load(const char *nam);
extern void fnt_loadslot(const char *nam, int i);
extern void fnt_use(const char *nam);
extern void fnt_useslot(int i);
extern void fnt_color3f(float r, float g, float b);
extern void fnt_color4f(float r, float g, float b, float a);
extern void fnt_color3fv(float c[3]);
extern void fnt_color4fv(float c[4]);
extern void fnt_getcolor(float c[4]);
extern void fnt_scale(float s);
extern void fnt_printf(int posx, int posy, int alignment, const char *fmt, ...);

extern void cs_printf(const char *fmt, ...);

extern void cmd_addtxt(char *t);
extern int cmd_geti(int *i);
extern int cmd_getf(float *f);
extern int cmd_gets(char *s);
extern void cmd_execnow(char *cmd);
extern void cmd_exec();
extern void cmd_add(const char *nam, void (*func)());

extern void vu_use(const char *nam);
extern void vu_force_bgn(char *nam);
extern void vu_end();

extern void vu_bgn(char *nam);
extern void vu_kcmd(int k, char *cmd);
extern void vu_cmd(char *cmd);




typedef struct cvec_s{
	char *nam;
	float v[16];
	int id;
} cvec_t;

#define cvec(aa) \
	static cvec_t cvec__N_a_M_e_##aa __attribute__((aligned (4) )) ={ .nam= #aa }; \
	static float *aa __attribute__((aligned (4) ))=cvec__N_a_M_e_##aa.v;

#define cvec_end(filename) \
	static cvec_t cvec_end_##filename __attribute__((aligned (4) ))= {.nam= 0 }; \
 	void cvec_init_##filename(cvec_t **first, int *n) {	\
		cvec_t *begin=&cvec_bgn_##filename; \
		cvec_t *end=&cvec_end_##filename;\
		*first=(float**)(begin+1)+1;				\
		*n = ((int)end- (int)(*first))/(sizeof(cvec_t)+sizeof(float*))   ; \
	}

#define cvec_bgn(filename) \
	static cvec_t cvec_bgn_##filename __attribute__((aligned (4) ))={ .nam= 0, .id=__LINE__} ; \
	static float *cvec_placeholder_##filename __attribute__((aligned (4) ))=cvec_bgn_##filename.v;


enum {TOK_IF = 256, TOK_ELSE, TOK_WHILE, TOK_BREAK, TOK_CONTINUE, TOK_RETURN, TOK_INT,
      TOK_FLOAT, TOK_FOR, TOK_VOID, TOK_GE, TOK_LE, TOK_EQ,
      TOK_NEQ, TOK_INT_CONST, TOK_FLOAT_CONST, TOK_STRING_CONST, TOK_END, TOK_ID,
      TOK_L_OR, TOK_L_AND, TOK_MAX};

#ifdef USE_MEH_TYPE
typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned char uchar;
#endif

extern void cam3_update(char *nam, float du, float dt);
extern void cam3_movetg(char *nam, float x, float y, float z);
extern void cam3_use(char *nam);
extern void cam3_add(char *nam, v3_t target, v3_t pos, float Kd, float angKd);

#endif
