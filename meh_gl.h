#ifndef MEH_GL_H
#define MEH_GL_H

// meh graphics lib

extern int vfs_load(const char *nam);
extern void vfs_use(int id);
extern void vfs_reload(const char *nam);

extern void uni1f(const char *nam, float v0);
extern void uni2f(const char *nam, float v0, float v1);
extern void uni3f(const char *nam, float v0, float v1, float v2);
extern void uni4f(const char *nam, float v0, float v1, float v2, float v3);
extern void unifv(const char *nam, const float *dat);
extern void samp2d(const char *nam, unsigned tex);

extern int vb_alloc(int n, int sz, void *dat);
extern int eb_alloc(int n, int elemsz, void *dat);
extern void vb_free(int vb);
extern void eb_free(int eb);
extern void vb_update(int vb, int n, int sz, void *dat);
extern void eb_update(int eb, int nbytes, void *dat);
extern int vb_getcap(int vb);

extern void rend(const char *fmt, int eb, ...);
extern void rendq(const char *fmt, int eb, ...);
extern void rendrect(const char *fmt, ...);
extern void rnd_tg(const char *fmt, ...);
extern void rnd_scrncoord();
extern void rnd_scrncoordsz(int w, int h);
extern void rnd_quad(float rgba[4], float qbnd[4]);
extern void rnd_quadsz(int w, int h, float rgba[4], float qbnd[4]);
extern void rnd_texquad(unsigned tex, float tcbnd[4], float qbnd[4]);
extern void rnd_texquadsz(unsigned tex, int w, int h, float tcbnd[4], float qbnd[4]);

enum {RT_4B, RT_4H, RT_DEPTH, RT_MAX};
extern unsigned rt_alloc(int type, int level);
extern void rt_free(unsigned tex);

extern unsigned tex1b(int w, int h, int filter, void *data);
extern unsigned tex3b(int w, int h, int filter, void *data);
extern unsigned tex4b(int w, int h, int filter, void *data);
extern unsigned tex1h(int w, int h, int filter, void *data);
extern unsigned tex3h(int w, int h, int filter, void *data);
extern unsigned tex4h(int w, int h, int filter, void *data);
extern unsigned tex1f(int w, int h, int filter, void *data);
extern unsigned tex3f(int w, int h, int filter, void *data);
extern unsigned tex4f(int w, int h, int filter, void *data);
extern unsigned tex3bdxt1(int w, int h, int filter, void *data);
extern unsigned dtex(int w, int h, int filter, void *data);
extern unsigned dtex2b(int w, int h, int filter, void *data);
extern unsigned dtex3b(int w, int h, int filter, void *data);
extern unsigned dstex(int w, int h, int filter, void *data);


#endif
