#ifndef CLUSOL_H_
#define CLUSOL_H_

#include <stdint.h>

void clu1fac(
  int64_t* m,
  int64_t* n,
  int64_t* nelem,
  int64_t* lena,
  int64_t* ap,
  int64_t* aq, 
  int64_t* frank,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* iploc,
  int64_t* iqloc,
  int64_t* ipinv,
  int64_t* iqinv,
  double* w,
  int64_t* inform);

void clu6sol(
  int64_t* mode,
  int64_t* m,
  int64_t* n,
  double* v,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform);

void clu8rpc(
  int64_t* mode1,
  int64_t* mode2,
  int64_t* m,
  int64_t* n,
  int64_t* jrep,
  double* v,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform,
  double* diag,
  double* vnorm);

void clu6mul(
  int64_t* mode,
  int64_t* m,
  int64_t* n,
  double* v,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr);

void clu8adc(
  int64_t* mode,
  int64_t* m,
  int64_t* n,
  double* v,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform,
  double* diag,
  double* vnorm);

void clu8adr(
  int64_t* m,
  int64_t* n,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform,
  double* diag);

void clu8dlc(
  int64_t* m,
  int64_t* n,
  int64_t* jdel,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform);

void clu8dlr(
  int64_t* mode,
  int64_t* m,
  int64_t* n,
  int64_t* idel,
  double* v,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform);

void clu9mod(
  int64_t* m,
  int64_t* n,
  double* beta,
  double* v,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc, 
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform);

void clu8rpr(
  int64_t* mode1,
  int64_t* mode2,
  int64_t* m,
  int64_t* n,
  int64_t* irep,
  double* v,
  double* w,
  double* wnew,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform);

void clu9rpc(
  int64_t* mode1,
  int64_t* mode2,
  int64_t* m,
  int64_t* n,
  int64_t* irep,
  double* v,
  double* w,
  double* wnew,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform);

void clu9clr(
  int64_t* m,
  int64_t* n,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr, 
  int64_t* p, 
  int64_t* q, 
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  double* c, 
  int64_t* inform); 

void clu1or2(
  int64_t* n, 
  int64_t* numa,
  int64_t* lena,
  double* a, 
  int64_t* inum,
  int64_t* jnum,
  int64_t* lenc,
  int64_t* locc);

void clu9maxs(
  int64_t* m,
  int64_t* n,
  int64_t* nrank,
  int64_t* lenlu,
  int64_t* luparm,
  double* parmlu,
  double* lua,
  int64_t* luindc,
  int64_t* luindr,
  int64_t* lup,
  int64_t* luq,
  int64_t* lulenc,
  int64_t* lulenr,
  int64_t* lulocc,
  int64_t* lulocr,
  int64_t* annz,
  double* av,
  int64_t* ai,
  int64_t* aj,
  double* u,
  double* uS,
  double* v,
  double* w,
  int64_t* s_r,
  int64_t* s_c, 
  double* alpha
  );

void clu9swp(
  int64_t* m,
  int64_t* n,
  int64_t* a_r,
  int64_t* a_c,
  int64_t* s_r,
  int64_t* s_c,
  int64_t* annz,
  double* av,
  int64_t* ai,
  int64_t* aj,
  double* v1,
  double* c1,
  double* w1,
  double* v2,
  double* c2, 
  double* w2,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* ip,
  int64_t* iq,
  int64_t* ap, 
  int64_t* aq, 
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform);
#endif // CLUSOL_H_
