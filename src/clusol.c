#include "clusol.h"

// declarations for fortran function calls
void __lusol_MOD_lu1fac(
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

void __lusol_MOD_lu6sol(
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

void __lusol_MOD_lu8rpc(
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

void __lusol_MOD_lu1or2(
  int64_t* n, 
  int64_t* numa,
  int64_t* lena,
  double* a, 
  int64_t* inum,
  int64_t* jnum,
  int64_t* lenc,
  int64_t* locc);

void lu6mul_(
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

void lu8adc_(
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

void lu8adr_(
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

void lu8dlc_(
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

void lu8dlr_(
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

void lu8mod_(
  int64_t* mode,
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

void lu8rpr_(
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

void lu9rpc_(
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

void lu9clr_(
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

void lu9maxs_(
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

void lu9swp_(
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

// c interface function definitions
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
  int64_t* inform) {
  __lusol_MOD_lu1fac(m,n,nelem,lena,ap, aq, frank,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,iploc,iqloc,ipinv,iqinv,w,inform);
}

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
  int64_t* inform) {
  __lusol_MOD_lu6sol(mode,m,n,v,w,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,inform);
}

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
  double* vnorm) {
  __lusol_MOD_lu8rpc(mode1,mode2,m,n,jrep,v,w,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,inform,diag,vnorm);
}

void clu1or2(
  int64_t* n, 
  int64_t* numa,
  int64_t* lena,
  double* a, 
  int64_t* inum,
  int64_t* jnum,
  int64_t* lenc,
  int64_t* locc) {
  __lusol_MOD_lu1or2(n,numa,lena,a,inum,jnum,lenc,locc);
}

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
  int64_t* locr) {
  lu6mul_(mode,m,n,v,w,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr);
}

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
  double* vnorm) {
  lu8adc_(mode,m,n,v,w,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,inform,diag,vnorm);
}

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
  double* diag) {
  lu8adr_(m,n,w,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,inform,diag);
}

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
  int64_t* inform) {
  lu8dlc_(m,n,jdel,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,inform);
}

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
  int64_t* inform) {
  lu8dlr_(mode,m,n,idel,v,w,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,inform);
}

void clu9diagu(
  int64_t* m,
  int64_t* n, 
  int64_t* lena,
  int64_t* luparm,
  double* diagU,
  double* a, 
  int64_t* indc,
  int64_t* indr,
  int64_t* ip,
  int64_t* iq,
  int64_t* lenr,
  int64_t* locr,
  int64_t* inform) {
  lu9diagu_(m,n,lena,luparm,diagU,a,indc,indr,ip,iq,lenr,locr,inform);
}

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
  int64_t* inform) {
  lu9mod_(m,n,beta,v,w,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,inform);
}

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
  int64_t* inform) {
  lu8rpr_(mode1,mode2,m,n,irep,v,w,wnew,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,inform);
}

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
  int64_t* inform) {
  lu9rpc_(mode1,mode2,m,n,irep,v,w,wnew,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,inform);
}

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
  int64_t* inform) {
  lu9clr_(m, n, lena, luparm, parmlu, a, indc, indr, p, q, lenr, locc, locr, c, inform); 
}


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
  ) {
  lu9maxs_(m, n, nrank, lenlu, luparm, parmlu, lua, luindc, luindr, lup, luq, lulenc, lulenr, lulocc, lulocr, annz, av, ai, aj, u, uS, v, w, s_r, s_c, alpha);
}

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
  int64_t* inform){
    lu9swp_(m,n,a_r,a_c,s_r,s_c,annz,av,ai,aj,v1,c1, w1,v2,c2,w2,lena,luparm,parmlu,a,indc,indr,ip,iq,ap,aq,lenc,lenr,locc,locr,inform);
  }


