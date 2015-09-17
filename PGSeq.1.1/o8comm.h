/* **************************************************************************** */
/*                                  o8comm.h                                    */
/* **************************************************************************** */
                                      
#include "o8fuco.h"
    
X RREAL      runtim;
X RREAL      optite;
X DDOUBLE    **accinf;
X IINTEGER   itstep,phase;

X DDOUBLE    upsi,upsi0,upsi1,upsist,psi,psi0,
            psi1,psist,psimin,
            phi,phi0,phi1,phimin,fx,fx0,fx1,
            fxst,donlp2_fmin,b2n,b2n0,xnorm,x0norm,sig0,dscal,ddnorm,d0norm;
X DDOUBLE    sig,sigmin,dirder,cosphi,upsim;
X DDOUBLE    *donlp2_x,*x0,*x1,*xmin,*d,*d0,
            *dd,*difx,*resmin;

X DDOUBLE    *gradf,gfn,*qgf,*gphi0,*gphi1,
            **gres,*gresn;

X IINTEGER   *perm,*perm1,*colno,rank;
X DDOUBLE    **qr,*betaq,*diag,
            *cscal,
            *colle;

/* colno also used o8qpso with double length ! */

X DDOUBLE    **a,scalm,scalm2,*diag0,matsc;

X IINTEGER   *violis,*alist,*o8bind,
            *o8bind0, *aalist,*clist;
                        
X DDOUBLE    *u,*u0,
            *w,*w1,*res,
            *res0,*res1,
            *resst,scf,scf0,
            *yu,*slack,infeas,*work;

X IINTEGER   iterma;
X DDOUBLE    del,del0,del01,delmin,tau0,tau,ny;
X DDOUBLE    smalld,smallw,rho,rho1,eta,epsx,c1d,
            scfmax,updmy0,tauqp,taufac,taumax;

X DDOUBLE    alpha,bbeta,theta,sigsm,sigla,delta,stptrm;
X DDOUBLE    delta1,stmaxl;

X DDOUBLE    level;
X IINTEGER   clow,lastdw,lastup,lastch;

X FILE      *prou,*meu;

X DDOUBLE    *ug,*og;

X DDOUBLE    *low,*up,big;

X IINTEGER   nreset;

X DDOUBLE    *xst;
