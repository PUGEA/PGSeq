/* **************************************************************************** */
/*                                  o8fuco.h                                    */
/* **************************************************************************** */
    
X LLOGICAL   *val,*llow,*lup;

X IINTEGER   n,nr,nres,nlin,nonlin, nstep, ndualm, mdualm;

X DDOUBLE    epsmac,tolmac,deldif;

X char      name[41];

X DDOUBLE    epsdif;

X LLOGICAL   intakt,te0,te1,te2,te3,singul;
X LLOGICAL   ident,eqres,silent,analyt,cold;

X IINTEGER   icf,icgf,cfincr,*cres,*cgres;

X LLOGICAL   ffuerr,*confuerr;

/*  special variables for the interface to old fashioned function */
/*  specification                                                 */
/*  can be removed is only problems in the new formulation are to be used */
X IINTEGER nh,ng;

X IINTEGER *nonlinlist;

X IINTEGER **gunit; 

X LLOGICAL *gconst; 

X LLOGICAL *cfuerr; 
/* this is necessary because of the different index positions used */
/*  in the old and new versions   */
