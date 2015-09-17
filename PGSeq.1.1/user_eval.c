#include "o8para.h"              

void user_eval(DDOUBLE xvar[],IINTEGER mode) {
/* ************************************************************************ */
/* called only if bloc = TTRUE                                               */
/* interface to a users evaluation code of functions                        */
/* which define the optimization problem                                    */
/* mode = -1 : evaluate function values at xvar only, store them in fu      */
/* mode = 1 : evaluate f , grad(f), con[], grad(con[]) and store            */
/*             them in  fu and fugrad                                       */
/* mode = 2  : evaluate gradients only at xvar and store them in fugrad     */
/* The users evaluation routine                                             */
/*             eval_extern(mode)                                            */
/* is called with two settings of its mode only:                            */
/* on one variable only (gunit[1][i] = 1).                                  */
/* mode = 1: set function value                                             */
/* fu using xtr, mode = 2 set function values fu and gradients fugrad using */
/* xtr                                                                      */
/* Since donlp2 uses internally a scaled value of donlp2_x, the external value xtr */
/* is obtained by multiplying with xsc. In evaluating the gradients,        */
/* the user must not take care of this scaling. Scaling of the gradients    */
/* is done by donlp2 internally.                                            */
/* that means xtr and the function values returned below to the unscaled    */
/* original problem                                                         */
/* If the user has set analyt = TTRUE, he is responsible for computing       */
/* the gradients himself correctly                                          */
/* otherwise, if analyt = FFALSE, then this routine does numerical           */
/* differentiation according to the method defined by difftype:             */
/* difftype = 1   ordinary forward difference quotient, using               */
/*                min(0.1*sqrt(epsfcn),1.e-2) as stepsize                   */
/*                epsfcn is the assumed precision of function values        */
/*                and must be set by the user in setup0                     */
/* difftype = 2   central difference quotient with stepsize                 */
/*                min(0.1*pow(epsfcn,2/3),1.e-2) as stepsize                */
/* difftype = 3   6-th order Richardson extrapolation method using          */
/*                6 function values for each component of grad,             */
/*                with stepsize min(0.1*pow(epsfcn,1/7),0.01)               */
/* This is done by several calls of eval_extern                             */
/* xtr is the current design to be used by the external routine             */
/* Values must be returned in the arrays fu (and fugrad).                   */
/*  The order of information is as follows:                                 */
/* fu[0]  objective function value                                          */
/* fugrad[i][0] : objective function gradient, i = 1,...,n the components   */
/* fu[1],...,fu[nonlin]  the values of c  (if any)                          */
/* fugrad[i][1],.., fugrad[i][nonlin] gradients of  constraint functions c  */
/* i=1,..,n                                                                 */
/* bound constraints and linear constraints must not be evaluated here      */
/* The method uses the parameters                                           */
/*      epsfcn  =  relative accuracy of function values                     */
/*      taubnd  =  amount by which bound constraints may be violated        */
/*                 during finite differencing                               */
/*      difftype   (see above)                                              */
/*      analyt     (see above)                                              */
/* These must be set prior in user_init                                     */
/* ************************************************************************ */
    
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #undef   X
    #include "o8cons.h"
    
    void    eval_extern(IINTEGER mode);
    DDOUBLE* d1_malloc(IINTEGER size1, IINTEGER init); 

    static DDOUBLE   *fusave,xhelp,xincr,sd1,sd2,sd3,d1,d2,d3;
    static IINTEGER  i,j;
    static LLOGICAL  changed;
    static IINTEGER  firstTime = 1;
    
    if (firstTime) {
        firstTime = 0;
	fusave = d1_malloc(nlin+nonlin+1, 1);
    }	
	
    
    if ( mode < -1 || mode  > 2 || mode == 0 ) { 
        fprintf(stderr,"donlp2: call of user_eval with undefined mode\n");
        exit(1);
    }
    if ( mode == -1 ) {
        changed = FFALSE;
        for (i = 1 ; i <= n ; i++) {
            if ( xtr[i] != xsc[i]*xvar[i] ) changed = TTRUE;
        }
        if ( changed || ! valid ) {
            for (i = 1 ; i <= n ; i++) {
                xtr[i] = xsc[i]*xvar[i];
            }
            eval_extern(1);
            
            valid = TTRUE;
        }
        return;
    }
    if ( mode >= 1 ) {
    
        /* evaluate functions and its gradients */
        
        for (i = 1 ; i <= n ; i++) {
            xtr[i] = xsc[i]*xvar[i];
        }
        /* even if xtr did not change we must evaluate externally, since now */
        /* gradients are required                                            */
        
        if ( analyt ) {
        
            eval_extern(2);
            
            valid = TTRUE;
            
            return;
            
        } else {
            if ( mode == 1 ) eval_extern(1);
        }
        /* continue with computing the gradients : mode = 1 or mode = 2 */
        /* require gradients                                            */
        
        if ( difftype == 1 ) {
            deldif = min(tm1*sqrt(epsfcn),tm2);
            
            if ( mode == 2 ) eval_extern(1);
            
            for (j = 0 ; j <= nonlin ; j++) {
                fusave[j] = fu[j];
            }
            for (i = 1 ; i <= n ; i++) {
                xhelp = xtr[i];
                xincr = min(min(tm2,deldif*fabs(xhelp)+deldif),taubnd);
                if ( xhelp >= zero ) {
                    xtr[i] = xhelp+xincr;
                } else {
                    xtr[i] = xhelp-xincr;
                }
                eval_extern(1);
                
                for (j = 0 ; j <= nonlin ; j++) {
                    fugrad[i][j] = (fu[j]-fusave[j])/(xtr[i]-xhelp);
                }
                xtr[i] = xhelp;
            }
            for (j = 0 ; j <= nonlin ; j++) {
                fu[j] = fusave[j];
            }
            valid = TTRUE;
            
            return;
            
        } else if ( difftype == 2 ) {
            deldif = min(tm1*pow(epsfcn,one/three),tm2);
            for (j = 0 ; j <= nonlin ; j++) {
                fusave[j] = fu[j];
            }
            for (i = 1 ; i <= n ; i++) {
                xhelp  = xtr[i];
                xincr  = min(min(tm2,deldif*fabs(xhelp)+deldif),taubnd);
                xtr[i] = xhelp+xincr;
                
                eval_extern(1);
                
                for (j = 0 ; j <= nonlin ; j++) {
                    fud[j][1] = fu[j];
                }
                xtr[i] = xhelp-xincr;
                
                eval_extern(1);
                
                for (j = 0 ; j <= nonlin ; j++) {
                    fud[j][2] = fu[j];
                }
                for (j = 0 ; j <= nonlin ; j++) {
                    fugrad[i][j] = (fud[j][1]-fud[j][2])/(xincr+xincr);
                }
                xtr[i] = xhelp;
            }
            for (j = 0 ; j <= nonlin ; j++) {
                fu[j] = fusave[j];
            }
            valid = TTRUE;
            
            return;
            
        } else {
            deldif = min(tm1*pow(epsfcn,one/seven),tm2);
            for (j = 0 ; j <= nonlin ; j++) {
                fusave[j] = fu[j];
            }
            for (i = 1 ; i <= n ; i++) {
                xhelp  = xtr[i];
                xincr  = min(min(tm2,deldif*fabs(xhelp)+deldif),taubnd/four);
                xtr[i] = xhelp-xincr;
                
                eval_extern(1);
        
                for (j = 0 ; j <= nonlin ; j++) {
                    fud[j][1] = fu[j];
                }
                xtr[i] = xhelp+xincr;
                
                eval_extern(1);
                
                for (j = 0 ; j <= nonlin ; j++) {
                    fud[j][2] = fu[j];
                }
                xincr  = xincr+xincr;
                d1     = xincr;
                xtr[i] = xhelp-xincr;
                
                eval_extern(1);
                
                for (j = 0 ; j <= nonlin ; j++) {
                    fud[j][3] = fu[j];
                }
                xtr[i] = xhelp+xincr;
                
                eval_extern(1);
                
                for (j = 0 ; j <= nonlin ; j++) {
                    fud[j][4] = fu[j];
                }
                xincr  = xincr+xincr;
                d2     = xincr;
                xtr[i] = xhelp-xincr;
                
                eval_extern(1);
                
                for (j = 0 ; j <= nonlin ; j++) {
                    fud[j][5] = fu[j];
                }
                xtr[i] = xhelp+xincr;
                
                eval_extern(1);
                
                for (j = 0 ; j <= nonlin ; j++) {
                    fud[j][6] = fu[j];
                }
                xtr[i] = xhelp;
                d3     = xincr+xincr;
                for (j = 0 ; j <= nonlin ; j++) {
                    sd1 = (fud[j][2]-fud[j][1])/d1;
                    sd2 = (fud[j][4]-fud[j][3])/d2;
                    sd3 = (fud[j][6]-fud[j][5])/d3;
                    sd3 = sd2-sd3;
                    sd2 = sd1-sd2;
                    sd3 = sd2-sd3;
                    fugrad[i][j] = sd1+p4*sd2+sd3/c45;
                }
            }
            for (j = 0 ; j <= nonlin ; j++) {
                fu[j] = fusave[j];
            }
            valid = TTRUE;
            
            return;
        }
    }
}


