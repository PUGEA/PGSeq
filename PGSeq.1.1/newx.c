#include "o8para.h"
/**********************************************************************/
/***    enable user inspection of iteration data                *******/
/**********************************************************************/
void newx ( DDOUBLE x[] , DDOUBLE u[] , IINTEGER itstep , DDOUBLE **accinf, 
            LLOGICAL *cont )
{
/*   x = current solution                                                 */
/*   u = current multiplier estimates                                     */
/*   itstep = last step performed                                         */
/*   accinf : accumulated short information   defined from accinf[1][] to */
/*            accinf[itstep][]                                            */
/*   continue must be set TTRUE . otherwise donlp2 will terminate          */
/*   irrespective of the state of the optimization                        */
/*   for description of accinf see below                                  */
/*  accinf : a c c u m u l a t e d   i n f o r m a t i o n                */
/*  on iteration sequence : components 0--32 for one step                 */
/*  0: not used                                                           */
/*  1: step-nr                                                            */
/*  2: f(x-k) current value of objective (zero in feasibility improvement */
/*            phase (-1)  )                                               */
/*  3: scf    internal scaling of objective (zero in phase -1)            */
/*  4: psi    the weighted penalty-term                                   */
/*  5: upsi   the unweighted penalty-term (l1-norm of constraint vector)  */
/*  6: del_k_1 bound for currently active constraints                     */
/*  7: b2n0   weigthed l2-norm of projected gradient, based on inequality */
/*            constraints <= delmin and equality constraints              */
/*  8: b2n    l2-norm of projected gradient based on del_k_1              */
/*  9: nr     number of binding constraints                               */
/*  10: sing  if 1, the binding constraints don't satisfy the regularity  */
/*            condition                                                   */
/*  11: umin   infinity norm of negative part of multiplier               */
/*  12: presently not used                                                */
/*  13: cond_r condition number of diagonal part of qr-decomposition      */
/*             of normalized gradients of binding constraints             */
/*  14: cond_h condition number of diagonal of cholesky-factor            */
/*             of updated full hessian                                    */
/*  15: scf0   the relative damping of tangential component if upsi>tau0/2*/
/*  16: xnorm  l2-norm of x                                               */
/*  17: ddnorm  l2-norm of d (correction from qp -subproblem, unscaled)    */
/*  18: phase  -1:infeasibility improvement phase, 0:initial optimization */
/*            1: binding constraints unchanged ,                          */
/*            2: d small, maratos correction in use                       */
/*  19: c_k    number of decreases of penalty weights                     */
/*  20: wmax   infinity norm of current penalty weights                   */
/*  21: sig_k  stepsize from unidimensional minimization (backtracking)   */
/*  22: cfincr number of objective evaluations for stepsize-algorithm     */
/*  23: dirder directional derivative of penalty-function along d (scaled)*/
/*  24: dscal  scaling factor for d                                       */
/*  25: cosphi cos of arc between d and d_previous. if cosphi>=theta      */
/*             then  sig larger  than     one (up to sigla) is tried      */
/*  26: violis[0] number of constraints not binding at x but hit during   */
/*             line search                                                */
/*  27:        type of update for hessian:                                */
/*             1: normal p&m-bfgs-update,                                 */
/*             0: update suppress$                                        */
/*            -1 restart with scaled unit matrix                          */
/*             2 standard bfgs (unconstrained)                            */
/*             3 bfgs modified by powells device                          */
/*  28: ny_k resp. t_k  modification factor for damping the projector in  */
/*             the bfgs  resp. the pantoja-mayne update                   */
/*  29: 1-my_k/xsi_k modification factor for damping the                  */
/*             quasi-newton-condition in bfgs resp the regularization term*/
/*            for unmodified bfgs ny_k should be larger than updmy0       */
/*            (near one)                                                  */
/*            and 1-my_k equal one.                                       */
/*  30: qpterm 0, if sing=-1                                              */
/*                termination indicator of qp-solver otherwise            */
/*             1: successful,                                             */
/*             -1: tau becomes larger than tauqp without slack-           */
/*                 variables becoming sufficiently small .                */
/*             -2: infeasible qp-problem (theoretically impossible)       */
/*  31: tauqp: weight of slack-variables in qp-solver                     */
/*  32: infeas l1-norm of slack-variables in qp-solver                    */


     *cont = TTRUE ;
}
 
 

