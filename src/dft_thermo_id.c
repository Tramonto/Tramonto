/* ---------------------------------------------------------
Calculate the thermodynamic properties of an ideal gas fluid.
------------------------------------------------------------*/
#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"

/******************************************************************************/
/* calc_ideal_gas: This routine computes ideal gas properties at a density of interest */
double calc_ideal_gas(double *rho,double *betamu)
{
   double betap=0.0;
   int i;
   
   for (i=0;i<Ncomp;i++){
       betamu[i] = log(rho[i]);
                  /*           - 3.0*log(Sigma_ff[icomp][icomp]) -
                               1.5*log(Mass[icomp]*Temp); */
       betap+=rho[i];
   }
   return(betap);
}
/******************************************************************************/
