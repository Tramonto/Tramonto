/* ---------------------------------------------------------
Calculate the thermodynamic properties for attractions in the bulk fluid.
------------------------------------------------------------*/
#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"
/*************************************************************
calc_att_properties: In this routine calculate the strict mean field
                     attraction contribution to the pressure and
                     chemical potential */
double calc_att_properties(double *betamu_att, double *rho)
{
  int icomp,jcomp;
  double betap_att;

  betap_att = 0.0; 
  for (icomp=0; icomp<Ncomp; icomp++){
      betamu_att[icomp] = 0.0;
      for (jcomp=0; jcomp<Ncomp; jcomp++) Avdw[icomp][jcomp]=0.0;
  }

  for (icomp=0; icomp<Ncomp; icomp++) {
     for (jcomp=0; jcomp<Ncomp;jcomp++){
       int_stencil_bulk(U_ATTRACT,icomp,jcomp,NULL);
       Avdw[icomp][jcomp] = Temporary_sum;
       betamu_att[icomp] += rho[jcomp]*Avdw[icomp][jcomp];
       betap_att += 0.5*Avdw[icomp][jcomp]*rho[icomp]*rho[jcomp];
     }
  }
  return(betap_att);
}
/*************************************************************************
dp_drho_att: the derivative of the attractive part of the pressure 
             with respect to rho*/
double dp_drho_att(double *rho)
{
 double dp_drho;

 dp_drho = Avdw[0][0]*rho[0];
 return (dp_drho);
}
/*************************************************************************
dmu_drho_att: the derivative of the attractive part of the chemical potential 
               with respect to rho*/
double dmu_drho_att(double *rho)
{
   double dmu_drho;

   dmu_drho = Avdw[0][0];
   return (dmu_drho);
}
/******************************************************************************/
