/*
//@HEADER
// ********************************************************************
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/* ---------------------------------------------------------
Calculate the relevant thermodynamic properties of the bulk fluid.
We use the Percus-Yevick compressibilty equation of state.

Input:  bulk densities of each component: Rho_b[icomp]

Output:	1) hard sphere pressure:  p sigma_ff[1]^3 / kT
        2) excess hard sphere chemical potentials: mu_hs_ex/kT
        3)       hard sphere chemical potentials: mu_hs/kT
------------------------------------------------------------*/
#include "dft_thermo.h"
/*#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"*/

/*double   coexistence();*/
/*******************************************************************************/
/* left over coexistence logic in the main thermodynamics routine..... note that 
   the coexistence routines need to be completely rewritten - if possible before release 2.1.
   note that this code is not currently compiled or called in the source code */
{

  /*   if (Ncomp == 1 && Ipot_ff_n==2 && Iliq_vap != -2){
      p_coex = coexistence();
      if (Proc==0 &&Iwrite!=NO_SCREEN) {
         if (print_flag) {
           printf("For kT/epsilon=:%9.6f  the coexistence densities are:\n",1.0/Eps_ff[0][0]);
           printf("\t\t\tliquid: %9.6f  vapor: %9.6f\n",Rho_coex[0],Rho_coex[1]);
           printf("and the pressure at coexistence is:%9.6f\n",p_coex);
           printf("\n\n your density is %9.6f\n",Rho_b[0]);
         }
         else {
           printf("\n\tThermodynamics called for For kT/epsilon=:%9.6f and bulk density %9.6f\n",
                   1.0/Eps_ff[0][0], Rho_b[0]);
         }
      }

      if (Iliq_vap == 1 ){
         if (Proc==0 && print_flag&&Iwrite!=NO_SCREEN) printf("\n RESETTING RHO_BULK = RHO_VAPOR\n");
         Rho_b[0] = Rho_coex[1];
      }
      else if (Iliq_vap == 2 ||Iliq_vap==3){
         if (Proc==0 && print_flag&&Iwrite!=NO_SCREEN) printf("\n RESETTING RHO_BULK = RHO_LIQUID\n");
         Rho_b[0] = Rho_coex[0];
      }
      if (Rho_b[0] < Rho_coex[0] && Rho_b[0] > Rho_coex[1] &&
          Iliq_vap == -1 && Loca.method == -1){
           printf("PROBLEMS: The density you've selected is\n");
           printf("          in the coexistence region\n ");
           printf("          setting density to rho_vapor    \n ");
           Rho_b[0] = Rho_coex[1];
      }
   }*/


/*   if (Ncomp == 1 && Ipot_ff_n==2 && Iliq_vap != -2) P_over_po=Betap/p_coex;*/


}
/*******************************************************************************/
/* coexistence: This routine finds the coexisting densities of a pure fluid
                given a temperature */
double coexistence()
{
  double resid[2],Jac[2][2],Jac_inv[2][2];
  double betamu_hs_l[1],betamu_hs_g[1];
  double betamu_att_l[1],betamu_att_g[1];
  double rho_coex_old[2],p_coex;
  int i,j,niter;

  /* for one component only */

  /* Initial guesses for the densitites */
   Rho_coex[0] = 0.7;
   Rho_coex[1] = 0.01;

   /* calculate residuals at these densities */

   resid[0] = ( Rho_coex[0]+calc_hs_properties(betamu_hs_l,&(Rho_coex[0])) 
                      + calc_att_properties(betamu_att_l,&(Rho_coex[0])) ) -
              ( Rho_coex[1]+calc_hs_properties(betamu_hs_g,&(Rho_coex[1])) 
                      + calc_att_properties(betamu_att_g,&(Rho_coex[1])) );

   resid[1] = ( log(Rho_coex[0])+betamu_hs_l[0] + betamu_att_l[0]) - 
              ( log(Rho_coex[1])+betamu_hs_g[0] + betamu_att_g[0]);

   niter = 0;
if (niter==0) printf("Rho_coex=%9.6f  %9.6f\n",Rho_coex[0],Rho_coex[1]);
if (niter==0) printf("betamu_hs_l=%9.6f  betamu_hs_g=%9.6f  betamu_att_l=%9.6f  betamu_att_g=%9.6f\n",
                      betamu_hs_l[0],betamu_hs_g[0],betamu_att_l[0],betamu_att_g[0]);
if (niter==0) printf("resid= %9.6f  %9.6f\n", resid[0],resid[1]);


   /* perform Newton iterations to find the root */

   while (fabs(resid[0]) > 1e-10 || fabs(resid[1]) > 1e-10){

        rho_coex_old[0] = Rho_coex[0];
        rho_coex_old[1] = Rho_coex[1];

        Jac[0][0] = dp_drho_hs (&(Rho_coex[0])) + dp_drho_att(&(Rho_coex[0]));
        Jac[0][1] = dp_drho_hs (&(Rho_coex[1])) + dp_drho_att(&(Rho_coex[1]));
        Jac[1][0] = dmu_drho_hs(&(Rho_coex[0])) + dmu_drho_att(&(Rho_coex[0]));
        Jac[1][1] = dmu_drho_hs(&(Rho_coex[1])) + dmu_drho_att(&(Rho_coex[1]));

if (niter==0) printf("Jac: %9.6f  %9.6f  %9.6f  %9.6f\n",Jac[0][0],Jac[0][1],Jac[1][0],Jac[1][1]);


        /* matrix inversion is simple for this two by two */

        Jac_inv[1][0] = 1.0/(Jac[0][1]-Jac[0][0]*Jac[1][1]/Jac[1][0]);
        Jac_inv[1][1] = 1.0/(Jac[1][1]-Jac[1][0]*Jac[0][1]/Jac[0][0]);
        Jac_inv[0][0] = -(Jac[1][1]/Jac[1][0])*Jac_inv[1][0];
        Jac_inv[0][1] = -(Jac[0][1]/Jac[0][0])*Jac_inv[1][1];
if (niter==0) printf("Jac: %9.6f  %9.6f  %9.6f  %9.6f\n",Jac_inv[0][0],Jac_inv[0][1],Jac_inv[1][0],Jac_inv[1][1]);

        for (i=0; i< 2; i++) {
          for (j=0; j< 2; j++) {
             Rho_coex[i] += Jac_inv[i][j]*resid[j];
          }
        }

        if (Rho_coex[1] < 0.0) Rho_coex[1] = 0.0001 ;
        if (Rho_coex[0] > 1.0) Rho_coex[0] = rho_coex_old[0] ;


        resid[0] = ( Rho_coex[0]+calc_hs_properties(betamu_hs_l,&(Rho_coex[0])) 
                      + calc_att_properties(betamu_att_l,&(Rho_coex[0])) ) -
                   ( Rho_coex[1]+calc_hs_properties(betamu_hs_g,&(Rho_coex[1])) 
                      + calc_att_properties(betamu_att_g,&(Rho_coex[1])) );

        resid[1] = ( log(Rho_coex[0])+betamu_hs_l[0] + betamu_att_l[0]) - 
                   ( log(Rho_coex[1])+betamu_hs_g[0] + betamu_att_g[0]);
        niter++;
        if (niter > 100) {
           printf("problems locating coexistence....more than 100 iterations");
           exit(-1);
        }
   }
   p_coex = calc_hs_properties(betamu_hs_l,&(Rho_coex[0])) 
                      + calc_att_properties(betamu_att_l,&(Rho_coex[0])) ;
   return(p_coex);
}
/******************************************************************************
