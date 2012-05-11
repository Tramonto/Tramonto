/*
//@HEADER
// ********************************************************************
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/* -----------------------------------------------------------------------------------------
dft_thermo_elec_MSA.c: Calculate the thermodynamic properties for a charged system where MSA
corrections are turned on.
--------------------------------------------------------------------------------------------*/
#include "dft_thermo_elec_MSA.h"

/********************************************************************/
/*chempot_ELEC_MSA_RPM: Here we compute the chemical potential contribution due
      to cross correlations between the hard sphere and coulomb parts of 
      the potential for the restricted primitive model based using the 
      mean spherical approximation (following the work of Tang and Davis)*/
void chempot_ELEC_MSA_RPM(double *rho)
{
   int icomp,jcomp;

   Deltac_b = (double *) array_alloc (1, Ncomp, sizeof(double));
   for (icomp=0; icomp<Ncomp; icomp++) Deltac_b[icomp] = 0.0;

   for (icomp=0; icomp<Ncomp; icomp++) 
      for (jcomp=0; jcomp<Ncomp;jcomp++){
          Deltac_b[icomp]+= rho[jcomp]*int_stencil_bulk(THETA_CR_RPM_MSA,icomp,jcomp,NULL);
      }
   return;
}
/*******************************************************************************/
/* deltaC_MSA:  calculate deltac_MSA */

double deltaC_MSA(double r,int i, int j)
{
  double deltac,B,kappa,kappa_sq;
  int icomp;

  kappa_sq = 0.0;
  for(icomp = 0; icomp<Ncomp; icomp++)
     kappa_sq += (4.0*PI/Temp_elec)*Rho_b[icomp]*Charge_f[icomp]*Charge_f[icomp];
  kappa = sqrt(kappa_sq);

  /* note that kappa is in units of (kappa*sigma_ref) so we need to explicitly include
     hard sphere diameters (units of HS_diam/sigma_ref) */

  /******* NOTE ----- epsilon in the Davis and Hanson papers should be interpreted 
   as 4*pi*epsilon_r*epsilon_o to be consistent with SI units used in the code here.
   Also note that Sigma has been replaced with HS_diam below (4/2010 LF) */

  B = (kappa*HS_diam[i] + 1.0 - sqrt(1.0+2.0*kappa*HS_diam[i]))/(kappa*HS_diam[i]);

  if (r == 0.0) printf("trouble with deltaC term .... r=0");
  if (r <= HS_diam[i] && r>0) {                   

     deltac = -Charge_f[i]*Charge_f[j]/Temp_elec*
              (  2*B/HS_diam[i] - 1.0/r
               - POW_DOUBLE_INT(B/HS_diam[i],2)*r );
  }
  else deltac = 0.0;

  return deltac;
}
/*******************************************************************************/
/* deltaC_MSA_int:  given range of integrattion, r, calculate the 
           integral of deltac_MSA over all space */

double deltaC_MSA_int(double r,int i, int j)
{
  double deltac_int,B,kappa,kappa_sq;
  int icomp;

  kappa_sq = 0.0;
  for(icomp = 0; icomp<Ncomp; icomp++)
     kappa_sq += (4.0*PI/Temp_elec)*Rho_b[icomp]*
                  Charge_f[icomp]*Charge_f[icomp];
  kappa = sqrt(kappa_sq);
  B = (kappa*HS_diam[i] + 1.0 - sqrt(1.0+2.0*kappa*HS_diam[i]))/(kappa*HS_diam[i]);


   /********* NOTE --- does 4PI r^2 here is the integration over the spherical surface at a distance r, and is not
       related to the plasma parameters (T_elec).*/

   deltac_int = -(4*PI*r*r)*(Charge_f[i]*Charge_f[j]/Temp_elec)*
               (  2*B*r/(3.0*HS_diam[i]) - 0.5
               - 0.25*POW_DOUBLE_INT(B/HS_diam[i],2)*r*r );

/*printf("deltaC_int:::: r=%g i=%d  j=%d  deltac_int=%g  \n",r,i,j,deltac_int);
printf("\t prefac=%g  term1(r)=%g  term2(no r)=%g  term3(rsq)=%g\n",
      -(4*PI*r*r)*(Charge_f[i]*Charge_f[j]/Temp_elec),
      2*B*r/(3.0*HS_diam[i]), -0.5,- 0.25*POW_DOUBLE_INT(B/HS_diam[i],2)*r*r);*/

  return deltac_int;
}
/*******************************************************************************/
