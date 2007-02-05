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

/********************************************************************
chempot_ELEC_MSA: Here we compute the chemical potential contribution due
      to cross correlations between the hard sphere and coulomb parts of 
      the potential for the restricted primitive model based using the 
      mean spherical approximation (following the work of Tang and Davis)*/
void chempot_ELEC_MSA(double *rho)
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
/* deltaC_MSA:  given r12, calculate the attractive part of a cut and
           shifted 12-6 LJ potential. */

double deltaC_MSA(double r,int i, int j)
{
  double deltac,B,kappa,kappa_sq;
  int icomp;

  kappa_sq = 0.0;
  for(icomp = 0; icomp<Ncomp; icomp++)
     kappa_sq += (4.0*PI/Temp_elec)*Rho_b[icomp]*
                  Charge_f[icomp]*Charge_f[icomp];
  kappa = sqrt(kappa_sq);
  B = (kappa + 1.0 - sqrt(1.0+2.0*kappa))/kappa;

/*  printf("\t r: %9.6f icomp: %d  jcomp: %d  kappa: %9.6f  B: %9.6f  Sigma_ff: %9.6f\n",
          r,i,j,kappa,B,Sigma_ff[i][j]);*/

  if (r == 0.0) printf("trouble with deltaC term .... r=0");
  if (r <= Sigma_ff[i][j] && r>0) {

     deltac = -Charge_f[i]*Charge_f[j]/Temp_elec*
              (  2*B/Sigma_ff[i][j] - 1.0/r
               - POW_DOUBLE_INT(B/Sigma_ff[i][j],2)*r );
  }
  else deltac = 0.0;

  return deltac;
}
/*******************************************************************************/
/* deltaC_MSA_int:  given range of integrattion, r, calculate the definite
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
  B = (kappa + 1.0 - sqrt(1.0+2.0*kappa))/kappa;

  deltac_int = -(4*PI*Charge_f[i]*Charge_f[j]/Temp_elec)*
                r*r*
               (  2*B*r/(3.0*Sigma_ff[i][j]) - 0.5
               - 0.25*POW_DOUBLE_INT(B/Sigma_ff[i][j],2)*r*r );

  return deltac_int;
}
/*******************************************************************************/
