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

/*  file dft_stencil_theta_RPMmsa.c: contains functions that set up properties
    of a theta function (multiplied by an attractive pair potential) 
    stencil with a range of the pair potential cutoff between two atoms. */

#include "dft_stencil_theta_RPMmsa.h"

/*********************************************************************/
double StenTheta_RPMmsa_sten_rad(int icomp)
{ return(Sigma_ff[icomp][icomp]); }

/*********************************************************************/
double StenTheta_RPMmsa_sten_vol(int i,int j)
{

   double r_max,vol_sten;

   r_max = Sigma_ff[i][j];
   vol_sten = deltaC_MSA_int(r_max,i,j);

   return(vol_sten);
}
/*********************************************************************/
int StenTheta_RPMmsa_Njcomp()
{ return(Ncomp); }

/*********************************************************************/
int StenTheta_RPMmsa_NquadPtsBoundary()
{
   int npt;
   switch(Ndim){
      case 3: 
      case 2:
      case 1: npt=20; break;
   }
   return(npt); 
}
/*********************************************************************/
int StenTheta_RPMmsa_NquadPtsGauss(double r)
{
   int npt;
   if (r <= 4.000000001) npt=6;
   else                  npt=3;
   return(npt);
}
/*********************************************************************/
int StenTheta_RPMmsa_NquadPtsGaussIntegrand(double r)
{
   int npt;
   if (r <= 4.000000001)     npt=40;
   else if (r<= 16.00000001) npt=20;
   else                      npt=12;

   npt=20;   /* setting from old stencil code */ 
   return(npt);
}
/*********************************************************************/
double StenTheta_RPMmsa_GetWeightFromSten(int icomp, int jcomp, double rsq,
                                         int ngpu, double *gpu, double *gwu)
{
  double weight, zmax, z, rho;
  int i;

  weight = 0.0;
  zmax = sqrt(1 - rsq);
  switch(Ndim){
  case 1:
     for (i=0; i < ngpu; i++) {
         z = zmax * gpu[i];
         /* perhaps rho doesn't need Sigma_ff multiplier */
         rho = sqrt(rsq + z*z) * Sigma_ff[icomp][jcomp];
         weight += gwu[i] * z *  deltaC_MSA(rho, icomp, jcomp);
     }
     return(2.0 * PI * weight * Sigma_ff[icomp][jcomp]
                            * Sigma_ff[icomp][jcomp] * zmax);
     break;

  case 2:
     for (i=0; i < ngpu; i++) {
         z = zmax * gpu[i];
         /* perhaps rho doesn't need Sigma_ff multiplier */
         rho = sqrt(rsq + z*z) * Sigma_ff[icomp][jcomp];
         weight += gwu[i] *  deltaC_MSA(rho, icomp, jcomp);
     }
     return(2.0 * weight * Sigma_ff[icomp][jcomp] * zmax);
     break;

  case 3:
     rho = sqrt(rsq) * Sigma_ff[icomp][jcomp];
     weight = deltaC_MSA(rho, icomp, jcomp);
     return(weight);
     break;
  }
}
/*********************************************************************/
