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

/*  file dft_stencil_theta_GENmsa.c: contains functions that set up properties
    of a theta function multiplied by a direct correlation function correction
    that represents the difference between a charged and a hard sphere fluid 
    at a constant bulk density.  Analytical MSA results are used for delta c */ 

#include "dft_stencil_theta_GENmsa.h"

/*********************************************************************/
double StenTheta_GENmsa_sten_rad(int icomp,int jcomp)
{ return((HS_diam[icomp]+HS_diam[jcomp])/2.0); }

/*********************************************************************/
double StenTheta_GENmsa_sten_vol(int i,int j)
{

   double r_max,vol_sten;

   r_max = (HS_diam[i]+HS_diam[j])/2.0;
   vol_sten = deltaC_GENERAL_MSA_int(r_max,i,j);

   return(vol_sten);
}
/*********************************************************************/
int StenTheta_GENmsa_Njcomp()
{ return(Ncomp); }

/*********************************************************************/
int StenTheta_GENmsa_NquadPtsBoundary()
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
int StenTheta_GENmsa_NquadPtsGauss(double r)
{
   int npt;
   if (r <= 4.000000001) npt=6;
   else                  npt=3;
   return(npt);
}
/*********************************************************************/
int StenTheta_GENmsa_NquadPtsGaussIntegrand(double r)
{
   int npt;
   if (r <= 4.000000001)     npt=40;
   else if (r<= 16.00000001) npt=20;
   else                      npt=12;

   /*npt=20;*/   /* setting from old stencil code */ 
   return(npt);
}
/*********************************************************************/
double StenTheta_GENmsa_GetWeightFromSten(int icomp, int jcomp, double rsq,
                                         int ngpu, double *gpu, double *gwu)
{
  double weight, zmax, z, rho, hs;
  int i;
    
    hs = (HS_diam[icomp]+HS_diam[jcomp])/2.;
    
  weight = 0.0;
  zmax = sqrt(1 - rsq);
  switch(Ndim){
  case 1:
     for (i=0; i < ngpu; i++) {
         z = zmax * gpu[i];
         rho = sqrt(rsq + z*z) * hs;
         weight += gwu[i] * z *  deltaC_GENERAL_MSA(rho, icomp, jcomp);
     }
     return(2.0 * PI * weight * hs * hs * zmax);
     break;

  case 2:
     for (i=0; i < ngpu; i++) {
         z = zmax * gpu[i];
         rho = sqrt(rsq + z*z) * hs;
         weight += gwu[i] *  deltaC_GENERAL_MSA(rho, icomp, jcomp);
     }
     return(2.0 * weight * hs * zmax);
     break;

  case 3:
     rho = sqrt(rsq) * hs;
     weight = deltaC_GENERAL_MSA(rho, icomp, jcomp);
     return(weight);
     break;

  default:
     return(0.0);
     break;
  }
}
/*********************************************************************/
