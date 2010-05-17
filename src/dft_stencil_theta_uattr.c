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

/*  file dft_stencil_theta_uattr.c: contains functions that set up properties
    of a theta function (multiplied by an attractive pair potential) 
    stencil with a range of the pair potential cutoff between two atoms. */

#include "dft_stencil_theta_uattr.h"

/*********************************************************************/
double StenTheta_uattr_sten_rad(int icomp,int jcomp)
{ return(Cut_ff[icomp][jcomp]); }

/*********************************************************************/
double StenTheta_uattr_sten_vol(int i,int j)
{

   double r_min,r_cut,vol_sten,rCore_left,rCore_right,epsCore;
   int    LCore;

   r_cut = Cut_ff[i][j];
   pairPot_InnerCore_switch(i,j,Type_pairPot,&rCore_left,&rCore_right,&epsCore);


   if (fabs(epsCore)<1.e-8) LCore=FALSE;
   else                     LCore=TRUE;

   if (LCore==TRUE){
       if (fabs(rCore_left)<Esize_x[0]){   /* WCA-like case where there is only one constant in the core */
          r_min=rCore_right;
          vol_sten =  (4.0/3.0)*PI*pow(r_min,3.0)*pairPot_ATT_noCS_switch(r_min,i,j,Type_pairPot)
                    - (4.0/3.0)*PI*pow(r_cut,3.0)*pairPot_ATT_noCS_switch(r_cut,i,j,Type_pairPot)
                     + pairPot_integral_switch(r_cut,i,j,Type_pairPot) 
                     - pairPot_integral_switch(r_min,i,j,Type_pairPot);

       }
       else{                                  /* Jain-Dominik-Chapman case two constants in the core.  Note that
                                                 the first constant (from r=0 to r=rCore_left is assumed to be 
                                                 zero (u=0 for r<Sigma in the JDC case). */
          r_min=rCore_right;
          vol_sten =  (4.0/3.0)*PI*pow(r_min,3.0)*pairPot_ATT_noCS_switch(r_min,i,j,Type_pairPot)
              - (4.0/3.0)*PI*pow(r_cut,3.0)*pairPot_ATT_noCS_switch(r_cut,i,j,Type_pairPot)
              - (4.0/3.0)*PI*pow(rCore_left,3.0)*pairPot_ATT_noCS_switch(r_min,i,j,Type_pairPot)
              + (4.0/3.0)*PI*pow(rCore_left,3.0)*pairPot_ATT_noCS_switch(r_cut,i,j,Type_pairPot)
              + pairPot_integral_switch(r_cut,i,j,Type_pairPot) - pairPot_integral_switch(r_min,i,j,Type_pairPot);

       }
   }
   else{
       if (Type_pairPot==PAIR_SW){   /* this one case is different because there is a strict discontinuity at r=rcut */
          r_min = rCore_right;
          vol_sten = pairPot_integral_switch(r_cut,i,j,Type_pairPot) - pairPot_integral_switch(r_min,i,j,Type_pairPot);
       }
       else{ 
          r_min = rCore_right;
          /* volume of stencil from r_min to r_cut -- adjusted for cut and shift*/
          vol_sten = (pairPot_integral_switch(r_cut,i,j,Type_pairPot) 
                      - (4.0/3.0)*PI*pow(r_cut,3.0)*pairPot_ATT_noCS_switch(r_cut,i,j,Type_pairPot))
                - (pairPot_integral_switch(r_min,i,j,Type_pairPot) 
                      - (4.0/3.0)*PI*pow(r_min,3.0)*pairPot_ATT_noCS_switch(r_cut,i,j,Type_pairPot));
       }
   }
   return(vol_sten);
}
/*********************************************************************/
int StenTheta_uattr_Njcomp()
{ return(Ncomp); }

/*********************************************************************/
int StenTheta_uattr_NquadPtsBoundary()
{
   int npt;
   switch(Ndim){
      case 3: npt = 1; break;
      case 2:
      case 1: npt=20; break;
   }
return(npt);
}
/*********************************************************************/
int StenTheta_uattr_NquadPtsGauss(double r)
{
   int npt;

   npt=1;  /* independent of distance because this stencil can be very large */

   return(npt);
}
/*********************************************************************/
int StenTheta_uattr_NquadPtsGaussIntegrand(double r)
{
   int npt;
   if (r <= 4.000000001)     npt=40;
   else if (r<= 16.00000001) npt=20;
   else                      npt=12;
   return(npt);
}

/*********************************************************************/
double StenTheta_uattr_GetWeightFromSten(int icomp, int jcomp, double rsq,
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
            rho = sqrt(rsq + z*z) * Cut_ff[icomp][jcomp];
            weight += gwu[i] * z * pairPot_ATT_CS_switch(rho, icomp, jcomp,Type_pairPot);
         }
         return(2.0 * PI * weight * Cut_ff[icomp][jcomp] * Cut_ff[icomp][jcomp] * zmax);
         break;

      case 2:
         for (i=0; i < ngpu; i++) {
             z = zmax * gpu[i];
             rho = sqrt(rsq + z*z) * Cut_ff[icomp][jcomp];
             weight += gwu[i] * pairPot_ATT_CS_switch(rho, icomp, jcomp,Type_pairPot);
         }
         return(2.0 * weight * Cut_ff[icomp][jcomp] * zmax);
         break;

      case 3:
         rho = sqrt(rsq) * Cut_ff[icomp][jcomp];
         weight = pairPot_ATT_CS_switch(rho, icomp, jcomp,Type_pairPot);
         return(weight);
         break;
  }
}
/*********************************************************************/
