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

/*  file dft_stencil_theta_CrCMS.c: contains functions that set up properties
    of a theta function (multiplied with a direct correlation function for 
    CMS DFT) stencil with a range of the radius of the direct correlation function. */

#include "dft_stencil_theta_CrCMS.h"

/*********************************************************************/
double StenTheta_CrCMS_sten_rad(int icomp,int jcomp)
{ 
  return(Cr_rad[icomp][jcomp]); 
}

/*********************************************************************/
double StenTheta_CrCMS_sten_vol(int icomp,int jcomp)
{ 
   return((double)NO_RENORMALIZATION_FLAG);  /* note return a flag to skip renormalization of the
                                                stencil.  Since the direct correlation function is
                                                read in (not an analytical function) it is not
                                                possible to compute an "exact" integral value */
}

/*********************************************************************/
int StenTheta_CrCMS_Njcomp()
{ return(Ncomp); }

/*********************************************************************/
int StenTheta_CrCMS_NquadPtsBoundary()
{
   int npt;
   switch(Ndim){
      case 3: npt = 20; break;
      case 2:
      case 1: npt=20; break;
   }
return(npt);
}
/*********************************************************************/
int StenTheta_CrCMS_NquadPtsGauss(double r)
{
   int npt;
   if (r <= 4.000000001) npt=6;
   else                  npt=3;
   return(npt);
}
/*********************************************************************/
double StenTheta_CrCMS_GetWeightFromSten(int icomp, int jcomp, double rsq, double R)
{
   double rmin, rlast_nz, r_upp,slope_dr,r_low,zsq,rx_low,weight;
   int irmin;                       

   switch(Ndim){

   case 3:
     rmin = sqrt(rsq*R*R);  /* need real distance for this case */
     irmin = (int) (rmin/Deltar_cr + 0.00000001);

     if (irmin < Last_nz_cr)
         slope_dr = (Rism_cr[icomp][jcomp][irmin+1] -
                       Rism_cr[icomp][jcomp][irmin]);
     else
         slope_dr = (Rism_cr[icomp][jcomp][irmin] -
                       Rism_cr[icomp][jcomp][irmin-1]);

         weight = Rism_cr[icomp][jcomp][irmin] +
                            slope_dr*(rmin - irmin*Deltar_cr)/Deltar_cr;
     break;

   case 2:
   case 1:
      rmin = sqrt(rsq*R*R);  /* need real distance for this case */
      zsq = rmin*rmin;
      irmin = (int) (rmin/Deltar_cr + 0.00000001);
/*    last_nz_cr = (int) (R/Deltar_cr + 0.00000001);     */
      rlast_nz   = Deltar_cr*Last_nz_cr;
      weight = 0.;
      if (irmin < Last_nz_cr) {
         r_upp = Deltar_cr*(irmin+1);
         if (r_upp > R) r_upp = R;
         slope_dr = (Rism_cr[icomp][jcomp][irmin+1] -
                       Rism_cr[icomp][jcomp][irmin]);
         rx_low = 0.;
         weight = int_cr(rmin,r_upp,slope_dr,icomp,jcomp,irmin,zsq,&rx_low);
         for (irmin = irmin+1; irmin < Last_nz_cr; irmin++) {
             if (r_upp < R){
              r_low = r_upp;
              r_upp = Deltar_cr*(irmin+1);
              if (r_upp > R) r_upp = R;
              slope_dr = (Rism_cr[icomp][jcomp][irmin+1] -
                       Rism_cr[icomp][jcomp][irmin]);
             weight += int_cr(r_low,r_upp,slope_dr,icomp,jcomp,irmin,zsq,&rx_low);
             }
         }
      }
      if (rlast_nz < R) {
        if (rlast_nz > rmin){
           r_low = rlast_nz;
           rx_low = sqrt(r_low*r_low - zsq);
         }
         else{
              r_low = rmin;
              rx_low = 0.;
         }
      r_upp = R;
      slope_dr = (Rism_cr[icomp][jcomp][Last_nz_cr] -
                      Rism_cr[icomp][jcomp][Last_nz_cr-1]);
      weight += int_cr(r_low,r_upp,slope_dr,icomp,jcomp,Last_nz_cr,zsq,&rx_low);
      }
      break;

    }
    return(weight);
}
/********************************************************************************/
/*  int_cr:  integrate the direct correlation function numerically              */
double int_cr(double r_low,double r_upp,double slope_dr,int icomp,int jcomp,
                  int irmin, double zsq, double *rx_low)
{
  double temp,rusq,rlsq,rx_upp;
/*
if (icomp==0 && jcomp==0 && fabs(r_low-3.411270)<1.e-4 && fabs(r_upp- 3.42500)<1.e-4){
   printf("Deltar_cr=%9.6f\n",Deltar_cr);
}
*/
  rusq = r_upp*r_upp;
  rlsq = r_low*r_low;
  if (Ndim == 1)   {
    temp = 2.*PI*(
      ( rusq*(Rism_cr[icomp][jcomp][irmin] - slope_dr*irmin)/2. +
        slope_dr*rusq*r_upp/(3.*Deltar_cr) ) -
      ( rlsq*(Rism_cr[icomp][jcomp][irmin] - slope_dr*irmin)/2. +
        slope_dr*rlsq*r_low/(3.*Deltar_cr) )   );
  }
  else{
    rx_upp = sqrt(r_upp*r_upp - zsq);
    temp = 2.*(rx_upp - *rx_low)*(Rism_cr[icomp][jcomp][irmin] - slope_dr*irmin)
           + slope_dr*( (r_upp*rx_upp - r_low* *rx_low) +
           zsq*log((r_upp + rx_upp)/(r_low + *rx_low)) )/Deltar_cr;
    *rx_low = rx_upp;
  }
  return temp;
}
/*********************************************************************/

