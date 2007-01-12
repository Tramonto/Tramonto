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

/*  file dft_stencil_theta_R.c: contains functions that set up properties
    of a theta function stencil with a range of the radius of an atom */

#include "dft_stencil_theta_R.h"

/*********************************************************************/
double StenTheta_R_sten_rad(int icomp)
{ return(HS_diam[icomp]/2.0); }

/*********************************************************************/
double StenTheta_R_sten_vol(int i)
{ return(PI * POW_DOUBLE_INT(HS_diam[i],3)/6.0); }

/*********************************************************************/
int StenTheta_R_Njcomp()
{ return(1); }

/*********************************************************************/
int StenTheta_R_NquadPtsBoundary()
{ 
   int npt;
   switch(Ndim){
      case 3: npt = 20; break;
      case 2:  
      case 1: npt=40; break;
   }
return(npt); 
}
/*********************************************************************/
int StenTheta_R_NquadPtsGauss(double r)
{
   int npt;
   if (r <= 4.000000001) npt=6;
   else                  npt=3;
   return(npt);
}
/*********************************************************************/
double StenTheta_R_GetWeightFromSten(double rsq, double R)
{
   double weight;

    switch(Ndim){
      case 1:  weight=(1.0-rsq) * PI*R*R;        break;
      case 2:  weight=sqrt(1.0-rsq) * 2.0*R; break;
      case 3:  weight=1.0;                 break;
    }
    return(weight);
}
/*********************************************************************/
