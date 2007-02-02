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

/*  file dft_stencil_delta_Bond.c: contains functions that set up properties
    of a delta function stencil with a range of the bond length of the polymers
    in a CMS DFT calculations */

#include "dft_stencil_delta_Bond.h"

/*********************************************************************/
double StenDelta_Bond_sten_rad(int icomp, int jcomp)
{ return(Bond_ff[icomp][jcomp]); }

/*********************************************************************/
double StenDelta_Bond_sten_vol(int icomp,int jcomp)
{ 

/*        return (4.0 * PI * POW_DOUBLE_INT(Bond_ff[i][j],2) );*/
/*        in order to avoid having to carry prefactors, we will set the
          stencil volume to 1.0.  This will effectively perform the multiplication
          of 1/(4 PI Bond[icomp][jcomp]^2) times the native stencil.*/

   return(1.0); 
}
/*********************************************************************/
int StenDelta_Bond_Njcomp()
{ return(Ncomp); }

/*********************************************************************/
int StenDelta_Bond_NquadPtsBoundary()
{
   int npt;
   switch(Ndim){
      case 3: npt = 40; break;
      case 2:
      case 1: npt=150; break;
   }
return(npt);
}
/*********************************************************************/
int StenDelta_Bond_NquadPtsGauss(double r)
{
   int npt;
   if (r <= 4.000000001) npt=6; 
   else                  npt=3;
   return(npt);
}
/*********************************************************************/
double StenDelta_Bond_GetWeightFromSten(double rsq,double R)
{
   double weight;

    switch(Ndim){
      case 1:  weight=2.0 * PI * R;        break;
      case 2:  weight=2.0 / sqrt(1.0-rsq); break;
      case 3:  weight=1.0;                 break;
    }

    return(weight);
}
/*********************************************************************/
