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

/*  file dft_stencil_delta_BondCMS.c: contains functions that set up properties
    of a delta function stencil with a range of the bond length of the polymers
    in a CMS DFT calculations */

#include "dft_stencil_delta_BondCMS.h"
/*********************************************************************/
double StenDelta_BondCMS_sten_rad(int icomp,int jcomp)
{ return(Bond_ff[icomp][jcomp]); }

/*********************************************************************/
double StenDelta_BondCMS_sten_vol(int icomp,int jcomp)
{ return(1.0); }

/*********************************************************************/
int StenDelta_BondCMS_Njcomp()
{ return(Ncomp); }

/*********************************************************************/
int StenDelta_BondCMS_NquadPtsBoundary()
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
int StenDelta_BondCMS_NquadPtsGauss(double r)
{
   int npt;
   if (r <= 4.000000001) npt=6;
   else                  npt=3;
   return(npt);
}
/*********************************************************************/
double StenDelta_BondCMS_GetWeightFromSten(int icomp,int jcomp,double rsq,double R)
{
   double weight,sigsq;

    sigsq = Bond_ff[icomp][jcomp]*Bond_ff[icomp][jcomp];
/*  sigsq = Sigma_ff[icomp][jcomp]*Sigma_ff[icomp][jcomp];*/
         /* note need to check which of these should be used - if Bond_ff, note
            that Bond_ff[i][j]=sten_rad=R in this function so the equations
            below can be cast in terms of R and rsq only....the icomp,jcomp
            should not be needed here. */
    switch(Ndim){
      case 1:  weight=R / (2.0*sigsq);        break;
      case 2:  weight=0.5 / (sqrt(1.0-rsq) * PI*sigsq); break;
      case 3:  weight=0.25 / (PI*sigsq);                 break;
    }
    return(weight);
}
/*********************************************************************/
