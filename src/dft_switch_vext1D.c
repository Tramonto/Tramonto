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

/*
 *  FILE: dft_vext1D.c
 *
 *  This file contains routines selects and computes the correct external
 *  field potential (and its derivative) to use when setting up a 1D external field
 *
 */

#include "dft_switch_vext1D.h"

/******************************************************************************/
/* Vext_1D:  given a wall-fluid point separation, calculate the 
           9-3 Lennard Jones potential (note that prefactors are
           calculated in calling routine)                           */
double Vext_1D(double x,int icomp, int iwall_type)
{
  double vext;

  switch(Vext_PotentialID[iwall_type]){
      case LJ9_3_CS:
        vext = Vext_LJ9_3_CS(x,icomp,iwall_type);
        break;
      case LJ9_3_v2_CS:
        vext = Vext_LJ9_3_v2_CS(x,icomp,iwall_type);
        break;
      case LJ9_3_noCS:
        vext = Vext_LJ9_3_noCS(x,icomp,iwall_type);
        break;
      case LJ9_3_shiftX_CS:
        vext = Vext_LJ9_3_shiftX_CS(x,icomp,iwall_type);
        break;
      case REPULSIVE9_noCS:
        vext = Vext_REPULSIVE9_noCS(x,icomp,iwall_type);
        break;
      case EXP_ATT_noCS:
        vext = Vext_EXP_ATT_noCS(x,icomp,iwall_type);
        break;
      case R7_YUKAWA_SUM_CS:
        vext = Vext_LJ7YukawaSum(x,icomp,iwall_type);
        break;
      case LINEAR_noCS:
        vext = Vext_LINEAR_noCS(x,icomp,iwall_type);
        break;
      default:
         if (Iwrite_screen !=SCREEN_NONE) printf("problems with your selection of Vext_PotentialID[iwall_type=%d]\n",iwall_type);
         exit(-1);
         break;
  }
  return vext;
}
/******************************************************************************/
/* Vext_1D_dash:  given a wall-fluid point separation, calculate the derivative of
               the 9-3 Lennard Jones potential (note that prefactors are
               calculated in calling routine)                           */
double Vext_1D_dash(double x,int icomp, int iwall_type)
{
  double vdash;

  switch(Vext_PotentialID[iwall_type]){
      case LJ9_3_CS:
      case LJ9_3_noCS:
        vdash = Vextderiv_LJ9_3(x,icomp,iwall_type);
        break;
      case LJ9_3_v2_CS:
        vdash = Vextderiv_LJ9_3_v2(x,icomp,iwall_type);
        break;
      case REPULSIVE9_noCS:
        vdash = Vextderiv_REPULSIVE9(x,icomp,iwall_type);
        break;
      case EXP_ATT_noCS:
        vdash = Vextderiv_EXP_ATT(x,icomp,iwall_type);
        break;
      case R7_YUKAWA_SUM_CS:
        vdash = Vextderiv_LJ7YukawaSum(x,icomp,iwall_type);
        break;
      case LINEAR_noCS:
        vdash = Vextderiv_LINEAR(x,icomp,iwall_type);
        break;
      case LJ9_3_shiftX_CS:
         vdash=0.0;  /* vdash code not yet written */
         break;
      default:
         if (Iwrite_screen !=SCREEN_NONE) printf("problems with your selection of Vext_PotentialID[iwall_type=%d]\n",iwall_type);
         exit(-1);
         break;
  }

  return vdash;
}
/******************************************************************************/
