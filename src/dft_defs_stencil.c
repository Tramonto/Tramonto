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
 *  FILE: dft_defs_stencil.c
 *
 *  This file contains routines to control setup of the integration stencils for Tramonto.
 *
 */

#include "dft_defs_stencil.h"
/*********************************************************************************/
void setup_stencil_logicals()
{
  int isten;
  /* Read in Functional Switches  and set up the Sten_Type array.*/
  for (isten=0; isten<NSTEN; ++isten) Sten_Type[isten]=FALSE;

  if (Type_poly == CMS || Type_poly == CMS_SCFT){
      Sten_Type[DELTA_FN_BOND]=TRUE;
      if (Type_poly == CMS) Sten_Type[THETA_CR_DATA]=TRUE;
  }
  else{
     if (Type_func !=NONE){
        Sten_Type[DELTA_FN_R]=TRUE;
        Sten_Type[THETA_FN_R]=TRUE;
      }
      if (Type_attr != NONE) Sten_Type[THETA_PAIRPOT_RCUT]=TRUE;
      if (Type_coul == DELTAC_RPM) Sten_Type[THETA_CR_RPM_MSA]=TRUE;
      if (Type_coul == DELTAC_GENERAL) Sten_Type[THETA_CR_GENERAL_MSA]=TRUE;

      if (Type_poly == WTC || Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3) {
        Sten_Type[THETA_FN_SIG]=TRUE;
        Sten_Type[DELTA_FN_BOND]=TRUE;
      }
  }
  return;
}
/***********************************************************************************/
/*stencil_deltaLogical: return a logical (mapped to an int) that
   identifies all stencils based on delta functions. */
int stencil_deltaLogical(int sten)
{

  switch(sten){
     case DELTA_FN_R:      return TRUE;
     case THETA_FN_R:      return FALSE;
     case THETA_FN_SIG:  return FALSE;
     case THETA_PAIRPOT_RCUT:     return FALSE;
     case THETA_CR_RPM_MSA:  return FALSE;
     case THETA_CR_GENERAL_MSA:  return FALSE;
     case THETA_CR_DATA:    return FALSE;
     case DELTA_FN_BOND: return TRUE;
     default:
         printf("problem with stencil definitions: stencil_deltaLogical");
         exit(-1); break;
  }
}
/****************************************************************************/


