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
      Sten_Type[DELTA_FN]=TRUE;
      Sten_Type[U_ATTRACT]=FALSE;   /* attractions handled differently for polymers */
      Sten_Type[POLYMER_CR]=TRUE;
      if (Type_poly==CMS_SCFT){
        printf ("To do SCFT with CMS theory, we need to test and debug all code !\n");
        exit(-1);
      }
  }
  else{
     if (Type_func !=NONE){
        Sten_Type[DELTA_FN]=TRUE;
        Sten_Type[THETA_FN]=TRUE;
      }
      if (Type_attr != NONE) Sten_Type[U_ATTRACT]=TRUE;
      if (Type_coul == DELTAC) Sten_Type[THETA_CHARGE]=TRUE;
      if (Type_poly == WTC) {
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
     case DELTA_FN:      return TRUE;
     case THETA_FN:      return FALSE;
     case THETA_FN_SIG:  return FALSE;
     case U_ATTRACT:     return FALSE;
     case THETA_CHARGE:  return FALSE;
     case POLYMER_CR:    return FALSE;
     case DELTA_FN_BOND: return TRUE;
     default:
         printf("problem with stencil definitions: stencil_deltaLogical");
         exit(-1); break;
  }
}
/****************************************************************************/


