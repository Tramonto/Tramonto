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
 *  FILE: dft_guess_noinfo_WJDC.c
 *
 *  This file contains routines that set up an initial guess for the
 *  WJDC polymer DFT.
 *
 */

#include "dft_guess_WJDC.h"
 
/*********************************************************/
/*setup_polymer_field: in this routine sets up the initial guess for the WJDC field variable */
void setup_polymer_field_wjdc(double **xOwned)
{
  int loc_inode,itype_mer,irho, iunk,i,Nloop;
  double field;

  Nloop=Ncomp;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     for (i=0; i<Nloop; i++){
         iunk=Phys2Unk_first[WJDC_FIELD]+i;
        /* if (Field_WJDC_b[i]<10.) */ xOwned[iunk][loc_inode]=Field_WJDC_b[i];
/*         else xOwned[iunk][loc_inode]=1.;*/
     }
   }
   return;
}
/*********************************************************/
/*setup_polymer_G_wjdc: in this routine sets up the initial guess for the chain variable
in the wjdc functional */
void setup_polymer_G_wjdc(double **xOwned)
{
  int loc_inode,itype_mer,irho, iunk,i,Nloop;
  double field;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     for (i=0; i<Nbonds; i++){
         iunk=Phys2Unk_first[G_CHAIN]+i;
        /* if (G_WJDC_b[i]<10.)*/ xOwned[iunk][loc_inode]=G_WJDC_b[i];
/*         else xOwned[iunk][loc_inode]=1.;*/
     }
   }
   return;
}
/*********************************************************/
