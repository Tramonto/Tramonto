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
 *  FILE: dft_guess_MFATT.c
 *
 *  This file contains routines to set up an initial guess for the
 *  variables associated with mean field attractions.  
 *
 */

#include "dft_guess_MFATT.h"
 
/************************************************************/
/*setup_mf_attract: set up variables for mean field attractions. */
void setup_mf_attract(double **xOwned)
{
  int loc_inode,inode_box,inode,ijk[3],iunk,icomp,jcomp;
  double vol,area,x_dist,sum;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     for (icomp = 0; icomp < Ncomp; icomp++){
       iunk = Phys2Unk_first[MF_EQ] + icomp;
       if (Lsteady_state || (Nwall == 0 && Iliq_vap == 3)){
           printf("attractions not set up for steady state yet. \n");
           exit(-1);
       }
       else {
          sum=0.0;
          for (jcomp=0; jcomp<Ncomp;jcomp++){
            sum += Avdw[icomp][jcomp]*Rho_b[jcomp];
          }
          xOwned[iunk][loc_inode] = sum;
       }
     }
  }
  return;
}
/************************************************************/
