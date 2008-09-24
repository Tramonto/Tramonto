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
 *  FILE: dft_guess_HSFMT.c
 *
 *  This file contains routines to set up an initial guess for the
 *  nonlocal density variables associated with FMT based perturbation DFTs.
 *
 */

#include "dft_guess_HSFMT.h"
 
/************************************************************/
/*setup_rho_bar: set up the rhobar initial guesses.  For now
                use Rho_b to set initial guess.  Later calculate
                based on rho initial guess. */
void setup_rho_bar(double **xInBox)
{
  int loc_inode,inode_box,inode,ijk[3],iunk,irb;
  double vol,area,x_dist;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     for (irb = 0; irb < Nrho_bar; irb++){
       iunk = Phys2Unk_first[HSRHOBAR] + irb;

       if (Lsteady_state!=UNIFORM_INTERFACE){
           inode     = B2G_node[inode_box];
           node_to_ijk(inode,ijk); 
           x_dist = Esize_x[Grad_dim]*ijk[Grad_dim];

           xInBox[iunk][inode_box] = Rhobar_b_LBB[irb] + 
                      (Rhobar_b_RTF[irb]-Rhobar_b_LBB[irb])*
                           x_dist/Size_x[Grad_dim];
       }
       else {
          xInBox[iunk][inode_box] = Rhobar_b[irb];
       }
     }
  }
  return;
}
/************************************************************/
/*calc_init_rho_bar: set up the rhobar initial guesses based on a known
                density profile. Set these variables up on local nodes only.*/
void calc_init_rho_bar(double **xInBox)
{
  int loc_inode,inode_box,inode,ijk[3],iunk,irb;
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box=L2B_node[loc_inode];
     for (irb = 0; irb < Nrho_bar; irb++){
       iunk = Phys2Unk_first[HSRHOBAR] + irb;
       xInBox[iunk][inode_box]=int_stencil_HSFMT(xInBox,inode_box,iunk);
       if (irb==0 && xInBox[iunk][inode_box] >1.0) xInBox[iunk][inode_box]=0.98;
     }
  }
  return;
}
/************************************************************/
