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
 *  FILE: dft_guess_WTC.c
 *
 *  This file contains routines that to set up an initial guess for 
 *  bonding and cavity function variables associated with the 
 *  Tripathi-Chapman functionals (based on Wertheim's theory).
 *
 */

#include "dft_guess_WTC.h"
 
/************************************************************/
/*setup_Xi_cavWTC: set up the cavity function initial guesses for 
                 Wertheim-Tripathi-Chapman functionals.  For now
                use bulk segment densities initial guess.  Later calculate
                based on rho initial guess. */
void setup_Xi_cavWTC(double **xInBox)
{
  int loc_inode,iunk,icav,inode_box;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box=L2B_node[loc_inode];
     for (icav = 0; icav < Phys2Nunk[CAVWTC]; icav++){
       iunk = Phys2Unk_first[CAVWTC] + icav;
       xInBox[iunk][inode_box] = Xi_cav_b[icav+2];
     }
  }
  return;
}
/************************************************************/
/*setup_BondWTC: set up the bond function initial guesses for 
                 Wertheim-Tripathi-Chapman functionals.  For now
                use bulk segment densities initial guess.  Later calculate
                based on rho initial guess. */
void setup_BondWTC(double **xInBox)
{
  int loc_inode,iunk,ibond,inode_box;;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box=L2B_node[loc_inode];
     for (ibond = 0; ibond < Nbonds; ibond++){
       iunk = Phys2Unk_first[BONDWTC] + ibond;
       xInBox[iunk][inode_box] = BondWTC_b[ibond];
     }
  }
  return;
}
/************************************************************/
/*calc_init_Xi_cavWTC: compute the cavity function initial guesses 
                 based on known density profiles.  */
void calc_init_Xi_cavWTC(double **xInBox)
{ 
  int loc_inode,iunk,icav,inode_box;
  
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box=L2B_node[loc_inode];
     for (icav = 0; icav < Phys2Nunk[CAVWTC]; icav++){
       iunk = Phys2Unk_first[CAVWTC] + icav;
       xInBox[iunk][inode_box] = int_stencil_CAV(xInBox,inode_box,iunk);
     }
  }
  return;
}
/************************************************************/
/*calc_init_BondWTC: set up the bond function initial guesses for 
                 Wertheim-Tripathi-Chapman functionals.  For now 
                use bulk segment densities initial guess.  Later calculate
                based on rho initial guess. */
void calc_init_BondWTC(double **xInBox)
{ 
  int loc_inode,iunk,ibond,inode_box;;
  
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box=L2B_node[loc_inode];
     for (ibond = 0; ibond < Nbonds; ibond++){
       iunk = Phys2Unk_first[BONDWTC] + ibond;
       xInBox[iunk][inode_box] = int_stencil_BondWTC(xInBox,inode_box,iunk);
     }
  }
  return;
}
/************************************************************/


