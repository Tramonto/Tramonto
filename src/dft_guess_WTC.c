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

/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

/*
 *  FILE: dft_guess.c
 *
 *  This file contains routines that allocate and pick an initial
 *  guess for the solution vector.
 *
 */

#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"
#include <string.h>
 
/************************************************************/
/*setup_Xi_cavWTC: set up the cavity function initial guesses for 
                 Wertheim-Tripathi-Chapman functionals.  For now
                use bulk segment densities initial guess.  Later calculate
                based on rho initial guess. */
void setup_Xi_cavWTC(double **xOwned)
{
  int loc_inode,inode_box,inode,ijk[3],iunk,icav;
  double vol,area,x_dist;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     for (icav = 0; icav < Phys2Nunk[CAVWTC]; icav++){
       iunk = Phys2Unk_first[CAVWTC] + icav;
       xOwned[iunk][loc_inode] = Xi_cav_b[icav+2];
     }
  }
  return;
}
/************************************************************/
/*setup_BondWTC: set up the bond function initial guesses for 
                 Wertheim-Tripathi-Chapman functionals.  For now
                use bulk segment densities initial guess.  Later calculate
                based on rho initial guess. */
void setup_BondWTC(double **xOwned)
{
  int loc_inode,inode_box,inode,ijk[3],iunk,ibond;
  double vol,area,x_dist;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     for (ibond = 0; ibond < Nbonds; ibond++){
       iunk = Phys2Unk_first[BONDWTC] + ibond;
       xOwned[iunk][loc_inode] = BondWTC_b[ibond];
     }
  }
  return;
}
/************************************************************/
