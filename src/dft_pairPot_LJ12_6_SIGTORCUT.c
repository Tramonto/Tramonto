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
 *  FILE: dft_pairPot_LJ12_6_SIGTORCUT.c
 *
 *  This file contains routines specific to a strict mean field implementation of
 *  a Lennard-Jones fluid with a Weeks-Chandler-Anderson split of the LJ potential
 *  into an infinite hard core and an attractive piece.  Note that the hard core
 *  diameters can be modified by the Barker-Henderson approach if desired.
 */

#include "dft_pairPot_LJ12_6_SIGTORCUT.h"

/******************************************************************************/
/* uLJ12_6_ATT_SIGTORCUT_CS: the pair potential (based on a 12-6 LJ fluid) that will be used 
                  as the attractive perturbation to a hard sphere reference fluid
                  in strict mean field DFT calculations */
double uLJ12_6_ATT_SIGTORCUT_CS(double r,int i, int j)
{
  double uatt,sigma2,sigma6;
  double r_inv,r2_inv,r6_inv,r12_inv,r_min;
  double rc_inv,rc2_inv,rc6_inv,rc12_inv;

  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;

  rc_inv   = 1.0/Cut_ff[i][j];
  rc2_inv  = rc_inv*rc_inv;
  rc6_inv  = rc2_inv*rc2_inv*rc2_inv;
  rc12_inv = rc6_inv*rc6_inv;


  if (r<Sigma_ff[i][j]){
     uatt=0.0;
  }
  else if (r <= Cut_ff[i][j]) {
     r_min = Sigma_ff[i][j] * pow(2.0,1.0/6.0);
     if (r < r_min) r = r_min;

     r_inv = 1.0/r;
     r2_inv  = r_inv*r_inv;
     r6_inv  = r2_inv*r2_inv*r2_inv;
     r12_inv = r6_inv*r6_inv;

     uatt = 4.0 * Eps_ff[i][j]* sigma6 * (
               sigma6*(r12_inv - rc12_inv)
                    - (r6_inv  - rc6_inv ) );
  }
  else uatt = 0.0;

  return uatt;
}
/****************************************************************************/
/* uLJ12_6_ATT_SIGTORCUT_noCS:  calculate the attractive part of
                  12-6 LJ potential at the minimum. */

double uLJ12_6_ATT_SIGTORCUT_noCS(double r,int i, int j)
{
  double uatt,r_min,sigma2,sigma6;
  double r_inv,r2_inv,r6_inv,r12_inv;

  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;

  if (r<Sigma_ff[i][j]) uatt=0.0;
  else{
     r_min = Sigma_ff[i][j] * pow(2.0,1.0/6.0);
     if (r < r_min) r = r_min;

     r_inv = 1.0/r;
     r2_inv  = r_inv*r_inv;
     r6_inv  = r2_inv*r2_inv*r2_inv;
     r12_inv = r6_inv*r6_inv;
     uatt = 4.0 * Eps_ff[i][j]* sigma6 * ( sigma6*r12_inv  - r6_inv);
  }

  return uatt;
}
/****************************************************************************/

