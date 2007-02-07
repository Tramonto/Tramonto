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
 *  FILE: dft_pairPot_LJ12_6.c
 *
 *  This file contains routines specific to a strict mean field implementation of
 *  a Lennard-Jones fluid with a Weeks-Chandler-Anderson split of the LJ potential
 *  into an infinite hard core and an attractive piece.  Note that the hard core
 *  diameters can be modified by the Barker-Henderson approach if desired.
 */

#include "dft_pairPot_LJ12_6.h"

/******************************************************************************/
/* uLJ12_6_CS: The cut and shifted LJ 12-6 potential                         */

/*double uLJ12_6_cut(double r,double sigma, double eps, double rcut)*/
double uLJ12_6_CS(double r,double sigma, double eps, double rcut)
{
  double u;
  
  if (r <= rcut) {
     u = ( POW_DOUBLE_INT(sigma/r,12) - POW_DOUBLE_INT(sigma/rcut,12) ) -
            ( POW_DOUBLE_INT(sigma/r,6) - POW_DOUBLE_INT(sigma/rcut,6) );
  }
  else u = 0.0;
  return (4.0*eps*u);
}
/*******************************************************************************/
/* uLJ12_6_DERIV1D: The derivative of a 12-6 potential in the x (or y or z) direction */

/*double uderiv_LJ12_6(double r,double x,double sigma, double eps, double rcut)*/
double uLJ12_6_DERIV1D(double r,double x,double sigma, double eps, double rcut)
{
  double uderiv;
  
  if (r <= rcut) {
     uderiv = (1.0/(sigma*sigma)) * (
            -12.*x*POW_DOUBLE_INT(sigma/r,14) + 6.*x*POW_DOUBLE_INT(sigma/r,8) );
  }
  else uderiv = 0.0;
  return (4.0*eps*uderiv);
}
/******************************************************************************/
/* uLJ12_6_ATT_CS: the pair potential (based on a 12-6 LJ fluid) that will be used 
                  as the attractive perturbation to a hard sphere reference fluid
                  in strict mean field DFT calculations */
/*double uLJatt_n(double r,int i, int j)*/
double uLJ12_6_ATT_CS(double r,int i, int j)
{
  double uatt,r_min,sigma2,sigma6;
  double r_inv,r2_inv,r6_inv,r12_inv;
  double rc_inv,rc2_inv,rc6_inv,rc12_inv;

  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;

  rc_inv   = 1.0/Cut_ff[i][j];
  rc2_inv  = rc_inv*rc_inv;
  rc6_inv  = rc2_inv*rc2_inv*rc2_inv;
  rc12_inv = rc6_inv*rc6_inv;


  if (r <= Cut_ff[i][j]) {
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
/* uLJ12_6_ATT_noCS:  calculate the attractive part of
                  12-6 LJ potential at the minimum. */

/*double uLJatt_n_noshift(double r,int i, int j)*/
double uLJ12_6_ATT_noCS(double r,int i, int j)
{
  double uatt,r_min,sigma2,sigma6;
  double r_inv,r2_inv,r6_inv,r12_inv;

  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;

  r_min = Sigma_ff[i][j] * pow(2.0,1.0/6.0);
  if (r < r_min) r = r_min;

  r_inv = 1.0/r;

  r2_inv  = r_inv*r_inv;
  r6_inv  = r2_inv*r2_inv*r2_inv;
  r12_inv = r6_inv*r6_inv;

  uatt = 4.0 * Eps_ff[i][j]* sigma6 * ( sigma6*r12_inv  - r6_inv);

  return uatt;
}
/****************************************************************************/
/* uLJ12_6_IntStencil:  the integral of the 12-6 potential that is used
                        to define the magnitude of the DFTMFT UATTRACT
                        integration stencil. */

/*double uLJatt_n_int(double r,int i, int j)*/
double uLJ12_6_Integral(double r,int i, int j)
{
  double uatt_int, sigma2,sigma6;
  double r_inv,r3_inv,r9_inv;

  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;

  r_inv = 1.0/r;

  r3_inv  = r_inv*r_inv*r_inv;
  r9_inv  = r3_inv*r3_inv*r3_inv;

  uatt_int = 16 * PI * Eps_ff[i][j]* sigma6 * ( - sigma6*r9_inv/9.0  + r3_inv/3.0 );

  return uatt_int;
}
/****************************************************************************/
