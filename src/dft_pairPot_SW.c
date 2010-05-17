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

#include "dft_pairPot_SW.h"

/******************************************************************************/
/* uSW: The Square-Well potential                         */

double uSW(double r,double sigma, double eps, double rcut)
{
  double u;
  
  if (r <= rcut) {
	  u = -eps;
  }
  else u = 0.0;
  return (u);
}
/*******************************************************************************/
/* uSW_setparams: The parameters for the Square-Well potential */
void uSW_setparams(int context, int i, int j, double *param1,double *param2, double *param3)
{
  switch (context){
     case FLUID_FLUID:
        *param1 = Sigma_ff[i][j];
        *param2 = Eps_ff[i][j];
        *param3 = Cut_ff[i][j];
        break;
     case WALL_FLUID:
        *param1 = Sigma_wf[i][WallType[j]];
        *param2 = Eps_wf[i][WallType[j]];
        *param3 = Cut_wf[i][WallType[j]];
        break;
     case WALL_WALL:
        *param1 = Sigma_ww[WallType[i]][WallType[j]];
        *param2 = Eps_ww[WallType[i]][WallType[j]];
        *param3 = Cut_ww[WallType[i]][WallType[j]];
        break;
     default:
        printf("problem with potential context uSW_setparams\n");
        exit(-1);
   }
   return;
}
/*******************************************************************************/
/* uSW_DERIV1D: The derivative of a SW potential in the x (or y or z) direction
 Note that we don't take account of the discontinuity correctly, so this contribution to the 
 force won't be right */

double uSW_DERIV1D(double r,double x,double sigma, double eps, double rcut)
{
  double uderiv;
  
  if (r <= rcut) {
	  uderiv = 0.0;
  }
  else uderiv = 0.0;
  return (uderiv);
}
/******************************************************************************/
/* uSW_InnerCore : define the properties of the inner core of the potential based on
                  input parameters */
void uSW_InnerCore(int i, int j,double *rCore_left, double *rCore_right, double *epsCore)
{
   *rCore_left=0.0;
   *rCore_right=Sigma_ff[i][j];
   *epsCore=0.0;

   return;
}
/******************************************************************************/
/* uSW_ATT_CS: the attractive part of the potential for a cut and shifted 
			SW system */
double uSW_ATT_CS(double r,int i, int j)
{
	double uatt,r_min,rcut,sigma,alpha,eps;
	sigma=Sigma_ff[i][j];
	rcut=Cut_ff[i][j];
	eps=Eps_ff[i][j];
	
	if(r<sigma) uatt=0.0;
	else if (r<=rcut){
		uatt = -eps;
	}
	else uatt=0.0;
	return uatt;
}
/****************************************************************************/
/* uSW_ATT_noCS:  calculate the attractive part of
                  SW potential at the minimum. */

double uSW_ATT_noCS(double r,int i, int j)
{
	double uatt,sigma,alpha,r_min,eps;
	sigma=Sigma_ff[i][j];
	eps=Eps_ff[i][j];
	
	if (r<sigma) uatt = 0.0;
	else
		uatt = -eps;
	
	return uatt;
}
/****************************************************************************/
/* uSW_IntStencil:  the integral of the SW potential. Note that at the moment we don't actually use this at all */

double uSW_Integral(double r,int i, int j)
{
	double uatt_int,eps;
	
	eps=Eps_ff[i][j];	
	uatt_int = -(4.0/3.0)*PI*Eps_ff[i][j]*pow(r,3);
	
	return uatt_int;
}
/****************************************************************************/


