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
 *  FILE: dft_pairPot_Exp.c
 *
 *  This file contains routines specific to a strict mean field implementation of
 *  an exponential attraction.
 *
 */

#include "dft_pairPot_EXP.h"

/******************************************************************************/
/* uEXP_CS: The cut and shifted exponential potential                         */

double uEXP_CS(double r,double sigma, double eps, double rcut,double yukawaK)
{
  double u,alpha;
  alpha=yukawaK;
  
  if (r <= rcut) {
     u = -eps*exp(-(r-sigma)/alpha) + eps*exp(-(rcut-sigma)/alpha);
  }
  else u = 0.0;
  return (u);
}
/*******************************************************************************/
/* uEXP_CS_setparams: The parameters for the cut and shifted exponential potential */
void uEXP_CS_setparams(int context, int i, int j, double *param1,double *param2,double *param3,double *param4)
{    
  switch (context){
     case FLUID_FLUID:
        *param1 = Sigma_ff[i][j];
        *param2 = Eps_ff[i][j];
        *param3 = Cut_ff[i][j];
        *param4 = YukawaK_ff[i][j];
        break;
     case WALL_FLUID:
        *param1 = Sigma_wf[i][WallType[j]];
        *param2 = Eps_wf[i][WallType[j]];
        *param3 = Cut_wf[i][WallType[j]];
        *param4 = YukawaK_wf[i][WallType[j]];
        break;
     case WALL_WALL:
        *param1 = Sigma_ww[WallType[i]][WallType[j]];
        *param2 = Eps_ww[WallType[i]][WallType[j]];
        *param3 = Cut_ww[WallType[i]][WallType[j]];
        *param4 = YukawaK_ww[WallType[i]][WallType[j]];
        break;
     default:
        printf("problem with potential context uYUKAWA_CS_setparams\n");
        exit(-1);
   }
   return;
}
/*******************************************************************************/
/* uEXP_DERIV1D: The derivative of the exponential potential in the x (or y or z) direction */

double uEXP_DERIV1D(double r,double x,double sigma, double eps, double rcut,double yukawaK)
{
  double uderiv,alpha;
  alpha=yukawaK;
  
  if (r <= rcut) {
     uderiv = (eps/alpha)*exp(-(r-sigma)/alpha);
  }
  else uderiv = 0.0;
  return (uderiv);
}
/******************************************************************************/
/* uEXP_InnerCore : define the properties of the inner core of the potential based on
                  input parameters */
void uEXP_InnerCore(int i, int j,double *rCore_left, double *rCore_right, double *epsCore)
{
   *rCore_left=0.0;
   *rCore_right=Sigma_ff[i][j];   /* monotonic potential - UMIN & UZERO options not available */

   switch(Type_CoreATT_CONST){
      case CORECONST_UCONST:   *epsCore=uEXP_ATT_noCS(*rCore_right,i,j); break;
      case CORECONST_ZERO:     *epsCore=0.0; break;
      default:
        printf("Problem with Type_CoreATT_CONST - set to %d\n",Type_CoreATT_CONST);
        exit(-1);
   }
   return;
}
/******************************************************************************/
/* uEXP_ATT_CS: the attractive part of the potential for a cut and shifted 
                   exponential system */
double uEXP_ATT_CS(double r,int i, int j)
{
  double uatt,rcut,sigma,alpha,eps;
  sigma=Sigma_ff[i][j];
  rcut=Cut_ff[i][j];
  eps=Eps_ff[i][j];
  alpha=YukawaK_ff[i][j];

  if(r<sigma){
     if (Type_CoreATT_CONST==CORECONST_ZERO) uatt = 0.0;
     else                                    uatt = eps*(-1.0+exp(-(rcut-sigma)/alpha));
  }
  else if (r<=rcut){
     uatt = -eps*exp(-(r-sigma)/alpha) + eps*exp(-(rcut-sigma)/alpha);
  }
  else uatt=0.0;
  return uatt;
}
/******************************************************************************/
/* uEXP_ATT_noCS: the attractive part of the potential for an exponential potential */
double uEXP_ATT_noCS(double r,int i, int j)
{
  double uatt,sigma,alpha,eps;
  sigma=Sigma_ff[i][j];
  eps=Eps_ff[i][j];
  alpha=YukawaK_ff[i][j];

  if (r<sigma){
       if (Type_CoreATT_CONST==CORECONST_ZERO) uatt = 0.0;
       else                                    uatt = -eps;
  }
  else   uatt = -eps*exp(-(r-sigma)/alpha);

  return uatt;
}
/****************************************************************************/
/* uEXP_IntStencil:  the integral of the exponential potential that is used
                        to define the magnitude of the DFTMFT UATTRACT
                        integration stencil. */

double uEXP_Integral(double r,int i, int j)
{
  double uatt_int, sigma,a3,alpha,eps;

  sigma=Sigma_ff[i][j];
  eps=Eps_ff[i][j];
  alpha=YukawaK_ff[i][j];
  a3 = alpha*alpha*alpha;

  uatt_int = -4.0*PI*eps*exp(-(r-sigma)/alpha)*(-alpha*r*r + 2*a3*(-r/alpha - 1.));

  return uatt_int;
}
/****************************************************************************/

