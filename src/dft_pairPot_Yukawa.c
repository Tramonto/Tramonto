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
 *  FILE: dft_pairPot_Yukawa.c
 *
 *  This file contains routines specific to a strict mean field implementation of
 *  a Yukawa fluid .
 *
 */

#include "dft_pairPot_Yukawa.h"

/******************************************************************************/
/* uYUKAWA_CS: The cut and shifted YUKAWA potential                         */

double uYUKAWA_CS(double r,double sigma, double eps, double rcut,double yukawaK)
{
  double u,alpha;
  alpha=yukawaK*sigma;

  /* note that writing the potential this way implies that the Yukawa parameters in the
     input file will be entered as K*sigma_ij */
  
  if (r <= rcut) {
     u = eps*(exp(-alpha*(r/sigma-1.))/(r/sigma)- exp(-alpha*(rcut/sigma-1.))/(rcut/sigma));
  }
  else u = 0.0;
  return (u);
}
/*******************************************************************************/
/* uYUKAWA_CS_setparams: The parameters for the cut and shifted LJ12_6 potential */
void uYUKAWA_CS_setparams(int context, int i, int j, double *param1,double *param2,double *param3,double *param4)
{    
  switch (context){
     case FLUID_FLUID:
        *param1 = Sigma_ff[i][j];
        *param2 = EpsYukawa_ff[i][j];
        *param3 = Cut_ff[i][j];
        *param4 = YukawaK_ff[i][j];
        break;
     case WALL_FLUID:
        *param1 = Sigma_wf[i][WallType[j]];
        *param2 = EpsYukawa_wf[i][WallType[j]];
        *param3 = Cut_wf[i][WallType[j]];
        *param4 = YukawaK_wf[i][WallType[j]];
        break;
     case WALL_WALL:
        *param1 = Sigma_ww[WallType[i]][WallType[j]];
        *param2 = EpsYukawa_ww[WallType[i]][WallType[j]];
        *param3 = Cut_ww[WallType[i]][WallType[j]];
        *param4 = YukawaK_ww[WallType[i]][WallType[j]];
        break;
     default:
        if (Iwrite_screen != NONE) printf("problem with potential context uYUKAWA_CS_setparams\n");
        exit(-1);
   }
   return;
}
/*******************************************************************************/
/* uYUKAWA_DERIV1D: The derivative of the Yukawa potential in the x (or y or z) direction */

double uYUKAWA_DERIV1D(double r,double x,double sigma, double eps, double rcut,double yukawaK)
{
  double uderiv,alpha;
  alpha=yukawaK*sigma;
  
  if (r <= rcut) {
     uderiv = -eps*sigma*x*exp(-alpha*(r/sigma-1.))*((1./r)+(alpha/(sigma)))/(r*r);
  }
  else uderiv = 0.0;
  return (uderiv);
}
/******************************************************************************/
/* uYUKAWA_InnerCore : define the properties of the inner core of the potential based on
                  input parameters */
void uYUKAWA_InnerCore(int i, int j,double *rCore_left, double *rCore_right, double *epsCore)
{
   *rCore_left=0.0;
   *rCore_right=Sigma_ff[i][j];  /* monotonic potential UMIN and UZERO options not available */

   switch(Type_CoreATT_CONST){
      case CORECONST_UCONST:   *epsCore=uYUKAWA_ATT_noCS(*rCore_right,i,j); break;
      case CORECONST_ZERO:     *epsCore=0.0; break;
      default:
        if (Iwrite_screen != NONE) printf("Problem with Type_CoreATT_CONST - set to %d\n",Type_CoreATT_CONST);
        exit(-1);
   }
   return;
}
/******************************************************************************/
/* uYUKAWA_ATT_CS: the attractive part of the potential for a cut and shifted 
                   yukawa system */
double uYUKAWA_ATT_CS(double r,int i, int j)
{
  double uatt,r_min,rcut,sigma,alpha,eps;
  sigma=Sigma_ff[i][j];
  rcut=Cut_ff[i][j];
  eps=Eps_ff[i][j];
  alpha=YukawaK_ff[i][j]*sigma;
  
  r_min=sigma;
  if (r<r_min){
     if (Type_CoreATT_CONST==CORECONST_ZERO) uatt=0.0;
     else uatt=eps*(1.0-exp(-alpha*(rcut/sigma-1.))/(rcut/sigma));
  }
  else if (r<=rcut){
     uatt=eps*exp(-alpha*(r/sigma-1.))/(r/sigma)
        - eps*exp(-alpha*(rcut/sigma-1.))/(rcut/sigma);
  }
  else uatt=0.0;

  return uatt;
}
/******************************************************************************/
/* uYUKAWA_ATT_noCS: the attractive part of the potential for a yukawa potential */
double uYUKAWA_ATT_noCS(double r,int i, int j)
{
  double uatt,sigma,alpha,r_min,eps;
  sigma=Sigma_ff[i][j];
  eps=Eps_ff[i][j];
  alpha=YukawaK_ff[i][j]*sigma;

  r_min=sigma;
  if (r<r_min){
    if (Type_CoreATT_CONST==CORECONST_ZERO) uatt=0.0;
    else                                    uatt=eps;
  } 
  else   uatt=eps*exp(-alpha*(r/sigma-1.))/(r/sigma);

  return uatt;
}
/****************************************************************************/
/* uYUKAWA_IntStencil:  the integral of the Yukawa potential that is used
                        to define the magnitude of the DFTMFT UATTRACT
                        integration stencil. */

double uYUKAWA_Integral(double r,int i, int j)
{
  double uatt_int, sigma,c,alpha,eps;

  sigma=Sigma_ff[i][j];
  eps=Eps_ff[i][j];
  alpha=YukawaK_ff[i][j]*sigma;
  c=alpha/sigma;

  uatt_int = 4*PI*eps*sigma*exp(alpha)*((exp(-c*r)/(c*c))*(-c*r-1.0));

  return uatt_int;
}
/****************************************************************************/

