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
 *  FILE: dft_pairPot_COULOMB.c
 *
 *  This file contains routines specific to the implementation of Coulombic
 *  systems using a strict mean field approach.  
 */

#include "dft_pairPot_COULOMB.h"

/******************************************************************************/
/* uCOULOMB_CS: The Coulomb potential                   */

double uCOULOMB_CS(double r,double z1,double z2,double rcut)
{
  double u;

  u = z1*z2/(Temp_elec*r)-z1*z2/(Temp_elec*rcut);
  return (u);
}
/******************************************************************************/
/* uCOULOMB: The Coulomb potential                   */

double uCOULOMB(double r,double z1,double z2)
{
  double u;

  u = z1*z2/(Temp_elec*r);
  return (u);
}
/*******************************************************************************/
/* uCOULOMB_CS_setparams: The parameters for the cut and shifted Coulomb potential */
void uCOULOMB_CS_setparams(int context, int i, int j, double *param1,double *param2, double *param3)
{
  switch (context){
     case FLUID_FLUID:
        *param1 = Charge_f[i];
        *param2 = Charge_f[j];
        *param3 = Cut_ff[i][j];
        break;
     case WALL_FLUID:
        *param1 = Charge_f[i];
        *param2 = Elec_param_w[j];
        *param3 = Cut_wf[i][WallType[j]];
        break;
     case WALL_WALL:
        *param1 = Elec_param_w[i];
        *param2 = Elec_param_w[j];
        *param3 = Cut_ww[WallType[i]][WallType[j]];
        break;
     default:
        printf("problem with potential context uCOULOMB_CS_setparams\n");
        exit(-1);
   }
   return;
}
/*******************************************************************************/
/* uCOULOMB_setparams: The parameters for the full Coulomb potential */
void uCOULOMB_setparams(int context, int i, int j, double *param1,double *param2, double *param3)
{


  switch (context){
     case FLUID_FLUID:
        *param1 = Charge_f[i];
        *param2 = Charge_f[j];
        *param3 = 1.e6;
        break;
     case WALL_FLUID:
        *param1 = Charge_f[i];
        *param2 = Elec_param_w[j];
        *param3 = 1.e6;   
        break;
     case WALL_WALL:
        *param1 = Elec_param_w[i];
        *param2 = Elec_param_w[j];
        *param3 = 1.e6;
        break;
     default:
        printf("problem with potential context uCOULOMB_setparams\n");
        exit(-1);
   }
   return;
}
/******************************************************************************/
/* uCOULOMB_CS_DERIV1D: The derivative of the Coulomb potential in the x (or y or z)
                   direction....                                            */

double uCOULOMB_CS_DERIV1D(double r,double x,double z1,double z2,double rcut)
{
  double uderiv;

  if (r<rcut) uderiv = -z1*z2*x/(Temp_elec*r*r*r);
  else uderiv=0.0;

  return (uderiv);
}
/****************************************************************************/
/* uCOULOMB_InnerCore : define the properties of the inner core of the potential based on
                  input parameters */
void uCOULOMB_InnerCore(int i, int j,double *rCore_left, double *rCore_right, double *epsCore)
{
   *rCore_left=0.0;
   *rCore_right=Sigma_ff[i][j]; /* monotonic potential - UMIN & UZERO options not available */
   switch(Type_CoreATT_CONST){
      case CORECONST_UCONST:   *epsCore=uCOULOMB_ATT_noCS(*rCore_right,i,j);break;
      case CORECONST_ZERO:     *epsCore=0.0; break;
      default:
        printf("Problem with Type_CoreATT_CONST - set to %d\n",Type_CoreATT_CONST);
        exit(-1);
   }
   return;
}
/******************************************************************************/

/* uCOULOMB_CS_ATT_CS:  the pair potential (based on a Coulomb potential model fluid) 
                  that will be used as the attractive perturbation to a hard sphere 
                  reference fluid in strict mean field DFT calculations */

double uCOULOMB_ATT_CS(double r,int i, int j)
{
  double uatt,r_min,rcut;
  rcut=Cut_ff[i][j];
  r_min = Sigma_ff[i][j];

  if (r<r_min && Type_CoreATT_CONST==CORECONST_ZERO) uatt=0.0;
  else{
     if (r <= rcut) {
        if (r < r_min) r = r_min;
        uatt = Charge_f[i]*Charge_f[j]*(1./(r*Temp_elec)-1./(rcut*Temp_elec));
     }
  else uatt = 0.0;
  }

  return uatt;
}
/****************************************************************************/
/* uCOULOMB_ATT_noCS:  the pair potential (based on a Coulomb potential model fluid) 
                  that will be used as the attractive perturbation to a hard sphere 
                  reference fluid in strict mean field DFT calculations */

double uCOULOMB_ATT_noCS(double r,int i, int j)
{
  double uatt,r_min;

  r_min = Sigma_ff[i][j];

  if (r<r_min && Type_CoreATT_CONST==CORECONST_ZERO) uatt=0.0;
  else{
     if (r < r_min) r = r_min;
     uatt = Charge_f[i]*Charge_f[j]/(r*Temp_elec);
  }
  return uatt;
}
/******************************************************************************/
/* uCOULOMB_DERIV1D: The derivative of the Coulomb potential in the x (or y or z)
                   direction....                                            */

double uCOULOMB_DERIV1D(double r,double x,double z1,double z2)
{
  double uderiv;

  uderiv = -z1*z2*x/(Temp_elec*r*r*r);
  return (uderiv);
}
/****************************************************************************/
/* uCOULOMB_ATT_CnoS:  the pair potential (based on a Coulomb potential model fluid) 
                  that will be used as the attractive perturbation to a hard sphere 
                  reference fluid in strict mean field DFT calculations */

double uCOULOMB_ATT_CnoS(double r,int i, int j)
{
  double uatt,r_min;
  r_min = Sigma_ff[i][j];
  

  if (r<r_min && Type_CoreATT_CONST==CORECONST_ZERO) uatt=0.0;
  else{
    if (r <= Cut_ff[i][j]) {
       if (r < r_min) r = r_min;
       uatt = Charge_f[i]*Charge_f[j]/(r*Temp_elec);
    }
  else uatt = 0.0;
  }
  return uatt;
}
/****************************************************************************/
/* uCOULOMB_Integral:  the integral of the COULOMB potential to a distance r.
                         note that this is a non-converging function.  This lack of
                         convergence makes a strict mean field approach incorrect for
                         electrostatics.         */

double uCOULOMB_Integral(double r,int i, int j)
{
  double uCOUL_int;

  uCOUL_int = 4.*PI*Charge_f[i]*Charge_f[j]*r*r/Temp_elec;

  return uCOUL_int;
}
/****************************************************************************/
