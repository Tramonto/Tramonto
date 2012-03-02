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
 *  FILE: dft_switch_pairPot.c
 *
 *  This file contains routines selects the correct pair interaction potential
 *  functions to call for computing external fields based on pair potentials and
 *  for performing mean field DFT calculations.
 *
 */

#include "dft_switch_pairPot.h"

/******************************************************************************/
/* pairPot_switch:  switch to choose the correct pair potential to use in an
           external field calculation when using integrated potential surfaces
           or atomistic surface. This routine is also used when computing
           Barker-Henderson diameters */
double pairPot_switch(double r,double param1, double param2, double param3,double param4,double param5,double param6,int typePairPot)
{
  double u;

  switch(typePairPot){
      case PAIR_LJ12_6_CS:
        u= uLJ12_6_CS(r,param1,param2,param3);
        break;
      case PAIR_COULOMB_CS:
        u = uCOULOMB_CS(r,param1,param2,param3);
        break;
      case PAIR_COULOMB:
        u = uCOULOMB(r,param1,param2);
        break;
      case PAIR_YUKAWA_CS:
        u = uYUKAWA_CS(r,param1,param2,param3,param4);
        break;
      case PAIR_LJandYUKAWA_CS:
        u = uLJandYUKAWA_CS(r,param1,param2,param3,param4,param5);
        break;
      case PAIR_r12andYUKAWA_CS:
        u = ur12andYUKAWA_CS(r,param1,param2,param3,param4,param5);
        break;
      case PAIR_r18andYUKAWA_CS:
        u = ur18andYUKAWA_CS(r,param1,param2,param3,param4,param5);
        break;
      case PAIR_rNandYUKAWA_CS:
        u = urNandYUKAWA_CS(r,param1,param2,param3,param4,param5,param6);
        break;
      case PAIR_EXP_CS:
        u = uEXP_CS(r,param1,param2,param3,param4);
        break;
      case PAIR_SW:
        u = uSW(r,param1,param2,param3);
        break;
      default:
         if (Iwrite_screen !=SCREEN_NONE)printf("problems with your selection of typePairPotin pairPot_switch: typePairPot=%d\n",typePairPot);
         exit(-1);
         break;
  }
  return u;
}
/******************************************************************************/
/* pairPotparams_switch:  switch to set the correct parameters for a given choice of potential.
           Note that these parameters must map correctly to the potential functions (i.e. uLJ12_6_CS). */
void pairPotparams_switch(int typePairPot,int context, int i, int j,
      double *param1, double *param2, double *param3,double *param4,double *param5,double *param6)

{
  switch(typePairPot){
      case PAIR_LJ12_6_CS:
        uLJ12_6_CS_setparams(context,i,j,param1,param2,param3);
        break;
      case PAIR_COULOMB_CS:
        uCOULOMB_CS_setparams(context,i,j,param1,param2,param3);
        break;
      case PAIR_COULOMB:
        uCOULOMB_setparams(context,i,j,param1,param2,param3);
        break;
      case PAIR_YUKAWA_CS:
        uYUKAWA_CS_setparams(context,i,j,param1,param2,param3,param4);
        break;
      case PAIR_LJandYUKAWA_CS:
        uLJandYUKAWA_CS_setparams(context,i,j,param1,param2,param3,param4,param5);
        break;
      case PAIR_r12andYUKAWA_CS:
        ur12andYUKAWA_CS_setparams(context,i,j,param1,param2,param3,param4,param5);
        break;
      case PAIR_r18andYUKAWA_CS:
        ur18andYUKAWA_CS_setparams(context,i,j,param1,param2,param3,param4,param5);
        break;
      case PAIR_rNandYUKAWA_CS:
        urNandYUKAWA_CS_setparams(context,i,j,param1,param2,param3,param4,param5,param6);
        break;
      case PAIR_EXP_CS:
        uEXP_CS_setparams(context,i,j,param1,param2,param3,param4);
        break;
      case PAIR_SW:
        uSW_setparams(context,i,j,param1,param2,param3);
        break;
      default:
        if (Iwrite_screen !=SCREEN_NONE) printf("problems with your selection of typePairPot in pairPotparams_switch typePairPot=%d\n",typePairPot);
        exit(-1);
        break;
  }
  return;
}
/******************************************************************************/
/* pairPot_deriv_switch:  switch to choose the correct pair potential derivative
           needed for force calculations */
double pairPot_deriv_switch(double r, double x, double param1, double param2, double param3,
                            double param4,double param5,double param6,int typePairPot)
{
  double uderiv;

  switch(typePairPot){
      case PAIR_LJ12_6_CS:
        uderiv= uLJ12_6_DERIV1D(r,x,param1,param2,param3);
        break;
      case PAIR_COULOMB_CS:
        uderiv = uCOULOMB_CS_DERIV1D(r,x,param1,param2,param3);
        break;
      case PAIR_COULOMB:
        uderiv = uCOULOMB_DERIV1D(r,x,param1,param2);
        break;
      case PAIR_YUKAWA_CS:
        uderiv = uYUKAWA_DERIV1D(r,x,param1,param2,param3,param4);
        break;
      case PAIR_LJandYUKAWA_CS:
        uderiv = uLJandYUKAWA_DERIV1D(r,x,param1,param2,param3,param4,param5);
        break;
      case PAIR_r12andYUKAWA_CS:
        uderiv = ur12andYUKAWA_DERIV1D(r,x,param1,param2,param3,param4,param5);
        break;
      case PAIR_r18andYUKAWA_CS:
        uderiv = ur18andYUKAWA_DERIV1D(r,x,param1,param2,param3,param4,param5);
        break;
      case PAIR_rNandYUKAWA_CS:
        uderiv = urNandYUKAWA_DERIV1D(r,x,param1,param2,param3,param4,param5,param6);
        break;
      case PAIR_EXP_CS:
        uderiv = uEXP_DERIV1D(r,x,param1,param2,param3,param4);
        break;
      case PAIR_SW:
        uderiv = uSW_DERIV1D(r,x,param1,param2,param3);
        break;
      default:
         if (Iwrite_screen !=SCREEN_NONE) printf("problems with your selection of typePairPot in pairPot_deriv_switch typePairPot=%d\n",typePairPot);
         exit(-1);
         break;
  }
  return uderiv;
}
/******************************************************************************/
/* pairPot_InnerCore_switch:  switch to choose the correct properties of the 
           inner core of attractions as used in the DFT calculation.       */
void pairPot_InnerCore_switch(int icomp, int jcomp,int typePairPot, 
              double *rCore_left, double *rCore_right, double *epsCore)
{
  switch(typePairPot){
      case PAIR_LJ12_6_CS:
        uLJ12_6_InnerCore(icomp,jcomp,rCore_left,rCore_right,epsCore);
        break;
      case PAIR_COULOMB:
      case PAIR_COULOMB_CS:
        uCOULOMB_InnerCore(icomp,jcomp,rCore_left,rCore_right,epsCore);
        break;
      case PAIR_YUKAWA_CS:
         uYUKAWA_InnerCore(icomp,jcomp,rCore_left,rCore_right,epsCore);
         break;
      case PAIR_LJandYUKAWA_CS:
         uLJandYUKAWA_InnerCore(icomp,jcomp,rCore_left,rCore_right,epsCore);
         break;
      case PAIR_r12andYUKAWA_CS:
         ur12andYUKAWA_InnerCore(icomp,jcomp,rCore_left,rCore_right,epsCore);
         break;
      case PAIR_r18andYUKAWA_CS:
         ur18andYUKAWA_InnerCore(icomp,jcomp,rCore_left,rCore_right,epsCore);
         break;
      case PAIR_rNandYUKAWA_CS:
         urNandYUKAWA_InnerCore(icomp,jcomp,rCore_left,rCore_right,epsCore);
         break;
      case PAIR_EXP_CS:
         uEXP_InnerCore(icomp,jcomp,rCore_left,rCore_right,epsCore);
         break;
      case PAIR_SW:
         uSW_InnerCore(icomp,jcomp,rCore_left,rCore_right,epsCore);
         break;
      default:
         if (Iwrite_screen !=SCREEN_NONE) printf("problems with your selection of typePairPot in pairPot_innerCore_switch typePairPot=%d\n",typePairPot);
         exit(-1);
         break;
  }
  return;
}
/******************************************************************************/
/* pairPot_ATT_CS_switch:  switch to choose the correct pair potential for
           mean field perturbation DFT calculations with a hard sphere
           reference fluid.  This is the core function used to compute integration
           stencils.  It is also used in computing bulk thermodynamics for these terms.*/
double pairPot_ATT_CS_switch(double r, int icomp, int jcomp,int typePairPot)
{
  double u;

  switch(typePairPot){
      case PAIR_LJ12_6_CS:
        u= uLJ12_6_ATT_CS(r,icomp,jcomp);
        break;
      case PAIR_COULOMB_CS:
        u = uCOULOMB_ATT_CS(r,icomp,jcomp);  
        break;
      case PAIR_COULOMB:
        u = uCOULOMB_ATT_CnoS(r,icomp,jcomp);  
        break;
      case PAIR_YUKAWA_CS:
        u = uYUKAWA_ATT_CS(r,icomp,jcomp);  
        break;
      case PAIR_LJandYUKAWA_CS:
        u = uLJandYUKAWA_ATT_CS(r,icomp,jcomp);  
        break;
      case PAIR_r12andYUKAWA_CS:
        u = ur12andYUKAWA_ATT_CS(r,icomp,jcomp);  
        break;
      case PAIR_r18andYUKAWA_CS:
        u = ur18andYUKAWA_ATT_CS(r,icomp,jcomp);  
        break;
      case PAIR_rNandYUKAWA_CS:
        u = urNandYUKAWA_ATT_CS(r,icomp,jcomp);  
        break;
      case PAIR_EXP_CS:
        u = uEXP_ATT_CS(r,icomp,jcomp);  
        break;
      case PAIR_SW:
        u = uSW_ATT_CS(r,icomp,jcomp);  
        break;
      default:
         if (Iwrite_screen !=SCREEN_NONE) printf("problems with your selection of typePairPot in pairPot_ATT_CS_switch typePairPot=%d\n",typePairPot);
         exit(-1);
         break;
  }
  return u;
}
/******************************************************************************/
/* pairPot_ATT_noCS_switch:  switch to choose the correct full pair 
          potential (no cut and shift) used in setting up integration stencils 
          for mean field DFT calculations. */
double pairPot_ATT_noCS_switch(double r, int icomp, int jcomp,int typePairPot)
{
  double u;

  switch(typePairPot){
      case PAIR_LJ12_6_CS:
        u= uLJ12_6_ATT_noCS(r,icomp,jcomp);
        break;
      case PAIR_COULOMB:
      case PAIR_COULOMB_CS:
        u = uCOULOMB_ATT_noCS(r,icomp,jcomp);
        break;
      case PAIR_YUKAWA_CS:
        u = uYUKAWA_ATT_noCS(r,icomp,jcomp);
        break;
      case PAIR_LJandYUKAWA_CS:
        u = uLJandYUKAWA_ATT_noCS(r,icomp,jcomp);
        break;
      case PAIR_r12andYUKAWA_CS:
        u = ur12andYUKAWA_ATT_noCS(r,icomp,jcomp);
        break;
      case PAIR_r18andYUKAWA_CS:
        u = ur18andYUKAWA_ATT_noCS(r,icomp,jcomp);
        break;
      case PAIR_rNandYUKAWA_CS:
        u = urNandYUKAWA_ATT_noCS(r,icomp,jcomp);
        break;
      case PAIR_EXP_CS:
        u = uEXP_ATT_noCS(r,icomp,jcomp);
        break;
      case PAIR_SW:
        u = uSW_ATT_noCS(r,icomp,jcomp);
        break;	  
      default:
         if (Iwrite_screen !=SCREEN_NONE) printf("problems with your selection of typePairPot in pairPot_ATT_noCS_switch typePairPot=%d\n",typePairPot);
         exit(-1);
         break;
  }
  return u;
}
/******************************************************************************/
/* pairPot_integral_switch:  switch to choose the correct integrated MFDFT 
           pair potential (e.g. only the attractive part of the 12-6 potential)
           for setting up the integration stencils for DFT calculations. */
double pairPot_integral_switch(double r, int icomp, int jcomp,int typePairPot)
{
  double u;

  switch(typePairPot){
      case PAIR_LJ12_6_CS:
        u= uLJ12_6_Integral(r,icomp,jcomp);
        break;
      case PAIR_COULOMB:
      case PAIR_COULOMB_CS:
        u = uCOULOMB_Integral(r,icomp,jcomp);
        break;
      case PAIR_YUKAWA_CS:
         u = uYUKAWA_Integral(r,icomp,jcomp);
         break;
      case PAIR_LJandYUKAWA_CS:
         u = uLJandYUKAWA_Integral(r,icomp,jcomp);
         break;
      case PAIR_r12andYUKAWA_CS:
         u = ur12andYUKAWA_Integral(r,icomp,jcomp);
         break;
      case PAIR_r18andYUKAWA_CS:
         u = ur18andYUKAWA_Integral(r,icomp,jcomp);
         break;
      case PAIR_rNandYUKAWA_CS:
         u = urNandYUKAWA_Integral(r,icomp,jcomp);
         break;
      case PAIR_EXP_CS:
         u = uEXP_Integral(r,icomp,jcomp);
         break;
      case PAIR_SW:
         u = uSW_Integral(r,icomp,jcomp);
         break;
      default:
         if (Iwrite_screen !=SCREEN_NONE) printf("problems with your selection of typePairPot in pairPot_integral_switch typePairPot=%d\n",typePairPot);
         exit(-1);
         break;
  }
  return u;
}
/******************************************************************************/
