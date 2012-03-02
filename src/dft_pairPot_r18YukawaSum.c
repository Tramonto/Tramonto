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
 *  FILE: dft_pairPot_r18YukawaSum.c
 *
 *  This file contains routines specific to a strict mean field implementation of
 *  a summed repulsive and Yukawa fluid.  Note that the repulsive r^18 form was 
 *  chosen for a specific study of polyplex systems (Sachs/Braun UMN)
 *
 */

#include "dft_pairPot_r18YukawaSum.h"

/******************************************************************************/
/* ur18andYUKAWA_CS: The cut and shifted r18+YUKAWA potential                         */

double ur18andYUKAWA_CS(double r,double sigma, double eps, double rcut,double yukawaK,double AYukawa)
{
  double u,alpha;
  alpha=yukawaK*sigma;

  /* note that writing the potential this way implies that the Yukawa parameters in the
     input file will be entered as K*sigma_ij */
  
  if (r <= rcut) {
     u = eps*( POW_DOUBLE_INT(sigma/r,18) - POW_DOUBLE_INT(sigma/rcut,18) ) + 
         AYukawa*(exp(-alpha*(r/sigma-1.0))/(r/sigma)- exp(-alpha*(rcut/sigma-1.0))/(rcut/sigma));
  }
  else u = 0.0;
  return (u);
}
/*******************************************************************************/
/* ur18andYUKAWA_CS_setparams: The parameters for a potential that sums a r18 repulsive
core and a Yukawa potential. */
void ur18andYUKAWA_CS_setparams(int context, int i, int j, double *param1,double *param2,double *param3,double *param4,double *param5)
{    
  switch (context){
     case FLUID_FLUID:
        *param1 = Sigma_ff[i][j];
        *param2 = Eps_ff[i][j];
        *param3 = Cut_ff[i][j];
        *param4 = YukawaK_ff[i][j];
        *param5 = EpsYukawa_ff[i][j];
        break;
     case WALL_FLUID:
        *param1 = Sigma_wf[i][WallType[j]];
        *param2 = Eps_wf[i][WallType[j]];
        *param3 = Cut_wf[i][WallType[j]];
        *param4 = YukawaK_wf[i][WallType[j]];
        *param5 = EpsYukawa_wf[i][WallType[j]];
        break;
     case WALL_WALL:
        *param1 = Sigma_ww[WallType[i]][WallType[j]];
        *param2 = Eps_ww[WallType[i]][WallType[j]];
        *param3 = Cut_ww[WallType[i]][WallType[j]];
        *param4 = YukawaK_ww[WallType[i]][WallType[j]];
        *param5 = EpsYukawa_ww[WallType[i]][WallType[j]];
        break;
     default:
        if (Iwrite_screen != NONE) printf("problem with potential context ur18andYUKAWA_CS_setparams\n");
        exit(-1);
   }
   return;
}
/*******************************************************************************/
/* ur18andYUKAWA_DERIV1D: The derivative of the summed r18 and Yukawa potential in one cartesian direction x,y,z*/

double ur18andYUKAWA_DERIV1D(double r,double x,double sigma, double eps, double rcut,double yukawaK,double AYukawa)
{
  double uderiv,alpha;
  alpha=yukawaK*sigma;
  
  if (r <= rcut) {
     uderiv = (eps/(sigma*sigma)) * (-18.*x*POW_DOUBLE_INT(sigma/r,20) )
            -AYukawa*sigma*x*exp(-alpha*(r/sigma-1.0))*((1./r)+(alpha/sigma))/(r*r);
  }
  else uderiv = 0.0;
  return (uderiv);
}
/******************************************************************************/
/* ur18andYUKAWA_InnerCore : define the properties of the inner core of the potential based on
                  input parameters */
void ur18andYUKAWA_InnerCore(int i, int j,double *rCore_left, double *rCore_right, double *epsCore)
{
   switch(Type_CoreATT_R){
         /* note we have full flexibility --- if the potential is purely repulsive (monotonic),
            both Rmin_ff and Rzero_ff will be set to Sigma_ff */
      case ATTCORE_SIGMA:      *rCore_right=Sigma_ff[i][j]; *rCore_left=0.0; break;
      case ATTCORE_UMIN:       *rCore_right=Rmin_ff[i][j]; *rCore_left=0.0; break;
      case ATTCORE_UCSZERO:    *rCore_right=Rzero_ff[i][j]; *rCore_left=0.0; break;
      case ATTCORE_SIGTOUMIN:   *rCore_right=Rmin_ff[i][j]; *rCore_left=Sigma_ff[i][j]; break;
      default:
        if (Iwrite_screen != NONE) printf("Problem with Type_CoreATT_R - set to %d\n",Type_CoreATT_R);
        exit(-1);
   } 
   switch(Type_CoreATT_CONST){
      case CORECONST_UCONST:   *epsCore=ur18andYUKAWA_ATT_noCS(*rCore_right,i,j); break;
      case CORECONST_ZERO:     *epsCore=0.0; break;
      default:
        if (Iwrite_screen != NONE) printf("Problem with Type_CoreATT_CONST - set to %d\n",Type_CoreATT_CONST);
        exit(-1);
   }
   return;
}
/******************************************************************************/
/* ur18andYUKAWA_ATT_CS: the attractive part of the potential for a summed r18 and 
                          yukawa potential. */
double ur18andYUKAWA_ATT_CS(double r,int i, int j)
{
  double uatt,r_min,rcut,sigma,alpha,eps,Ayukawa;
  double sigma2,sigma6;
  double r_inv,r2_inv,r6_inv,r12_inv,r18_inv;
  double rc_inv,rc2_inv,rc6_inv,rc12_inv,rc18_inv;

  sigma=Sigma_ff[i][j];
  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;

  rc_inv   = 1.0/Cut_ff[i][j];
  rc2_inv  = rc_inv*rc_inv;
  rc6_inv  = rc2_inv*rc2_inv*rc2_inv;
  rc12_inv = rc6_inv*rc6_inv;
  rc18_inv = rc6_inv*rc12_inv;

  rcut=Cut_ff[i][j];
  eps=Eps_ff[i][j];
  alpha=YukawaK_ff[i][j]*sigma;
  Ayukawa=EpsYukawa_ff[i][j];

  /* note that Rmin and Rzero will both be Sigma_ff[i][j] for a monotonic potential */
  switch(Type_CoreATT_R){
     case ATTCORE_SIGMA:      r_min=sigma; break;
     case ATTCORE_SIGTOUMIN:
     case ATTCORE_UMIN:       r_min=Rmin_ff[i][j]; break; 
     case ATTCORE_UCSZERO:    r_min=Rzero_ff[i][j]; break; 
  }

  if ((r<r_min && Type_CoreATT_CONST==CORECONST_ZERO) ||
      (r<sigma && Type_CoreATT_R==ATTCORE_SIGTOUMIN))       uatt=0.0;
  else{
     if (r<=rcut){
        if (r<r_min) r=r_min;
   
        r_inv = 1.0/r;
        r2_inv  = r_inv*r_inv;
        r6_inv  = r2_inv*r2_inv*r2_inv;
        r12_inv = r6_inv*r6_inv;
        r18_inv = r6_inv*r12_inv;

        uatt= eps*sigma6*sigma6*sigma6*(r18_inv - rc18_inv)+
             Ayukawa*exp(-alpha*(r/sigma-1.0))/(r/sigma)
           - Ayukawa*exp(-alpha*(rcut/sigma-1.0))/(rcut/sigma);
     }
     else uatt=0.0;
  }
  return uatt;
}
/******************************************************************************/
/* ur18andYUKAWA_ATT_noCS: the attractive part of the potential for a 
       combined r18 and yukawa potential */
double ur18andYUKAWA_ATT_noCS(double r,int i, int j)
{
  double uatt,sigma,alpha,r_min,eps,Ayukawa;
  double sigma2,sigma6;
  double r_inv,r2_inv,r6_inv,r12_inv,r18_inv;

  sigma=Sigma_ff[i][j];
  eps=Eps_ff[i][j];
  alpha=YukawaK_ff[i][j]*sigma;
  Ayukawa=EpsYukawa_ff[i][j];
  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;

  /* note that Rmin and Rzero will both be Sigma_ff[i][j] for a monotonic potential */
  switch(Type_CoreATT_R){
     case ATTCORE_SIGMA:      r_min=sigma; break;
     case ATTCORE_SIGTOUMIN:
     case ATTCORE_UMIN:       r_min=Rmin_ff[i][j]; break; 
     case ATTCORE_UCSZERO:    r_min=Rzero_ff[i][j]; break; 
  }

  if ((r<r_min && Type_CoreATT_CONST==CORECONST_ZERO) ||
      (r<sigma && Type_CoreATT_R==ATTCORE_SIGTOUMIN))       uatt=0.0;
  else{
     if (r<r_min) r=r_min;
     r_inv = 1.0/r;

     r2_inv  = r_inv*r_inv;
     r6_inv  = r2_inv*r2_inv*r2_inv;
     r12_inv = r6_inv*r6_inv;
     r18_inv = r6_inv*r12_inv;

     uatt= eps*sigma6*sigma6*sigma6*r18_inv + Ayukawa*exp(-alpha*(r/sigma-1.0))/(r/sigma);
  }

  return uatt;
}
/****************************************************************************/
/* ur18andYUKAWA_IntStencil:  the integral of the summed r18 and Yukawa potential 
                        that is used to define the magnitude of the DFTMFT UATTRACT
                        integration stencil. */

double ur18andYUKAWA_Integral(double r,int i, int j)
{
  double uatt_int, sigma,c,alpha,eps,Ayukawa;
  double sigma2,sigma6;
  double r_inv,r3_inv,r9_inv,r15_inv;

  sigma=Sigma_ff[i][j];
  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;
  eps=Eps_ff[i][j];
  Ayukawa=EpsYukawa_ff[i][j];
  alpha=YukawaK_ff[i][j]*sigma;
  c=alpha/sigma;

  r_inv = 1.0/r;

  r3_inv  = r_inv*r_inv*r_inv;
  r9_inv  = r3_inv*r3_inv*r3_inv;
  r15_inv  = r9_inv*r3_inv*r3_inv;


  uatt_int = -4.*PI*eps*sigma6*sigma6*sigma6*r15_inv/15.0 +
             4*PI*Ayukawa*sigma*exp(alpha)*((exp(-c*r)/(c*c))*(-c*r-1.0));

  return uatt_int;
}
/****************************************************************************/

