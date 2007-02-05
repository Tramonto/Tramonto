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
 *  FILE: dft_switch_stencil.c
 *
 *  This file contains switches to select among varous stencils.  These switches
 *  are used in the course of setting up the various integration stencils.
 *
 */

#include "dft_switch_stencil.h"
/****************************************************************************/
/*stencil_Njcomp_switch: Logic to select correct type of stencil.  Typically,
  stencils depend on either one component or two.  Higher order stencils are 
  not currently implemented in Tramonto */
int stencil_Njcomp_switch(int sten)
{ 
  int njcomp;

  switch(sten){
     case DELTA_FN_R:
         njcomp=StenDelta_R_Njcomp();
         break;
     case THETA_FN_R:
         njcomp=StenTheta_R_Njcomp();
         break;
     case THETA_FN_SIG:
         njcomp=StenTheta_Sigma_Njcomp();
         break;
     case THETA_PAIRPOT_RCUT:
         njcomp=StenTheta_uattr_Njcomp();
         break;
     case THETA_CR_RPM_MSA:
         njcomp=StenTheta_RPMmsa_Njcomp();
         break;
     case THETA_CR_DATA:
         njcomp=StenTheta_CrCMS_Njcomp();
         break;
     case DELTA_FN_BOND:
         njcomp=StenDelta_Bond_Njcomp();
         break;
     default:
         printf("problem with stencil definitions: stencil_Njcomp_switch ");
         exit(-1); break;
  }
  return(njcomp);
}
/****************************************************************************/
/*stencil_radius_switch: Logic to set the range of the stencil of interest. */
double stencil_radius_switch(int sten,int icomp,int jcomp)
{ 
  double sten_rad;

  switch(sten){
     case DELTA_FN_R: 
         sten_rad=StenDelta_R_sten_rad(icomp);
         break;
     case THETA_FN_R:
         sten_rad=StenTheta_R_sten_rad(icomp);
         break;
     case THETA_FN_SIG:
         sten_rad=StenTheta_Sigma_sten_rad(icomp);
         break;
     case THETA_PAIRPOT_RCUT:
         sten_rad=StenTheta_uattr_sten_rad(icomp,jcomp);
         break;
     case THETA_CR_RPM_MSA:
         sten_rad=StenTheta_RPMmsa_sten_rad(icomp); /* approximation is for a one component system only */
         break;
     case THETA_CR_DATA:
         sten_rad=StenTheta_CrCMS_sten_rad(icomp,jcomp);
         break;
     case DELTA_FN_BOND:
         sten_rad=StenDelta_Bond_sten_rad(icomp,jcomp);
         break;
     default:
         printf("problem with stencil definitions: stencil_radius_switch ");
         exit(-1); break;
  }
  return(sten_rad);
}
/****************************************************************************/
/*stencil_volume_switch: Logic to set the analytical volume of the stencil.  This
  value will be used for renormalization of the numerically generated stencil.  */
double stencil_volume_switch(int sten,int icomp,int jcomp)
{ 
  double sten_vol;

  switch(sten){
     case DELTA_FN_R:
         sten_vol=StenDelta_R_sten_vol(icomp);
         break;
     case THETA_FN_R:
         sten_vol=StenTheta_R_sten_vol(icomp);
         break;
     case THETA_FN_SIG:
         sten_vol=StenTheta_Sigma_sten_vol(icomp);
         break;
     case THETA_PAIRPOT_RCUT:
         sten_vol=StenTheta_uattr_sten_vol(icomp,jcomp);
         break;
     case THETA_CR_RPM_MSA:
         sten_vol=StenTheta_RPMmsa_sten_vol(icomp,jcomp); /* approximation is for a one component system only */
         break;
     case THETA_CR_DATA:
         sten_vol=StenTheta_CrCMS_sten_vol(icomp,jcomp);
         break;
     case DELTA_FN_BOND:
         sten_vol=StenDelta_Bond_sten_vol(icomp,jcomp);
         break;
     default:
         printf("problem with stencil definitions: stencil_volume_switch ");
         exit(-1); break;
  }
  return(sten_vol);
}
/****************************************************************************/
/*stencil_GetWeight_switch: Logic to set stencil weights at various distances from
  the origin.  Note that for some stencils, this weight is computed analytically, for others
  a gaussian quadrature is applied with analytical integrands, for others a numerical integration
  is required along with the gaussian quadrature to compute the weight.  This accounts for the
  different variables passed through to the various GetWeightFromSten functions.*/
double stencil_GetWeight_switch(int sten, int icomp, int jcomp, double rsq,
                                 double sten_rad, int ngpu, double *gpu, double *gwu)
{ 
  double weight;

  switch(sten){
     case DELTA_FN_R:
         weight=StenDelta_R_GetWeightFromSten(rsq,sten_rad);
         break;
     case THETA_FN_R:
         weight=StenTheta_R_GetWeightFromSten(rsq,sten_rad);
         break;
     case THETA_FN_SIG:
         weight=StenTheta_Sigma_GetWeightFromSten(rsq,sten_rad);
         break;
     case THETA_PAIRPOT_RCUT:
         weight=StenTheta_uattr_GetWeightFromSten(icomp,jcomp,rsq,ngpu,gpu,gwu);
         break;
     case THETA_CR_RPM_MSA:
         weight=StenTheta_RPMmsa_GetWeightFromSten(icomp,jcomp,rsq,ngpu,gpu,gwu);
         break;
     case THETA_CR_DATA:
         weight=StenTheta_CrCMS_GetWeightFromSten(icomp,jcomp,rsq,sten_rad);
         break;
     case DELTA_FN_BOND:
         weight=StenDelta_Bond_GetWeightFromSten(rsq,sten_rad);
         break;
     default:
         printf("problem with stencil definitions: stencil_volume_switch ");
         exit(-1); break;
  }
  return(weight);
}
/****************************************************************************/
/*stencil_QuadBoundaryEl_switch: call stencil specific routines to set up the number
   of quadrature points to be used for computing the weight of a given element in
   a given stencil for the case where the element stradles the stencil cutoff.*/
int stencil_quadBoundaryEl_switch(int sten)
{ 
 int num_quad_pts;

  switch(sten){
     case DELTA_FN_R:
         num_quad_pts=StenDelta_R_NquadPtsBoundary();
         break;
     case THETA_FN_R:
         num_quad_pts=StenTheta_R_NquadPtsBoundary();
         break;
     case THETA_FN_SIG:
         num_quad_pts=StenTheta_Sigma_NquadPtsBoundary();
         break;
     case THETA_PAIRPOT_RCUT:
         num_quad_pts=StenTheta_uattr_NquadPtsBoundary();
         break;
     case THETA_CR_RPM_MSA:
         num_quad_pts=StenTheta_RPMmsa_NquadPtsBoundary();
         break;
     case THETA_CR_DATA:
         num_quad_pts=StenTheta_CrCMS_NquadPtsBoundary();
         break;
     case DELTA_FN_BOND:
         num_quad_pts=StenDelta_Bond_NquadPtsBoundary();
         break;
     default:
         printf("problem with stencil definitions: stencil_quadBoundaryEl_switch ");
         exit(-1); break;
  }
  return(num_quad_pts);
}
/****************************************************************************/
/*stencil_quadGauss_switch: call stencil specific routines to set up the number
   of quadrature points to be used for computing the weight of a given element in
   a given stencil for the case where the element stradles the stencil cutoff.*/
int stencil_quadGauss_switch(int sten,double r)
{ 
 int num_quad_pts;

  switch(sten){
     case DELTA_FN_R:
         num_quad_pts=StenDelta_R_NquadPtsGauss(r);
         break;
     case THETA_FN_R:
         num_quad_pts=StenTheta_R_NquadPtsGauss(r);
         break;
     case THETA_FN_SIG:
         num_quad_pts=StenTheta_Sigma_NquadPtsGauss(r);
         break;
     case THETA_PAIRPOT_RCUT:
         num_quad_pts=StenTheta_uattr_NquadPtsGauss(r);
         break;
     case THETA_CR_RPM_MSA:
         num_quad_pts=StenTheta_RPMmsa_NquadPtsGauss(r);
         break;
     case THETA_CR_DATA:
         num_quad_pts=StenTheta_CrCMS_NquadPtsGauss(r);
         break;
     case DELTA_FN_BOND:
         num_quad_pts=StenDelta_Bond_NquadPtsGauss(r);
         break;
     default:
         printf("problem with stencil definitions: stencil_quadGauss_switch ");
         exit(-1); break;
  }
  return(num_quad_pts);
}
/****************************************************************************/
/*stencil_quadGaussIntegrand_switch: call stencil specific routines to set up the number
   of quadrature points to be used for computing the value of the integrand for 
   stencil functions that require integration */
int stencil_quadGaussIntegrand_switch(int sten,double r)
{ 
 int num_quad_pts=0;

  switch(sten){
     case DELTA_FN_R: break;
     case THETA_FN_R: break;
     case THETA_FN_SIG: break;
     case THETA_PAIRPOT_RCUT:
         num_quad_pts=StenTheta_uattr_NquadPtsGaussIntegrand(r); break;
     case THETA_CR_RPM_MSA:
         num_quad_pts=StenTheta_RPMmsa_NquadPtsGaussIntegrand(r); break;
     case THETA_CR_DATA: break;
     case DELTA_FN_BOND: break;
     default:
         printf("problem with stencil definitions: stencil_quadGaussIntegrand_switch ");
         exit(-1); break;
  }
  return(num_quad_pts);
}
/****************************************************************************/
