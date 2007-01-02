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
 *  FILE: dft_physics_FMT2.c
 *
 *  This file contains the essential physics of the Rosenfeld et.al. 
 *  functionals that correct the zero-dimensional crossover behavior of
 *  the original functionals.
 *
 */

#include "dft_physics_FMT2.h"

/*****************************************************************************/
double FMT2_energy_density(double *n)
{
    double dot_12,dot_22,phi_1,phi_2,phi_3;
    int iv1,iv2,idim;

    dot_12=0.0; dot_22=0.0;
    for (idim = 0; idim<Ndim; idim++){
         iv1=Nrho_bar_s+idim;
         iv2=Nrho_bar_s+Ndim+idim;
         dot_12 += n[iv1]*n[iv2];
         dot_22 += n[iv2]*n[iv2];
    }

     phi_1 = -n[0]*log(1.0-n[3]);
     phi_2 = (n[1]*n[2] - dot_12)/(1.0-n[3]);
     phi_3 = POW_DOUBLE_INT(n[2]-dot_22/n[2],3)/(24.0*PI*(1.0-n[3])*(1.0-n[3]));

     return(phi_1+phi_2+phi_3);
}
/****************************************************************************/
void FMT2_1stderiv(double *n,double DOT_12,double DOT_22,double *inv_n3, double *dphi_drb_loc)
{
   int idim;
   double alpha, gamma[NDIM_MAX],alpha_sq,alpha_cb,beta;

   /* same as FMT1 contributions */
   dphi_drb_loc[0] =  log(inv_n3[1]);
   dphi_drb_loc[1] =  n[2]*inv_n3[1];
   dphi_drb_loc[2] =  n[1]*inv_n3[1];
   dphi_drb_loc[3] = n[0]*inv_n3[1] + (n[1]*n[2]-DOT_12) * inv_n3[2];

   for (idim = 0; idim < Ndim; idim++) {
       dphi_drb_loc[Nrho_bar_s+idim] = -n[Nrho_bar_s+Ndim+idim]*inv_n3[1];
       dphi_drb_loc[Nrho_bar_s+Ndim+idim] = -n[Nrho_bar_s+idim]*inv_n3[1];
   }

   /* new contributions */
   if (n[2] > 1.e-15){
         alpha=n[2]-DOT_22/n[2];
         beta = 1.0+DOT_22/(n[2]*n[2]);
         for (idim = 0; idim < Ndim; idim++){
             gamma[idim] = n[Nrho_bar_s+Ndim+idim]/n[2];
         }
   }
   else{
       alpha=n[2];
       beta = 1.0;
       for (idim = 0; idim < Ndim; idim++) gamma[idim] = 0.0;
   }
   alpha_sq=alpha*alpha;
   alpha_cb=alpha_sq*alpha;

   dphi_drb_loc[2] += alpha_sq*beta*inv_n3[2]/(8.0*PI);
   dphi_drb_loc[3] += alpha_cb*inv_n3[3]/(12.0*PI);

   for (idim = 0; idim < Ndim; idim++) {
       dphi_drb_loc[Nrho_bar_s+Ndim+idim] -= inv_n3[2]*alpha_sq*gamma[idim]/(4.0*PI);
   }
   return;
}
/*******************************************************************************************/
/* d2phi_drb2_delta_rb_FMT2:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                 for the dphi_drb that use Delta_Fn Stencils (all but S3) */

struct RB_Struct d2phi_drb2_delta_rb_FMT2(double *n, int *offset, double *sign,int icomp)
{
  struct RB_Struct tmp;
  double n1v[NDIM_MAX], n2v[NDIM_MAX];
  double inv_one_m_n3, inv_one_m_n3_sq, inv_one_m_n3_3rd, inv_one_m_n3_4th;
  double vector[NDIM_MAX];
  int idim;
  double DOT_rho22,alpha,alpha_sq,alpha_cb,beta,gamma[3],DOT_gamma,eps;
  
  for (idim = 0; idim<Ndim; idim++) {
    n2v[idim] = sign[idim]*n[Nrho_bar_s+Ndim+idim];
    n1v[idim] = sign[idim]*n[Nrho_bar_s+idim];
  }

  DOT_rho22 = 0.0;
  for (idim = 0; idim < Ndim; idim++) {
      DOT_rho22 += n2v[idim] * n2v[idim];
  }
  if (n[2] > 1.e-15){
     alpha=n[2]-DOT_rho22/n[2];
     beta = 1.0+DOT_rho22/(n[2]*n[2]);
     DOT_gamma=0.0;
     for (idim = 0; idim < Ndim; idim++){
          gamma[idim] = n2v[idim]/n[2];
          DOT_gamma+=gamma[idim]*gamma[idim];
     }
     eps=alpha/n[2];
  }
  else{
      alpha=n[2];
      beta = 1.0;
      for (idim = 0; idim < Ndim; idim++) gamma[idim] = 0.0;
      DOT_gamma=0.0;
      eps=1.0;
  }
  alpha_sq=alpha*alpha;
  alpha_cb=alpha_sq*alpha;

  inv_one_m_n3 = 1.0 / (1.0 - n[3]);
  inv_one_m_n3_sq = inv_one_m_n3*inv_one_m_n3;
  inv_one_m_n3_3rd = inv_one_m_n3_sq*inv_one_m_n3;
  inv_one_m_n3_4th = inv_one_m_n3_sq*inv_one_m_n3_sq;
  for (idim = 0; idim<Ndim; idim++) 
    vector[idim] = offset[idim] * Esize_x[idim] * Inv_rad[icomp];

  tmp.S2 = inv_one_m_n3*Inv_4pir[icomp];      /*old rosen n2-n1*/

  tmp.S2 += Inv_4pi*inv_one_m_n3_sq*     /* new rosen n2-n2*/
            alpha*(beta*beta-eps*DOT_gamma);

  for (idim = 0; idim<Ndim; idim++)              /*new rosen n2-n2v*/
      tmp.S2 -= Inv_4pi*inv_one_m_n3_sq*
                alpha*(eps -2*beta)*gamma[idim]*vector[idim];

  tmp.S3 = (Inv_4pirsq[icomp]*inv_one_m_n3   /* old rosen*/
         + n[2]*Inv_4pir[icomp]*inv_one_m_n3_sq
         + n[1]*inv_one_m_n3_sq);

  for (idim = 0; idim<Ndim; idim++)             /* old rosen n3-nv1/v2*/
     tmp.S3 += inv_one_m_n3_sq *
               vector[idim] * (n2v[idim]*Inv_4pir[icomp] + n1v[idim]);
  
  tmp.S3 += Inv_4pi*inv_one_m_n3_3rd*alpha_sq*beta; /*new rosen n3-n2*/

  for (idim = 0; idim<Ndim; idim++)
     tmp.S3 += 2.0*Inv_4pi*inv_one_m_n3_3rd *    /*new rosen n3-n2v*/
               alpha_sq*gamma[idim]*vector[idim];

  tmp.S0 = 0.0;                     /*old rosen n0-(n0-n2,nv1,nv2)*/
  tmp.S1 = inv_one_m_n3;    /*old rosen n1-n2*/

  for (idim = 0; idim<Ndim; idim++){                               /*old rosen*/
    tmp.V1[idim] = sign[idim]*(inv_one_m_n3)*vector[idim];     /*n1v-n2v*/
    tmp.V2[idim] = sign[idim]*(Inv_4pir[icomp]*inv_one_m_n3)*vector[idim]; /*n2v-n1v*/
  }

  for (idim = 0; idim<Ndim; idim++){                         
    tmp.V2[idim] += sign[idim]*Inv_4pi*inv_one_m_n3_sq*  /*new n2v-n2*/
                    alpha*(eps-2.0*beta)*gamma[idim];
                                                                
    tmp.V2[idim] -= sign[idim]*Inv_4pi*inv_one_m_n3_sq*  /*new n2v-n2v*/
                    alpha*(4.0*DOT_gamma-eps)*vector[idim];
  }
  return(tmp);
}
/****************************************************************************/
/* d2phi_drb2_theta_rb_FMT2:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                    for the dphi_drb that use Theta_Fn Stencils (S3)      */

struct RB_Struct d2phi_drb2_theta_rb_FMT2(double *n)
{
  struct RB_Struct tmp;
  double n1v[NDIM_MAX], n2v[NDIM_MAX];
  double inv_one_m_n3, inv_one_m_n3_sq, inv_one_m_n3_3rd, inv_one_m_n3_4th;
  int idim;
  double DOT_rho22,DOT_rho12,alpha,alpha_sq,alpha_cb,beta,gamma[3];

  for (idim = 0; idim<Ndim; idim++) {
    n2v[idim] = n[Nrho_bar_s+Ndim+idim];
    n1v[idim] = n[Nrho_bar_s+idim];
  }

  inv_one_m_n3 = 1.0 / (1.0 - n[3]);
  inv_one_m_n3_sq = inv_one_m_n3*inv_one_m_n3;
  inv_one_m_n3_3rd = inv_one_m_n3_sq*inv_one_m_n3;
  inv_one_m_n3_4th = inv_one_m_n3_3rd*inv_one_m_n3;

  DOT_rho22 = 0.0;
  DOT_rho12 = 0.0;
  for (idim = 0; idim < Ndim; idim++) {
      DOT_rho22 += n2v[idim] * n2v[idim];
      DOT_rho12 += n1v[idim] * n2v[idim];
  }

  if (n[2] > 1.e-15){
        alpha=n[2]-DOT_rho22/n[2];
        beta = 1.0+DOT_rho22/(n[2]*n[2]);
        for (idim = 0; idim < Ndim; idim++){
            gamma[idim] = n2v[idim]/n[2];
        }
  }
  else{
      alpha=n[2];
      beta = 1.0;
      for (idim = 0; idim < Ndim; idim++) gamma[idim] = 0.0;
  }
  alpha_sq=alpha*alpha;
  alpha_cb=alpha_sq*alpha;

  tmp.S2 = n[1]*inv_one_m_n3_sq;                  /*old n2-n3*/
  tmp.S2 += Inv_4pi*inv_one_m_n3_3rd*alpha_sq*beta;      /*new n2-n3*/

  tmp.S3 = (n[0]*inv_one_m_n3_sq                   /*old n3-n3*/
                   + 2.0*(n[1]*n[2] - DOT_rho12)*inv_one_m_n3_3rd);
  tmp.S3 += Inv_4pi*inv_one_m_n3_4th*alpha_cb;    /*new n3-n3*/

  tmp.S0 = inv_one_m_n3;            /*old n0-n3*/
  tmp.S1 = n[2] * inv_one_m_n3_sq;   /*old n1-n3*/

  for (idim = 0; idim<Ndim; idim++){         
      tmp.V1[idim] = - n2v[idim]*inv_one_m_n3_sq;  /*old nv1-n3*/
      tmp.V2[idim] = - n1v[idim]*inv_one_m_n3_sq;  /*old nv2-n3*/
  }

  for (idim = 0; idim<Ndim; idim++)                            /*new nv2-n3*/
      tmp.V2[idim] -= 2.0*Inv_4pi*inv_one_m_n3_3rd*alpha_sq*gamma[idim];

  return (tmp);
}
/*******************************************************************************************/
