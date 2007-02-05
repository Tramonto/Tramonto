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
 *  FILE: dft_physics_FMT1.c
 *
 *  This file contains the essential physics of the original Rosenfeld 
 *  fundamental measures theory DFT
 *
 */

#include "dft_physics_FMT1.h"

/*****************************************************************************/
double FMT1_energy_density(double *n)
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
     phi_3 = (n[2]*n[2]*n[2] - 3.0*n[2]*dot_22)/(24.0*PI*(1.0-n[3])*(1.0-n[3]));

     return(phi_1+phi_2+phi_3);
}
/*****************************************************************************/
void FMT1_1stderiv(double *n,double DOT_12,double DOT_22,double *inv_n3, double *dphi_drb_loc)
{
   int idim;

   dphi_drb_loc[0] = log(inv_n3[1]);
   dphi_drb_loc[1] = n[2]*inv_n3[1];
   dphi_drb_loc[2] = n[1]*inv_n3[1] + (n[2]*n[2]-DOT_22)*inv_n3[2] / (8.0*PI);
   dphi_drb_loc[3] = n[0]*inv_n3[1] + (n[1]*n[2]-DOT_12)*inv_n3[2] +
                            (n[2]*n[2]*n[2]-3.*n[2]*DOT_22 )*inv_n3[3]/(12.0*PI);

   for (idim=0;idim<Ndim;idim++){
       dphi_drb_loc[Nrho_bar_s+idim]= -n[Nrho_bar_s+Ndim+idim]*inv_n3[1];
       dphi_drb_loc[Nrho_bar_s+idim+Ndim] = -n[Nrho_bar_s+idim]*inv_n3[1] -
                           n[2]*n[Nrho_bar_s+Ndim+idim]*inv_n3[2]/(4.0*PI);
   }
   return;
}
/*******************************************************************************************/
/* d2phi_drb2_delta_rb_FMT1:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                 for the dphi_drb that use Delta_Fn Stencils (all but S3) */
struct RB_Struct d2phi_drb2_delta_rb_FMT1(double *n, int *offset, double *sign, int icomp)
{
  struct RB_Struct tmp;
  double n1v[NDIM_MAX], n2v[NDIM_MAX];
  double inv_one_m_n3, inv_one_m_n3_sq, inv_one_m_n3_3rd;
  double vector[NDIM_MAX];
  int idim;

  for (idim = 0; idim<Ndim; idim++) {   
     n2v[idim] = sign[idim]*n[Nrho_bar_s+Ndim+idim];  
     n1v[idim] = sign[idim]*n[Nrho_bar_s+idim]; 
  }

  inv_one_m_n3 = 1.0 / (1.0 - n[3]);
  inv_one_m_n3_sq = inv_one_m_n3*inv_one_m_n3;
  inv_one_m_n3_3rd = inv_one_m_n3_sq*inv_one_m_n3;
  for (idim = 0; idim<Ndim; idim++) 
    vector[idim] = offset[idim] * Esize_x[idim] * Inv_rad[icomp];

  tmp.S2 = (inv_one_m_n3*Inv_4pir[icomp]
         + n[2]*Inv_4pi*inv_one_m_n3_sq);

  for (idim = 0; idim<Ndim; idim++)
     tmp.S2 += vector[idim]
              * n2v[idim]*Inv_4pi * inv_one_m_n3_sq;

  tmp.S3 = (Inv_4pirsq[icomp]*inv_one_m_n3 
         + n[2]*Inv_4pir[icomp]*inv_one_m_n3_sq
         + n[1]*inv_one_m_n3_sq+n[2]*n[2]*Inv_4pi*inv_one_m_n3_3rd);

  for (idim = 0; idim<Ndim; idim++)
     tmp.S3 += inv_one_m_n3_sq *
           ( -n2v[idim]*n2v[idim]*Inv_4pi*inv_one_m_n3
           + vector[idim]
           * ( n2v[idim]* Inv_4pir[icomp] + n1v[idim]
              + n[2]*n2v[idim]*inv_one_m_n3/(2*PI) ) );

  tmp.S0 = 0.0;  
  tmp.S1 = inv_one_m_n3;  

  for (idim = 0; idim<Ndim; idim++)
    tmp.V1[idim] = sign[idim]*inv_one_m_n3*vector[idim];
  for (idim = 0; idim<Ndim; idim++)
    tmp.V2[idim] = sign[idim]*(- n2v[idim]*Inv_4pi*inv_one_m_n3_sq
                    - (-Inv_4pir[icomp]*inv_one_m_n3 
                    - n[2]* Inv_4pi*inv_one_m_n3_sq) * vector[idim]);
  return(tmp);
}
/****************************************************************************/
/* d2phi_drb2_theta_rb_FMT1:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                    for the dphi_drb that use Theta_Fn Stencils (S3)      */

struct RB_Struct d2phi_drb2_theta_rb_FMT1(double *n)
{
  struct RB_Struct tmp;
  double n1v[NDIM_MAX], n2v[NDIM_MAX];
  double inv_one_m_n3, inv_one_m_n3_sq, inv_one_m_n3_3rd, inv_one_m_n3_4th;
  int idim;

  for (idim = 0; idim<Ndim; idim++) {
     n2v[idim] = n[Nrho_bar_s+Ndim+idim];  
     n1v[idim] = n[Nrho_bar_s+idim]; 
  }

  inv_one_m_n3 = 1.0 / (1.0 - n[3]);
  inv_one_m_n3_sq = inv_one_m_n3*inv_one_m_n3;
  inv_one_m_n3_3rd = inv_one_m_n3_sq*inv_one_m_n3;
  inv_one_m_n3_4th = inv_one_m_n3_3rd*inv_one_m_n3;

  tmp.S2 = (n[1]*inv_one_m_n3_sq + n[2]*n[2]*Inv_4pi*inv_one_m_n3_3rd);
  for (idim = 0; idim<Ndim; idim++)
     tmp.S2 += - n2v[idim]*n2v[idim]
                         *Inv_4pi*inv_one_m_n3_3rd;
  tmp.S3 = n[0]*inv_one_m_n3_sq 
           + 2*n[1]*n[2]*inv_one_m_n3_3rd
           + n[2]*n[2]*n[2]*Inv_4pi*inv_one_m_n3_4th;
  for (idim = 0; idim<Ndim; idim++)
     tmp.S3 += -(2*n1v[idim]*n2v[idim]*inv_one_m_n3_3rd
            +3*n[2]*n2v[idim]*n2v[idim]*Inv_4pi*inv_one_m_n3_4th);

  tmp.S0 = inv_one_m_n3;             
  tmp.S1 = n[2] * inv_one_m_n3_sq;    

  for (idim = 0; idim<Ndim; idim++)
      tmp.V1[idim] = - n2v[idim]*inv_one_m_n3_sq;
  for (idim = 0; idim<Ndim; idim++)
      tmp.V2[idim] = -(n1v[idim]*inv_one_m_n3_sq
                     + n[2]*n2v[idim]*inv_one_m_n3_3rd/(2.0*PI) );
  return (tmp);
}
/****************************************************************************/
