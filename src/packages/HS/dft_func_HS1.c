/*
//@HEADER
// ********************************************************************
// Copyright (2006) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000, there is a non-exclusive license for use of this
// work by or on behalf of the U.S. Government. Export of this program
// may require a license from the United States Government.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// ********************************************************************
//@HEADER
*/


/*
 *  FILE: dft_func_HS1.c
 *
 *  This file contains first and second derivatives of the free energy densities
 *  of Rosenfeld's Fundamental Measures Theory
 *
 */

/*********************************************************************/
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
/*****************************************************************************/
/* d2phi_drb2_delta_rb_FMT1:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                 for the dphi_drb that use Delta_Fn Stencils (all but S3) */

static struct RB_Struct d2phi_drb2_delta_rb_FMT1(int junk, int jnode_box,double **x,
                                            double weight, int *offset, double *sign,
                                            double inv_rad, double inv_4pir,
                                            double inv_4pirsq)

{
  struct RB_Struct tmp;
  double rb0, rb1, rb2, rb3, rb1v[NDIM_MAX], rb2v[NDIM_MAX];
  double inv_one_m_rb3, inv_one_m_rb3_sq, inv_one_m_rb3_3rd;
  double vector[NDIM_MAX];
  int idim,loc_js,loc_jv;

  rb3=x[junk][jnode_box];
  rb2=x[junk+1][jnode_box];
  rb1=x[junk+2][jnode_box];
  rb0=x[junk+3][jnode_box];

  for (idim = 0; idim<Ndim; idim++) {
     rb2v[idim] = sign[idim]*x[junk+Nrho_bar_s+idim][jnode_box];
     rb1v[idim] = sign[idim]*x[junk+Nrho_bar_s+Ndim+idim][jnode_box];
  }

  inv_one_m_rb3 = 1.0 / (1.0 - rb3);
  inv_one_m_rb3_sq = inv_one_m_rb3*inv_one_m_rb3;
  inv_one_m_rb3_3rd = inv_one_m_rb3_sq*inv_one_m_rb3;
  for (idim = 0; idim<Ndim; idim++)
    vector[idim] = offset[idim] * Esize_x[idim] * inv_rad;

  tmp.S2 = (inv_one_m_rb3*inv_4pir
         + rb2*Inv_4pi*inv_one_m_rb3_sq)*weight;

  for (idim = 0; idim<Ndim; idim++)
     tmp.S2 += weight * vector[idim]
              * rb2v[idim]*Inv_4pi * inv_one_m_rb3_sq;

  tmp.S3 = weight * (inv_4pirsq*inv_one_m_rb3
         + rb2*inv_4pir*inv_one_m_rb3_sq
         + rb1*inv_one_m_rb3_sq+rb2*rb2*Inv_4pi*inv_one_m_rb3_3rd);

  for (idim = 0; idim<Ndim; idim++)
     tmp.S3 += weight * inv_one_m_rb3_sq *
           ( -rb2v[idim]*rb2v[idim]*Inv_4pi*inv_one_m_rb3
           + vector[idim]
           * ( rb2v[idim]* inv_4pir + rb1v[idim]
              + rb2*rb2v[idim]*inv_one_m_rb3/(2*PI) ) );

  tmp.S0 = 0.0;
  tmp.S1 = inv_one_m_rb3*weight;

  for (idim = 0; idim<Ndim; idim++)
    tmp.V1[idim] = sign[idim]*(weight * inv_one_m_rb3) * vector[idim];
  for (idim = 0; idim<Ndim; idim++)
    tmp.V2[idim] = sign[idim]*weight *(- rb2v[idim]*Inv_4pi*inv_one_m_rb3_sq
                    - (-inv_4pir*inv_one_m_rb3
                    - rb2* Inv_4pi*inv_one_m_rb3_sq) * vector[idim]);
  return(tmp);
}
/****************************************************************************/
/* d2phi_drb2_theta_rb_FMT1:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                    for the dphi_drb that use Theta_Fn Stencils (S3)      */

static struct RB_Struct d2phi_drb2_theta_rb_FMT1(int junk, int jnode_box,double **x,double weight,
                                            int *offset)
{
  struct RB_Struct tmp;
  double rb0, rb1, rb2, rb3, rb1v[NDIM_MAX], rb2v[NDIM_MAX];
  double inv_one_m_rb3, inv_one_m_rb3_sq, inv_one_m_rb3_3rd, inv_one_m_rb3_4th;
  int idim,loc_js,loc_jv;


  rb3 = x[junk][jnode_box];
  rb2 = x[junk+1][jnode_box];
  rb1 = x[junk+2][jnode_box];
  rb0 = x[junk+3][jnode_box];
 
  for (idim = 0; idim<Ndim; idim++) {
    rb2v[idim] = x[junk+Nrho_bar_s+idim][jnode_box];
    rb1v[idim] = x[junk+Nrho_bar_s+Ndim+idim][jnode_box];
  }

  inv_one_m_rb3 = 1.0 / (1.0 - rb3);
  inv_one_m_rb3_sq = inv_one_m_rb3*inv_one_m_rb3;
  inv_one_m_rb3_3rd = inv_one_m_rb3_sq*inv_one_m_rb3;
  inv_one_m_rb3_4th = inv_one_m_rb3_3rd*inv_one_m_rb3;

  tmp.S2 = weight * (rb1*inv_one_m_rb3_sq
         + rb2*rb2*Inv_4pi*inv_one_m_rb3_3rd);
  for (idim = 0; idim<Ndim; idim++)
     tmp.S2 += - weight * rb2v[idim]*rb2v[idim]
                         *Inv_4pi*inv_one_m_rb3_3rd;
  tmp.S3 = weight * (rb0*inv_one_m_rb3_sq
                     + 2*rb1*rb2*inv_one_m_rb3_3rd
                     + rb2*rb2*rb2*Inv_4pi*inv_one_m_rb3_4th);
  for (idim = 0; idim<Ndim; idim++)
     tmp.S3 += -weight * (2*rb1v[idim]*rb2v[idim]*inv_one_m_rb3_3rd
            +3*rb2*rb2v[idim]*rb2v[idim]*Inv_4pi*inv_one_m_rb3_4th);

  tmp.S0 = weight * inv_one_m_rb3;
  tmp.S1 = weight * rb2 * inv_one_m_rb3_sq;

  for (idim = 0; idim<Ndim; idim++)
      tmp.V1[idim] = - weight * rb2v[idim]*inv_one_m_rb3_sq;
  for (idim = 0; idim<Ndim; idim++)
      tmp.V2[idim] = - weight * (rb1v[idim]*inv_one_m_rb3_sq
                     + rb2*rb2v[idim]*inv_one_m_rb3_3rd/(2.0*PI) );
  return (tmp);
}
