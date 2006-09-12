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
 *  FILE: dft_func_HS2.c
 *
 *  This file contains first and second derivatives of the free energy densities
 *  of a modified Fundamental Measures Theory that correct for crossover to zero-dimensions
 *  i.e. freezing.
 *
 */


/*******************************************************************************************/
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
/****************************************************************************/
/* d2phi_drb2_delta_rb_FMT2:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                 for the dphi_drb that use Delta_Fn Stencils (all but S3) */

static struct RB_Struct d2phi_drb2_delta_rb_FMT2(int junk, int jnode_box,double **x,
                                             double weight, int *offset, double *sign,
                                             double inv_rad, double inv_4pir,
                                             double inv_4pirsq)
{
  struct RB_Struct tmp;
  double rb0, rb1, rb2, rb3, rb1v[NDIM_MAX], rb2v[NDIM_MAX];
  double inv_one_m_rb3, inv_one_m_rb3_sq, inv_one_m_rb3_3rd, inv_one_m_rb3_4th;
  double vector[NDIM_MAX];
  int idim,loc_js,loc_jv;
  double DOT_rho22,alpha,alpha_sq,alpha_cb,beta,gamma[3],DOT_gamma,eps;
 
  rb3 = x[junk][jnode_box];
  rb2 = x[junk+1][jnode_box];
  rb1 = x[junk+2][jnode_box];
  rb0 = x[junk+3][jnode_box];

  for (idim = 0; idim<Ndim; idim++) {
    rb2v[idim] = sign[idim]*x[junk+Nrho_bar_s+idim][jnode_box];
    rb1v[idim] = sign[idim]*x[junk+Nrho_bar_s+Ndim+idim][jnode_box];
  }

  DOT_rho22 = 0.0;
  for (idim = 0; idim < Ndim; idim++) {
      DOT_rho22 += rb2v[idim] * rb2v[idim];
  }
  if (rb2 > 1.e-15){
     alpha=rb2-DOT_rho22/rb2;
     beta = 1.0+DOT_rho22/(rb2*rb2);
     DOT_gamma=0.0;
     for (idim = 0; idim < Ndim; idim++){
          gamma[idim] = rb2v[idim]/rb2;
          DOT_gamma+=gamma[idim]*gamma[idim];
     }
     eps=alpha/rb2;
  }
  else{
      alpha=rb2;
      beta = 1.0;
      for (idim = 0; idim < Ndim; idim++) gamma[idim] = 0.0;
      DOT_gamma=0.0;
      eps=1.0;
  }
  alpha_sq=alpha*alpha;
  alpha_cb=alpha_sq*alpha;

  inv_one_m_rb3 = 1.0 / (1.0 - rb3);
  inv_one_m_rb3_sq = inv_one_m_rb3*inv_one_m_rb3;
  inv_one_m_rb3_3rd = inv_one_m_rb3_sq*inv_one_m_rb3;
  inv_one_m_rb3_4th = inv_one_m_rb3_sq*inv_one_m_rb3_sq;
  for (idim = 0; idim<Ndim; idim++)
    vector[idim] = offset[idim] * Esize_x[idim] * inv_rad;

  tmp.S2 = inv_one_m_rb3*inv_4pir*weight;      /*old rosen rb2-rb1*/

  tmp.S2 += weight*Inv_4pi*inv_one_m_rb3_sq*     /* new rosen rb2-rb2*/
            alpha*(beta*beta-eps*DOT_gamma);

  for (idim = 0; idim<Ndim; idim++)              /*new rosen rb2-rb2v*/
      tmp.S2 -= weight*Inv_4pi*inv_one_m_rb3_sq*
                alpha*(eps -2*beta)*gamma[idim]*vector[idim];

  tmp.S3 = weight * (inv_4pirsq*inv_one_m_rb3   /* old rosen*/
         + rb2*inv_4pir*inv_one_m_rb3_sq
         + rb1*inv_one_m_rb3_sq);

  for (idim = 0; idim<Ndim; idim++)             /* old rosen rb3-rbv1/v2*/
     tmp.S3 += weight * inv_one_m_rb3_sq *
               vector[idim] * (rb2v[idim]*inv_4pir + rb1v[idim]);

  tmp.S3 += weight*Inv_4pi*inv_one_m_rb3_3rd*alpha_sq*beta; /*new rosen rb3-rb2*/

  for (idim = 0; idim<Ndim; idim++)
     tmp.S3 += weight * 2.0*Inv_4pi*inv_one_m_rb3_3rd *    /*new rosen rb3-rb2v*/
               alpha_sq*gamma[idim]*vector[idim];

  tmp.S0 = 0.0;                     /*old rosen rb0-(rb0-rb2,rbv1,rbv2)*/
  tmp.S1 = inv_one_m_rb3*weight;    /*old rosen rb1-rb2*/

  for (idim = 0; idim<Ndim; idim++){                               /*old rosen*/
    tmp.V1[idim] = sign[idim]*(weight*inv_one_m_rb3)*vector[idim];     /*rb1v-rb2v*/
    tmp.V2[idim] = sign[idim]*(weight*inv_4pir*inv_one_m_rb3)*vector[idim]; /*rb2v-rb1v*/
  }

  for (idim = 0; idim<Ndim; idim++){
    tmp.V2[idim] += sign[idim]*weight*Inv_4pi*inv_one_m_rb3_sq*  /*new rb2v-rb2*/
                    alpha*(eps-2.0*beta)*gamma[idim];

    tmp.V2[idim] -= sign[idim]*weight*Inv_4pi*inv_one_m_rb3_sq*  /*new rb2v-rb2v*/
                    alpha*(4.0*DOT_gamma-eps)*vector[idim];
  }
  return(tmp);
}
/****************************************************************************/
/* d2phi_drb2_theta_rb_FMT2:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                    for the dphi_drb that use Theta_Fn Stencils (S3)      */

static struct RB_Struct d2phi_drb2_theta_rb_FMT2(int junk, int jnode_box,double **x,double weight,
                                             int *offset)
{
  struct RB_Struct tmp;
  double rb0, rb1, rb2, rb3, rb1v[NDIM_MAX], rb2v[NDIM_MAX];
  double inv_one_m_rb3, inv_one_m_rb3_sq, inv_one_m_rb3_3rd, inv_one_m_rb3_4th;
  int idim,loc_js,loc_jv;
  double DOT_rho22,DOT_rho12,alpha,alpha_sq,alpha_cb,beta,gamma[3];

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

  DOT_rho22 = 0.0;
  DOT_rho12 = 0.0;
  for (idim = 0; idim < Ndim; idim++) {
      DOT_rho22 += rb2v[idim] * rb2v[idim];
      DOT_rho12 += rb1v[idim] * rb2v[idim];
  }

  if (rb2 > 1.e-15){
        alpha=rb2-DOT_rho22/rb2;
        beta = 1.0+DOT_rho22/(rb2*rb2);
        for (idim = 0; idim < Ndim; idim++){
            gamma[idim] = rb2v[idim]/rb2;
        }
  }
  else{
      alpha=rb2;
      beta = 1.0;
      for (idim = 0; idim < Ndim; idim++) gamma[idim] = 0.0;
  }
  alpha_sq=alpha*alpha;
  alpha_cb=alpha_sq*alpha;

  tmp.S2 = weight *rb1*inv_one_m_rb3_sq;                  /*old rb2-rb3*/
  tmp.S2 += weight * Inv_4pi*inv_one_m_rb3_3rd*alpha_sq*beta;      /*new rb2-rb3*/

  tmp.S3 = weight *(rb0*inv_one_m_rb3_sq                   /*old rb3-rb3*/
                   + 2.0*(rb1*rb2 - DOT_rho12)*inv_one_m_rb3_3rd);
  tmp.S3 += weight*Inv_4pi*inv_one_m_rb3_4th*alpha_cb;    /*new rb3-rb3*/

  tmp.S0 = weight * inv_one_m_rb3;            /*old rb0-rb3*/
  tmp.S1 = weight * rb2 * inv_one_m_rb3_sq;   /*old rb1-rb3*/

  for (idim = 0; idim<Ndim; idim++){
      tmp.V1[idim] = - weight * rb2v[idim]*inv_one_m_rb3_sq;  /*old rbv1-rb3*/
      tmp.V2[idim] = - weight * rb1v[idim]*inv_one_m_rb3_sq;  /*old rbv2-rb3*/
  }

  for (idim = 0; idim<Ndim; idim++)                            /*new rbv2-rb3*/
      tmp.V2[idim] -= weight * 2.0*Inv_4pi*inv_one_m_rb3_3rd*alpha_sq*gamma[idim];

  return (tmp);
}

