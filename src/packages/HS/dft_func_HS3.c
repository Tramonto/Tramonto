/*****************************************************************************/
/* d2phi_drb2_delta_rb_FMT3:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                 for the dphi_drb that use Delta_Fn Stencils (all but S3) */

static struct RB_Struct d2phi_drb2_delta_rb_FMT3(int junk, int jnode_box,double **x,
                                            double weight, int *offset, double *sign,
                                            double inv_rad, double inv_4pir,
                                            double inv_4pirsq)

{
  struct RB_Struct tmp;
  double n[4+NDIM_MAX],n3sq,n3cb,n3_4th,n2sq,n2cb;
  double inv_n3[5],DOT_22,DOT_12,fac1,fac2,fac3,vector[NDIM_MAX];
  int idim,i2v,i1v;


  n[3] = x[junk][jnode_box];
  n[2] = x[junk+1][jnode_box];
  n[1] = x[junk+2][jnode_box];
  n[0] = x[junk+3][jnode_box];
 
  DOT_22 = 0.0;
  DOT_12 = 0.0;
  for (idim = 0; idim<Ndim; idim++) {
    i2v=Nrho_bar_s+Ndim+idim;
    i1v=Nrho_bar_s+idim;
    n[i2v] = x[junk+Nrho_bar_s+idim][jnode_box];
    n[i1v]= x[junk+Nrho_bar_s+Ndim+idim][jnode_box];
    DOT_22 += n[i2v] * n[i2v];
    DOT_12 += n[i1v] * n[i2v];
  }

  inv_n3[0] = (1.-n[3]);
  inv_n3[1] = 1.0 / (1.0 - n[3]);
  inv_n3[2] = inv_n3[1]*inv_n3[1];
  inv_n3[3] = inv_n3[2]*inv_n3[1];
  inv_n3[4] = inv_n3[3]*inv_n3[1];

  n3sq=n[3]*n[3];
  n3cb=n3sq*n[3];
  n3_4th=n3cb*n[3];

  n2sq=n[2]*n[2];
  n2cb=n2sq*n[2];

  if (n[3]>1.e-10){
     fac1 = n[3]-2*(1.-n[3])*log(1.-n[3]);
     fac2 = n[3]+(1.-n[3])*(1.-n[3])*log(1.-n[3]);
     fac3 = -(inv_n3[2]/n3cb) + (inv_n3[3]/n3sq);

     for (idim = 0; idim<Ndim; idim++)
       vector[idim] = offset[idim] * Esize_x[idim] * inv_rad;

     tmp.S2 = weight*(inv_n3[1]*inv_4pir
         + n[2]*fac2*Inv_4pi*inv_n3[2]/(6.*PI*n3sq));

     for (idim = 0; idim<Ndim; idim++){
        i2v=Nrho_bar_s+Ndim+idim;
        tmp.S2 += weight * vector[idim]
              *n[i2v]*fac2*Inv_4pi*inv_n3[2]/(6.*PI*n3sq);
     }

     tmp.S3 = weight * (inv_4pirsq*inv_n3[1]
         + inv_4pir*n[2]*inv_n3[2]
         + Inv_4pi*n[1]*inv_n3[2]+(n2sq-DOT_22)*(
      (-2.*fac2*fac3) + (fac1*inv_n3[2]/n3sq))/(12.*PI)  );

     for (idim = 0; idim<Ndim; idim++){
        i2v=Nrho_bar_s+Ndim+idim;
        i1v=Nrho_bar_s+idim;
        tmp.S3 += weight * inv_n3[2] * vector[idim] * (
               n[i2v]*inv_4pir + n[i1v]
              + n[2]*n[i2v]*(
              (-2.*fac2*fac3) + (fac1*inv_n3[2]/n3sq))/(6.*PI) ) ;
     }

     tmp.S0 = 0.0; 
     tmp.S1 = weight*inv_n3[1];

     for (idim = 0; idim<Ndim; idim++){
       tmp.V1[idim] = sign[idim]*(weight * inv_n3[1]) * vector[idim];
     }

     for (idim = 0; idim<Ndim; idim++){
       i2v=Nrho_bar_s+Ndim+idim;
       i1v=Nrho_bar_s+idim;
       tmp.V2[idim] = sign[idim]*weight *
                    ( 2.*n[i2v]*Inv_4pi*inv_n3[2]*fac2/(12.*PI*n3sq)
                    - (-inv_4pir*inv_n3[1] - n[2]*Inv_4pi*fac2*inv_n3[2]/(6.*PI*n3sq))*vector[idim]);
     }
  }
  else{
    tmp.S2 = 0.0;
    tmp.S3 = 0.0;
    tmp.S0 = 0.0;
    tmp.S1 = 0.0;

    for (idim = 0; idim<Ndim; idim++){
      tmp.V1[idim] = 0.0;
      tmp.V2[idim] = 0.0;
    }
  }
  return(tmp);
}
/****************************************************************************/
/* d2phi_drb2_theta_rb_FMT3:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                    for the dphi_drb that use Theta_Fn Stencils (S3)      */

static struct RB_Struct d2phi_drb2_theta_rb_FMT3(int junk, int jnode_box,double **x,double weight,
                                            int *offset)
{
  struct RB_Struct tmp;
  double n[4+NDIM_MAX],n3sq,n3cb,n3_4th,n2sq,n2cb;
  double inv_n3[5],DOT_22,DOT_12,fac1,fac2,fac3;
  int idim,i2v,i1v;


  n[3] = x[junk][jnode_box];
  n[2] = x[junk+1][jnode_box];
  n[1] = x[junk+2][jnode_box];
  n[0] = x[junk+3][jnode_box];
 
  DOT_22 = 0.0;
  DOT_12 = 0.0;
  for (idim = 0; idim<Ndim; idim++) {
    i2v=Nrho_bar_s+Ndim+idim;
    i1v=Nrho_bar_s+idim;
    n[i2v] = x[junk+Nrho_bar_s+idim][jnode_box];
    n[i1v]= x[junk+Nrho_bar_s+Ndim+idim][jnode_box];
    DOT_22 += n[i2v] * n[i2v];
    DOT_12 += n[i1v] * n[i2v];
  }

  inv_n3[0] = (1.-n[3]);
  inv_n3[1] = 1.0 / (1.0 - n[3]);
  inv_n3[2] = inv_n3[1]*inv_n3[1];
  inv_n3[3] = inv_n3[2]*inv_n3[1];
  inv_n3[4] = inv_n3[3]*inv_n3[1];

  n3sq=n[3]*n[3];
  n3cb=n3sq*n[3];
  n3_4th=n3cb*n[3];

  n2sq=n[2]*n[2];
  n2cb=n2sq*n[2];

  if (n[3]>1.e-10){
    fac1 = n[3]-2*(1.-n[3])*log(1.-n[3]);
    fac2 = n[3]+(1.-n[3])*(1.-n[3])*log(1.-n[3]);
    fac3 = -(inv_n3[2]/n3cb) + (inv_n3[3]/n3sq);

    tmp.S2 = weight * (
           n[1]*inv_n3[2]
         + (n[2]*n[2]-DOT_22)*(fac1*inv_n3[2]/n3sq + 2.*fac2*fac3)/(12*PI));

    tmp.S3 = weight * (
           n[1]*inv_n3[2]
           + 2*(n[1]*n[2]-DOT_12)*inv_n3[3]
           + ((n2cb-3.*n[2]*DOT_22)/(36.*PI))*(
           + (((2.*log(1.-n[3])+1)*inv_n3[2])/n3sq)
           + (4.*fac1*fac3) + (6.*fac2*inv_n3[2]/n3_4th)
           + (inv_n3[4]/n3sq) )  );

    tmp.S0 = weight * inv_n3[1];
    tmp.S1 = weight * n[2] * inv_n3[2];

    for (idim = 0; idim<Ndim; idim++){
      i2v=Nrho_bar_s+idim+Ndim;
      tmp.V1[idim] = - weight * n[i2v]*inv_n3[2];
    }

    for (idim = 0; idim<Ndim; idim++){
      i1v=Nrho_bar_s+idim;
      i2v=Nrho_bar_s+idim+Ndim;
      tmp.V2[idim] = - weight * (n[i1v]*inv_n3[2]
        + n[2]*n[i2v]*((fac1*inv_n3[2]/n3sq)+(2.*fac2*fac3))/(6.*PI) );
    }
  }
  else{
    tmp.S2 = 0.0;
    tmp.S3 = 0.0;
    tmp.S0 = 0.0;
    tmp.S1 = 0.0;

    for (idim = 0; idim<Ndim; idim++){
      tmp.V1[idim] = 0.0;
      tmp.V2[idim] = 0.0;
    }
  }
  return (tmp);
}
/*******************************************************************************************/
void FMT3_1stderiv(double *n,double DOT_12,double DOT_22,double *inv_n3, double *dphi_drb_loc)
{
   int idim,iv1,iv2;
;
   dphi_drb_loc[0] =  log(inv_n3[1]);
   dphi_drb_loc[1] =  n[2]*inv_n3[1];

   if (n[3]>1.e-10){
   dphi_drb_loc[2] =  n[1]*inv_n3[1] +
                    (n[2]*n[2]-DOT_22)*inv_n3[2]*
                    (n[3]+(1.-n[3])*(1.-n[3])*log(1.-n[3]))/(12.0*PI*n[3]*n[3]);


   dphi_drb_loc[3] = n[0]*inv_n3[1] + (n[1]*n[2] - DOT_12)*inv_n3[2]
                + (((n[2]*n[2]*n[2]-3.0*n[2]*DOT_22)*inv_n3[2])/(36.*PI*n[3]*n[3]))*
           (   + 2.*(n[3]+(1.-n[3])*(1.-n[3])*log(1.-n[3]))*inv_n3[1]
               + (-2.*(1.-n[3])*log(1.-n[3])+n[3])
               -2.*((n[3]+(1.-n[3])*(1.-n[3])*log(1.-n[3]))/n[3]) );
  }
  else{
     dphi_drb_loc[2]=0.0;
     dphi_drb_loc[3]=0.0;
  }
 
   for (idim=0;idim<Ndim;idim++){
      iv1=Nrho_bar_s+idim;
      iv2=Nrho_bar_s+Ndim+idim;
      dphi_drb_loc[iv1] = -n[iv2]*inv_n3[1];
      if (n[3]>1.e-10)
           dphi_drb_loc[iv2] = -n[iv1]*inv_n3[1] -
                     n[2]*n[iv2]*inv_n3[2]*(n[3]+(1.0-n[3])*(1.0-n[3])*log(1.0-n[3]))/
                                                                    (6.0*PI*n[3]*n[3]);
      else dphi_drb_loc[iv2]=0.0;
   }
   return;
}
/*******************************************************************************************/

