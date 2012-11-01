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
 *  FILE: dft_physics_FMT4.c
 *
 *  This file contains the essential physics of the FMT functional of Oleksy and 
 *  Hansen, Mol Phys 104, 2871 (2006) (see appendix in that paper)
 */

#include "dft_physics_FMT4.h"

/*****************************************************************************/
double FMT4_energy_density(double *n)
{
    double dot_12,dot_22,phi_1,phi_2,phi_3,n2sq;
    int iv1,iv2,idim;

    dot_12=0.0; dot_22=0.0;
    for (idim = 0; idim<Ndim; idim++){
         iv1=Nrho_bar_s+idim;
         iv2=Nrho_bar_s+Ndim+idim;
         dot_12 += n[iv1]*n[iv2];
         dot_22 += n[iv2]*n[iv2];
    }
    n2sq = n[2]*n[2];

     phi_1 = -n[0]*log(1.0-n[3]);
     phi_2 = (n[1]*n[2] - dot_12)/(1.0-n[3]);
     phi_3 = n[2]*n2sq*(1-dot_22/n2sq)*(1-dot_22/n2sq)*(1-dot_22/n2sq)*
             (n[3]+(1.0-n[3])*(1.0-n[3])*log(1.0-n[3]))/
                 (36.0*PI*n[3]*n[3]*(1.0-n[3])*(1.0-n[3]));

     return(phi_1+phi_2+phi_3);
}
/*******************************************************************************************/
void FMT4_1stderiv(double *n,double DOT_12,double DOT_22,double *inv_n3, double *dphi_drb_loc)
{
   int idim,iv1,iv2;
   double fac, fac1, fac2;

   dphi_drb_loc[0] =  log(inv_n3[1]);
   dphi_drb_loc[1] =  n[2]*inv_n3[1];

   if (n[3]>1.e-10){
   fac = 1-DOT_22/(n[2]*n[2]);
   fac1 = n[3] - 2.*(1.-n[3])*log(1.-n[3]);
   fac2 = (n[3]+(1.-n[3])*(1.-n[3])*log(1.-n[3]));
   dphi_drb_loc[2] =  n[1]*inv_n3[1] +
                    (n[2]*n[2]*fac*fac*fac+2*fac*fac*DOT_22)*inv_n3[2]*
                    fac2/(12.0*PI*n[3]*n[3]);


   dphi_drb_loc[3] = n[0]*inv_n3[1] + (n[1]*n[2] - DOT_12)*inv_n3[2] 
                + fac1*n[2]*n[2]*n[2]*fac*fac*fac*inv_n3[2]/(36.0*PI*n[3]*n[3])
                + fac2*n[2]*n[2]*n[2]*fac*fac*fac*(-inv_n3[2]/(n[3]*n[3]*n[3])+inv_n3[3]/(n[3]*n[3]))/(18.0*PI);
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
                     n[2]*n[iv2]*fac2*inv_n3[2]*fac*fac/(6.0*PI*n[3]*n[3]);
      else dphi_drb_loc[iv2]=0.0;
   }
   return;
}
/*****************************************************************************/
/* d2phi_drb2_delta_rb_FMT4:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                 for the dphi_drb that use Delta_Fn Stencils (all but S3) */

struct RB_Struct d2phi_drb2_delta_rb_FMT4(double *n, int *offset, double *sign, int icomp)
{
  struct RB_Struct tmp;
  double n3sq,n3cb,n2sq,n2cb;
  double inv_n3[5],DOT_22,DOT_12,DOT_nw,fac,fac1,fac2,fac3,vector[NDIM_MAX];
  int idim,i2v,i1v;

  DOT_22 = 0.0;
  DOT_12 = 0.0;
  for (idim = 0; idim<Ndim; idim++) {
    i2v=Nrho_bar_s+Ndim+idim;
    i1v=Nrho_bar_s+idim;
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

  n2sq=n[2]*n[2];
  n2cb=n2sq*n[2];

  if (n[3]>1.e-10){
     fac1 = n[3]-2*(1.-n[3])*log(1.-n[3]);
     fac2 = n[3]+(1.-n[3])*(1.-n[3])*log(1.-n[3]);
     fac3 = -(inv_n3[2]/n3cb) + (inv_n3[3]/n3sq);
     fac = 1-DOT_22/(n[2]*n[2]);

     for (idim = 0; idim<Ndim; idim++) 
       vector[idim] = offset[idim] * Esize_x[idim] * Inv_rad[icomp];

       tmp.S2 = (inv_n3[1]*Inv_4pir[icomp]
         + fac2*inv_n3[2]*(n[2]*fac*fac*fac+3*DOT_22*fac*fac/n[2]+4*DOT_22*DOT_22*fac/(n2cb))/(6.*PI*n3sq));

     for (idim = 0; idim<Ndim; idim++){
        i2v=Nrho_bar_s+Ndim+idim;
	tmp.S2 += vector[idim]
              *n[i2v]*(fac2*inv_n3[2]/(6.*PI*n3sq))*(fac*fac+4*DOT_22*fac/(n2sq));
     }

       tmp.S3 = (Inv_4pirsq[icomp]*inv_n3[1] 
         + Inv_4pir[icomp]*n[2]*inv_n3[2]
         + n[1]*inv_n3[2]+(n2sq*fac*fac*fac+2*DOT_22*fac*fac)*(
      (2.*fac2*fac3) + (fac1*inv_n3[2]/n3sq))/(12.*PI)  );

     for (idim = 0; idim<Ndim; idim++){
        i2v=Nrho_bar_s+Ndim+idim;
        i1v=Nrho_bar_s+idim;
	tmp.S3 +=  vector[idim] * ( 
	      n[i2v]*Inv_4pir[icomp]*inv_n3[2] + n[i1v]*inv_n3[2] 
              + n[2]*n[i2v]*fac*fac*( 
              (2.*fac2*fac3) + (fac1*inv_n3[2]/n3sq))/(6.*PI) ) ;
     }

     tmp.S0 = 0.0;  
     tmp.S1 = inv_n3[1];  

     for (idim = 0; idim<Ndim; idim++){
       tmp.V1[idim] = sign[idim]*(inv_n3[1]) * vector[idim];
     }

    DOT_nw=0.0;
    for(idim=0; idim<Ndim; idim++) {
        i2v=Nrho_bar_s+Ndim+idim;
        DOT_nw += n[i2v]*vector[idim];
    }
    
     for (idim = 0; idim<Ndim; idim++){
       i2v=Nrho_bar_s+Ndim+idim;
/*       i1v=Nrho_bar_s+idim;*/

       tmp.V2[idim] = sign[idim]*
                    ( -(n[i2v]*inv_n3[2]*fac2/(6.*PI*n3sq))*(fac*fac-4*DOT_22*fac/n[2]+4*DOT_nw*fac/n[2])
                    - (-Inv_4pir[icomp]*inv_n3[1] - n[2]*fac*fac*fac2*inv_n3[2]/(6.*PI*n3sq))*vector[idim]);
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
/* d2phi_drb2_theta_rb_FMT4:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                    for the dphi_drb that use Theta_Fn Stencils (S3)      */

struct RB_Struct d2phi_drb2_theta_rb_FMT4(double *n)
{
  struct RB_Struct tmp;
  double n3sq,n3cb,n3_4th,n2sq,n2cb;
  double inv_n3[5],DOT_22,DOT_12,fac,fac1,fac2,fac3;
  int idim,i2v,i1v;

  DOT_22 = 0.0;
  DOT_12 = 0.0;
  for (idim = 0; idim<Ndim; idim++) {
    i2v=Nrho_bar_s+Ndim+idim;
    i1v=Nrho_bar_s+idim;
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
    fac = 1-DOT_22/(n[2]*n[2]);

    tmp.S2 = (
           n[1]*inv_n3[2] 
         + (n2sq*fac*fac*fac+2*fac*fac*DOT_22)*(fac1*inv_n3[2]/n3sq + 2.*fac2*fac3)/(12*PI));

      tmp.S3 = (
           n[0]*inv_n3[2] 
           + 2*(n[1]*n[2]-DOT_12)*inv_n3[3]
           + ((n2cb*fac*fac*fac)/(36.*PI))*(
           + (((2.*log(1.-n[3])+3)*inv_n3[2])/n3sq) 
           + (4.*fac1*fac3) + (fac2*(6.*inv_n3[2]/n3_4th 
           + 6.*inv_n3[4]/n3sq - 8.*inv_n3[3]/n3cb)) )  );

    tmp.S0 = inv_n3[1];             
    tmp.S1 = n[2] * inv_n3[2];    

    for (idim = 0; idim<Ndim; idim++){
      i2v=Nrho_bar_s+idim+Ndim;
      tmp.V1[idim] = - n[i2v]*inv_n3[2];
    }

    for (idim = 0; idim<Ndim; idim++){
      i1v=Nrho_bar_s+idim;
      i2v=Nrho_bar_s+idim+Ndim;
      tmp.V2[idim] = - (n[i1v]*inv_n3[2]
        + n[2]*n[i2v]*fac*fac*((fac1*inv_n3[2]/n3sq)+(2.*fac2*fac3))/(6.*PI) );
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
