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

#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

/****************************************************************************/
double integrand_hs_freen(int iunk, int inode_box,double **x)
{
  int i;
  double integrand,rho_bar[10];

  for(i=0; i<Phys2Nunk[HSRHOBAR]; i++)
      rho_bar[i] = x[i+Phys2Unk_first[HSRHOBAR]][inode_box];

  integrand = phispt(rho_bar);
  return(integrand);
}
/*******************************************************************************/
double integrand_hs_freen_bulk(int iunk, int inode_box,double **x)
{
  double integrand;

  if (Lsteady_state) integrand = phispt(Rhobar_b_RTF);
  else               integrand = phispt(Rhobar_b);

  return(integrand);
}
/****************************************************************************/
/*phispt: Calculate the hard sphere free energy contribution from
           scaled particle theory at a given node i.                      */
double phispt(double *rho_bar)
{  
  int idim,iv1,iv2;
  double n[4+2*NDIM_MAX];
  double phi_1=0.0,phi_2=0.0,phi_3,dot_12,dot_22;
   printf("in phistp\n"); 
  n[0] = rho_bar[3];
  n[1] = rho_bar[2]; 
  n[2] = rho_bar[1];
  n[3] = rho_bar[0];
  for (idim=0; idim<Ndim; idim++){
    iv1=Nrho_bar_s+idim;
    iv2=Nrho_bar_s+Ndim+idim;
    n[iv1] = rho_bar[iv2];
    n[iv2] = rho_bar[iv1];
  }
  dot_12 = 0.0;
  dot_22 = 0.0;
  for (idim = 0; idim<Ndim; idim++){
    iv1=Nrho_bar_s+idim;
    iv2=Nrho_bar_s+Ndim+idim;
    dot_12 += n[iv1]*n[iv2];
    dot_22 += n[iv2]*n[iv2];
  }

   printf("point 2 in phistp\n"); 
  if (n[3] < 1.0 && n[2] > 0.0 && n[3]>1.e-10){

     phi_1 = -n[0]*log(1.0-n[3]);
     phi_2 = (n[1]*n[2] - dot_12)/(1.0-n[3]);
     if (Type_func==FMT1){
        phi_3 = (n[2]*n[2]*n[2] - 3.0*n[2]*dot_22)/(24.0*PI*(1.0-n[3])*(1.0-n[3]));
     }
     else if (Type_func==FMT2){
        phi_3 = POW_DOUBLE_INT(n[2]-dot_22/n[2],3)/(24.0*PI*(1.0-n[3])*(1.0-n[3]));
     }
     else if (Type_func==FMT3){
        phi_3 = (n[2]*n[2]*n[2] - 3.0*n[2]*dot_22)*
                 (n[3]+(1.0-n[3])*(1.0-n[3])*log(1.0-n[3]))/
                    (36.0*PI*n[3]*n[3]*(1.0-n[3])*(1.0-n[3]));
     }
     else{
       printf("problem with type of HS FUNCTIONAL");
       exit(-1);
     }
/*   printf("point 3 in phistp  %9.6f %9.6f %9.6f \n",phi_1,phi_2,phi_3); */
     return(phi_1 + phi_2 + phi_3);
  }
  else {
/*   printf("point 4 in phistp\n"); */
 return(0.0);
  }
}
/****************************************************************************/

