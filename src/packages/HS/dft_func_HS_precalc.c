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
 *  FILE: dft_func_HS_precalc.c
 *
 *  This file implements calculations that are needed before the matrix fill 
 *  associated with hard sphere terms can be completed.
 *
 */

#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"
#include "HSpkg.h"

void precalc_HS(struct RB_Struct *dphi_drb)
{
  if (Type_func==FMT1){
     for (inode_box=0;inode_box<Nnodes_box; inode_box++){
       calc_FMT_derivatives(&FMT1_1stderiv,inode_box,x,dphi_drb);
     }
  }
  else if (Type_func==FMT2){
     for (inode_box=0;inode_box<Nnodes_box; inode_box++)
       calc_FMT_derivatives(&FMT2_1stderiv,inode_box,x,dphi_drb);
  }
  else if (Type_func==FMT3){
    for (inode_box=0;inode_box<Nnodes_box; inode_box++)
      calc_FMT_derivatives(&FMT3_1stderiv,inode_box,x,dphi_drb);
  }
   return;
}
/****************************************************************************/
void calc_FMT_derivatives(void(*fp_FMTderiv)(double *,double,double,double *,double *),
                     int inode_box,double **x, struct RB_Struct *dphi_drb)
{
  double n[4+2*NDIM_MAX], inv_n3[5],dphi_drb_loc[4+2*NDIM_MAX];
  double DOT_22,DOT_12;
  int iunk,idim;

  inv_n3[0]=inv_n3[1]=inv_n3[2]=inv_n3[3]=inv_n3[4]=0.0;

  iunk = Phys2Unk_first[RHOBAR_ROSEN];      
  n[3] = x[iunk][inode_box];                
  n[2] = x[iunk+1][inode_box];
  n[1] = x[iunk+2][inode_box];
  n[0] = x[iunk+3][inode_box];        
  
  for (idim = 0; idim<Ndim; idim++) {
    n[Nrho_bar_s+idim+Ndim] = x[iunk+Nrho_bar_s+idim][inode_box];   
    n[Nrho_bar_s+idim] = x[iunk+Nrho_bar_s+Ndim+idim][inode_box];   
  }

  inv_n3[0]= (1.0 - n[3]);
  inv_n3[1] = 1.0 / inv_n3[0];
  inv_n3[2] = inv_n3[1]*inv_n3[1];
  inv_n3[3] = inv_n3[2]*inv_n3[1];
  inv_n3[4] = inv_n3[3]*inv_n3[1];
                       
  DOT_22 = 0.0;        
  DOT_12 = 0.0;        
  for (idim = 0; idim < Ndim; idim++) {
      DOT_22 += n[Nrho_bar_s+Ndim+idim] * n[Nrho_bar_s+Ndim+idim];
      DOT_12 += n[Nrho_bar_s+idim] * n[Nrho_bar_s+Ndim+idim];
  }
  
  (*fp_FMTderiv)(n,DOT_12,DOT_22,inv_n3,dphi_drb_loc);
  dphi_drb[inode_box].S0=dphi_drb_loc[0];
  dphi_drb[inode_box].S1=dphi_drb_loc[1];
  dphi_drb[inode_box].S2=dphi_drb_loc[2];
  dphi_drb[inode_box].S3=dphi_drb_loc[3];
  for (idim=0;idim<Ndim;idim++){
        dphi_drb[inode_box].V1[idim]=dphi_drb_loc[4+idim];
        dphi_drb[inode_box].V2[idim]=dphi_drb_loc[4+Ndim+idim];
  }
  return;
}

