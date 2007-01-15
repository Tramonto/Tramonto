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

/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

/*
 *  FILE: dft_fill_pde.c
 *
 *  This file contains routines to set up finite element weights for
 *  various PDEs that may be solved.
 */

#include "dft_fill_pde.h"

/****************************************************************************/
void set_fem_weights(double **wt_laplace_ptr, double **wt_source_ptr)
/* This routine sets the fem weights for the laplace and source terms */

{
  int num_nodes;
  double dy_x, dx_y, evol, dxy_z, dxz_y, dyz_x;

  num_nodes = POW_INT(3, Ndim);

  if (*wt_laplace_ptr == NULL) {
    *wt_laplace_ptr = (double *) array_alloc(1, num_nodes, sizeof(double));
    *wt_source_ptr  = (double *) array_alloc(1, num_nodes, sizeof(double)); 
  }

  if (Ndim == 1) {
    /*Finite element weights for Laplace terms on rectangular grid */
    (*wt_laplace_ptr)[0] = -1.0/Esize_x[0];
    (*wt_laplace_ptr)[1] =  2.0/Esize_x[0];
    (*wt_laplace_ptr)[2] = -1.0/Esize_x[0];

    /*Finite element weights for source terms on rectangular grid */
    (*wt_source_ptr)[0]  = Esize_x[0]/6.0;
    (*wt_source_ptr)[1]  = Esize_x[0]*2.0/3.0;
    (*wt_source_ptr)[2]  = Esize_x[0]/6.0;
  }

  else  if (Ndim == 2) {
    dx_y = Esize_x[0]/Esize_x[1];
    dy_x = Esize_x[1]/Esize_x[0];

    /*Finite element weights for Laplace terms on rectangular grid */
    /* corners */
    (*wt_laplace_ptr)[0] =  (*wt_laplace_ptr)[2] = (*wt_laplace_ptr)[6] =
      (*wt_laplace_ptr)[8] = - (dy_x + dx_y)/6.0;
    /* top and bottom */
    (*wt_laplace_ptr)[1] = (*wt_laplace_ptr)[7] = (dy_x - 2.0*dx_y)/3.0;
    /* left and right */
    (*wt_laplace_ptr)[3] = (*wt_laplace_ptr)[5] = (-2.0*dy_x + dx_y)/3.0;
    /* center */
    (*wt_laplace_ptr)[4] = 4.0*(dy_x + dx_y)/3.0;;

    /*Finite element weights for source terms on rectangular grid */
    evol = Esize_x[0]*Esize_x[1];
    /* corners */
    (*wt_source_ptr)[0] =  (*wt_source_ptr)[2] = (*wt_source_ptr)[6] =
      (*wt_source_ptr)[8] = evol / 36.0;
    /* top and bottom, left and right */
    (*wt_source_ptr)[1] = (*wt_source_ptr)[7] = (*wt_source_ptr)[3] = 
      (*wt_source_ptr)[5] = evol / 9.0;
    /* center */
    (*wt_source_ptr)[4] = 4.0 * evol / 9.0;

  }

  else {

    /*Finite element weights for Laplace terms on rectangular grid */
    dxy_z = Esize_x[0]*Esize_x[1]/Esize_x[2];
    dyz_x = Esize_x[1]*Esize_x[2]/Esize_x[0];
    dxz_y = Esize_x[2]*Esize_x[0]/Esize_x[1];
    /* center */
    (*wt_laplace_ptr)[13] = 8.0*(dyz_x + dxz_y + dxy_z)/9.0;
    /* +- ix */
    (*wt_laplace_ptr)[12] = (*wt_laplace_ptr)[14] = 2.0*(-2*dyz_x + dxz_y + dxy_z)/9.0;
    /* +- iy */
    (*wt_laplace_ptr)[10] = (*wt_laplace_ptr)[16] = 2.0*(dyz_x - 2.0*dxz_y +dxy_z)/9.0;
    /* +- iz */
    (*wt_laplace_ptr)[4] = (*wt_laplace_ptr)[22] = 2.0*(dyz_x + dxz_y -2.0*dxy_z)/9.0;
    /* +- ix and +- iy */
    (*wt_laplace_ptr)[9] = (*wt_laplace_ptr)[11] = (*wt_laplace_ptr)[15]
        =  (*wt_laplace_ptr)[17] = (-dyz_x - dxz_y + 0.5*dxy_z)/9.0;
    /* +- ix and +- iz */
    (*wt_laplace_ptr)[3] = (*wt_laplace_ptr)[5] = (*wt_laplace_ptr)[21]
        =  (*wt_laplace_ptr)[23] = (-dyz_x + 0.5*dxz_y - dxy_z)/9.0;
    /* +- iy and +- iz */
    (*wt_laplace_ptr)[1] = (*wt_laplace_ptr)[7] = (*wt_laplace_ptr)[19]
        =  (*wt_laplace_ptr)[25] = (0.5*dyz_x - dxz_y - dxy_z)/9.0;
    /* (corners) +- ix and +- iy and +- iz */
    (*wt_laplace_ptr)[ 0] = (*wt_laplace_ptr)[ 2] = (*wt_laplace_ptr)[6]
        =  (*wt_laplace_ptr)[ 8] = (*wt_laplace_ptr)[18] = (*wt_laplace_ptr)[20]
        =  (*wt_laplace_ptr)[24] = (*wt_laplace_ptr)[26] = -(dyz_x + dxz_y + dxy_z)/36.0;

    /*Finite element weights for source terms on rectangular grid */
    evol = Esize_x[0]*Esize_x[1]*Esize_x[2];

    /* center */
    (*wt_source_ptr)[13] = 64.0 * evol/216.0;

    /* One away from center */
    (*wt_laplace_ptr)[12] = (*wt_laplace_ptr)[14] = 
      (*wt_laplace_ptr)[10] = (*wt_laplace_ptr)[16] = 
      (*wt_laplace_ptr)[ 4] = (*wt_laplace_ptr)[22] = 16.0 * evol/216.0;

    /* two away from center */
    (*wt_source_ptr)[ 9] = (*wt_source_ptr)[11] = (*wt_source_ptr)[15]
        =  (*wt_source_ptr)[17] 
        =  (*wt_source_ptr)[ 3] = (*wt_source_ptr)[ 5] = (*wt_source_ptr)[21]
        =  (*wt_source_ptr)[23]
        =  (*wt_source_ptr)[ 1] = (*wt_source_ptr)[ 7] = (*wt_source_ptr)[19]
        =  (*wt_source_ptr)[25] = 4.0 * evol/216.0;

    /* (corners) +- ix and +- iy and +- iz */
    (*wt_source_ptr)[ 0] = (*wt_source_ptr)[ 2] = (*wt_source_ptr)[ 6]
        =  (*wt_source_ptr)[ 8] = (*wt_source_ptr)[18] = (*wt_source_ptr)[20]
        =  (*wt_source_ptr)[24] =  (*wt_source_ptr)[26] = evol/216.0;
  }

}
/****************************************************************************/
void set_fem_1el_weights(double **wt_lp_1el_ptr, double **wt_s_1el_ptr,
                         int ***elem_permute)
/* This routine sets the fem weights for the laplace and source terms */

{
  int iln;
  double dy_x, dx_y, evol, dxy_z, dxz_y, dyz_x;

  if (*wt_lp_1el_ptr == NULL) {
    *wt_lp_1el_ptr = (double *) array_alloc(1, Nnodes_per_el_V,sizeof(double));
    *wt_s_1el_ptr  = (double *) array_alloc(1, Nnodes_per_el_V,sizeof(double));
    *elem_permute  = (int **) array_alloc(2, Nnodes_per_el_V,Nnodes_per_el_V, sizeof(int));
  }

  if (Ndim == 1) {
    /*Finite element weights for Laplace terms on rectangular grid */
    (*wt_lp_1el_ptr)[0] =  1.0/Esize_x[0];
    (*wt_lp_1el_ptr)[1] = -1.0/Esize_x[0];

    /*Finite element weights for source terms on rectangular grid */
    (*wt_s_1el_ptr)[0]  = Esize_x[0]/3.0;
    (*wt_s_1el_ptr)[1]  = Esize_x[0]/6.0;

    (*elem_permute)[0][0]=0;
    (*elem_permute)[1][0]=1;
    (*elem_permute)[0][1]=1;
    (*elem_permute)[1][1]=0;

  }

  else  if (Ndim == 2) {
    dx_y = Esize_x[0]/Esize_x[1];
    dy_x = Esize_x[1]/Esize_x[0];

    /*Finite element weights for Laplace terms on rectangular grid */
    (*wt_lp_1el_ptr)[0] = (dy_x + dx_y) /3.0;
    (*wt_lp_1el_ptr)[1] = (-dy_x + dx_y/2.0) /3.0;
    (*wt_lp_1el_ptr)[2] = ( dy_x/2.0 - dx_y) /3.0;
    (*wt_lp_1el_ptr)[3] = -(dy_x + dx_y) /6.0;

    /*Finite element weights for source terms on rectangular grid */
    evol = Esize_x[0]*Esize_x[1];
    (*wt_s_1el_ptr)[0] = evol/9.0;
    (*wt_s_1el_ptr)[1] = evol/18.0;
    (*wt_s_1el_ptr)[2] = evol/18.0;
    (*wt_s_1el_ptr)[3] = evol/36.0;

    for (iln=0; iln<Nnodes_per_el_V; iln++) {
       (*elem_permute)[iln][0] = iln;
       (*elem_permute)[iln][1] = (13-iln)%4;
       (*elem_permute)[iln][2] = (iln+6)%4;
       (*elem_permute)[iln][3] = 3 - iln;
    }
  }

  else {

    /*Finite element weights for Laplace terms on rectangular grid */
    dxy_z = Esize_x[0]*Esize_x[1]/Esize_x[2];
    dyz_x = Esize_x[1]*Esize_x[2]/Esize_x[0];
    dxz_y = Esize_x[2]*Esize_x[0]/Esize_x[1];

    (*wt_lp_1el_ptr)[0] = (dyz_x + dxz_y + dxy_z) /9.0;
    (*wt_lp_1el_ptr)[1] = (-dyz_x + dxz_y/2.0 + dxy_z/2.0) /9.0;
    (*wt_lp_1el_ptr)[2] = (dyz_x/2.0 - dxz_y + dxy_z/2.0) /9.0;
    (*wt_lp_1el_ptr)[3] = (-dyz_x/2.0 - dxz_y/2.0 + dxy_z/4.0) /9.0;
    (*wt_lp_1el_ptr)[4] = (dyz_x/2.0 + dxz_y/2.0 - dxy_z) /9.0;
    (*wt_lp_1el_ptr)[5] = (-dyz_x/2.0 + dxz_y/4.0 - dxy_z/2.0) /9.0;
    (*wt_lp_1el_ptr)[6] = (dyz_x/4.0 - dxz_y/2.0 - dxy_z/2.0) /9.0;
    (*wt_lp_1el_ptr)[7] = -(dyz_x + dxz_y + dxy_z) /36.0;
    
    /*Finite element weights for source terms on rectangular grid */
    evol = Esize_x[0]*Esize_x[1]*Esize_x[2];
    (*wt_s_1el_ptr)[0] = evol/27.0;
    (*wt_s_1el_ptr)[1] = evol/54.0;
    (*wt_s_1el_ptr)[2] = evol/54.0;
    (*wt_s_1el_ptr)[3] = evol/108.0;
    (*wt_s_1el_ptr)[4] = evol/54.0;
    (*wt_s_1el_ptr)[5] = evol/108.0;
    (*wt_s_1el_ptr)[6] = evol/108.0;
    (*wt_s_1el_ptr)[7] = evol/216.0;

    for (iln=0; iln<4; iln++) {
       (*elem_permute)[iln][0] = (*elem_permute)[iln+4][0+4] = iln;
       (*elem_permute)[iln][1] = (*elem_permute)[iln+4][1+4] = (13-iln)%4;
       (*elem_permute)[iln][2] = (*elem_permute)[iln+4][2+4] = (iln+6)%4;
       (*elem_permute)[iln][3] = (*elem_permute)[iln+4][3+4] = 3-iln; 
    }
    for (iln=4; iln<8; iln++) {
       (*elem_permute)[iln][0] = (*elem_permute)[iln-4][0+4] = iln;
       (*elem_permute)[iln][1] = (*elem_permute)[iln-4][1+4] = (13-iln)%4 + 4;
       (*elem_permute)[iln][2] = (*elem_permute)[iln-4][2+4] = (iln+6)%4 + 4;
       (*elem_permute)[iln][3] = (*elem_permute)[iln-4][3+4] = 7- (iln-4); 
    }
  }
}
/*****************************************************************************/
/* basis_fn_calc: calcs phi and grad_phi for regular hex element */
void basis_fn_calc(double **phi, double ***grad_phi, double *evol)
{
   int iln, igp;
   double gp[2];

   gp[0] = (1.0 - 5.773502691896258e-01) / 2.0;
   gp[1] = (1.0 + 5.773502691896258e-01) / 2.0;

   if (Ndim == 3) {

      *evol = Esize_x[0] * Esize_x[1] * Esize_x[2];
      
      for (iln=0; iln<8; iln++) {
        for (igp=0; igp<8; igp++) {
          phi[iln][igp] = gp[(iln+igp)%2] * gp[((iln/2) + (igp/2))%2]
                                          * gp[((iln/4) + (igp/4))%2];
          grad_phi[iln][igp][0] = gp[((iln/2) + (igp/2))%2]
                                * gp[((iln/4) + (igp/4))%2] / Esize_x[0];
          if (((iln+igp)%2)==0) grad_phi[iln][igp][0] *= -1.0;
          grad_phi[iln][igp][1] = gp[(iln+igp)%2]  
                                * gp[((iln/4) + (igp/4))%2] / Esize_x[1];
          if (((iln/2) + (igp/2))%2==0) grad_phi[iln][igp][1] *= -1.0;
          grad_phi[iln][igp][2] = gp[(iln+igp)%2]  
                                * gp[((iln/2) + (igp/2))%2] / Esize_x[2];
          if (((iln/4) + (igp/4))%2==0) grad_phi[iln][igp][2] *= -1.0;
        }
      }
   }
   else {
     if (Proc==0) printf("ERROR basis_fn_calc: Only for 3D\n");
     exit(-1);
   }
}
/************************************************************************/
