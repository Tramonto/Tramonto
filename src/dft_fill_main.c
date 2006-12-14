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
 *  FILE: dft_fill.c
 *
 *  This file contains the fill of the residual equations and Jacobian
 *  matrix.
 */

/*#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"*/
#include "dft_fill_main.h"

void print_residuals(int,int,double *);
void load_standard_node(int,int,int *,int,double **,struct  RB_Struct *, double *,int,int);
/****************************************************************************/
void fill_resid_and_matrix (double **x, int iter, int resid_only_flag,int unk_flag)
{
 /*
  * Local variable declarations
  */

  char   *yo = "fill_resid_and_matrix";
  int     loc_inode, inode_box,ijk_box[3],iunk,junk,iunk_start,iunk_end;
  int     mesh_coarsen_flag_i;
  struct  RB_Struct *dphi_drb=NULL;
  double *resid_unk;

  if (Proc == 0 && !resid_only_flag && Iwrite != NO_SCREEN) printf("\n\t%s: Doing fill of residual and matrix\n",yo);
  resid_unk = (double *) array_alloc (1, Nunk_per_node, sizeof(double));

  /* pre calculations required for the Hard sphere (FMT) functionals only*/

  if (Type_func !=NONE){
     dphi_drb = (struct RB_Struct *) array_alloc
                    (1, Nnodes_box, sizeof(struct RB_Struct));

     if (Type_func==FMT1){
        for (inode_box=0;inode_box<Nnodes_box; inode_box++){
          calc_FMT_derivatives(&FMT1_1stderiv,inode_box,x,dphi_drb);
        }
     }
     else if (Type_func==FMT2)
        for (inode_box=0;inode_box<Nnodes_box; inode_box++)
          calc_FMT_derivatives(&FMT2_1stderiv,inode_box,x,dphi_drb);
     else if (Type_func==FMT3)
        for (inode_box=0;inode_box<Nnodes_box; inode_box++)
          calc_FMT_derivatives(&FMT3_1stderiv,inode_box,x,dphi_drb);
  }

  /* for debugging print out profiles on each iteration */
  if (Iwrite==VERBOSE) print_profile_box(x, "dens_iter.dat");

  if (unk_flag == NODAL_FLAG){
      iunk_start = 0;
      iunk_end = Nunk_per_node;
  } 
  else{
      iunk_start = unk_flag;
      iunk_end = unk_flag+1;
  }

  /* Load residuals and matrix */
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {

    /* convert local node to global */

    inode_box = L2B_node[loc_inode];
    node_box_to_ijk_box(inode_box, ijk_box);

    if ( ((Mesh_coarsening != FALSE) && (Nwall_type >0)) || L1D_bc) mesh_coarsen_flag_i = Mesh_coarsen_flag[inode_box];
    else mesh_coarsen_flag_i = 0;

    for (iunk=iunk_start; iunk<iunk_end; iunk++) {

      resid_unk[iunk]=0.0;

      if (mesh_coarsen_flag_i == FLAG_1DBC) load_coarse_node_1dim(loc_inode,inode_box,ijk_box,iunk,x);

      else if (mesh_coarsen_flag_i < 0 && 
               mesh_coarsen_flag_i != FLAG_BULK && 
               mesh_coarsen_flag_i != FLAG_PBELEC) load_coarse_node_Ndim(loc_inode,inode_box,iunk,x);

      else load_standard_node(loc_inode,inode_box,ijk_box,iunk,x,dphi_drb,
                               resid_unk,mesh_coarsen_flag_i,resid_only_flag);

     /* print for debugging purposes call this print routine */ 
     /*print_residuals(loc_inode,iunk,resid_unk);*/

    } /* end of loop over # of unknowns per node */
  } /* end of loop over local nodes */

  if (Type_func != NONE) safe_free((void *) &dphi_drb);
  safe_free((void *) &resid_unk);
  return;
}
/*************************************************************************************************/
/* load_standard_node:  this routine controls the invocation of different physics routines that load the 
                        matrix problem of interest.  */
void load_standard_node(int loc_inode,int inode_box, int *ijk_box, int iunk, double **x,
                        struct  RB_Struct *dphi_drb, double *resid_unk,
                        int mesh_coarsen_flag_i,int resid_only_flag)
{
   int izone,i;
                                           
   /* IZONE: izone is the zone number that controls the quadrature scheme to be used. */

   izone = mesh_coarsen_flag_i;

   switch(Unk2Phys[iunk]){
       case DENSITY: 
             if (Type_poly==NONE || Type_poly==WTC) 
                resid_unk[iunk]=load_euler_lagrange(iunk,loc_inode,inode_box,ijk_box,izone,x,
                                                    dphi_drb,mesh_coarsen_flag_i,resid_only_flag);
             else                              
                resid_unk[iunk]=load_CMS_density(iunk,loc_inode,inode_box,x,resid_only_flag);
             break;

       case HSRHOBAR: 
          if (iunk == Phys2Unk_first[HSRHOBAR]){
             resid_unk[iunk]=load_rho_bar_s(THETA_FN,x,iunk,loc_inode,inode_box,izone,ijk_box,
                            resid_only_flag);
          }
          else if (iunk < Phys2Unk_first[HSRHOBAR]+Nrho_bar_s){
             resid_unk[iunk]=load_rho_bar_s(DELTA_FN,x,iunk,loc_inode,inode_box,izone,ijk_box,
                            resid_only_flag);
          }
          else if (iunk >= Phys2Unk_first[HSRHOBAR]+Nrho_bar_s){
              resid_unk[iunk]=load_rho_bar_v(x,iunk,loc_inode,inode_box,izone,ijk_box, resid_only_flag);
          }
        break;

       case POISSON: resid_unk[iunk]=load_poisson_control(iunk,loc_inode,inode_box,ijk_box,x); break;

       case DIFFUSION:
            if (Linear_transport)
                resid_unk[iunk]=load_linear_transport_eqn(iunk,loc_inode,inode_box,ijk_box,x);
            else
                resid_unk[iunk]=load_nonlinear_transport_eqn(iunk,loc_inode,inode_box,ijk_box,x);
            break;

       case CAVWTC: 
            resid_unk[iunk]=load_cavity_wtc(iunk,loc_inode,inode_box,ijk_box,izone,x,resid_only_flag);
            break;

       case BONDWTC:
            resid_unk[iunk]=load_bond_wtc(iunk,loc_inode,inode_box,ijk_box,izone,x,resid_only_flag);
            break;

       case CMS_FIELD: 
            resid_unk[iunk]=load_CMS_field(iunk,loc_inode,inode_box,ijk_box,izone,x,resid_only_flag);
            break;

       case CMS_G:
            resid_unk[iunk]=load_CMS_Geqns(iunk,loc_inode,inode_box,ijk_box,izone,x,resid_only_flag);
            break;

   }  /* end of physics switch */
   return;
}
/*************************************************************************************************/
/* print_residuals: tool to assist in debugging residual fills.  This routine prints the
   residuals at a given node for each unknown in the problem. Note that these tools are 
   most often used for debugging physics on a single processor, and may need modification
   for multi-processor debugging. */
   void print_residuals(int loc_inode,int iunk,double *resid_unk)
{
               /* note: translate local coordinates to global coordinates with L2G_node[loc_inode]*/
               /* note: translate local coordinates to box coordinates with L2B_node[loc_inode]*/
               /* note: if you want to print to a file you need to modify print statements below */
   FILE *ifp; 
   char filename[20]="resid.out";
               /* also note: to separate parts of the physics constributions (e.g. different parts 
                  of the euler-lagrange equation you will need to modify the return parameters from the
                  various load_ functions. Note that this
                  kind of analysis may require multiple runs and so output to a file is recommended. */

    /* PRINT STATEMENTS FOR PHYSICS DEBUGGING .... CHECK RESIDUALS INDEPENDENTLY  */
    switch(Unk2Phys[iunk]){
       case DENSITY:  printf("loc_inode=%d iunk_rho=%d ", loc_inode,iunk); break;
       case HSRHOBAR: printf("loc_inode=%d iunk_rbar=%d resid=%9.6f ", loc_inode,iunk); break;
       case POISSON:  printf("loc_inode=%d iunk_poisson=%d ", loc_inode,iunk); break;
       case DIFFUSION: printf(" loc_inode=%d  iunk_diffusion=%d ",loc_inode,iunk); break;
       case CAVWTC: printf(" loc_inode=%d  iunk_cavity=%d ",loc_inode,iunk); break;
       case CMS_FIELD: printf(" loc_inode=%d  iunk_cmsfield=%d ",loc_inode,iunk); break;
       case CMS_G: printf(" loc_inode=%d  iunk_cmsG=%d ",loc_inode,iunk); break;
    }
    printf(" resid=%9.6f \n",resid_unk[iunk]); 

    return;
}
/*************************************************************************************************/
