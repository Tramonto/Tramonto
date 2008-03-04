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
 *  FILE: dft_fill.c
 *
 *  This file contains the fill of the residual equations and Jacobian
 *  matrix.
 */

#include "dft_fill_main.h"

/****************************************************************************/
void fill_resid_and_matrix (double **x, int iter, int resid_only_flag,int unk_flag)
{
 /*
  * Local variable declarations
  */

  char   *yo = "fill_resid_and_matrix";
  int     loc_inode, inode_box,ijk_box[3],iunk,junk,iunk_start,iunk_end;
  int     mesh_coarsen_flag_i,switch_constmatrix;
  struct  RB_Struct *dphi_drb=NULL;
  double *resid_unk;

  if (Proc == 0 && !resid_only_flag && Iwrite != NO_SCREEN) printf("\n\t%s: Doing fill of residual and matrix\n",yo);
  resid_unk = (double *) array_alloc (1, Nunk_per_node, sizeof(double));

  /* pre calculations required for the Hard sphere (FMT) functionals only*/

  if (Type_func !=NONE){
     dphi_drb = (struct RB_Struct *) array_alloc
                    (1, Nnodes_box, sizeof(struct RB_Struct));
     FMT1stDeriv_switch(inode_box,x,dphi_drb);
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

      if (mesh_coarsen_flag_i == FLAG_1DBC) load_coarse_node_1dim(loc_inode,inode_box,ijk_box,iunk,x,resid_only_flag);

      else if (mesh_coarsen_flag_i < 0 && 
               mesh_coarsen_flag_i != FLAG_BULK && 
               mesh_coarsen_flag_i != FLAG_PBELEC) load_coarse_node_Ndim(loc_inode,inode_box,iunk,x,resid_only_flag);

      else{
          /*switch_constmatrix=FALSE;
          if (iter>1 && resid_only_flag==FALSE && Constant_row_flag[Unk2Phys[iunk]]==TRUE) {
             resid_only_flag=TRUE; switch_constmatrix=TRUE;
          }*/
          load_standard_node(loc_inode,inode_box,ijk_box,iunk,x,dphi_drb,
                               resid_unk,mesh_coarsen_flag_i,resid_only_flag);
/*          if (switch_constmatrix) resid_only_flag=FALSE;*/
      }

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
   int izone,i,icomp;
                                           
   /* IZONE: izone is the zone number that controls the quadrature scheme to be used. */

   izone = mesh_coarsen_flag_i;

   switch(Unk2Phys[iunk]){
       case DENSITY: 
             if (L_HSperturbation && Type_poly != WJDC)
                resid_unk[iunk]=load_euler_lagrange(iunk,loc_inode,inode_box,ijk_box,izone,x,
                                                    dphi_drb,mesh_coarsen_flag_i,resid_only_flag);
             else if(Type_poly==CMS)                              
                resid_unk[iunk]=load_CMS_density(iunk,loc_inode,inode_box,x,resid_only_flag);
             else if(Type_poly==WJDC)
                resid_unk[iunk]=load_WJDC_density(iunk,loc_inode,inode_box,x,resid_only_flag);
             break;

       case HSRHOBAR: 
          if (iunk == Phys2Unk_first[HSRHOBAR]){
             resid_unk[iunk]=load_rho_bar_s(THETA_FN_R,x,iunk,loc_inode,inode_box,izone,ijk_box, resid_only_flag);
          }
          else if (iunk < Phys2Unk_first[HSRHOBAR]+Nrho_bar_s){
             resid_unk[iunk]=load_rho_bar_s(DELTA_FN_R,x,iunk,loc_inode,inode_box,izone,ijk_box, resid_only_flag);
          }
          else if (iunk >= Phys2Unk_first[HSRHOBAR]+Nrho_bar_s){
              resid_unk[iunk]=load_rho_bar_v(x,iunk,loc_inode,inode_box,izone,ijk_box, resid_only_flag);
          }
        break;

       case MF_EQ: 
           icomp=iunk-Phys2Unk_first[MF_EQ];
           resid_unk[iunk]=load_mean_field(THETA_PAIRPOT_RCUT,iunk,loc_inode,
                                                   icomp,izone,ijk_box, x, resid_only_flag);
           break;

       case POISSON: resid_unk[iunk]=load_poisson_control(iunk,loc_inode,inode_box,ijk_box,x,resid_only_flag); break;

       case DIFFUSION:
            if (Linear_transport)
                resid_unk[iunk]=load_linear_transport_eqn(iunk,loc_inode,inode_box,ijk_box,x,resid_only_flag);
            else
                resid_unk[iunk]=load_nonlinear_transport_eqn(iunk,loc_inode,inode_box,ijk_box,x,resid_only_flag);
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

       case WJDC_FIELD: 
            resid_unk[iunk]=load_WJDC_field(iunk,loc_inode,inode_box,ijk_box,izone,x,                                                    dphi_drb,mesh_coarsen_flag_i,resid_only_flag);
            break;

       case G_CHAIN:
            if (Type_poly==CMS){
               resid_unk[iunk]=load_CMS_Geqns(iunk,loc_inode,inode_box,ijk_box,izone,x,resid_only_flag);
            }
            else if (Type_poly==WJDC){
               resid_unk[iunk]=load_WJDC_Geqns(iunk,loc_inode,inode_box,ijk_box,izone,x,resid_only_flag);
            }
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
       case DENSITY:  printf("Proc=%d: loc_inode=%d of %d (Global val=%d) iunk_rho=%d ", Proc,loc_inode,Nnodes_per_proc,L2G_node[loc_inode],iunk); break;
       case HSRHOBAR: printf("Proc=%d: loc_inode=%d iunk_rbar=%d ", Proc,loc_inode,iunk); break;
       case POISSON:  printf("Proc=%d: loc_inode=%d iunk_poisson=%d ", Proc,loc_inode,iunk); break;
       case DIFFUSION: printf("Proc=%d: loc_inode=%d  iunk_diffusion=%d ",Proc,loc_inode,iunk); break;
       case CAVWTC: printf("Proc=%d: loc_inode=%d  iunk_cavity=%d ",Proc,loc_inode,iunk); break;
       case BONDWTC: printf("Proc=%d: loc_inode=%d  iunk_bondwtc=%d ",Proc,loc_inode,iunk); break;
       case CMS_FIELD: printf("Proc=%d: loc_inode=%d  iunk_cmsfield=%d ",Proc,loc_inode,iunk); break;
       case WJDC_FIELD: printf("Proc=%d: loc_inode=%d  iunk_wjdc_field=%d ",Proc,loc_inode,iunk); break;
       case G_CHAIN: printf("Proc=%d: loc_inode=%d  iunk_Gchain=%d ",Proc,loc_inode,iunk); break;
    }
    printf(" resid=%11.8f \n",resid_unk[iunk]); 

    return;
}
/*************************************************************************************************/
