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
 *  FILE: dft_rhobar.c
 *
 *  This file contains calculations of rhobar on the mesh
 *  this is a pre-calculation step for the fill.
 */
#include <stdio.h>
#include <assert.h>
#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"
#define HIT_FLAG 999

/****************************************************************************/
/* pre_calc_rho_bar:  all rho_bars are calculated at each wall & fluid node*/

void pre_calc_rho_bar(struct RB_Struct *rho_bar, double *x,
                      int fill_flag,int iter,
                      double ***jac_weights_hs, int ***jac_columns_hs,
                      double *fill_time)

{
  int isten, idim, icomp,iunk;
  int izone=0, inode_end, ilist=0,izone_mesh;
  int inode_box,ijk_box[3],ijk[3],inode_sten_box;
  int *offset, inode_tmp;
  double weight, rho=0;
  struct RB_Struct *rb_inode=NULL;
  struct Stencil_Struct *sten;
  int **sten_offset, reflect_flag[NDIM_MAX];
  double *sten_weight;
  int jsten_delta=0, jsten_theta=0;
  int mesh_coarsen_flag_i;
  double *rhobar3_old;

  /* calculate rho_bars at all fluid and wall nodes for preprocessing of fill
     calculate one rho_bar for post processing - calculation of omega. */

  if (fill_flag == MSR_PREPROCESS) inode_end = Nnodes_per_proc;
  else  inode_end = Nnodes_box;

  for (inode_tmp=0; inode_tmp < inode_end; inode_tmp++) {

    if (fill_time != NULL && B2L_node[inode_tmp] !=-1) 
           fill_time[B2L_node[inode_tmp]] = - MPI_Wtime();

    if (fill_flag == MSR_PREPROCESS)
      inode_box = L2B_node[inode_tmp];
    else
      inode_box = inode_tmp;

/*    if (Mesh_coarsening != FALSE && Nwall_type >0 || L1D_bc)
	mesh_coarsen_flag_i = Mesh_coarsen_flag[inode_box];
    else
       mesh_coarsen_flag_i = 0;


    if (fill_flag == MSR_PREPROCESS ||  mesh_coarsen_flag_i >= 0){ */
                                    /* proceed with usual calculation if this node
                                       has not been coarsened out of the mesh */


    node_box_to_ijk_box(inode_box, ijk_box);
    ijk_box_to_ijk(ijk_box,ijk);

    if (fill_flag ==MSR_PREPROCESS || B2L_1stencil[inode_box]  != -1) {

    /* initialize counters */
    if (fill_flag == RHOBAR_JAC_SAVE) {
       jsten_delta = 1;
       jsten_theta = 1;
    }

    /* Set pointer to rho_bars at this node, initialize all rho_bars to 0.0 */

    if (fill_flag != MSR_PREPROCESS && fill_flag != RHOBAR_JAC_SAVE) {
      rb_inode = &(rho_bar[B2L_1stencil[inode_box]]);

      rb_inode->S0 = 0.0;
      rb_inode->S1 = 0.0;
      rb_inode->S2 = 0.0;
      rb_inode->S3 = 0.0;
      for (idim=0; idim<Ndim; idim++) {
         rb_inode->V1[idim] = 0.0;
         rb_inode->V2[idim] = 0.0;
      }
    }

    /* Identify the quadrature zone where this node is located.  */
    /* Use most accurate stencil, unless computing coarsened     */
    /* Jacobian or Mesh (i.e. residual) entries                  */

    if (fill_flag != RHOBAR_JAC_SAVE && Mesh_coarsening==FALSE)
       izone = 0;
    else {
       izone_mesh = Nodes_to_zone[inode_box];


       if (fill_flag == RHOBAR_JAC_SAVE){
          if (Coarser_jac == 0 || izone_mesh == Nzone-1) izone = izone_mesh;
          else{
             if (Coarser_jac == 1){
                 if (izone_mesh ==0 ) izone = izone_mesh + 1;  
                 else                 izone = izone_mesh;
             }
             else if (Coarser_jac == 2) izone = izone_mesh + 1;
             else if (Coarser_jac == 3) izone = Nzone-1;
             else if (Coarser_jac == 4) izone = Nzone-2;
             else if (Coarser_jac == 5) izone = Nzone-1;
          }
       }
       else izone = Nodes_to_zone[inode_box];
    }
    
    /*
     * Loop over Delta_Fn Stencil to assemble the following rho_bars: 
     * S0, S1, S2, V1[Ndim], V2[Ndim]
     */
 
    for (icomp=0; icomp<Ncomp; icomp++) {
      if (Nlists_HW == 1 || Nlists_HW == 2) ilist = 0;
      else if (Nlists_HW > 2)  ilist = icomp;

      sten = &(Stencil[DELTA_FN][izone][icomp]);
      sten_offset = sten->Offset;
      sten_weight = sten->Weight;

      for (isten = 0; isten < sten->Length; isten++) {

        /* get offset and weight for this stencil point */

        offset = sten_offset[isten];
        weight = sten_weight[isten];

        /* Find in the Stencil position on overall mesh */

        inode_sten_box = offset_to_node_box(ijk_box, offset, reflect_flag); 

        /*
         * Flag value of -1 indicates bulk boundary conditiom
         * Otherwise, translate stencil position to reduced mesh (wall & fluid)
         * If this point isn't on reduced mesh (flag value -1), or is in the
         * wall, set rho=0.0.
         * Otherwise, get the rho value from the x array of unknowns
         */

        if      (inode_sten_box == -1) rho = Rho_b[icomp];
        else if (inode_sten_box == -3) {
              if (Lsteady_state==TRUE) rho = Rho_b_LBB[icomp];
              else                     rho = Rho_coex[1];
        }
        else if (inode_sten_box == -4) {
              if (Lsteady_state==TRUE) rho = Rho_b_RTF[icomp];
              else                     rho = Rho_coex[0];
        }
        else if (inode_sten_box == -2) rho = 0.0;
        else {

          /*
           * If we are doing Hard Walls, and this is a boundary node, 
           * modify the weight to only include contributions from fluid elements
           * and not wall elements (where the density is 0 but, through linear
           * interpolation, would otherwise give a nonzero contribution).
           */

          if (fill_flag != MSR_PREPROCESS){

             if (Lhard_surf) 
               if (Nodes_2_boundary_wall[ilist][inode_sten_box]!=-1) {
                   weight = HW_boundary_weight
                        (icomp,ilist,sten->HW_Weight[isten], inode_sten_box,reflect_flag);
               }

             if (Zero_density_TF[inode_sten_box][icomp]) {
                rho = 0.0;
             }
             else {
               iunk=Unk_start_eq[DENSITY]+icomp;
               if (fill_flag != RHOBAR_JAC_SAVE) {
		 assert( B2L_unknowns[loc_find(iunk,inode_sten_box,BOX)] != -1 );
		 rho = x[B2L_unknowns[loc_find(iunk,inode_sten_box,BOX)]];
               }
               else {
                  jac_weights_hs[B2L_1stencil[inode_box]][0][jsten_delta] = weight;
                  jac_columns_hs[B2L_1stencil[inode_box]][0][jsten_delta++] = loc_find(iunk,inode_sten_box,BOX);
                }
             }
          }
        }

        /* For each component, add contributions to each rho_bar */

         if (fill_flag != MSR_PREPROCESS && fill_flag != RHOBAR_JAC_SAVE) {
           rb_inode->S0 += weight * rho  * Inv_4pirsq[icomp];
           rb_inode->S1 += weight * rho  * Inv_4pir[icomp];
           rb_inode->S2 += weight * rho;

           /*
            * The quantity (offset[idim] * Esize_x[idim] / radius[icomp])
            * is the component of the radial unit vector at the stencil
            * node in the idim direction
            */
  
           for (idim=0; idim<Ndim; idim++) {
             rb_inode->V1[idim] += weight * rho
                   * (offset[idim] * Esize_x[idim] * Inv_rad[icomp])
                   * Inv_4pir[icomp];
             rb_inode->V2[idim] += weight * rho
                   * (offset[idim] * Esize_x[idim] * Inv_rad[icomp]);
           }
         }
         else if (fill_flag == MSR_PREPROCESS) {
           if (inode_sten_box >=0){
                if  ( !Zero_density_TF[inode_box][icomp] )
                    B2L_1stencil[inode_sten_box] = HIT_FLAG;
                if ( !Zero_density_TF[inode_sten_box][icomp])
                    B2L_1stencil[inode_box] = HIT_FLAG;

/*                    B2L_1stencil[inode_sten_box] = HIT_FLAG;*/
           }
         }

      } /* End of loop over Delta_Fn Stencil Length*/

      /* Now do the same for the lone rho_bar that needs the Theta_Fn  */

      sten = &(Stencil[THETA_FN][izone][icomp]);
      sten_offset = sten->Offset;
      sten_weight = sten->Weight;

      for (isten = 0; isten < sten->Length; isten++) {

        offset = sten_offset[isten];
        weight = sten_weight[isten];

        inode_sten_box = offset_to_node_box(ijk_box, offset, reflect_flag);
        if (inode_sten_box == -1 || inode_sten_box == -3 || inode_sten_box == -4)
            rho = constant_boundary(Unk_start_eq[DENSITY]+icomp,inode_sten_box);
        else if (inode_sten_box == -2) rho = 0.0;
        else {

          /*
           * If we are doing Hard Walls, and this is a bondary node, 
           * modify the weight to only include contributions from fluid elements
           * and not wall elements (where the density is 0 but, through linear
           * interpolation, would otherwise give a nonzero contribution).
           */

          if (fill_flag != MSR_PREPROCESS){

             if (Lhard_surf) 
               if (Nodes_2_boundary_wall[ilist][inode_sten_box]!=-1) {
                   weight = HW_boundary_weight
                       (icomp,ilist,sten->HW_Weight[isten], inode_sten_box,reflect_flag);
               }

             if (Zero_density_TF[inode_sten_box][icomp]) {
                rho = 0.0;
             }
             else { 
               iunk = Unk_start_eq[DENSITY]+icomp;
               if (fill_flag != RHOBAR_JAC_SAVE) {
		 assert( B2L_unknowns[loc_find(iunk,inode_sten_box,BOX)] != -1 );
		 rho = x[B2L_unknowns[loc_find(iunk,inode_sten_box,BOX)]];
               }
               else {
                  jac_weights_hs[B2L_1stencil[inode_box]][1][jsten_theta] = weight;
                  jac_columns_hs[B2L_1stencil[inode_box]][1][jsten_theta++] = loc_find(iunk,inode_sten_box,BOX);
               }
             } 
          }
        }

        if (fill_flag != MSR_PREPROCESS && fill_flag != RHOBAR_JAC_SAVE){
           rb_inode->S3 += weight * rho;
        }
        else if (fill_flag == MSR_PREPROCESS) {
           if (inode_sten_box >=0){ 
                if  ( !Zero_density_TF[inode_box][icomp] )
                    B2L_1stencil[inode_sten_box] = HIT_FLAG;
                if ( !Zero_density_TF[inode_sten_box][icomp])
                    B2L_1stencil[inode_box] = HIT_FLAG;

/*                    B2L_1stencil[inode_sten_box] = HIT_FLAG;*/
           }
        }

      } /* End of loop over Theta_Fn Stencil Length*/
    }  /* End of loop over icomp */

    /* in first position, put the length of the jac_column entries */

    if (fill_flag == RHOBAR_JAC_SAVE) {
       jac_columns_hs[B2L_1stencil[inode_box]][0][0] = jsten_delta;
       jac_columns_hs[B2L_1stencil[inode_box]][1][0] = jsten_theta;
    }

    /* don't allow rho_bar_3 >= 1.0 */
    if (fill_flag !=MSR_PREPROCESS && fill_flag != RHOBAR_JAC_SAVE) {
       if (iter >1 &&  rb_inode->S3 >=1.0) 
            rb_inode->S3 = 0.5*(1.0-Rhobar3_old[inode_tmp]);
       if (iter >= 1)
          Rhobar3_old[inode_tmp] = rb_inode->S3;
    }
    

  }  /* End check for this node in 1stencil */
  /*}*/ /* End check for Mesh coarsening flag */

    if (fill_time != NULL && B2L_node[inode_tmp] !=-1) 
      fill_time[B2L_node[inode_tmp]] += MPI_Wtime();
  }  /* End of loop over inode */
}

/****************************************************************************/
/* pre_calc_coarse_rho_bar:  here just average appropriately when
     we are doing mesh coarsening so that we do have a value for rho_bar
     on the entire mesh.*/

void pre_calc_coarse_rho_bar(struct RB_Struct *rho_bar)

{
  int inode_box,jnode_box,mesh_coarsen_flag_i,ijk_box[3],ijk_tmp[3];
  int idim,icomp,reflect_flag[3];
  double rb0,rb1,rb2,rb3;
  struct RB_Struct *rb_inode;
  struct RB_Struct *rb_jnode;
  /* the 6 offset patterns for nearest neighbors */
  int offset_idim_pm[18] = {1,0,0,  0,1,0,  0,0,1,  -1,0,0,  0,-1,0,  0,0,-1};
  int *offset_ptr; /* pointer into above */


  for (idim=0; idim<Ndim; idim++) reflect_flag[idim]=FALSE;
 
  for (inode_box=0; inode_box < Nnodes_box; inode_box++) {
    mesh_coarsen_flag_i = Mesh_coarsen_flag[inode_box];

    if (mesh_coarsen_flag_i == FLAG_1DBC) {
        rb_inode = &(rho_bar[B2L_1stencil[inode_box]]);
        node_box_to_ijk_box(inode_box,ijk_box);
        for (idim=0; idim<Ndim; idim++) ijk_tmp[idim]=0;
        ijk_tmp[Grad_dim]=ijk_box[Grad_dim];
        jnode_box = ijk_box_to_node_box(ijk_tmp);
        rb_jnode = &(rho_bar[B2L_1stencil[jnode_box]]);
        rb_inode->S0 = rb_jnode->S0;
        rb_inode->S1 = rb_jnode->S1;
        rb_inode->S2 = rb_jnode->S2;
        rb_inode->S3 = rb_jnode->S3;
        for (idim=0; idim<Ndim; idim++) {
              rb_inode->V1[idim] = rb_jnode->V1[idim];
              rb_inode->V2[idim] = rb_jnode->V2[idim];
        }
    }
    else if (mesh_coarsen_flag_i == FLAG_BULK){
       rb_inode = &(rho_bar[B2L_1stencil[inode_box]]);
       rb_inode->S0 = 0.0;
       rb_inode->S1 = 0.0;
       rb_inode->S2 = 0.0;
       rb_inode->S3 = 0.0;
       for (idim=0; idim<Ndim; idim++) {
          rb_inode->V1[idim] = 0.0;
          rb_inode->V2[idim] = 0.0;
       }
       for (icomp=0; icomp<Ncomp; icomp++) {
            rb_inode->S0 += Rho_b[icomp];
            rb_inode->S1 += Rho_b[icomp]* Sigma_ff[icomp][icomp]/2.0;
            rb_inode->S2 += Rho_b[icomp]* PI* Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
            rb_inode->S3 += Rho_b[icomp]* (PI / 6.0) * Sigma_ff[icomp][icomp]
                    * Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
         }
    }
    else if (mesh_coarsen_flag_i < 0) {  /* compute coarsened rho_bar at this node */


       rb_inode = &(rho_bar[B2L_1stencil[inode_box]]);
       rb_inode->S0 = 0.0;
       rb_inode->S1 = 0.0;
       rb_inode->S2 = 0.0;
       rb_inode->S3 = 0.0;
       for (idim=0; idim<Ndim; idim++) {
          rb_inode->V1[idim] = 0.0;
          rb_inode->V2[idim] = 0.0;
       }

        offset_ptr = &offset_idim_pm[3*(-mesh_coarsen_flag_i - 1)];
        jnode_box = offset_to_node_box(ijk_box, offset_ptr, reflect_flag);

        if (jnode_box >= 0) {
           rb_jnode = &(rho_bar[B2L_1stencil[jnode_box]]);
           rb_inode->S0 += 0.5*rb_jnode->S0;
           rb_inode->S1 += 0.5*rb_jnode->S1;
           rb_inode->S2 += 0.5*rb_jnode->S2;
           rb_inode->S3 += 0.5*rb_jnode->S3;
           for (idim=0; idim<Ndim; idim++) {
              rb_inode->V1[idim] += 0.5*rb_jnode->V1[idim];
              rb_inode->V2[idim] += 0.5*rb_jnode->V2[idim];
           }
        }
        else {
          if(jnode_box == -1){   /* bulk boundary */
            rb0 = rb1 = rb2 = rb3 = 0.0;
            for (icomp=0; icomp<Ncomp; icomp++) {
               rb0 += Rho_b[icomp];
               rb1 += Rho_b[icomp]* Sigma_ff[icomp][icomp]/2.0;
               rb2 += Rho_b[icomp]* PI* Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
               rb3 += Rho_b[icomp]* (PI / 6.0) * Sigma_ff[icomp][icomp]
                       * Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
            }
          }
          else if(jnode_box == -3){    /* LBB boundary */
            rb0 = rb1 = rb2 = rb3 = 0.0;
            for (icomp=0; icomp<Ncomp; icomp++) {
               rb0 += Rho_b_LBB[icomp];
               rb1 += Rho_b_LBB[icomp]* Sigma_ff[icomp][icomp]/2.0;
               rb2 += Rho_b_LBB[icomp]* PI* Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
               rb3 += Rho_b_LBB[icomp]* (PI / 6.0) * Sigma_ff[icomp][icomp]
                       * Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
            }
          }
          else if(jnode_box == -4){  /* RTF boundary */
            rb0 = rb1 = rb2 = rb3 = 0.0;
            for (icomp=0; icomp<Ncomp; icomp++) {
               rb0 += Rho_b_RTF[icomp];
               rb1 += Rho_b_RTF[icomp]* Sigma_ff[icomp][icomp]/2.0;
               rb2 += Rho_b_RTF[icomp]* PI* Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
               rb3 += Rho_b_RTF[icomp]* (PI / 6.0) * Sigma_ff[icomp][icomp]
                       * Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
            }
          }
          rb_inode->S0 += 0.5*rb0;
          rb_inode->S1 += 0.5*rb1;
          rb_inode->S2 += 0.5*rb2;
          rb_inode->S3 += 0.5*rb3;
        }


        /* Go to node 1 lower in the appropriate ijk direction */

        offset_ptr += 9;  /* gets negative of offset used above */
        jnode_box = offset_to_node_box(ijk_box, offset_ptr, reflect_flag);

        if (jnode_box >= 0) {
           rb_jnode = &(rho_bar[B2L_1stencil[jnode_box]]);
           rb_inode->S0 += 0.5*rb_jnode->S0;
           rb_inode->S1 += 0.5*rb_jnode->S1;
           rb_inode->S2 += 0.5*rb_jnode->S2;
           rb_inode->S3 += 0.5*rb_jnode->S3;
           for (idim=0; idim<Ndim; idim++) {
              rb_inode->V1[idim] += 0.5*rb_jnode->V1[idim];
              rb_inode->V2[idim] += 0.5*rb_jnode->V2[idim];
           }
        }
        else {
          if(jnode_box == -1){   /* bulk boundary */
            rb0 = rb1 = rb2 = rb3 = 0.0;
            for (icomp=0; icomp<Ncomp; icomp++) {
               rb0 += Rho_b[icomp];
               rb1 += Rho_b[icomp]* Sigma_ff[icomp][icomp]/2.0;
               rb2 += Rho_b[icomp]* PI* Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
               rb3 += Rho_b[icomp]* (PI / 6.0) * Sigma_ff[icomp][icomp]
                       * Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
            }
          }
          else if(jnode_box == -3){    /* LBB boundary */
            rb0 = rb1 = rb2 = rb3 = 0.0;
            for (icomp=0; icomp<Ncomp; icomp++) {
               rb0 += Rho_b_LBB[icomp];
               rb1 += Rho_b_LBB[icomp]* Sigma_ff[icomp][icomp]/2.0;
               rb2 += Rho_b_LBB[icomp]* PI* Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
               rb3 += Rho_b_LBB[icomp]* (PI / 6.0) * Sigma_ff[icomp][icomp]
                       * Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
            }
          }
          else if(jnode_box == -4){  /* RTF boundary */
            rb0 = rb1 = rb2 = rb3 = 0.0;
            for (icomp=0; icomp<Ncomp; icomp++) {
               rb0 += Rho_b_RTF[icomp];
               rb1 += Rho_b_RTF[icomp]* Sigma_ff[icomp][icomp]/2.0;
               rb2 += Rho_b_RTF[icomp]* PI* Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
               rb3 += Rho_b_RTF[icomp]* (PI / 6.0) * Sigma_ff[icomp][icomp]
                       * Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
            }
          }
          rb_inode->S0 += 0.5*rb0;
          rb_inode->S1 += 0.5*rb1;
          rb_inode->S2 += 0.5*rb2;
          rb_inode->S3 += 0.5*rb3;
        }

    }
  }
}
/****************************************************************************/
double  HW_boundary_weight(int icomp,int ilist, double *hw_weight,
                           int inode_box, int *reflect_flag)
/* This routine assembles the weight of boundary nodes in the Hard Wall model*/

{
  int local_node,iel_box;
  double weight=0.0;

  /*
   * Wall_elems is equal to the wall number for wall elements & -1 for fluid elements
   * node_to_elem gives the element that has node inode as
   *   its local node local_node
   */

  for (local_node=0; local_node<Nnodes_per_el_V; local_node++){

    iel_box = node_box_to_elem_box_reflect(inode_box, local_node, reflect_flag);
    if (iel_box >=0)
      if (Wall_elems[ilist][iel_box] == -1 || Lsemiperm[WallType[Wall_elems[ilist][iel_box]]][icomp])
         weight += hw_weight[local_node];
  }

  return(weight);
}
/****************************************************************************/
/* find_jzone:  this little subroutine sets jzone for any izone given the options in the
code.  Previously this little piece of code appeared in many places making it rather
bug prone. */
int find_jzone(int izone)
{
  int jzone;
  if (Coarser_jac > 0 && izone<Nzone-1){
      if (Coarser_jac == 1){
          if (izone ==0 ) jzone = izone + 1;
          else                 jzone = izone;
      }
      else if (Coarser_jac == 2) jzone = izone + 1;
      else if (Coarser_jac == 3) jzone = Nzone-1;
      else if (Coarser_jac == 4) jzone = Nzone-2;
      else if (Coarser_jac == 5) jzone = Nzone-1;
  }
  else                         jzone = izone;
  return jzone;
}
/****************************************************************************/
/*constant_boundary:  This routine just returns a boundary condition 
given an unknown number and an indication of what kind of boundary we have.*/
double constant_boundary(int iunk,int jnode_box)
{
    double bcval;
    switch(Unk_to_eq_type[iunk]){
       case DENSITY:
          if (Type_poly_TC){
             printf("problem ... we haven't done the coarsening for TC polymers yet\n");
             exit(-1);
             if (jnode_box==-1){ }
             if (jnode_box==-2){ }
             else if (jnode_box==-3){ }
             else if (jnode_box==-4){ }
           }
           else{
              if (jnode_box==-1)     bcval = Rho_b[iunk-Unk_start_eq[DENSITY]];
              if (jnode_box==-2)     bcval = 0.0;
              else if (jnode_box==-3){
                  if (Lsteady_state) bcval = Rho_b_LBB[iunk-Unk_start_eq[DENSITY]];
                  else               bcval = Rho_coex[1];
              }
              else if (jnode_box==-4){
                  if (Lsteady_state) bcval = Rho_b_RTF[iunk-Unk_start_eq[DENSITY]];
                  else               bcval = Rho_coex[2];
              }
           }
           break;
       case RHOBAR_ROSEN:
           if (jnode_box==-1)       bcval = Rhobar_b[iunk-Unk_start_eq[RHOBAR_ROSEN]];
           if (jnode_box==-2)       bcval = 0.0;  /* assuming wall begins in domain and rhobars 
                                                     have decayed beyond the boundary */
           else if (jnode_box==-3)  bcval = Rhobar_b_LBB[iunk-Unk_start_eq[RHOBAR_ROSEN]];
           else if (jnode_box==-4)  bcval = Rhobar_b_RTF[iunk-Unk_start_eq[RHOBAR_ROSEN]];
           break;
       case POISSON:
           if (jnode_box==-1)       bcval = 0.;
           if (jnode_box==-2){      printf("can't define Electric Potential in a wall\n"); 
                                    exit(-1);
                              }
           else if (jnode_box==-3)  bcval = Elec_pot_LBB;
           else if (jnode_box==-4)  bcval = Elec_pot_RTF;
           break;
       case DIFFUSION:
           if (jnode_box==-1)       bcval = 0.;
           if (jnode_box==-2)       bcval = -VEXT_MAX;
           else if (jnode_box==-3)  bcval = Betamu_LBB[iunk - Unk_start_eq[DIFFUSION]];
           else if (jnode_box==-4)  bcval = Betamu_RTF[iunk - Unk_start_eq[DIFFUSION]];
           break;
       case DENSITY_SEG:
           if (jnode_box==-1){bcval=0.; }
           else if (jnode_box==-3){bcval=0.;}
           else if (jnode_box==-4){bcval=0.;}
           break;
       case CAVITY_WTC:
           if (jnode_box==-1){bcval=0.; }
           else if (jnode_box==-3){bcval=0.;}
           else if (jnode_box==-4){bcval=0.;}
           break;
       case BOND_WTC:
           if (jnode_box==-1){bcval=0.; }
           else if (jnode_box==-3){bcval=0.;}
           else if (jnode_box==-4){bcval=0.;}
           break;
       case CMS_FIELD:
           bcval=0.0; 
           break;
       case CMS_G:
           bcval=1.0; 
           break;
   }
   return(bcval);
}
/****************************************************************************/



