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
 *  FILE: dft_fill_mf.c
 *
 *  This file contains the fill of the residual equations and Jacobian
 *  matrix for mean field attractions and electrostatics.
 */

#include "dft_globals_const.h"
#include "rf_allo.h"

/*****************************************************************************/
/* load_mean_field: In this routine we load the residual and Jacobian for all
                    terms that have the mean field form.  The attractions and
                    the electrostatic correlation terms.                     */
void load_mean_field(int sten_type, int loc_i, int icomp, int izone, 
                     int *ijk_box, int fill_flag, double *resid,
                     double *mat_row, int *bindx_tmp, double *x, 
                     int resid_only_flag)
{
  int   **sten_offset, *offset, isten;
  double *sten_weight,  weight;
  struct Stencil_Struct *sten;

  int   **sten_offsetJ, *offsetJ;
  double *sten_weightJ,weightJ;
  struct Stencil_Struct *stenJ;

  double sign;
  int jcomp, jlist;
  int reflect_flag[NDIM_MAX];
  int j_box, jnode_box;
  int jzone=0;
  int jnode_boxJ;

  sign = 1.0;
  if (sten_type == THETA_CHARGE) sign = -1.0;

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

  for (jcomp=0; jcomp<Ncomp; jcomp++){
      if (Nlists_HW <= 2) jlist = 0;
      else                jlist = jcomp;
           
      sten = &(Stencil[sten_type][izone][icomp+Ncomp*jcomp]);
      sten_offset = sten->Offset;
      sten_weight = sten->Weight;

      stenJ = &(Stencil[sten_type][jzone][icomp+Ncomp*jcomp]);
      sten_offsetJ = stenJ->Offset;
      sten_weightJ = stenJ->Weight;

      for (isten = 0; isten < sten->Length; isten++) {
        offset = sten_offset[isten];
        weight = sten_weight[isten];

         jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);

         if (fill_flag != MSR_PREPROCESS){

            if (jnode_box >= 0 && !Zero_density_TF[jnode_box][jcomp]) {
               if (Lhard_surf) {
                   if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1)
                      weight = HW_boundary_weight
                       (jcomp,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
               }
               j_box=jnode_box*Nunk_per_node+jcomp;
               resid[loc_i] +=  sign*weight*x[B2L_unknowns[j_box]];
            }
            else if (jnode_box == -1 || jnode_box ==-3 || jnode_box == -4){
               if      (jnode_box == -1) resid[loc_i] +=  sign*weight*Rho_b[jcomp];
               else if (jnode_box == -3) {
                  if (Lsteady_state == TRUE)
                       resid[loc_i] +=  sign*weight*Rho_b_LBB[jcomp];
                  else resid[loc_i] +=  sign*weight*Rho_coex[1];
               }
               else if (jnode_box == -4) {
                  if (Lsteady_state == TRUE)
                       resid[loc_i] +=  sign*weight*Rho_b_RTF[jcomp];
                  else
                    resid[loc_i] +=  sign*weight*Rho_coex[0];
               }
            }
         }   

         if (!resid_only_flag) {
           if (isten < stenJ->Length){
           offsetJ = sten_offsetJ[isten];
           weightJ = sten_weightJ[isten];
           jnode_boxJ = offset_to_node_box(ijk_box, offsetJ, reflect_flag);
              if (jnode_boxJ >=0 && !Zero_density_TF[jnode_boxJ][jcomp]){
                 if (fill_flag != MSR_PREPROCESS){
                    if (Lhard_surf) {
                       if (Nodes_2_boundary_wall[jlist][jnode_boxJ]!=-1)
                          weightJ = HW_boundary_weight (jcomp, jlist,
                            stenJ->HW_Weight[isten], jnode_boxJ, reflect_flag);
                    }
                    j_box=jnode_boxJ*Nunk_per_node+jcomp;
                    mat_row[B2L_unknowns[j_box]] += sign*weight;
                 }
                 else {
                     j_box=jnode_boxJ*Nunk_per_node+jcomp;
                     bindx_tmp[j_box] = TRUE;
                 }
              }
           }
        }
      }
  }
  return;
}
/*****************************************************************************/
/* load_mean_field_old: In this routine we load the residual and Jacobian for all
                    terms that have the mean field form.  The attractions and
                    the electrostatic correlation terms.                     
void load_mean_field_old(int sten_type, int loc_i, int icomp, int izone, 
                int *ijk_box, int fill_flag,
                double *resid, double *mat_row, int *bindx_tmp, double *x)
{
  int   **sten_offset, *offset, isten;
  double *sten_weight,  weight;
  struct Stencil_Struct *sten;

  double sign;
  int jcomp, jlist;
  int reflect_flag[NDIM_MAX];
  int j_box, jnode_box;

  sign = 1.0;
  if (sten_type == THETA_CHARGE) sign = -1.0;

  for (jcomp=0; jcomp<Ncomp; jcomp++){
      if (Nlists_HW <= 2) jlist = 0;
      else                jlist = jcomp;
           
      sten = &(Stencil[sten_type][izone][icomp+Ncomp*jcomp]);
      sten_offset = sten->Offset;
      sten_weight = sten->Weight;

      for (isten = 0; isten < sten->Length; isten++) {
        offset = sten_offset[isten];
        weight = sten_weight[isten];

         jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);

         if (jnode_box >= 0 ) {
            if (Lhard_surf && fill_flag != MSR_PREPROCESS) {
                if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1) 
                   weight = HW_boundary_weight 
                    (jcomp,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
            }
             j_box=jnode_box*Nunk_per_node+jcomp;
             if (fill_flag != MSR_PREPROCESS){
                 resid[loc_i] +=  sign*weight*x[B2L_unknowns[j_box]];
                 mat_row[B2L_unknowns[j_box]] += sign*weight;
             }
             else 
                 bindx_tmp[j_box] = TRUE;
         }
         else if ( (jnode_box == -1 || jnode_box ==-3 || jnode_box == -4) 
                              && fill_flag != MSR_PREPROCESS){
              if      (jnode_box == -1) resid[loc_i] +=  sign*weight*Rho_b[jcomp];
              else if (jnode_box == -3) {
                 if (Lsteady_state == TRUE)
                      resid[loc_i] +=  sign*weight*Rho_b_LBB[jcomp];
                 else resid[loc_i] +=  sign*weight*Rho_coex[1];
              }
              else if (jnode_box == -4) {
                 if (Lsteady_state == TRUE)
                      resid[loc_i] +=  sign*weight*Rho_b_RTF[jcomp];
                 else
                   resid[loc_i] +=  sign*weight*Rho_coex[0];
              }
        }
      }
  }
  return;
}
*/
/****************************************************************************/
