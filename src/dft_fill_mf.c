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
  int jcomp, jlist,junk;
  int reflect_flag[NDIM_MAX];
  int j_box, jnode_box;
  int jzone=0;
  int jnode_boxJ;

  sign = 1.0;
  if (sten_type == THETA_CHARGE) sign = -1.0;

  jzone = find_jzone(izone);

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
               j_box=loc_find(jcomp+Phys2Unk_first[DENSITY],jnode_box,BOX);
               resid[loc_i] +=  sign*weight*x[B2L_unknowns[j_box]];
            }
            else if (jnode_box == -1 || jnode_box ==-3 || jnode_box == -4){
               junk=jcomp+Phys2Unk_first[DENSITY];
               resid[loc_i] += sign*weight*constant_boundary(junk,jnode_box);
            }
         }   

         if (!resid_only_flag) {
           if (isten < stenJ->Length){
           offsetJ = sten_offsetJ[isten];
           weightJ = sten_weightJ[isten];
           jnode_boxJ = offset_to_node_box(ijk_box, offsetJ, reflect_flag);
           j_box = loc_find(Phys2Unk_first[DENSITY]+jcomp,jnode_boxJ,BOX);
              if (jnode_boxJ >=0 && !Zero_density_TF[jnode_boxJ][jcomp]){
                 if (fill_flag != MSR_PREPROCESS){
                    if (Lhard_surf) {
                       if (Nodes_2_boundary_wall[jlist][jnode_boxJ]!=-1)
                          weightJ = HW_boundary_weight (jcomp, jlist,
                            stenJ->HW_Weight[isten], jnode_boxJ, reflect_flag);
                    }
                    mat_row[B2L_unknowns[j_box]] += sign*weight;
                 }
                 else  bindx_tmp[j_box] = TRUE;
              }
           }
        }
      }
  }
  return;
}
/*****************************************************************************/
