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

/* dft_scale_variables.c:  These routines are used for scaling variables during
   the course of arc length continuation calculations */

#include "dft_scale_variables.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* recalculate_stencils:  free the old Stencil structure and recompute stencil
   with the new parameters */
void recalculate_stencils()
{
   int izone, isten,jmax,icomp,jcomp;
   struct Stencil_Struct *sten;

   for (izone=0; izone<Nzone; izone++){
     for (isten=0; isten<NSTEN; isten++)

       if (Sten_Type[isten]) {
          jmax = stencil_Njcomp_switch(isten);

          for (icomp=0; icomp<Ncomp; icomp++) {
            for (jcomp=0; jcomp<jmax; jcomp++) {
               sten = &(Stencil[isten][izone][icomp+Ncomp*jcomp]);
               safe_free((void **) &sten->Offset);
               safe_free((void **) &sten->Weight);
               if (Lhard_surf) safe_free((void **) &sten->HW_Weight);
            }
          }
       }
  }
  safe_free((void **) &Stencil);
  calc_stencils();
  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* scale_external_field: when changing eps wall only a
   simple scaling is required. */
void scale_vext_temp(double ratio)
{
   int loc_inode,icomp,iwall,idim,iunk;

   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
      for (icomp=0; icomp<Ncomp; icomp++){
         if (Restart_Vext != READ_VEXT_STATIC){
            Vext[loc_inode][icomp] /= ratio;
         }
         else{
            Vext[loc_inode][icomp] -= Vext_static[loc_inode][icomp];
            Vext[loc_inode][icomp] /= ratio;
            Vext[loc_inode][icomp] += Vext_static[loc_inode][icomp];
            if (Vext_static[loc_inode][icomp] >= VEXT_MAX) Vext[loc_inode][icomp]=VEXT_MAX;
         }
         if (Lvext_dash && Restart_Vext == READ_VEXT_FALSE){
         for (iwall=0; iwall<Nwall; iwall++)
         for (idim=0; idim<Ndim; idim++){
            iunk = iwall*Ncomp + icomp;
            Vext_dash[loc_inode][iunk][idim] /= ratio;
         }
         }
      }
   }
   if (Iwrite==VERBOSE) {
      print_vext(Vext,"dft_vext_cont2.dat");
      if (Restart_Vext==READ_VEXT_STATIC) print_vext(Vext_static,"dft_vext_static.dat");
   }
   return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* scale_external_field: when changing eps wall-fluid for the
   first component, only a simple scaling is required. */

void scale_vext_epswf(double ratio, int icomp,int iwall)
{
   int loc_inode,idim,iunk;
   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
       if (Restart_Vext != READ_VEXT_STATIC){
            Vext[loc_inode][icomp] *= ratio;
       }
       else{
           Vext[loc_inode][icomp] -= Vext_static[loc_inode][icomp];
           Vext[loc_inode][icomp] *= ratio;
           Vext[loc_inode][icomp] += Vext_static[loc_inode][icomp];
       }
       if (Lvext_dash && Restart_Vext == READ_VEXT_FALSE){
       iunk = iwall*Ncomp+icomp;
       for (idim=0; idim<Ndim; idim++) Vext_dash[loc_inode][iunk][idim] *= ratio;
       }
   }
   if (Iwrite==VERBOSE) {
      print_vext(Vext,"dft_vext_cont1.dat");
      if (Restart_Vext==READ_VEXT_STATIC) print_vext(Vext_static,"dft_vext_static.dat");
   }
   return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* scale_elec_param: when changing eps wall only a simple scaling is required. */
void scale_elec_param(double ratio)
{
   int iwall,iel,loc_inode,idim;

   if (Surf_charge_flag){
      for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++)
           for (idim=0; idim<Ndim; idim++) Charge_w_sum_els[loc_inode][idim] *=ratio;
   }
   if (Vol_charge_flag) for (iel=0; iel<Nelements_box; iel++) Charge_vol_els[iel]*=ratio;
   for (iwall=0; iwall<Nwall; iwall++) Elec_param_w[iwall]*=ratio;

   if (Iwrite == VERBOSE){
       if (Proc==0) printf ("PRINTING CHARGE DISTRIBUTIONS: Scale_fac=%9.6f\n",ratio);
       if (Vol_charge_flag) print_charge_vol(Charge_vol_els,"dft_charge_vol.dat");
       if (Surf_charge_flag) print_charge_surf(Charge_w_sum_els,"dft_charge_surf.dat");
    }

   return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* set_new_membrane_potential: change the set point for the external field */
void set_new_membrane_potential(double param_old,double param_new,int icomp)
{
   int inode;

   for (inode=0;inode<Nnodes_per_proc;inode++){
        if (fabs(Vext[inode][icomp]-param_old)<1.e-10) Vext[inode][icomp]=param_new;
   }
   return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* scale_all_epsParams: change the set point for the external field */
void scale_all_epsParams(double ratio)
{
  int icomp,jcomp,i,j;

  if (Mix_type==0){
     if (Type_attr != NONE){
        for (icomp=0; icomp<Ncomp; icomp++) Eps_ff[icomp][icomp] /= ratio;
     }
     for (i=0; i<Nwall_type;i++){
        if(Ipot_wf_n[i] != VEXT_HARD) Eps_w[i] /= ratio;
     }
  }
  else if (Mix_type==1){
     for (icomp=0; icomp<Ncomp; icomp++){
        for(jcomp=0; jcomp<Ncomp; jcomp++) Eps_ff[icomp][jcomp] /= ratio;
        for (i=0; i<Nwall_type;i++) Eps_wf[icomp][i] /= ratio;
     }
     for (i=0; i<Nwall_type;i++) {
        for (j=0;j<Nwall_type;j++) Eps_ww[i][j] /= ratio;
     }
  }
  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/


