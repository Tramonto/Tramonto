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

/* ---------------------------------------------------------
dft_fill_CLSmf.c:

Here are the generic fill routines that use function pointers
for mean field terms in the density functional theory.
------------------------------------------------------------*/

#include "dft_fill_CLSmf.h"

/****************************************************************************/
double resid_and_Jac_sten_fill_sum_Ncomp (int sten_type, double **x, int iunk, 
  int loc_inode, int inode_box, int izone,
  int *ijk_box, int resid_only_flag, int jzone_flag,
  double (*fp_prefactor)(int,int,int *), double (*fp_resid)(int,int,double **), 
  double (*fp_jacobian)(int,int,double **))
{
  int   **sten_offset, *offset, isten,idim;
  int   **sten_offsetJ, *offsetJ;
  double *sten_weightJ,weightJ;
  double *sten_weight,  weight,fac;
  struct Stencil_Struct *sten;
  struct Stencil_Struct *stenJ;
  double resid=0.0,mat_val,resid_sum=0.0,weight_bulk,bulk_term;

  int icomp,jzone, jnode_box, jcomp,jlist,junk,loop_max,jloop,jseg,index;
  int jnode_boxJ,i;
  int reflect_flag[NDIM_MAX];

  for (idim=0;idim<Ndim;idim++) reflect_flag[idim]=FALSE;
  if (jzone_flag) jzone = find_jzone(izone,inode_box);
  else{
      if (izone >= 0) jzone = izone;
      else{
         if (izone=FLAG_BULK || izone==FLAG_PBELEC) jzone=0;
         else {
             printf("I'm confused about the zones see dft_utils.c\n");
             exit(-1);
         }
      }
  }    

  if (Lseg_densities) loop_max=Nseg_tot;
  else                loop_max=Ncomp;


  for (jloop=0; jloop<loop_max; jloop++){
      if (Lseg_densities) {
         jseg=jloop;
         jcomp=Unk2Comp[jseg];
         junk=jseg;
         if (Pol_Sym_Seg[jseg] != -1) junk=Pol_Sym_Seg[jseg];
         junk+=Phys2Unk_first[DENSITY];
      }
      else{
         jcomp=jloop;
         junk=Phys2Unk_first[DENSITY]+jcomp;
      }

      if (Nlists_HW <= 2) jlist = 0;
      else                jlist = jcomp;
      if (sten_type==THETA_PAIRPOT_RCUT) {

            if (Unk2Phys[iunk]==MF_EQ){
                  icomp=iunk-Phys2Unk_first[MF_EQ];
            }
            else{
               if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3) i=iunk-Phys2Unk_first[WJDC_FIELD];
               else                 i=iunk-Phys2Unk_first[DENSITY];

               if ((Type_poly==WTC || Type_poly==WJDC) && Unk2Phys[iunk] != MF_EQ)  icomp=Unk2Comp[i];
               else icomp=i;
            }
            index=icomp+Ncomp*jcomp;
      }
      else if (sten_type==THETA_CR_RPM_MSA || sten_type==THETA_CR_GENERAL_MSA){
            i=iunk-Phys2Unk_first[DENSITY];
            if (Type_poly==WTC || Type_poly==WJDC || Type_poly==WJDC2) icomp=Unk2Comp[i]; 
            else icomp=i;
            index=icomp+Ncomp*jcomp;
      }
      else if (sten_type==THETA_CR_DATA){
            if (Unk2Phys[iunk]==MF_EQ){
                  icomp=iunk-Phys2Unk_first[MF_EQ];
            }
            else  icomp=iunk-Phys2Unk_first[CMS_FIELD];
            index=icomp+Ncomp*jcomp;
      }
      else index=jcomp;


      sten = &(Stencil[sten_type][izone][index]);
      sten_offset = sten->Offset;
      sten_weight = sten->Weight;

      stenJ = &(Stencil[sten_type][jzone][index]);
      sten_offsetJ = stenJ->Offset;
      sten_weightJ = stenJ->Weight;

      for (isten = 0; isten < sten->Length; isten++) {
        offset = sten_offset[isten];
        weight = sten_weight[isten];
        weight_bulk=weight;
      
        if (fp_prefactor!=NULL) fac = (*fp_prefactor)(iunk,jcomp,offset);
        else fac=1.0;
        if (sten_type==THETA_CR_RPM_MSA || sten_type==THETA_CR_GENERAL_MSA || sten_type==THETA_CR_DATA) fac=-1.0;
        if (sten_type==THETA_CR_DATA) bulk_term= weight_bulk*Rho_b[jcomp];
        else bulk_term=0.0;

         jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);
         if (Lhard_surf && jnode_box >=0) {
             if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1)
                weight = HW_boundary_weight
                 (jcomp,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
         }
         resid =  fac*(weight*(*fp_resid)(junk,jnode_box,x) - bulk_term);
         resid_sum+=resid;
         if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY){
               dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
         }

         if (resid_only_flag==FALSE)
         if (isten < stenJ->Length){
            offsetJ = sten_offsetJ[isten];
            weightJ = sten_weightJ[isten];
            jnode_boxJ = offset_to_node_box(ijk_box, offsetJ, reflect_flag);
            if (jnode_boxJ >=0 && !Zero_density_TF[jnode_boxJ][jcomp]){
                  if (Lhard_surf) {
                     if (Nodes_2_boundary_wall[jlist][jnode_boxJ]!=-1)
                        weightJ = HW_boundary_weight
                         (jcomp,jlist,stenJ->HW_Weight[isten], jnode_boxJ, reflect_flag);
                   }
                   mat_val = weightJ*fac*(*fp_jacobian)(junk,jnode_box,x);
                   dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,jnode_boxJ,mat_val);
            }
         }
      }

  }
  return(resid_sum);
}
/****************************************************************************/
double resid_and_Jac_sten_fill (int sten_type, double **x, int iunk, int junk,
  int icomp, int jcomp, int loc_inode, int inode_box, int izone,
  int *ijk_box, int resid_only_flag, int jzone_flag,
  double(*fp_prefactor)(int,int,int *), double(*fp_resid)(int,int,double **), 
  double (*fp_jacobian)(int,int,double **))
{
  int   **sten_offset, *offset, isten,idim;
  int   **sten_offsetJ, *offsetJ;
  double *sten_weightJ,weightJ;
  double *sten_weight,  weight,fac;
  struct Stencil_Struct *sten;
  struct Stencil_Struct *stenJ;
  double resid,mat_val,resid_sum=0.0;

  int jzone, jnode_box, jlist,loop_max,jloop,jseg;
  int jnode_boxJ;
  int reflect_flag[NDIM_MAX];

  for (idim=0;idim<Ndim;idim++) reflect_flag[idim]=FALSE;

  if (jzone_flag) jzone = find_jzone(izone,inode_box);
  else            jzone = izone;


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

    if (fp_prefactor!=NULL) fac = (*fp_prefactor)(iunk,jcomp,offset);
    else fac=1.0;

     jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);
     if (Lhard_surf && jnode_box >=0) {
         if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1)
            weight = HW_boundary_weight
             (jcomp,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
     }
     resid =  weight*fac*(*fp_resid)(junk,jnode_box,x);
     if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
     resid_sum+=resid;

     if (resid_only_flag==FALSE)
     if (isten < stenJ->Length){
        offsetJ = sten_offsetJ[isten];
        weightJ = sten_weightJ[isten];
        jnode_boxJ = offset_to_node_box(ijk_box, offsetJ, reflect_flag);
        if (jnode_boxJ >=0 && !Zero_density_TF[jnode_boxJ][jcomp]){
              if (Lhard_surf) {
                 if (Nodes_2_boundary_wall[jlist][jnode_boxJ]!=-1)
                    weightJ = HW_boundary_weight
                     (jcomp,jlist,stenJ->HW_Weight[isten], jnode_boxJ, reflect_flag);
               }
               mat_val = weightJ*fac*(*fp_jacobian)(junk,jnode_box,x);
               dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,jnode_boxJ,mat_val);
        }
     }
  }
  return(resid_sum);
}
/****************************************************************************/
