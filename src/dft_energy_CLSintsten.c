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
dft_energy_CLSintsten.c:

Here are some routines that involve integration over
the integration stencils in the DFT code - these are used in 
energy calculations.
------------------------------------------------------------*/

#include "dft_energy_CLSintsten.h"

/***********************************************************************
int_stencil_bulk: this routine sums the appropriate stencil to get
                  the bulk contributions to various terms in the E-L
                  equation. Note that some terms (attractions, WCA electrostatics,
                  CMS polymers bury a multiplier function with the stencil weight function.*/
double int_stencil_bulk(int sten_type,int icomp,int jcomp,double(*fp_integrand)(double,int,int))
{
  int izone, isten,*offset,**sten_offset,idim;
  double sum, weight, *sten_weight,integrand,rsq;
  struct Stencil_Struct *sten;

  sum = 0.0;
  izone = 0;
  sten = &(Stencil[sten_type][izone][icomp+Ncomp*jcomp]);
  sten_weight = sten->Weight;
  sten_offset = sten->Offset;

  for (isten = 0; isten<sten->Length; isten++){
     
     if (fp_integrand!=NULL){
          offset = sten_offset[isten];
          rsq=0.;
          for (idim=0;idim<Ndim;idim++) rsq+=offset[idim]*offset[idim]*Esize_x[idim]*Esize_x[idim];
          if (rsq>1.e-8) integrand = (*fp_integrand)(rsq,icomp,jcomp);
          else integrand=1.0;
     }
     else integrand=1.0;
     weight = integrand*sten_weight[isten];
     sum += weight;
  }
  return(sum);
}
/*******************************************************************************/
/*int_stencil: Perform the integral sum(j)int rho_j(r')*weight[sten] for
  attraction and delta c mean field integrations.*/
 double int_stencil(double **x,int inode_box,int iunk,int sten_type)
{
  int isten,*offset,inode_sten,ijk_box[3],izone,idim;
  int j,jcomp,junk,icomp,jlist;
  double weight, sum;
  struct Stencil_Struct *current_sten;
  int **current_sten_offset, reflect_flag[NDIM_MAX];
  double *current_sten_weight;
  for (idim=0; idim<Ndim; idim++) reflect_flag[idim]=FALSE;

/*  izone = Nodes_to_zone_box[inode_box];*/
  izone = 0;

  sum = 0.0;
  node_box_to_ijk_box(inode_box,ijk_box);
  if (Lseg_densities) icomp=Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
  else                icomp=iunk-Phys2Unk_first[DENSITY];

  for (junk=Phys2Unk_first[DENSITY];junk<Phys2Unk_last[DENSITY];junk++){
     if (Lseg_densities) jcomp=Unk2Comp[junk-Phys2Unk_first[DENSITY]];
     else                jcomp=junk-Phys2Unk_first[DENSITY];

     if (Nlists_HW <= 2) jlist = 0;
     else                jlist = jcomp;

     current_sten = &(Stencil[sten_type][izone][icomp+Ncomp*jcomp]);
     current_sten_offset = current_sten->Offset;
     current_sten_weight = current_sten->Weight;

     for (isten = 0; isten < current_sten->Length; isten++) {
        offset = current_sten_offset[isten];
        weight = current_sten_weight[isten];

         /* Find in the Stencil position on overall mesh */
        inode_sten =offset_to_node_box(ijk_box, offset, reflect_flag);

        if (inode_sten >= 0 && !Zero_density_TF[inode_sten][jcomp]) {
           if (Lhard_surf) {
               if (Nodes_2_boundary_wall[jlist][inode_sten]!=-1)
                  weight = HW_boundary_weight
                    (jcomp,jlist,current_sten->HW_Weight[isten], inode_sten, reflect_flag);
           }
           if (inode_sten<Nnodes_box && inode_sten >=0){
               sum +=  weight*x[junk][inode_sten];
           }
        }
        else if (inode_sten<0){
             sum += weight*constant_boundary(junk,inode_sten);
        }

     }  /* end of loop over isten */
  }     /* end of loop over jcomp */
  return(sum);
}
/*******************************************************************************/
/*int_stencil_CMSField: Perform the integral sum(j)int (rho_j(r')-rho_b)*weight[sten] for
  attraction and delta c mean field integrations.*/
 double int_stencil_CMSField(double **x,int inode_box,int iunk,int sten_type)
{
  int isten,*offset,inode_sten,ijk_box[3],izone,idim;
  int j,jcomp,junk,icomp,jlist;
  double weight, weight_bulk,sum;
  struct Stencil_Struct *current_sten;
  int **current_sten_offset, reflect_flag[NDIM_MAX];
  double *current_sten_weight;
  for (idim=0; idim<Ndim; idim++) reflect_flag[idim]=FALSE;

/*  izone = Nodes_to_zone_box[inode_box];*/
  izone = 0;

  sum = 0.0;
  node_box_to_ijk_box(inode_box,ijk_box);
  if (Lseg_densities) icomp=Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
  else                icomp=iunk-Phys2Unk_first[DENSITY];

  for (junk=Phys2Unk_first[DENSITY];junk<Phys2Unk_last[DENSITY];junk++){
     if (Lseg_densities) jcomp=Unk2Comp[junk-Phys2Unk_first[DENSITY]];
     else                jcomp=junk-Phys2Unk_first[DENSITY];

     if (Nlists_HW <= 2) jlist = 0;
     else                jlist = jcomp;

     current_sten = &(Stencil[sten_type][izone][icomp+Ncomp*jcomp]);
     current_sten_offset = current_sten->Offset;
     current_sten_weight = current_sten->Weight;

     for (isten = 0; isten < current_sten->Length; isten++) {
        offset = current_sten_offset[isten];
        weight = current_sten_weight[isten];
   
        weight_bulk=weight; 

         /* Find in the Stencil position on overall mesh */
        inode_sten =offset_to_node_box(ijk_box, offset, reflect_flag);

        if (inode_sten >= 0 && !Zero_density_TF[inode_sten][jcomp]) {
           if (Lhard_surf) {
               if (Nodes_2_boundary_wall[jlist][inode_sten]!=-1)
                  weight = HW_boundary_weight
                    (jcomp,jlist,current_sten->HW_Weight[isten], inode_sten, reflect_flag);
           }
           if (inode_sten<Nnodes_box && inode_sten >=0){
               sum +=  (weight*x[junk][inode_sten]-weight_bulk*Rho_b[jcomp]);
           }
        }
        else if (inode_sten<0){ 
             sum += (weight*constant_boundary(junk,inode_sten)-weight_bulk*Rho_b[jcomp]);
        }
        else{
             sum += (-weight_bulk*Rho_b[jcomp]);
        }

     }  /* end of loop over isten */
  }     /* end of loop over jcomp */
  return(sum);
}
/*******************************************************************************/
/*int_stencil_HSFMT: Perform the integral sum(j)int rho_j(r')*prefac*weight[sten] 
for the FMT hard sphere functionals */
 double int_stencil_HSFMT(double **x,int inode_box,int iunk)
{
  int isten,*offset,inode_sten,ijk_box[3],izone,idim;
  int j,jcomp,junk,icomp,jlist,sten_type;
  double weight, sum,fac;
  struct Stencil_Struct *current_sten;
  int **current_sten_offset, reflect_flag[NDIM_MAX];
  double *current_sten_weight;
  for (idim=0; idim<Ndim; idim++) reflect_flag[idim]=FALSE;

/*  izone = Nodes_to_zone_box[inode_box];*/
  izone = 0;

  sum = 0.0;
  node_box_to_ijk_box(inode_box,ijk_box);

  if (iunk==Phys2Unk_first[HSRHOBAR]) sten_type=THETA_FN_R;
  else                                sten_type=DELTA_FN_R;

  for (junk=Phys2Unk_first[DENSITY];junk<Phys2Unk_last[DENSITY];junk++){

     if (Lseg_densities) jcomp=Unk2Comp[junk-Phys2Unk_first[DENSITY]];
     else                jcomp=junk-Phys2Unk_first[DENSITY];



     if (Nlists_HW <= 2) jlist = 0;
     else                jlist = jcomp;

     current_sten = &(Stencil[sten_type][izone][jcomp]);
     current_sten_offset = current_sten->Offset;
     current_sten_weight = current_sten->Weight;

     for (isten = 0; isten < current_sten->Length; isten++) {
        offset = current_sten_offset[isten];
        weight = current_sten_weight[isten];

        if (iunk < Phys2Unk_first[HSRHOBAR]+Nrho_bar_s) fac=prefactor_rho_bar_s(iunk,jcomp,offset);
        else                                            fac=prefactor_rho_bar_v(iunk,jcomp,offset);

         /* Find in the Stencil position on overall mesh */
        inode_sten =offset_to_node_box(ijk_box, offset, reflect_flag);

        if (inode_sten >= 0 && !Zero_density_TF[inode_sten][jcomp]) {
           if (Lhard_surf) {
               if (Nodes_2_boundary_wall[jlist][inode_sten]!=-1)
                  weight = HW_boundary_weight
                    (jcomp,jlist,current_sten->HW_Weight[isten], inode_sten, reflect_flag);
           }
           if (inode_sten<Nnodes_box && inode_sten >=0){
               sum +=  fac*weight*x[junk][inode_sten];
           }
        }
        else if (inode_sten<0){
             sum += fac*weight*constant_boundary(junk,inode_sten);
        }

     }  /* end of loop over isten */
  }     /* end of loop over jcomp */
  return(sum);
}
/****************************************************************************/
/*int_stencil_CAV: Perform the integral sum(j)int rho_j(r')*prefac*weight[sten] 
for the WTC and WJDC cavity functionals */
 double int_stencil_CAV(double **x,int inode_box,int iunk)
{
  int isten,*offset,inode_sten,ijk_box[3],izone,idim;
  int j,jcomp,junk,icomp,jlist,sten_type;
  double weight, sum,fac;
  struct Stencil_Struct *current_sten;
  int **current_sten_offset, reflect_flag[NDIM_MAX];
  double *current_sten_weight;
  for (idim=0; idim<Ndim; idim++) reflect_flag[idim]=FALSE;

/*  izone = Nodes_to_zone_box[inode_box];*/
  izone = 0;

  sum = 0.0;
  node_box_to_ijk_box(inode_box,ijk_box);

  sten_type=THETA_FN_SIG;

  for (junk=Phys2Unk_first[DENSITY];junk<Phys2Unk_last[DENSITY];junk++){

     if (Lseg_densities) jcomp=Unk2Comp[junk-Phys2Unk_first[DENSITY]];
     else                jcomp=junk-Phys2Unk_first[DENSITY];



     if (Nlists_HW <= 2) jlist = 0;
     else                jlist = jcomp;

     current_sten = &(Stencil[sten_type][izone][jcomp]);
     current_sten_offset = current_sten->Offset;
     current_sten_weight = current_sten->Weight;

     for (isten = 0; isten < current_sten->Length; isten++) {
        offset = current_sten_offset[isten];
        weight = current_sten_weight[isten];

        fac=prefactor_cavity_wtc(iunk,jcomp,offset);

         /* Find in the Stencil position on overall mesh */
        inode_sten =offset_to_node_box(ijk_box, offset, reflect_flag);

        if (inode_sten >= 0 && !Zero_density_TF[inode_sten][jcomp]) {
           if (Lhard_surf) {
               if (Nodes_2_boundary_wall[jlist][inode_sten]!=-1)
                  weight = HW_boundary_weight
                    (jcomp,jlist,current_sten->HW_Weight[isten], inode_sten, reflect_flag);
           }
           if (inode_sten<Nnodes_box && inode_sten >=0){
               sum +=  fac*weight*x[junk][inode_sten];
           }
        }
        else if (inode_sten<0){
             sum += fac*weight*constant_boundary(junk,inode_sten);
        }

     }  /* end of loop over isten */
  }     /* end of loop over jcomp */
  return(sum);
}
/****************************************************************************/
/*int_stencil_BondWTC: Perform the integral sum(j)int rho_j(r')*prefac*weight[sten] 
for the WTC and WJDC cavity functionals */
 double int_stencil_BondWTC(double **x,int inode_box,int iunk)
{
  int isten,*offset,inode_sten,ijk_box[3],izone,idim;
  int j,jcomp,junk,icomp,jlist,sten_type;
  int unk_bond,iseg,jseg,ibond;
  double weight, sum,fac;
  struct Stencil_Struct *current_sten;
  int **current_sten_offset, reflect_flag[NDIM_MAX];
  double *current_sten_weight;
  for (idim=0; idim<Ndim; idim++) reflect_flag[idim]=FALSE;

/*  izone = Nodes_to_zone_box[inode_box];*/
  izone = 0;

  sum = 0.0;
  node_box_to_ijk_box(inode_box,ijk_box);

  sten_type=DELTA_FN_BOND;

  unk_bond=iunk-Phys2Unk_first[BONDWTC];
  iseg=BondAll_to_isegAll[unk_bond];
  ibond=BondAll_to_ibond[unk_bond];
  jseg=Bonds_SegAll[iseg][ibond];

  jcomp = Unk2Comp[jseg];
  icomp = Unk2Comp[iseg];

  junk = jseg + Phys2Unk_first[DENSITY];

  if (Nlists_HW <= 2) jlist = 0;
     else            jlist = jcomp;

  current_sten = &(Stencil[sten_type][izone][icomp+Ncomp*jcomp]);
  current_sten_offset = current_sten->Offset;
  current_sten_weight = current_sten->Weight;

  for (isten = 0; isten < current_sten->Length; isten++) {
     offset = current_sten_offset[isten];
     weight = current_sten_weight[isten];

    /* Find in the Stencil position on overall mesh */
    inode_sten =offset_to_node_box(ijk_box, offset, reflect_flag);

    if (inode_sten >= 0 && !Zero_density_TF[inode_sten][jcomp]) {
        if (Lhard_surf) {
            if (Nodes_2_boundary_wall[jlist][inode_sten]!=-1)
               weight = HW_boundary_weight
                 (jcomp,jlist,current_sten->HW_Weight[isten], inode_sten, reflect_flag);
        }
        if (inode_sten<Nnodes_box && inode_sten >=0){
            sum +=  weight*x[junk][inode_sten];
        }
     }
     else if (inode_sten<0){
          sum += weight*constant_boundary(junk,inode_sten);
     }

  }  /* end of loop over isten */
  return(sum);
}
/****************************************************************************/
