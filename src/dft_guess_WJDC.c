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
 *  FILE: dft_guess_noinfo_WJDC.c
 *
 *  This file contains routines that set up an initial guess for the
 *  WJDC polymer DFT.
 *
 */

#include "dft_guess_WJDC.h"
 
/*********************************************************/
/*setup_polymer_field: in this routine sets up the initial guess for the WJDC field variable */
void setup_polymer_field_wjdc(double **xInBox)
{
  int loc_inode,itype_mer,irho, iunk,i,Nloop,inode_box,iref;
  double field;

  if (Type_poly==WJDC)                           Nloop=Nseg_tot; 
  else if (Type_poly==WJDC2 || Type_poly==WJDC3) Nloop=Ncomp;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box=L2B_node[loc_inode];
     for (i=0; i<Nloop; i++){
         iunk=Phys2Unk_first[WJDC_FIELD]+i;
         if (Type_poly==WJDC){
             if (!Zero_density_TF[inode_box][Unk2Comp[i]]) xInBox[iunk][inode_box]=Field_WJDC_b[Unk2Comp[i]];
             else                            xInBox[iunk][inode_box]=0.;
         }
         else if (Type_poly==WJDC2 || Type_poly==WJDC3){
             if (!Zero_density_TF[inode_box][i]) xInBox[iunk][inode_box]=Field_WJDC_b[i];
             else                            xInBox[iunk][inode_box]=0.;
         }
     }
   }
   return;
}
/*********************************************************/
/*calc_init_WJDC_field: in this routine sets up the initial guess for the WJDC field variable */
void calc_init_WJDC_field(double **xInBox)
{
  int loc_inode,inode_box,ijk_box[3],icomp,iunk,izone,mesh_coarsen_flag_i,i,Nloop;
  double resid_EL;
  struct  RB_Struct *dphi_drb=NULL;
  izone=0;
  mesh_coarsen_flag_i=0;

  if (Type_func !=NONE){
     dphi_drb = (struct RB_Struct *) array_alloc (1, Nnodes_box, sizeof(struct RB_Struct));
     FMT1stDeriv_switch(xInBox,dphi_drb);
  }

  if (Type_poly==WJDC) Nloop=Nseg_tot;
  if (Type_poly==WJDC2 || Type_poly==WJDC3) Nloop=Ncomp;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box=L2B_node[loc_inode];
     node_box_to_ijk_box(inode_box, ijk_box);

     for (i=0; i<Nloop; i++){  
        if (Type_poly==WJDC) icomp=Unk2Comp[i];
        else icomp=i;
        iunk=Phys2Unk_first[WJDC_FIELD]+i;

        if (!Zero_density_TF[inode_box][icomp]){
          resid_EL=load_euler_lagrange(iunk,loc_inode,inode_box,ijk_box,izone,
                          xInBox,dphi_drb,mesh_coarsen_flag_i,INIT_GUESS_FLAG);
          xInBox[iunk][inode_box]=resid_EL;
       }
       else xInBox[iunk][inode_box]=0.0;
     }
  }
  if (Type_func != NONE) safe_free((void *) &dphi_drb);
  return;
}
/*********************************************************/
/*setup_polymer_G_wjdc: in this routine sets up the initial guess for the chain variable
in the wjdc functional */
void setup_polymer_G_wjdc(double **xInBox)
{
  int loc_inode,itype_mer,irho, iunk,i,Nloop,inode_box,iseg,icomp_iseg;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box=L2B_node[loc_inode];
     for (i=0; i<Nbonds; i++){
         iseg=Unk_to_Seg[i];
         icomp_iseg=Unk2Comp[iseg];
         iunk=Phys2Unk_first[G_CHAIN]+i;
         if (!Zero_density_TF[inode_box][icomp_iseg]) xInBox[iunk][inode_box]=G_WJDC_b[i];
         else xInBox[iunk][inode_box]=0.0;
     }
   }
   return;
}
/*********************************************************/
/*calc_init_polymer_G_wjdc: in this routine sets up the initial guess for the chain variable
in the wjdc functional */
void calc_init_polymer_G_wjdc(double **xInBox)
{
  int loc_inode,itype_mer,irho, iunk,i,Nloop,inode_box,icomp_iseg;
  int ibond,jbond,index,iseg,jseg,pol_num,bond_num,test,ijk_box[3];
  double resid_G;
  int array_val[NMER_MAX*NBOND_MAX],array_fill,count_fill;
  double (*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **);
  double (*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **);

  fp_ResidG=&WJDC_Resid_GCHAIN;
  fp_ResidG_Bulk=&WJDC_Resid_Bulk_GCHAIN;

  /* need to be careful to generate the G's in the order dictated
     by the chain architecture.  Use same strategy as in dft_thermo_wjdc */

  for (ibond=0;ibond<Nbonds;ibond++) array_val[ibond]=FALSE;
  array_fill=FALSE;
  count_fill=0;

  while (array_fill==FALSE){
     for (ibond=0;ibond<Nbonds;ibond++){
        pol_num=Unk_to_Poly[ibond];
        iseg=Unk_to_Seg[ibond];
        icomp_iseg=Unk2Comp[iseg];
        bond_num=Unk_to_Bond[ibond];

        if (array_val[ibond]==FALSE){
           test=TRUE;  /* assume we will compute a bulk G */
           jseg=Bonds[pol_num][iseg][bond_num];
           if (jseg != -1 ){   /* may need to skip this G if we don't have all information yet  -
                                  always compute G for end segments flagged with -1 value */
              for (jbond=0;jbond<Nbond[pol_num][jseg];jbond++){
                 if (Bonds[pol_num][jseg][jbond] != iseg){ /* check all jbonds to see if we have necessary info */
                    index=Poly_to_Unk[pol_num][jseg][jbond]+Geqn_start[pol_num]-Geqn_start[0];
                    if (array_val[index]==FALSE) test=FALSE;
                 }
              }
           }
           if (test==TRUE){     /* compute a bulk G */
               for (loc_inode=0;loc_inode<Nnodes_per_proc;loc_inode++){
                     inode_box=L2B_node[loc_inode];
                     node_box_to_ijk_box(inode_box,ijk_box);
                     iunk=Phys2Unk_first[G_CHAIN]+ibond;
                     resid_G=load_Chain_Geqns(WJDC_FIELD,0,0,
                                             NULL,fp_ResidG,fp_ResidG_Bulk,
                                             iunk,loc_inode,inode_box,
                                             ijk_box,0,xInBox, INIT_GUESS_FLAG);
                
                     xInBox[iunk][inode_box]=resid_G;
               }
               communicate_to_fill_in_box_values(xInBox);  /* we need every G to be updated for all nodes in box as we
                                                             generate them so we will be able to perform the next integral */
              count_fill++;
              array_val[ibond]=TRUE;

           } /*end of test=TRUE condition */
        }
     }
     if (count_fill==Nbonds) array_fill=TRUE;
   }
   return;
}
/*********************************************************/
