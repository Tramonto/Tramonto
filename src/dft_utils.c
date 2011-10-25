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
dft_utils.c:

Here are some useful routines that are used many times in 
different places in the DFT code - 
------------------------------------------------------------*/


#include <stdio.h>
#include "dft_utils.h"

/***********************************************************************/
/* This routine assembles the weight of boundary nodes in the Hard Wall model*/
double  HW_boundary_weight(int icomp,int ilist, double *hw_weight,
                           int inode_box, int *reflect_flag)
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
int find_jzone(int izone,int inode_box)
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
  else if (Coarser_jac==0){
     if (Mesh_coarsening || Nwall==0) jzone = izone;
     else jzone = Mesh_coarsen_flag[inode_box];
  }
  else jzone=izone;
  return jzone;
}
/****************************************************************************/
/*constant_boundary:  This routine just returns a boundary condition
given an unknown number and an indication of what kind of boundary we have.*/
double constant_boundary(int iunk,int jnode_box)
{
    double bcval;
    switch(Unk2Phys[iunk]){
       case DENSITY:
          if (Lseg_densities){
             if (jnode_box==-1)      bcval=Rho_seg_b[iunk-Phys2Unk_first[DENSITY]];
             if (jnode_box==-2)      bcval=0.0;
             else if (jnode_box==-3) bcval=Rho_seg_LBB[iunk-Phys2Unk_first[DENSITY]];
             else if (jnode_box==-4) bcval=Rho_seg_RTF[iunk-Phys2Unk_first[DENSITY]];
           }
           else{
              if (jnode_box==-1)     bcval = Rho_b[iunk-Phys2Unk_first[DENSITY]];
              if (jnode_box==-2)     bcval = 0.0;
              else if (jnode_box==-3)bcval = Rho_b_LBB[iunk-Phys2Unk_first[DENSITY]];
              else if (jnode_box==-4)bcval = Rho_b_RTF[iunk-Phys2Unk_first[DENSITY]];
           }
           break;
       case HSRHOBAR:
           if (jnode_box==-1)       bcval = Rhobar_b[iunk-Phys2Unk_first[HSRHOBAR]];
           if (jnode_box==-2)       bcval = 0.0;  /* assuming wall begins in domain and rhobars
                                                     have decayed beyond the boundary */
           else if (jnode_box==-3)  bcval = Rhobar_b_LBB[iunk-Phys2Unk_first[HSRHOBAR]];
           else if (jnode_box==-4)  bcval = Rhobar_b_RTF[iunk-Phys2Unk_first[HSRHOBAR]];
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
           if   (jnode_box==-1)       bcval = 0.;
           else if (jnode_box==-2)       bcval = -VEXT_MAX;
           else{
              if (Type_poly==NONE){
                 if (jnode_box==-3)  bcval = Betamu_LBB[iunk - Phys2Unk_first[DIFFUSION]];
                 else if (jnode_box==-4)  bcval = Betamu_RTF[iunk - Phys2Unk_first[DIFFUSION]];
              }
              else{
                 if (jnode_box==-3)  bcval = Betamu_chain_LBB[iunk - Phys2Unk_first[DIFFUSION]];
                 else if (jnode_box==-4)  bcval = Betamu_chain_RTF[iunk - Phys2Unk_first[DIFFUSION]];
              }
           }
           break;
       case CAVWTC:
           if (jnode_box==-1)      bcval=Xi_cav_b[iunk-Phys2Unk_first[CAVWTC]+2];
	   else if (jnode_box==-2) bcval=0.0;
           else if (jnode_box==-3) bcval=Xi_cav_LBB[iunk-Phys2Unk_first[CAVWTC]+2];
           else if (jnode_box==-4) bcval=Xi_cav_RTF[iunk-Phys2Unk_first[CAVWTC]+2];
           break;
       case BONDWTC:
           if (jnode_box==-1)      bcval=BondWTC_b[iunk-Phys2Unk_first[BONDWTC]];
           else if (jnode_box==-3) bcval=BondWTC_LBB[iunk-Phys2Unk_first[BONDWTC]];
           else if (jnode_box==-4) bcval=BondWTC_RTF[iunk-Phys2Unk_first[BONDWTC]];
           break;
       case CMS_FIELD:
           if (jnode_box==-2) bcval=0.0;
		   else bcval=1.0; 
		   break;
       case WJDC_FIELD:
           if (jnode_box==-2) bcval=0.0;
           else if (jnode_box==-1) bcval=Field_WJDC_b[iunk-Phys2Unk_first[WJDC_FIELD]];
           else if (jnode_box==-3) bcval=Field_WJDC_LBB[iunk-Phys2Unk_first[WJDC_FIELD]];
           else if (jnode_box==-4) bcval=Field_WJDC_RTF[iunk-Phys2Unk_first[WJDC_FIELD]];
           break;
		case SCF_FIELD:
		   if (jnode_box==-2) bcval=0.0;
		   else if (jnode_box==-3 || jnode_box==-4){
               printf("diffusion boundaries not fully implemented for SCF fluid \n");
               exit(-1);
		   }
		   else bcval=1.0; 
		   break;
	   case SCF_CONSTR:
		   if (jnode_box==-2) bcval=0.0;
		   else if (jnode_box==-3 || jnode_box==-4){
               printf("diffusion boundaries not fully implemented for SCF fluid \n");
               exit(-1);
		   }
		   else bcval=1.0; 
		   break;

       case G_CHAIN:
           if (jnode_box==-2) bcval=0.0;
	   else{
               if (Type_poly == CMS || Type_poly == CMS_SCFT) bcval=1.0; 
               else if (Type_poly == WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){
                   if (jnode_box==-1) bcval = G_WJDC_b[iunk-Phys2Unk_first[G_CHAIN]];
                   if (jnode_box==-3) bcval = G_WJDC_LBB[iunk-Phys2Unk_first[G_CHAIN]];
                   if (jnode_box==-4) bcval = G_WJDC_RTF[iunk-Phys2Unk_first[G_CHAIN]];
               }
           }
			   break;
   }
   return(bcval);
}
/****************************************************************************/
int loc_find(int iunk,int inode,int flag)
{
  int loc_i;
  if (MATRIX_FILL_NODAL) loc_i = iunk + Nunk_per_node * inode;
  else{
     if (flag == LOCAL_N)
        loc_i = inode + Nnodes_per_proc*iunk;
     else if (flag == BOX)
        loc_i = inode + Nnodes_box*iunk;
     else if (flag == GLOBAL)
        loc_i = inode + Nnodes*iunk;
  }
  return loc_i;
}
/****************************************************************************/
void solutionVec_to_nOrdering(double *rhoBar_SVOrdering, double *n)
{
  int idim,iv1,iv2;

  n[0] = rhoBar_SVOrdering[3];
  n[1] = rhoBar_SVOrdering[2];
  n[2] = rhoBar_SVOrdering[1];
  n[3] = rhoBar_SVOrdering[0];

  for (idim=0; idim<Ndim; idim++){
    iv1=Nrho_bar_s+idim;
    iv2=Nrho_bar_s+Ndim+idim;
    n[iv1] = rhoBar_SVOrdering[iv2];
    n[iv2] = rhoBar_SVOrdering[iv1];
  }
  return;
}
/*****************************************************************************************************/
void pass_part_of_solnVector(double **xOwned, double **x,int iunk_start,int nunk_to_pass)
{
 int inode,iunk;
 double *nodal_owned, *nodal_box;
 nodal_owned = (double *) array_alloc (1, Nnodes_per_proc, sizeof(double));
 nodal_box = (double *) array_alloc (1, Nnodes_box, sizeof(double));

 for (iunk=iunk_start; iunk<iunk_start+nunk_to_pass; iunk++){
    for (inode=0;inode<Nnodes_per_proc;inode++) nodal_owned[inode]=xOwned[iunk][inode];
    for (inode=0;inode<Nnodes_box;inode++)      nodal_box[inode]=x[iunk][inode];
 
    (void) dft_linprobmgr_importnodalr2c(LinProbMgr_manager,nodal_owned,nodal_box);

    for (inode=0;inode<Nnodes_box;inode++)      x[iunk][inode]=nodal_box[inode];
 }

 safe_free((void *) &nodal_owned);
 safe_free((void *) &nodal_box);

 return;
}
/*****************************************************************************************************/
void print_to_screen(double val,char *var_label)
{
  printf("\t\t %s=%g\n",var_label,val); return;
}
/**************************************************************************************/
void print_to_screen_comp(int icomp,double val,char *var_label)
{
  printf("\t\t %s[icomp=%d]=%g\n",var_label,icomp,val); return;
}
/**************************************************************************************/
void print_to_file(FILE *fp,double val,char *var_label,int first)
{
  if (first != FALSE)  fprintf(fp,"%s  ",var_label); 
  if (first != TRUE)   fprintf(fp,"%g  ",val); 
  return;
}
/**************************************************************************************/
void print_to_file_comp(FILE *fp,int icomp,double val,char *var_label,int first)
{
  if (first != FALSE) fprintf(fp,"%s[%d]  ",var_label,icomp); 
  if (first != TRUE)  fprintf(fp,"%g  ",val); 
  return;
}
/**************************************************************************************/
