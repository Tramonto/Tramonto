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
    switch(Unk2Phys[iunk]){
       case DENSITY:
          if (Type_poly_TC){
             if (jnode_box==-1) bcval=Rho_seg_b[iunk-Phys2Unk_first[DENSITY]];
             if (jnode_box==-2) bcval=0.0;
             else if (jnode_box==-3) 
                 if(Lsteady_state) bcval=Rho_seg_LBB[iunk-Phys2Unk_first[DENSITY]];
             else if (jnode_box==-4) 
                 if(Lsteady_state) bcval=Rho_seg_RTF[iunk-Phys2Unk_first[DENSITY]];
           }
           else{
              if (jnode_box==-1)     bcval = Rho_b[iunk-Phys2Unk_first[DENSITY]];
              if (jnode_box==-2)     bcval = 0.0;
              else if (jnode_box==-3){
                  if (Lsteady_state) bcval = Rho_b_LBB[iunk-Phys2Unk_first[DENSITY]];
                  else               bcval = Rho_coex[1];
              }
              else if (jnode_box==-4){
                  if (Lsteady_state) bcval = Rho_b_RTF[iunk-Phys2Unk_first[DENSITY]];
                  else               bcval = Rho_coex[2];
              }
           }
           break;
       case RHOBAR_ROSEN:
           if (jnode_box==-1)       bcval = Rhobar_b[iunk-Phys2Unk_first[RHOBAR_ROSEN]];
           if (jnode_box==-2)       bcval = 0.0;  /* assuming wall begins in domain and rhobars 
                                                     have decayed beyond the boundary */
           else if (jnode_box==-3)  bcval = Rhobar_b_LBB[iunk-Phys2Unk_first[RHOBAR_ROSEN]];
           else if (jnode_box==-4)  bcval = Rhobar_b_RTF[iunk-Phys2Unk_first[RHOBAR_ROSEN]];
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
           else if (jnode_box==-3)  bcval = Betamu_LBB[iunk - Phys2Unk_first[DIFFUSION]];
           else if (jnode_box==-4)  bcval = Betamu_RTF[iunk - Phys2Unk_first[DIFFUSION]];
           break;
       case CAVITY_WTC:
           if (jnode_box==-1)      bcval=Xi_cav_b[iunk-Phys2Unk_first[CAVITY_WTC]+2]; 
           else if (jnode_box==-3) bcval=Xi_cav_LBB[iunk-Phys2Unk_first[CAVITY_WTC]+2];
           else if (jnode_box==-4) bcval=Xi_cav_RTF[iunk-Phys2Unk_first[CAVITY_WTC]+2];
           break;
       case BOND_WTC:
           if (jnode_box==-1)      bcval=BondWTC_b[iunk-Phys2Unk_first[BOND_WTC]]; 
           else if (jnode_box==-3) bcval=BondWTC_LBB[iunk-Phys2Unk_first[BOND_WTC]];
           else if (jnode_box==-4) bcval=BondWTC_RTF[iunk-Phys2Unk_first[BOND_WTC]];
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
int loc_find(int iunk,int inode,int flag)
{
  int loc_i;
  if (MATRIX_FILL_NODAL) loc_i = iunk + Nunk_per_node * inode;
  else{
     if (flag == LOCAL)
        loc_i = inode + Nnodes_per_proc*iunk;
     else if (flag == BOX)
        loc_i = inode + Nnodes_box*iunk;
     else if (flag == GLOBAL)
        loc_i = inode + Nnodes*iunk;
  }
  return loc_i;
}
/*****************************************************************************************************/
