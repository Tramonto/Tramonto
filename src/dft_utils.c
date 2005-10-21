/* ---------------------------------------------------------
Here are some useful routines that are used many times in 
different places in the DFT code - 
------------------------------------------------------------*/
#include <stdio.h>
#include <assert.h>
#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"
/***********************************************************************
int_stencil_bulk: this routine sums the appropriate stencil to get
                  the bulk contributions to various terms in the E-L
                  equation. Note that some terms (attractions, WCA electrostatics,
                  CMS polymers bury a multiplier function with the stencil weight function.*/
double int_stencil_bulk(int sten_type,int icomp,int jcomp)
{
  int izone, isten;
  double sum, weight, *sten_weight;
  struct Stencil_Struct *sten;

  sum = 0.0;
  izone = 0;
  sten = &(Stencil[sten_type][izone][icomp+Ncomp*jcomp]);
  sten_weight = sten->Weight;

  for (isten = 0; isten<sten->Length; isten++){
     weight = sten_weight[isten];
     sum += weight;
  }
  return(sum);
}
/*******************************************************************************/
/*int_stencil: Perform the integral sum(j)int rho_j(r')*weight[sten] */
 double int_stencil(double **x,int inode_box,int icomp,int sten_type)
{
  int isten,*offset,inode_sten,ijk_box[3],izone,idim;
  int j,jcomp,junk;
  double weight, sum;
  struct Stencil_Struct *current_sten;
  int **current_sten_offset, reflect_flag[NDIM_MAX];
  double *current_sten_weight;
  for (idim=0; idim<Ndim; idim++) reflect_flag[idim]=FALSE;

/*  izone = Nodes_to_zone_box[inode_box];*/
  izone = 0;

  sum = 0.0;
  node_box_to_ijk_box(inode_box,ijk_box);

  for (junk=Phys2Unk_first[DENSITY];junk<Phys2Unk_last[DENSITY];junk++){
     jcomp=Unk2Comp[junk];

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
               if (Nodes_2_boundary_wall[Nlists_HW-1][inode_sten]!=-1)
                  weight = HW_boundary_weight
                    (jcomp,Nlists_HW-1,current_sten->HW_Weight[isten],
                                      inode_sten, reflect_flag);
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
  return (sum);
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
    switch(Unk2Phys[iunk]){
       case DENSITY:
          if (Type_poly==WTC){
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
void integrateInSpace(double(*fp_integrand)(int,int,int,double **),int iunk,int i, double **x)
{

   double sum,sum_i;
   int loc_inode,inode_box,iwall;
   double area;

   sum_i=0.0,sum=0.0;

   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
      inode_box = L2B_node[loc_inode];
      sum_i += (*fp_integrand)(iunk,i,inode_box,x)*Vol_el/((double)Nnodes_per_el_V);
  }       /* end of loc_inode loop */

  sum=gsum_double(sum_i);

  if (Lper_area && Area>0.0) {
      sum /= Area;
      sum *= (Fac_vol/Fac_area);
  }
  else sum *= Fac_vol;

  Temporary_sum=sum;
  return;
}
/*****************************************************************************************************/
void integrateInSpace_SumInComp(double(*fp_integrand)(int,int,int,double**),int i,double **x)
{

   double sum,sum_i;
   int loc_inode,inode_box,iloop,nloop,iunk;

   sum_i=0.0,sum=0.0;

   nloop=Ncomp;
   if (Type_poly==WTC) nloop=Nseg_tot;

   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
    inode_box = L2B_node[loc_inode];
    for (iloop=0; iloop<nloop; iloop++){

         iunk = Phys2Unk_first[DENSITY]+iloop;
         sum_i += (*fp_integrand)(iunk,i,inode_box,x)*Vol_el/((double)Nnodes_per_el_V);

    }   
  }       /* end of loc_inode loop */

  sum=gsum_double(sum_i);

  if (Lper_area && Area>0.0) {
      sum /= Area;
      sum *= (Fac_vol/Fac_area);
  }
  else{
      sum*= Fac_vol;
  }
  Temporary_sum=sum;
  return;
}
/**************************************************************************************/
void print_to_screen(int i,double val,char *var_label)
{

  if (i==0)printf("\t=========================================\n");
  printf("\n \t Summary of %s RESULTS Results\n",var_label);

  if ((i==0 && !Lhard_surf) || (i==1 && Lhard_surf)){
       printf("\t===== Integrating over domain volume =====\n");
  }
  else{
       printf("\t===== Integrating over fluid volume only =====\n");
  }
  printf("\t\t %s=%9.6f\n",var_label,val);

  return;
}
/**************************************************************************************/
/**************************************************************************************/
void print_to_screen_comp(int i,int icomp,double val,char *var_label)
{

  if (icomp==0 && i==0)printf("\t=========================================\n");
  if (icomp==0) printf("\n \t Summary of %s RESULTS Results\n",var_label);

  if ((i==0 && !Lhard_surf) || (i==1 && Lhard_surf)){
        if (icomp==0) printf("\t===== Integrating over domain volume =====\n");
  }
  else{
       printf("\t===== Integrating over fluid volume only =====\n");
  }
  printf("\t\t %s[icomp=%d]=%9.6f\n",var_label,icomp,val);

  return;
}
/**************************************************************************************/
