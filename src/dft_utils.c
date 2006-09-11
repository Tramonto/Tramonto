/*
//@HEADER
// ********************************************************************
// Copyright (2006) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000, there is a non-exclusive license for use of this
// work by or on behalf of the U.S. Government. Export of this program
// may require a license from the United States Government.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// ********************************************************************
//@HEADER
*/

/* ---------------------------------------------------------
Here are some useful routines that are used many times in 
different places in the DFT code - 
------------------------------------------------------------*/
#include <stdio.h>
#include "dft_globals_const.h"
/***********************************************************************
int_stencil_bulk: this routine sums the appropriate stencil to get
                  the bulk contributions to various terms in the E-L
                  equation. Note that some terms (attractions, WCA electrostatics,
                  CMS polymers bury a multiplier function with the stencil weight function.*/
void int_stencil_bulk(int sten_type,int icomp,int jcomp,double(*fp_integrand)(double,int,int))
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
  Temporary_sum=sum;
  return;
}
/*******************************************************************************/
/*int_stencil: Perform the integral sum(j)int rho_j(r')*weight[sten] */
 void int_stencil(double **x,int inode_box,int iunk,int sten_type)
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
  if (Type_poly==WTC) icomp=Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
  else                icomp=iunk-Phys2Unk_first[DENSITY];

  for (junk=Phys2Unk_first[DENSITY];junk<Phys2Unk_last[DENSITY];junk++){
     if (Type_poly==WTC) jcomp=Unk2Comp[junk-Phys2Unk_first[DENSITY]];
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
  Temporary_sum=sum;
  return;
}
/****************************************************************************/
void resid_and_Jac_sten_fill_sum_Ncomp (int sten_type, double **x, int iunk, 
  int loc_inode, int inode_box, int izone,
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
  double resid,mat_val,resid_sum=0.0,weight_bulk,bulk_term;

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

  if (Type_poly==WTC) loop_max=Nseg_tot;
  else                loop_max=Ncomp;

  for (jloop=0; jloop<loop_max; jloop++){
      if (Type_poly==WTC) {
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


      if (sten_type==U_ATTRACT || sten_type==THETA_CHARGE) {
            i=iunk-Phys2Unk_first[DENSITY];
            if (Type_poly==WTC) icomp=Unk2Comp[i];
            else icomp=i;
            index=icomp+Ncomp*jcomp;
      }
      else if (sten_type==POLYMER_CR){
            icomp=iunk-Phys2Unk_first[CMS_FIELD];
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
        if (sten_type==THETA_CHARGE || sten_type==POLYMER_CR) fac=-1.0;
        if (sten_type==POLYMER_CR) bulk_term= weight_bulk*Rho_b[jcomp];
        else bulk_term=0.0;

         jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);
         if (Lhard_surf && jnode_box >=0) {
             if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1)
                weight = HW_boundary_weight
                 (jcomp,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
         }
         resid =  fac*(weight*((*fp_resid)(junk,jnode_box,x)) - bulk_term);
         dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
         resid_sum+=resid;

         if (!resid_only_flag)
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
  Temporary_sum=resid_sum;
  return;
}
/****************************************************************************/
void resid_and_Jac_sten_fill (int sten_type, double **x, int iunk, int junk,
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
     dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
     resid_sum+=resid;

     if (!resid_only_flag)
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
  Temporary_sum=resid_sum;
  return;
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
          if (Type_poly==WTC){
             if (jnode_box==-1) {
                 bcval=Rho_seg_b[iunk-Phys2Unk_first[DENSITY]];
             }
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
           if (jnode_box==-1)       bcval = 0.;
           if (jnode_box==-2)       bcval = -VEXT_MAX;
           else if (jnode_box==-3)  bcval = Betamu_LBB[iunk - Phys2Unk_first[DIFFUSION]];
           else if (jnode_box==-4)  bcval = Betamu_RTF[iunk - Phys2Unk_first[DIFFUSION]];
           break;
       case CAVWTC:
           if (jnode_box==-1)      bcval=Xi_cav_b[iunk-Phys2Unk_first[CAVWTC]+2];
           else if (jnode_box==-3) bcval=Xi_cav_LBB[iunk-Phys2Unk_first[CAVWTC]+2];
           else if (jnode_box==-4) bcval=Xi_cav_RTF[iunk-Phys2Unk_first[CAVWTC]+2];
           break;
       case BONDWTC:
           if (jnode_box==-1)      bcval=BondWTC_b[iunk-Phys2Unk_first[BONDWTC]];
           else if (jnode_box==-3) bcval=BondWTC_LBB[iunk-Phys2Unk_first[BONDWTC]];
           else if (jnode_box==-4) bcval=BondWTC_RTF[iunk-Phys2Unk_first[BONDWTC]];
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
void integrateInSpace(double(*fp_integrand)(int,int,double **),int iunk,
                                                int **nelhit,double **x,double *profile)
{

   double sum,sum_i,integrand;
   int loc_inode,inode_box,iwall;
   double area;

   sum_i=0.0,sum=0.0;
   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
      inode_box = L2B_node[loc_inode];
      integrand = (*fp_integrand)(iunk,inode_box,x);
      sum_i += integrand*nelhit[iunk-Phys2Unk_first[DENSITY]][inode_box]*Vol_el/((double)Nnodes_per_el_V);
      if (profile != NULL) profile[loc_inode]+=integrand;
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
void integrateInSpace_SumInComp(double(*fp_integrand)(int,int,double**),
                                                int **nelhit,double **x,double *profile)
{

   double sum,sum_i,integrand;
   int loc_inode,inode_box,iloop,nloop,iunk;

   sum_i=0.0,sum=0.0;

   nloop=Ncomp;
   if (Type_poly==WTC) nloop=Nseg_tot;

   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
    inode_box = L2B_node[loc_inode];
    for (iloop=0; iloop<nloop; iloop++){

         iunk = Phys2Unk_first[DENSITY]+iloop;
         integrand = (*fp_integrand)(iunk,inode_box,x);
         sum_i += integrand*
                   nelhit[iloop][inode_box]*Vol_el/((double)Nnodes_per_el_V);

         if (profile != NULL) profile[loc_inode]+=integrand;
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
/*****************************************************************************************************/
void integrateOverSurface(double(*fp_integrand)(int,int,int,double **),int iunk,double **x,double *profile)
{

  double sum,sum_i,integrand,fac;
  int loc_inode,inode_box,iwall,ijk[3],idim;

  sum_i=0.0,sum=0.0;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){

      inode_box = L2B_node[loc_inode];     
      iwall=Nodes_2_boundary_wall[Nlists_HW-1][inode_box];
      if  (iwall != -1){ /* identify surface nodes */

         integrand = (*fp_integrand)(iunk,inode_box,iwall,x);

         fac=1.0;
         node_to_ijk(L2G_node[loc_inode],ijk);
         for (idim=0;idim<Ndim;idim++){ 
            if ((ijk[idim]==0 || ijk[idim]==Nodes_x[idim]-1)  && Type_bc[idim][1] != PERIODIC ) fac*=0.5;
         } 
         sum_i += integrand*fac;
         if (profile != NULL) profile[loc_inode]+=integrand;
      }       /* end of loc_inode loop */
  }

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
/*setup_domain_multipliers: Here compute area and Fac_vol and Fac_area that
  are needed for all integrals */
void setup_domain_multipliers()
{
  int idim,iwall;

  /* compute Fac_vol and Fac_area */
   Fac_area = 1.0;
   Fac_vol = 1.0;
   for (idim = 0; idim<Ndim; idim++) {
       if (!(Type_bc[idim][0] == REFLECT && Type_bc[idim][1] == REFLECT) && Lcount_reflect){

       if (Type_bc[idim][0] == REFLECT){

          Fac_vol *= 2.0;
          if (WallPos[idim][0] == -0.5*Size_x[idim]) Fac_area *= 2.0;
       }
       else if (Type_bc[idim][1] == REFLECT) {

          Fac_vol *= 2.0;
          if (WallPos[idim][0] == 0.5*Size_x[idim]) Fac_area *= 2.0;
       }
       }
   }


  /* compute surface area */
  Area = 0.0;
  if (Nwall == 0) Area = 1.0;
  else{
     if (Nlink == Nwall) Area = S_area_tot[Nlists_HW-1][0];
     else
        for (iwall=0; iwall<Nwall; iwall++){
           if (Link[iwall]==0)
           Area += S_area_tot[Nlists_HW-1][iwall];
        }
  }

  return;
}
/***************************************************************************/
/* setup_integrals:  here we store arrays of elements hit per node and
                     the list propertiese for post processing integrated
                     parameters (adsorption, free energy etc.). */
void setup_integrals()
{
  int loc_inode, icomp, iel, iunk, inode_box, nel_hit,inode,iel_box,
      nel_hit2,ilist,idim,ielement,semiperm,iwall,jcomp;
  int reflect_flag[3],ijk[3],i,nloop,iloop;
  
  
  for (idim=0; idim<Ndim; idim++) reflect_flag[idim]=FALSE;
  
  if (Type_poly==WTC) nloop=Nseg_tot;
  else                nloop=Ncomp;

  Nel_hit = (int **) array_alloc (2,nloop, Nnodes_box, sizeof(int));
  Nel_hit2 = (int **) array_alloc (2,nloop, Nnodes_box, sizeof(int));
  
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
      inode = L2G_node[loc_inode];
      inode_box = L2B_node[loc_inode];
      node_to_ijk(inode,ijk);
      
      for (iloop=0; iloop<nloop; iloop++){
      
         if (Type_poly==WTC) icomp=Unk2Comp[iloop];
         else                icomp=iloop;
         
         if (Lhard_surf){
            if (Nlists_HW==1 || Nlists_HW==2) List[0] = 0;
            else                              List[0] = icomp;
            List[1] = Nlists_HW-1;            
         }  
         else{
             List[0]=0;
             List[1]=0;
         }

         if (Lhard_surf && Nwall>0) Imax = 2;
         else Imax = 1;

         /* define integration rules to density discontinuities */
         nel_hit = Nnodes_per_el_V;
         for (iel=0; iel<Nnodes_per_el_V; iel++){
             ielement = node_to_elem(inode,iel,reflect_flag);
             iel_box = el_to_el_box(ielement);
             if (ielement == -1) nel_hit--;
             else if (ielement != -2){
                iwall =  Wall_elems[List[1]][iel_box];
                semiperm=FALSE;
                for (jcomp=0; jcomp<Ncomp; jcomp++)
                  if (iwall>=0 && Lsemiperm[WallType[iwall]][jcomp]) semiperm=TRUE;
                if (iwall !=-1 && !semiperm) nel_hit--;
             }
         }

         /* define integration rules to core surfaces */
         nel_hit2 = Nnodes_per_el_V;
         for (iel=0; iel<Nnodes_per_el_V; iel++){
             ielement = node_to_elem(inode,iel,reflect_flag);
             iel_box = el_to_el_box(ielement);
             if (ielement == -1) nel_hit2--;
             else if (ielement != -2){
                iwall =  Wall_elems[List[0]][iel_box];
                semiperm=FALSE;
                for (jcomp=0; jcomp<Ncomp; jcomp++)
                  if (iwall>=0 && Lsemiperm[WallType[iwall]][jcomp]) semiperm=TRUE;
                if (iwall!=-1 && !semiperm) nel_hit2--;
             }
         }

         for (idim=0; idim<Ndim; idim++){
            if ( (Type_bc[idim][0] == REFLECT || Type_bc[idim][0] == IN_BULK || Type_bc[idim][0]==LAST_NODE)
                                                          &&  ijk[idim] == 0)  {
                  nel_hit /= 2;
                  nel_hit2 /= 2;
            }
            if ( (Type_bc[idim][1] == REFLECT || Type_bc[idim][1] == IN_BULK || Type_bc[idim][0]==LAST_NODE)
                                            &&  ijk[idim] == Nodes_x[idim]-1)  {
                  nel_hit /= 2;
                  nel_hit2 /= 2;
            }
         }

         Nel_hit[iloop][inode_box]=nel_hit;
         Nel_hit2[iloop][inode_box]=nel_hit2;

      } /* end iloop loop */
    }
}
/**************************************************************************************/
void print_to_screen(double val,char *var_label)
{
  printf("\t\t %s=%9.6f\n",var_label,val); return;
}
/**************************************************************************************/
void print_to_screen_comp(int icomp,double val,char *var_label)
{
  printf("\t\t %s[icomp=%d]=%9.6f\n",var_label,icomp,val); return;
}
/**************************************************************************************/
void print_to_file(FILE *fp,double val,char *var_label,int first)
{
  if (first){ 
              fprintf(fp,"%s=%9.6f  ",var_label,val); 
  }
  else{       fprintf(fp,"%9.6f  ",val); }
  return;
}
/**************************************************************************************/
void print_to_file_comp(FILE *fp,int icomp,double val,char *var_label,int first)
{
  if (first){
     fprintf(fp,"%s[%d]=%9.6f  ",var_label,icomp,val); return;
  }
  else{
     fprintf(fp,"%9.6f  ",val); return;
  }
  return;
}
/**************************************************************************************/
