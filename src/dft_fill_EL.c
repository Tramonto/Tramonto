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
 *  FILE: dft_fill_EL.c
 *
 *  This file contains the fill of the residual equations and corresponding
 *  Jacobian  matrix entries for the euler-lagrange equation.
 */

#include "dft_fill_EL.h"

/******************************************************************************************/
double load_euler_lagrange(int iunk,int loc_inode, int inode_box, int *ijk_box, int izone,
                    double **x,struct  RB_Struct *dphi_drb,int mesh_coarsen_flag_i, int resid_only_flag)
{
   int i,iseg,icomp,zero_TF,bulk_TF,sym_WTC_TF,iunk_att,first_unk,offset[3],reflect_flag[3],idim,jnode_box;
   double resid=0.0,resid_att,mat_val,resid_charge;


                  /* set icomp and iseg(WTC) */
   iseg=-1;
   if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3) { /* note that this routine fills the _field_variable rather than
                                                                   the density variable for a WJDC DFT problem */
             i=iunk-Phys2Unk_first[WJDC_FIELD];               
   } 
   else      i = iunk-Phys2Unk_first[DENSITY];

   if (Type_poly==WTC || Type_poly==WJDC){    /*be careful when treating WTC densities or WJDC fields on segment basis */
                iseg=i;
                icomp=Unk2Comp[iseg];
   }
   else         icomp=i;


                   /* note that there are several cases where the euler-lagrange fill
                      is preempted by a simpler residual equation.  We check for each
                      of these cases first.  If none are true, the full EL residual is filled */

   zero_TF=check_zero_density_EL(iunk,icomp,iseg,loc_inode,inode_box,x);
   if (zero_TF) {
       resid=fill_zero_value(iunk,loc_inode,inode_box,x,resid_only_flag);
       return(resid);
   }



   bulk_TF=FALSE;
   if (mesh_coarsen_flag_i == FLAG_BULK) bulk_TF=TRUE;
   if (bulk_TF){
       if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){
          resid=fill_bulk_field(iunk,icomp,iseg,loc_inode,inode_box,x,resid_only_flag);
       }
       else{
          resid=fill_bulk_density(iunk,icomp,iseg,loc_inode,inode_box,x,resid_only_flag);
       }
       if (resid_only_flag==INIT_GUESS_FLAG) return(-resid);
       else                                  return(resid);
   } 

   sym_WTC_TF=FALSE;
   if (Type_poly==WTC && Pol_Sym_Seg[iseg] != -1) sym_WTC_TF=TRUE;
   if (sym_WTC_TF){
       resid=fill_sym_WTC(iunk,iseg,loc_inode,inode_box,x,resid_only_flag);
       return(resid);
   }

/* pin one density (or field) at this point at mean between left and right interface*/
   first_unk=FALSE;
   if (Lconstrain_interface && Type_poly!=WJDC3 && Type_poly!=WJDC && Type_poly !=WJDC2){
      if (iunk==Phys2Unk_first[DENSITY]) first_unk=TRUE;
      if (Type_interface==PHASE_INTERFACE && B2G_node[inode_box]==(int)(0.5*Size_x[Grad_dim]/Esize_x[Grad_dim]) && first_unk){  
        if (Type_poly != WJDC && Type_poly != WJDC2 && Type_poly != WJDC3){
           resid=fill_constant_density(iunk,icomp,iseg,loc_inode,inode_box,x,resid_only_flag);
           return(resid);
        }
      }
      if (Type_interface==UNIFORM_INTERFACE && B2G_node[inode_box]==(int)(0.5*Size_x[Grad_dim]/Esize_x[Grad_dim]) && first_unk){  
            for (idim=0;idim<Ndim;idim++){
               if (idim==Grad_dim) offset[idim]=1;
               else                offset[idim]=0;
               reflect_flag[idim]=FALSE;
            }
            jnode_box=offset_to_node_box(ijk_box, offset, reflect_flag);
            resid=x[iunk][inode_box]-x[iunk][jnode_box];
            if (resid_only_flag !=CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
            if (resid_only_flag==FALSE){
                mat_val=1.0;
                if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[jnode_box]+Solver_Unk[iunk]*Nnodes]-=mat_val;
                dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);         
/*                dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,jnode_box,-mat_val);         */
            }
           return(resid);
      }
   }

   /* now fill EL physics dependent terms */ 
   resid=0.0; 
   resid+=fill_EL_ideal_gas(iunk,icomp,loc_inode,inode_box,x,resid_only_flag);
   if (Type_poly != WJDC && Type_poly != WJDC2 && Type_poly != WJDC3){
      resid+=fill_EL_chem_pot(iunk,icomp,iseg,loc_inode,inode_box,mesh_coarsen_flag_i,x,resid_only_flag);
   }

   resid+=fill_EL_ext_field(iunk,icomp,loc_inode,resid_only_flag);

   if (Type_coul != NONE){
         resid_charge=fill_EL_elec_field(iunk,icomp,loc_inode,inode_box,x,resid_only_flag);
         resid+=resid_charge;
   }

   if (mesh_coarsen_flag_i != FLAG_PBELEC){

   if (Type_func !=NONE) {
         resid +=load_nonlocal_hs_rosen_rb(DELTA_FN_R,iunk,loc_inode,inode_box,
                              icomp,izone, ijk_box,x,dphi_drb, resid_only_flag);


         resid +=load_nonlocal_hs_rosen_rb(THETA_FN_R,iunk,loc_inode,inode_box,
                               icomp,izone, ijk_box,x,dphi_drb, resid_only_flag);
   }

   if (Type_attr !=NONE) {
         if (ATTInA22Block==FALSE){
           iunk_att=Phys2Unk_first[MF_EQ]+icomp;
           resid_att = x[iunk_att][inode_box];
           resid += resid_att;
           if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) 
                dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid_att);
           if (resid_only_flag==FALSE){
              mat_val=1.0;
              if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk_att]*Nnodes]+=mat_val;
              dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk_att,inode_box,mat_val);
           }
         }
         else{
            resid_att=load_mean_field(THETA_PAIRPOT_RCUT,iunk,loc_inode,
                                icomp,izone,ijk_box, x, resid_only_flag);
            resid+=resid_att;
         }
     }
   }

   if (Type_coul==DELTAC_RPM) {   /* load electrostatics deltac correlations - RPM for now*/
         resid+=load_mean_field(THETA_CR_RPM_MSA,iunk,loc_inode,
                          icomp,izone,ijk_box,x, resid_only_flag);
   }
   else if (Type_coul==DELTAC_GENERAL) {   /* load electrostatics deltac correlations - RPM for now*/
         resid+=load_mean_field(THETA_CR_GENERAL_MSA,iunk,loc_inode,
                          icomp,izone,ijk_box,x, resid_only_flag);
   }


   if (Type_poly==WTC || Type_poly==WJDC || Type_poly==WJDC2){
       if (Type_poly==WTC){
       resid+=load_polyTC_diagEL(iunk,loc_inode,inode_box,icomp,
                                 izone,ijk_box,x,resid_only_flag);
       resid+=load_polyTC_bondEL(iunk,loc_inode,inode_box,icomp,
                                 izone,ijk_box,x,resid_only_flag);
       }
       resid+=load_polyTC_cavityEL(iunk,loc_inode,inode_box,icomp,
                                  izone,ijk_box,x,resid_only_flag);
   }
   else if (Type_poly==WJDC3){
       resid+=load_polyWJDC_cavityEL(iunk,loc_inode,inode_box,icomp,
                                  izone,ijk_box,x,resid_only_flag);
   }

/*if (loc_inode==30 && iunk==10) printf("resid after cavity term=%g\n",resid);*/

   if (resid_only_flag==INIT_GUESS_FLAG) return(exp(-resid));
   else                                  return(resid);

}

/******************************************************************************************/
/* check_zero_density_EL: this routine identifies a zero density node so that the full
                                EL fill may be skipped for this node. */
int check_zero_density_EL(int iunk, int icomp, int iseg, int loc_inode, int inode_box, double **x) 
{
   int zero_density_bond_check,ibond,unk_bond,zero_TF;
   double n;

   /* set a flag for zero density of bond parameters only applies when WTC functionals are present.*/
   zero_density_bond_check=FALSE;
   if (Type_poly==WTC){
        for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
             unk_bond = Poly_to_Unk_SegAll[iseg][ibond];
                if (Pol_Sym[unk_bond] != -1) unk_bond=Pol_Sym[Poly_to_Unk_SegAll[iseg][ibond]];
                unk_bond += Phys2Unk_first[BONDWTC];
                n=x[unk_bond][inode_box];
                if (fabs(n)<1.e-8) zero_density_bond_check=TRUE;
            }
   }

   zero_TF=FALSE;
   if (Zero_density_TF[inode_box][icomp] || Vext[loc_inode][icomp]==VEXT_MAX || zero_density_bond_check) zero_TF=TRUE;

   return zero_TF;
}
/******************************************************************************************/
double fill_zero_value(int iunk, int loc_inode, int inode_box, double **x,int resid_only_flag)
{
   double resid=0.0,mat_val;
  
   if (resid_only_flag != INIT_GUESS_FLAG){ 
      resid= x[iunk][inode_box];
      if (resid_only_flag !=CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

      if (resid_only_flag==FALSE){
         mat_val = 1.0;
         if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
         dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);         
      }
   }
    
   return(resid);
}
/******************************************************************************************/
double fill_sym_WTC(int iunk, int iseg, int loc_inode, int inode_box, double **x,int resid_only_flag)
{
   int junk,unkIndex[2],numEntries,jtmp;
   double values[2],resid=0.0;

    junk = Pol_Sym_Seg[iseg]+Phys2Unk_first[DENSITY];
    if (resid_only_flag != INIT_GUESS_FLAG) resid = x[iunk][inode_box];
    resid -=x[junk][inode_box];
    if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) 
        dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
    if (resid_only_flag==FALSE){
       unkIndex[0]=iunk; unkIndex[1]=junk;
       values[0]=1.0; values[1]=-1.0;
       numEntries=2;
       if (Iwrite_files==FILES_DEBUG_MATRIX){
           for (jtmp=0;jtmp<2;jtmp++) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[unkIndex[jtmp]]*Nnodes]+=values[jtmp];
       }
       dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                         unkIndex, inode_box, values, numEntries); 
    }
    return(resid);
}
/******************************************************************************************/
double fill_bulk_density(int iunk, int icomp, int iseg, int loc_inode, int inode_box, double **x,int resid_only_flag)
{
  double resid,mat_val;

  if (resid_only_flag != INIT_GUESS_FLAG){
     resid = log(x[iunk][inode_box]) ; 
     mat_val = 1.0/x[iunk][inode_box];
     if (resid_only_flag !=CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
     if (resid_only_flag==FALSE){
        if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
        dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);         
     }
   }

   if (Lseg_densities)   resid = -log(Rho_seg_b[iseg]);
   else                  resid = -log(Rho_b[icomp]);
   if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

   return resid;
}
/******************************************************************************************/
double fill_constant_density(int iunk, int icomp, int iseg, int loc_inode, int inode_box, double **x,int resid_only_flag)
{
  double resid,mat_val;

  if (resid_only_flag != INIT_GUESS_FLAG){
     resid = log(x[iunk][inode_box]) ; 
     mat_val = 1.0/x[iunk][inode_box];
     if (resid_only_flag !=CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
     if (resid_only_flag==FALSE){
        if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
        dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);         
     }
   }

   if (Lseg_densities)   resid = -log(0.5*(Rho_seg_LBB[iseg]+Rho_seg_RTF[iseg]));
   else                  resid = -log(0.5*(Rho_b_LBB[icomp]+Rho_b_RTF[icomp]));

   if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

   return resid;
}
/******************************************************************************************/
double fill_constant_field(int iunk, int icomp, int iseg, int loc_inode, int inode_box, double **x,int resid_only_flag)
{
  double resid,mat_val;

  if(resid_only_flag != INIT_GUESS_FLAG){
    resid = log(x[iunk][inode_box]) ; 
    mat_val = 1.0/x[iunk][inode_box];
    if (resid_only_flag !=CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
    if (resid_only_flag==FALSE){
       if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
       dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);         
    }
  }
  resid = -log((0.5*Field_WJDC_LBB[icomp]+0.5*Field_WJDC_RTF[icomp]));
  

  if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

  if (resid_only_flag==INIT_GUESS_FLAG) return(exp(-resid));
  else                                  return resid;
}
/******************************************************************************************/
double fill_bulk_field(int iunk, int icomp, int iseg, int loc_inode, int inode_box, double **x,int resid_only_flag)
{
  double resid,mat_val;

  if(resid_only_flag != INIT_GUESS_FLAG){
    resid = log(x[iunk][inode_box]) ; 
    mat_val = 1.0/x[iunk][inode_box];
    if (resid_only_flag !=CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
    if (resid_only_flag==FALSE){
       if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
       dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);         
    }
  }

  resid = -log(Field_WJDC_b[icomp]);
  if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

  if (resid_only_flag==INIT_GUESS_FLAG) return(exp(-resid));
  else                                  return resid;
}
/******************************************************************************************/
double fill_EL_ideal_gas(int iunk, int icomp, int loc_inode, int inode_box, double **x,int resid_only_flag)
{
   double resid,resid_ig=0.0,mat_val,scalefac;
   int pol_number;

   if (resid_only_flag != INIT_GUESS_FLAG){
      resid = log(x[iunk][inode_box]); 
      resid_ig=resid;
      if (resid_only_flag !=CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
      if (!resid_only_flag){
         mat_val = 1.0/x[iunk][inode_box];
         if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
         dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);         
      }

      if(LDeBroglie){
         resid = (- 3.0*log(Sigma_ff[icomp][icomp]) -1.5*log(Mass[icomp]*Temp));
         resid_ig += resid;
         dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
      }
   }

   if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){ 
      for (pol_number=0;pol_number<Npol_comp;pol_number++) {
       if (Nseg_type_pol[pol_number][icomp] !=0) scalefac=Scale_fac_WJDC[pol_number][icomp];
      }
      resid = -scalefac;
      resid_ig+=resid;
      if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag !=CALC_RESID_ONLY)
            dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
   }

   return(resid_ig);
}
/******************************************************************************************/
double fill_EL_chem_pot(int iunk, int icomp, int iseg, int loc_inode, int inode_box, 
                        int mesh_coarsen_flag_i, double **x,int resid_only_flag)
{
   double resid_mu,mat_val;
   int junk,usemu_test;

   usemu_test=FALSE;
   if (LBulk && Loca.cont_type1==CONT_BETAMU_I && icomp==Cont_ID[0][0]) usemu_test=TRUE;
   if (Loca.method==4 && LBulk && Loca.cont_type2==CONT_BETAMU_I && icomp==Cont_ID[1][0]) usemu_test=TRUE;


   resid_mu = 0.0;
   if (Type_interface != DIFFUSIVE_INTERFACE) {
      if (Type_interface ==UNIFORM_INTERFACE && !usemu_test){
        if (Lseg_densities) resid_mu -= log(Rho_seg_b[iseg]);
        else                resid_mu -= log(Rho_b[icomp]);

        if (mesh_coarsen_flag_i != FLAG_PBELEC){
            if (Type_func != NONE) resid_mu -= Betamu_hs_ex[icomp];
            if (Type_attr != NONE)    resid_mu -= Betamu_att[icomp];
            if (Type_poly == WTC)         resid_mu -= Betamu_wtc[iseg];
            if (Type_coul == DELTAC_RPM || Type_coul==DELTAC_GENERAL) resid_mu += Deltac_b[icomp];
        }
      }
      else if (Type_interface==PHASE_INTERFACE || usemu_test){
          if (Lseg_densities) resid_mu = -Betamu_seg[iseg];
          else {              resid_mu = -Betamu[icomp];
          }
      }
   }
   else if (Type_interface==DIFFUSIVE_INTERFACE){     
      if(Lseg_densities)  junk=Phys2Unk_first[DIFFUSION] + iseg;
      else                junk=Phys2Unk_first[DIFFUSION] + icomp;
      resid_mu = -x[junk][inode_box];
      if (resid_only_flag==FALSE){
         mat_val=-1.0;
         if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
         dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,inode_box,mat_val);
      }
   }

   if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY)
       dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid_mu);
 
   return(resid_mu);
}
/******************************************************************************************/
double fill_EL_ext_field(int iunk, int icomp, int loc_inode,int resid_only_flag)
{
   double resid_vext;
   if (Nwall > 0){
      resid_vext = Vext[loc_inode][icomp];
      if (resid_only_flag!=INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid_vext);
   }
   else resid_vext=0.0;
   return (resid_vext);
}
/******************************************************************************************/
double fill_EL_elec_field(int iunk, int icomp, int loc_inode, int inode_box, double **x,int resid_only_flag)
{
   int junk,numEntries, nodeIndices[3],jtmp;
   double resid_charge,mat_val,resid,fac_temp,gradphi,fac,values[3];

   junk = Phys2Unk_first[POISSON];
   resid_charge = Charge_f[icomp]*x[junk][inode_box];
   if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY){
      dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid_charge);
      if (!resid_only_flag){
         mat_val = Charge_f[icomp];
         if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[junk]*Nnodes]+=mat_val;
         dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,inode_box,mat_val);
      }
   }

   if (Lpolarize[icomp] && Ndim==1){

   fac_temp=Temp_elec/(4.0*PI*KAPPA_H2O);
   if (Nodes_2_boundary_wall[0][inode_box] != -1){
      if (Surf_normal[0][loc_inode][0] < 0){
          gradphi = 0.5*(3.0*x[junk][inode_box]-4.0*x[junk][inode_box-1]+x[junk][inode_box-2])/Esize_x[0];
          fac=Pol[icomp]*gradphi*fac_temp/Esize_x[0];

          numEntries=3;
          nodeIndices[0]=inode_box; nodeIndices[1]=inode_box-1; nodeIndices[2]=inode_box-2;
          values[0]=1.5*fac; values[1]=-2.0*fac; values[2]=0.5*fac;
      }
      else{
          gradphi = 0.5*(-3.0*x[junk][inode_box]+4.0*x[junk][inode_box+1]-x[junk][inode_box+2])/Esize_x[0];
          fac=Pol[icomp]*gradphi*fac_temp/Esize_x[0];

          numEntries=3;
          nodeIndices[0]=inode_box; nodeIndices[1]=inode_box+1; nodeIndices[2]=inode_box+2;
          values[0]=-1.5*fac; values[1]=2.0*fac; values[2]=-0.5*fac;
      }
   }
   else{
          gradphi = 0.5*(x[junk][inode_box+1]-x[junk][inode_box-1])/Esize_x[0];
          fac=Pol[icomp]*gradphi*fac_temp/Esize_x[0];
          numEntries=2;
          nodeIndices[0]=inode_box+1; nodeIndices[1]=inode_box-1;
          values[0]=0.5*fac; values[1]=-0.5*fac;
   }
   resid = 0.5*Pol[icomp]*gradphi*gradphi*fac_temp;
   resid_charge += resid;
   if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY){
      dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
      if (resid_only_flag==FALSE){
         if (Iwrite_files==FILES_DEBUG_MATRIX){
            for (jtmp=0;jtmp<numEntries;jtmp++) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[nodeIndices[jtmp]]+Solver_Unk[junk]*Nnodes]+=values[jtmp];
         }
         dft_linprobmgr_insertmultinodematrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                       junk,nodeIndices,values,numEntries);
      }
      }
   }

   return(resid_charge);
}
/******************************************************************************************/
