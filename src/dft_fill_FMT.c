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
 *  FILE: dft_fill_FMT.c
 *
 *  This file contains the fill for fundamental measures theory DFTs
 *
 */

#include "dft_fill_FMT.h"

/**********************************************************************/
/* load_nonlocal_hs_rosen_rb: Here we load all the dphi_drb terms for the 
                        Rosenfeld functional.                         */

double load_nonlocal_hs_rosen_rb(int sten_type, int iunk, int loc_inode, 
                       int inode_box, int icomp, int izone, int *ijk_box, 
                       double **x, struct RB_Struct *dphi_drb,
                       int resid_only_flag)
{
  int   **sten_offset, *offset, isten;
  int   **sten_offsetJ, *offsetJ;
  double *sten_weightJ,weightJ;
  double *sten_weight,  weight;
  struct Stencil_Struct *sten;
  struct Stencil_Struct *stenJ;

  int jzone=0, jnode_box, idim,j_box,unk_tmp,junk,i;
  int jnode_boxJ;
  int reflect_flag[NDIM_MAX];
  double  sign[3];
  struct  RB_Struct tmp;
  double resid=0.0,mat_val,resid_sum=0.0;
  int numEntries, indexUnks[4];
  double values[4];
  double n[4+2*NDIM_MAX];
  double rho_bar[4+2*NDIM_MAX];

  for (idim=0;idim<Ndim;idim++) reflect_flag[idim]=FALSE;
  jzone = find_jzone(izone,inode_box);

  sten = &(Stencil[sten_type][izone][icomp]);
  sten_offset = sten->Offset;
  sten_weight = sten->Weight;

  stenJ = &(Stencil[sten_type][jzone][icomp]);
  sten_offsetJ = stenJ->Offset;
  sten_weightJ = stenJ->Weight;

  for (isten = 0; isten < sten->Length; isten++) {

      offset = sten_offset[isten];
      weight = sten_weight[isten];

      jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);

      if (jnode_box >= 0) {

        if (sten_type == DELTA_FN_R) {
           resid = Fac_overlap_hs[icomp]*weight*
                  (dphi_drb[jnode_box].S0*Inv_4pirsq[icomp] +
                   dphi_drb[jnode_box].S1*Inv_4pir[icomp] +
                   dphi_drb[jnode_box].S2 );


           for (idim = 0; idim<Ndim; idim++){
              sign[idim]=1.0;
              if (reflect_flag[idim]) sign[idim]=-1.0;
              resid -= Fac_overlap_hs[icomp]*sign[idim]*weight * 
                      (  dphi_drb[jnode_box].V1[idim]*Inv_4pir[icomp]
                       + dphi_drb[jnode_box].V2[idim] ) *
                      (offset[idim] * Esize_x[idim]*Inv_rad[icomp]); 
           }
        }
        else if (sten_type == THETA_FN_R) resid = Fac_overlap_hs[icomp]*weight * dphi_drb[jnode_box].S3; 
      }
      else if (jnode_box == -1 || jnode_box == -3 || jnode_box == -4 ){
       if (jnode_box == -1) {
          if (sten_type == DELTA_FN_R) resid = Fac_overlap_hs[icomp]*weight*
                                      (Dphi_Drhobar_b[0]*Inv_4pirsq[icomp] +
                                       Dphi_Drhobar_b[1]*Inv_4pir[icomp] +
                                       Dphi_Drhobar_b[2] );
          else if (sten_type == THETA_FN_R) resid = Fac_overlap_hs[icomp]*weight*Dphi_Drhobar_b[3];
       }
       else if (jnode_box == -3) {
          if (sten_type == DELTA_FN_R) resid = Fac_overlap_hs[icomp]*weight*
                                      (Dphi_Drhobar_LBB[0]*Inv_4pirsq[icomp] +
                                       Dphi_Drhobar_LBB[1]*Inv_4pir[icomp] +
                                       Dphi_Drhobar_LBB[2] );
          else if (sten_type == THETA_FN_R) resid = Fac_overlap_hs[icomp]*weight*Dphi_Drhobar_LBB[3];
       }
       else if (jnode_box == -4) {
          if (sten_type == DELTA_FN_R) resid = Fac_overlap_hs[icomp]*weight*
                                       (Dphi_Drhobar_RTF[0]*Inv_4pirsq[icomp] +
                                        Dphi_Drhobar_RTF[1]*Inv_4pir[icomp] +
                                        Dphi_Drhobar_RTF[2] );
          else if (sten_type == THETA_FN_R) resid = Fac_overlap_hs[icomp]*weight*Dphi_Drhobar_RTF[3];
       }
      }
      else if (jnode_box==-2){ /* in the wall */
        resid=0.0;
      }
      if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
      resid_sum+=resid;
  
    if (resid_only_flag==FALSE)
    if (isten < stenJ->Length){
        if (jzone != izone){
           offsetJ = sten_offsetJ[isten];
    	     weightJ = sten_weightJ[isten];
           jnode_boxJ = offset_to_node_box(ijk_box, offsetJ, reflect_flag);
        }
        else{
            offsetJ = offset;
            weightJ = weight;
            jnode_boxJ = jnode_box;
        }
        if (jnode_boxJ >=0){
            junk=Phys2Unk_first[HSRHOBAR];

            for (idim = 0; idim<Ndim; idim++){
               if (reflect_flag[idim] == FALSE) sign[idim] = 1.0;
               else sign[idim] = -1.0;
            }

            for (i=0;i<Nrho_bar_s+2*Ndim;i++) rho_bar[i]=x[junk+i][jnode_boxJ];
            solutionVec_to_nOrdering(rho_bar,n);

           if (sten_type == DELTA_FN_R)       tmp = FMT2ndDerivDelta_switch(n,offsetJ,sign,icomp);
           else if (sten_type == THETA_FN_R)  tmp = FMT2ndDerivTheta_switch(n);

            numEntries=4;
            values[0]=Fac_overlap_hs[icomp]*weightJ*tmp.S3; values[1]=Fac_overlap_hs[icomp]*weightJ*tmp.S2; 
            values[2]=Fac_overlap_hs[icomp]*weightJ*tmp.S1; values[3]=Fac_overlap_hs[icomp]*weightJ*tmp.S0;

            indexUnks[0]=junk; indexUnks[1]=junk+1; indexUnks[2]=junk+2; indexUnks[3]=junk+3;
            dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                   indexUnks, jnode_boxJ, values, numEntries);

            for (idim = 0; idim<Ndim; idim++){
               numEntries=2;
               values[0]=Fac_overlap_hs[icomp]*weightJ*tmp.V2[idim]; 
               values[1]=Fac_overlap_hs[icomp]*weightJ*tmp.V1[idim];
               indexUnks[0]=junk+Nrho_bar_s+idim; 
               indexUnks[1]=indexUnks[0]+Ndim; 
               dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                    indexUnks, jnode_boxJ, values, numEntries);
            }
       }  
    }
  }
  return (resid_sum);
} 
/*****************************************************************************/
/* load_rho_bar_s:  Load scalar rho_bar definition equations  

   WITH JACOBIAN COARSENING .....*/

double load_rho_bar_s(int sten_type,double **x, int iunk,
                     int loc_inode, int inode_box, int izone,int *ijk_box,
                     int resid_only_flag)
{
  double resid_sum=0.0,resid,mat_val;
  int jzone_flag,junk;

  jzone_flag=FALSE;

  if (resid_only_flag != INIT_GUESS_FLAG){
     resid =-x[iunk][inode_box];
     resid_sum+=resid;
     if (resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
     if (resid_only_flag==FALSE){
        mat_val=-1.0;
        dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);
     }
  }
 
  if (iunk > Phys2Unk_first[HSRHOBAR]+1 && ((Lhard_surf && Nlists_HW == 2) ||
                                               (!Lhard_surf && Nwall>0 && Nlists_HW == 1))){
     junk=Phys2Unk_first[HSRHOBAR]+1;
     if (iunk == Phys2Unk_first[HSRHOBAR]+ 2){
        resid = x[junk][inode_box]*Inv_4pir[0];
        mat_val = Inv_4pir[0];
     }
     else{
        resid = x[junk][inode_box]*Inv_4pirsq[0];
        mat_val = Inv_4pirsq[0];
     }
     if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) 
                           dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
     if (resid_only_flag==FALSE) {
        dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,inode_box,mat_val);
     }
     resid_sum+=resid;
  }
  else{
    resid_sum+=resid_and_Jac_sten_fill_sum_Ncomp(sten_type,x,iunk,loc_inode,inode_box,izone,
                     ijk_box,resid_only_flag,jzone_flag,
                      &prefactor_rho_bar_s, &resid_rho_bar,&jac_rho_bar);
  }
  return(resid_sum);
}
/*****************************************************************************/
double prefactor_rho_bar_s(int iunk,int jcomp,int *offset)
{
  double fac;
  if      (iunk <= Phys2Unk_first[HSRHOBAR]+1) fac = 1.0;
  else if (iunk == Phys2Unk_first[HSRHOBAR]+2) fac = Inv_4pir[jcomp];
  else                                             fac = Inv_4pirsq[jcomp];
  fac *= Fac_overlap_hs[jcomp];

  return (fac);
}
/*****************************************************************************/
double resid_rho_bar(int junk,int jnode_box,double **x)
{
  int jcomp;
  double resid;


  if (Lseg_densities)  jcomp=Unk2Comp[junk-Phys2Unk_first[DENSITY]];
  else                 jcomp=junk-Phys2Unk_first[DENSITY];


  if (jnode_box >=0 && !Zero_density_TF[jnode_box][jcomp]){
       resid = x[junk][jnode_box];
  }
  else if (jnode_box == -1 || jnode_box ==-3 || jnode_box == -4){
       resid = constant_boundary(junk,jnode_box);
  }
  else resid=0.0;

  return (resid);
}
/*****************************************************************************/
double jac_rho_bar(int junk,int jnode_box,double **x)
{
  double jac;

  jac = 1.0;
  return (jac);
}
/*****************************************************************************/
/* load_rho_bar_v:  Load vector rho_bar definition equations  
      WITH JACOBIAN COARSENING ....*/

double load_rho_bar_v(double **x,int iunk, int loc_inode,int inode_box,
                    int izone,int *ijk_box, 
                    int resid_only_flag)
{
  double resid,resid_sum=0.0,mat_val;
  int junk,jzone_flag,idim;

  jzone_flag=FALSE;

  if (resid_only_flag != INIT_GUESS_FLAG){
     resid =-x[iunk][inode_box]; 
     resid_sum+=resid;
     if (resid_only_flag !=CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
     if (resid_only_flag==FALSE){
       mat_val=-1.0;
       dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);
     }
  }

  if (iunk >= Phys2Unk_first[HSRHOBAR]+Nrho_bar_s+Ndim && (
                                (Lhard_surf && Nlists_HW == 2) ||
                                (!Lhard_surf && Nwall>0 && Nlists_HW == 1))){
     idim = iunk - Phys2Unk_first[HSRHOBAR] - Nrho_bar_s - Ndim;
     junk = Phys2Unk_first[HSRHOBAR]+Nrho_bar_s+idim;

     resid = x[junk][inode_box]*Inv_4pir[0];
     resid_sum+=resid;
     if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) 
               dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
     if (!resid_only_flag){
        mat_val = Inv_4pir[0];
        dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,inode_box,mat_val);
     }
  }
  else { 
    resid_sum+=resid_and_Jac_sten_fill_sum_Ncomp(DELTA_FN_R,x,iunk,loc_inode,inode_box,izone,
                     ijk_box,resid_only_flag,jzone_flag,
                      &prefactor_rho_bar_v, &resid_rho_bar,&jac_rho_bar);
  }

  return(resid_sum);
}
/*****************************************************************************/
double prefactor_rho_bar_v(int iunk,int jcomp,int *offset)
{
  double fac,vector[3];
  int idim;

  if (iunk < Phys2Unk_first[HSRHOBAR] + Nrho_bar_s + Ndim) fac = 1.0;
  else                                                         fac = Inv_4pir[jcomp];

  if (iunk < Phys2Unk_first[HSRHOBAR]+Nrho_bar_s+Ndim)
          idim = iunk - Phys2Unk_first[HSRHOBAR] - Nrho_bar_s;
  else    idim = iunk - Phys2Unk_first[HSRHOBAR] - Nrho_bar_s - Ndim;

  vector[idim] = (double)offset[idim]*Esize_x[idim]*Inv_rad[jcomp];

  fac *= Fac_overlap_hs[jcomp]*vector[idim];

  return (fac);
}
/*****************************************************************************/
void calc_FMT_derivatives(void(*fp_FMTderiv)(double *,double,double,double *,double *),
                     int inode_box,double **x, struct RB_Struct *dphi_drb)
{
  double n[4+2*NDIM_MAX], rho_bar[4+2*NDIM_MAX];
  double inv_n3[5],dphi_drb_loc[4+2*NDIM_MAX];
  double DOT_22,DOT_12;
  int iunk,idim,i;

  inv_n3[0]=inv_n3[1]=inv_n3[2]=inv_n3[3]=inv_n3[4]=0.0;

  for (i=0;i<Nrho_bar_s+2*Ndim;i++) rho_bar[i]=x[Phys2Unk_first[HSRHOBAR]+i][inode_box];
  solutionVec_to_nOrdering(rho_bar,n);

  inv_n3[0]= (1.0 - n[3]);
  inv_n3[1] = 1.0 / inv_n3[0];
  inv_n3[2] = inv_n3[1]*inv_n3[1];
  inv_n3[3] = inv_n3[2]*inv_n3[1];
  inv_n3[4] = inv_n3[3]*inv_n3[1];

  DOT_22 = 0.0;
  DOT_12 = 0.0;
  for (idim = 0; idim < Ndim; idim++) {
      DOT_22 += n[Nrho_bar_s+Ndim+idim] * n[Nrho_bar_s+Ndim+idim];
      DOT_12 += n[Nrho_bar_s+idim] * n[Nrho_bar_s+Ndim+idim];
  }

  (*fp_FMTderiv)(n,DOT_12,DOT_22,inv_n3,dphi_drb_loc);
  dphi_drb[inode_box].S0=dphi_drb_loc[0];
  dphi_drb[inode_box].S1=dphi_drb_loc[1];
  dphi_drb[inode_box].S2=dphi_drb_loc[2];
  dphi_drb[inode_box].S3=dphi_drb_loc[3];
  for (idim=0;idim<Ndim;idim++){
        dphi_drb[inode_box].V1[idim]=dphi_drb_loc[4+idim];
        dphi_drb[inode_box].V2[idim]=dphi_drb_loc[4+Ndim+idim];
  }
  return;
}
/*******************************************************************************************/
