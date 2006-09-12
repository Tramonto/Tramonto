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
 *  FILE: dft_fill_EL_HSRHOBAR.c
 *
 *  This file implements the fill of the hard sphere terms in the 
 *  Euler-Lagrange Equations
 *
 */

#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"
#include "HSpkg.h"

double load_nonlocal_hs_rosen_rb(int, int, int,int,int, int,
                       int *,double **, struct RB_Struct *, int);

/**********************************************************************/
int fill_el_hsrhobar(int loc_inode,int iunk, int icomp, int inode_box, int izone, int *ijk_box,
                     int resid_only_flag,double **x,double resid)
{
    double resid_hs1,resid_hs2;

    resid_hs1 =load_nonlocal_hs_rosen_rb(DELTA_FN,iunk,loc_inode,inode_box,
                         icomp,izone, ijk_box,x,dphi_drb, resid_only_flag);
 
    resid_hs2=load_nonlocal_hs_rosen_rb(THETA_FN,iunk,loc_inode,inode_box,
                         icomp,izone, ijk_box,x,dphi_drb, resid_only_flag);
 
    *resid=resid_hs1+resid_hs2;

    return(FILL_BLOCK_FLAG);
}
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
  
  int jzone=0, jnode_box, idim,j_box,unk_tmp,junk;
  int jnode_boxJ;
  int reflect_flag[NDIM_MAX];
  double  sign[3];
  struct  RB_Struct tmp;
  double resid,mat_val,resid_sum=0.0;
  int numEntries, indexUnks[4];
  double values[4];

  
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
        
        if (sten_type == DELTA_FN) {
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
        else if (sten_type == THETA_FN) resid = Fac_overlap_hs[icomp]*weight * dphi_drb[jnode_box].S3;
      }
      else if (jnode_box == -1 || jnode_box == -3 || jnode_box == -4 ){
       if (jnode_box == -1) {
          if (sten_type == DELTA_FN) resid = Fac_overlap_hs[icomp]*weight* 
                                      (Dphi_Drhobar_b[0]*Inv_4pirsq[icomp] +
                                       Dphi_Drhobar_b[1]*Inv_4pir[icomp] +
                                       Dphi_Drhobar_b[2] );
          else if (sten_type == THETA_FN) resid = Fac_overlap_hs[icomp]*weight*Dphi_Drhobar_b[3];
       }
       else if (jnode_box == -3) {
          if (sten_type == DELTA_FN) resid = Fac_overlap_hs[icomp]*weight*
                                      (Dphi_Drhobar_LBB[0]*Inv_4pirsq[icomp] +
                                       Dphi_Drhobar_LBB[1]*Inv_4pir[icomp] +
                                       Dphi_Drhobar_LBB[2] );
          else if (sten_type == THETA_FN) resid = Fac_overlap_hs[icomp]*weight*Dphi_Drhobar_LBB[3];
       }
       }
       else if (jnode_box == -4) {
          if (sten_type == DELTA_FN) resid = Fac_overlap_hs[icomp]*weight*
                                       (Dphi_Drhobar_RTF[0]*Inv_4pirsq[icomp] +
                                        Dphi_Drhobar_RTF[1]*Inv_4pir[icomp] +
                                        Dphi_Drhobar_RTF[2] );
          else if (sten_type == THETA_FN) resid = Fac_overlap_hs[icomp]*weight*Dphi_Drhobar_RTF[3];
       }
      }
      else if (jnode_box==-2){ /* in the wall */
        resid=0.0;
      }
      dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
      resid_sum+=resid;
    
    if (!resid_only_flag)
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
            
            if (sten_type == DELTA_FN) {
               if (Type_func == FMT1) tmp = 
                          d2phi_drb2_delta_rb_FMT1(junk,jnode_boxJ,x,weightJ,offsetJ,
                          sign,Inv_rad[icomp],Inv_4pir[icomp],
                          Inv_4pirsq[icomp]);
               
               else if (Type_func ==FMT2)      tmp = 
                          d2phi_drb2_delta_rb_FMT2(junk,jnode_boxJ,x,weightJ,offsetJ,
                          sign,Inv_rad[icomp],Inv_4pir[icomp],
                          Inv_4pirsq[icomp]);
               else if (Type_func==FMT3) tmp=
                          d2phi_drb2_delta_rb_FMT1(junk,jnode_boxJ,x,weightJ,offsetJ,
                          sign,Inv_rad[icomp],Inv_4pir[icomp],
                          Inv_4pirsq[icomp]);
            }
            else if (sten_type == THETA_FN) {
               if (Type_func == FMT1) 
                       tmp = d2phi_drb2_theta_rb_FMT1(junk,jnode_boxJ,x,weightJ,offsetJ);
               else if (Type_func==FMT2)
                       tmp = d2phi_drb2_theta_rb_FMT2(junk,jnode_boxJ,x,weightJ,offsetJ);
               else if (Type_func==FMT3)
                       tmp = d2phi_drb2_theta_rb_FMT3(junk,jnode_boxJ,x,weightJ,offsetJ);
            }
            numEntries=4;            values[0]=Fac_overlap_hs[icomp]*tmp.S3; values[1]=Fac_overlap_hs[icomp]*tmp.S2; values[2]=Fac_overlap_hs[icomp]*tmp.S1; values[3]=Fac_overlap_hs[icomp]*tmp.S0;
            indexUnks[0]=junk; indexUnks[1]=junk+1; indexUnks[2]=junk+2; indexUnks[3]=junk+3;
            dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                   indexUnks, jnode_boxJ, values, numEntries);

            for (idim = 0; idim<Ndim; idim++){
               numEntries=2;
               values[0]=tmp.V2[idim]; values[1]=tmp.V1[idim];
               indexUnks[0]=junk+Nrho_bar_s+idim; indexUnks[1]=indexUnks[0]+Ndim;
               dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                    indexUnks, jnode_boxJ, values, numEntries);
            }
       }
    }
  }
  return (resid_sum);
}
/*****************************************************************************/

