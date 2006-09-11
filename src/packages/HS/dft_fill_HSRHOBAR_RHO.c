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
 *  FILE: dft_fill_rosen_rb.c
 *
 *  This file contains the fill for the rho and rhobar implementation
 *  of the rosenfeld functional.
 */

#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"
/*****************************************************************************/
int fill_HSRHOBAR_DENSITY(int inode,int iunk,int resid_only_flag)
{

   if (iunk == Phys2Unk_first[HSRHOBAR]){
        resid_rhobars+=load_rho_bar_s(THETA_FN,x,iunk,loc_inode,inode_box,izone,ijk_box,
                       resid_only_flag);
   }
   else if (iunk < Phys2Unk_first[HSRHOBAR]+Nrho_bar_s){
       resid_rhobars+=load_rho_bar_s(DELTA_FN,x,iunk,loc_inode,inode_box,izone,ijk_box,
                      resid_only_flag);
   }
   else if (iunk >= Phys2Unk_first[HSRHOBAR]+Nrho_bar_s){
        resid_rhobarv+=load_rho_bar_v(x,iunk,loc_inode,inode_box,izone,ijk_box, resid_only_flag);
   }

   return (FILL_BLOCK_FLAG);
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

  if (iunk > Phys2Unk_first[RHOBAR_ROSEN]+1 && ((Lhard_surf && Nlists_HW == 2) ||
                                               (!Lhard_surf && Nwall>0 && Nlists_HW == 1))){
     junk=Phys2Unk_first[RHOBAR_ROSEN]+1;
     if (iunk == Phys2Unk_first[RHOBAR_ROSEN]+ 2){
        resid = x[junk][inode_box]*Inv_4pir[0];
        mat_val = Inv_4pir[0];
     }
     else{
        resid = x[junk][inode_box]*Inv_4pirsq[0];
        mat_val = Inv_4pirsq[0];
     }
     dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
     dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,inode_box,mat_val);
     resid_sum+=resid;
  }
  else{
    resid_and_Jac_sten_fill_sum_Ncomp(sten_type,x,iunk,loc_inode,inode_box,izone,
                     ijk_box,resid_only_flag,jzone_flag,
                      &prefactor_rho_bar_s, &resid_rho_bar,&jac_rho_bar);
  }
  resid_sum+=Temporary_sum;
  return(resid_sum);
}
/*****************************************************************************/
double prefactor_rho_bar_s(int iunk,int jcomp,int *offset)
{
  double fac;
  if      (iunk <= Phys2Unk_first[RHOBAR_ROSEN]+1) fac = 1.0;
  else if (iunk == Phys2Unk_first[RHOBAR_ROSEN]+2) fac = Inv_4pir[jcomp];
  else                                             fac = Inv_4pirsq[jcomp];
  fac *= Fac_overlap_hs[jcomp];

  return (fac);
}
/*****************************************************************************/
double resid_rho_bar(int junk,int jnode_box,double **x)
{
  int jcomp;
  double resid;


  if (Type_poly==WTC)  jcomp=Unk2Comp[junk-Phys2Unk_first[DENSITY]];
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
  int jcomp;
  double jac;

  if (Type_poly==WTC)  jcomp=Unk2Comp[junk-Phys2Unk_first[DENSITY]];
  else                 jcomp=junk-Phys2Unk_first[DENSITY];

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

  if (iunk >= Phys2Unk_first[RHOBAR_ROSEN]+Nrho_bar_s+Ndim && (
                                (Lhard_surf && Nlists_HW == 2) ||
                                (!Lhard_surf && Nwall>0 && Nlists_HW == 1))){
     idim = iunk - Phys2Unk_first[RHOBAR_ROSEN] - Nrho_bar_s - Ndim;
     junk = Phys2Unk_first[RHOBAR_ROSEN]+Nrho_bar_s+idim;

     resid = x[junk][inode_box]*Inv_4pir[0];
     resid_sum+=resid;
     mat_val = Inv_4pir[0];
     dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
     dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,inode_box,mat_val);
  }
  else { 
    resid_and_Jac_sten_fill_sum_Ncomp(DELTA_FN,x,iunk,loc_inode,inode_box,izone,
                     ijk_box,resid_only_flag,jzone_flag,
                      &prefactor_rho_bar_v, &resid_rho_bar,&jac_rho_bar);
  }

  resid_sum+=Temporary_sum;
  return(resid_sum);
}
/*****************************************************************************/
double prefactor_rho_bar_v(int iunk,int jcomp,int *offset)
{
  double fac,vector[3];
  int idim;

  if (iunk < Phys2Unk_first[RHOBAR_ROSEN] + Nrho_bar_s + Ndim) fac = 1.0;
  else                                                         fac = Inv_4pir[jcomp];

  if (iunk < Phys2Unk_first[RHOBAR_ROSEN]+Nrho_bar_s+Ndim)
          idim = iunk - Phys2Unk_first[RHOBAR_ROSEN] - Nrho_bar_s;
  else    idim = iunk - Phys2Unk_first[RHOBAR_ROSEN] - Nrho_bar_s - Ndim;

  vector[idim] = (double)offset[idim]*Esize_x[idim]*Inv_rad[jcomp];

  fac *= Fac_overlap_hs[jcomp]*vector[idim];

  return (fac);
}
/*****************************************************************************/
