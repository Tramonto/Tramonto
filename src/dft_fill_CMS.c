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
 *  FILE: dft_fill_cms.c
 *
 *  This file contains the matrix fill routine associated with the CMS
 *  physics functionals.
 */

#include "dft_fill_CMS.h"

/****************************************************************************/
double load_CMS_field(int iunk, int loc_inode, int inode_box, int *ijk_box, int izone, double **x,int resid_only_flag)
{
   double resid_B,resid,mat_val;
   int itype_mer,junk,iunk_att;

   itype_mer = iunk - Phys2Unk_first[CMS_FIELD];

    if (Zero_density_TF[inode_box][itype_mer] || Vext[loc_inode][itype_mer] == VEXT_MAX){
         resid_B=fill_zero_value(iunk,loc_inode,inode_box,x,resid_only_flag);
    }
    else{
		/* note: mean-field part of Jacobian gets filled in resid_and_Jac_sten_fill_sum_Ncomp in dft_fill_CLSmf.c */
       if (ATTInA22Block==FALSE){
           iunk_att=Phys2Unk_first[MF_EQ]+itype_mer;
           resid_B=x[iunk_att][inode_box];
           if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY)
                dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid_B);
           if (resid_only_flag==FALSE){
              mat_val=1.0;
              dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk_att,inode_box,mat_val);
           }
       }
       else{
            resid_B = load_mean_field(THETA_CR_DATA,iunk,loc_inode,itype_mer,izone,ijk_box,x,resid_only_flag); 
       }
       resid = Vext[loc_inode][itype_mer]+log(x[iunk][inode_box]);
       resid_B+=resid;
       if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
       if(!resid_only_flag){
          mat_val = 1.0/x[iunk][inode_box];  
          dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);
       }

       /* add a mean field electrostatic contribution to the CMS field if there are charges in the problem */
       if (Type_coul != NONE){
          junk = Phys2Unk_first[POISSON];
          resid = Charge_f[itype_mer]*x[junk][inode_box];
          resid_B+=resid;
          if (resid_only_flag !=CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          if(!resid_only_flag){
               mat_val = Charge_f[itype_mer];
               dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,inode_box,mat_val);
          }
       }          
    }
    return(resid_B);
}
/****************************************************************************/
double load_CMS_density(int iunk, int loc_inode, int inode_box, double **x,int resid_only_flag)
{
   int itype_mer,unk_B;
   double resid_R;

   resid_R=0.0;

   itype_mer = iunk-Phys2Unk_first[DENSITY];

   if (Zero_density_TF[inode_box][itype_mer] || Vext[loc_inode][itype_mer] == VEXT_MAX){
         resid_R=fill_zero_value(iunk,loc_inode,inode_box,x,resid_only_flag);
   }
   else{
      unk_B=Phys2Unk_first[CMS_FIELD]+itype_mer; /* Boltzmann factor for this i */
      resid_R+=resid_and_Jac_ChainDensity (G_CHAIN,x,iunk,unk_B,loc_inode,inode_box,
                                           resid_only_flag, &prefactor_rho_cms);
   }
   return(resid_R);
}                                  
/****************************************************************************/
double prefactor_rho_cms(int itype_mer)
{
  int npol;
  double fac;

  npol = 0;
  while (Nmer_t[npol][itype_mer]==0) npol++;
	
  if(Grafted[npol])
	  fac = Rho_g[npol]/Gsum[npol];
  else
	  fac = Rho_b[itype_mer]/Nmer_t[npol][itype_mer]; 

  return (fac);
}
/****************************************************************************/
double load_CMS_Geqns(int iunk, int loc_inode, int inode_box, int *ijk_box, int izone, double **x,int resid_only_flag)
{
    int Njacobian_types;
    int Njacobian_sums;
	int field;
    double resid_G;
    void (*funcArray_Jac[3])(int,int,int,int,int,int,int,int,int *,double,double **);
    double (*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **);
    double (*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **);

    Njacobian_types=2;   /* derivatives with respect to field and G are needed */
    Njacobian_sums=1;    /* a sum(integral) must be computed for the field derivatives */

    funcArray_Jac[0]=&CMS_Jacobian_GCHAIN_derivG;
    funcArray_Jac[1]=&CMS_Jacobian_GCHAIN_derivFIELD;
    fp_ResidG=&CMS_Resid_GCHAIN;
    fp_ResidG_Bulk=&CMS_Resid_Bulk_GCHAIN;
	
	field=CMS_FIELD;

    resid_G=load_Chain_Geqns(field,Njacobian_types,Njacobian_sums,
                             funcArray_Jac,fp_ResidG,fp_ResidG_Bulk,
                             iunk, loc_inode,inode_box, 
                             ijk_box,izone,x, resid_only_flag);
    return(resid_G);
}
/****************************************************************************/
void CMS_Jacobian_GCHAIN_derivG(int iunk,int loc_inode,int pol_num,int jseg,int unk_B,
                                int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x)
{
    int i,j;
    double prefac,power,fac,mat_val;
    prefac= -1.0;   /*-x[unk_B][inode_box];*/
    power = -(Nbond[pol_num][jseg]-2);

    for (i=0;i < nunk-1; i++){
       fac=weight;              
       for (j=0; j<nunk-1; j++){                 
          if (j != i)  fac *= x[unk[j]][jnode_box];  /*Gs or Qs*/              
       }              
       mat_val = fac*prefac*POW_DOUBLE_INT(x[unk[nunk-1]][jnode_box],power); 
       dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk[i],jnode_box,mat_val);
    }
    return;
}
/****************************************************************************/
void CMS_Jacobian_GCHAIN_derivFIELD(int iunk,int loc_inode,int pol_num,int jseg,int unk_B,
                    int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x)
{
    int i;
    double prefac,power,fac,mat_val;

    prefac = (Nbond[pol_num][jseg]-2); /*x[unk_B][inode_box]*(Nbond[pol_num][jseg]-2);*/
    power = -(Nbond[pol_num][jseg]-1);

    fac=weight;
    for(i=0;i<nunk-1;i++) fac *=x[unk[i]][jnode_box];  /*Gs or Qs*/
    mat_val = fac*((double)prefac)*POW_DOUBLE_INT(x[unk[nunk-1]][jnode_box],power); /* Boltz Field Term */
    dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk[nunk-1],jnode_box,mat_val);
    return;
}
/****************************************************************************/
double CMS_Resid_GCHAIN(int iunk,int pol_num,int jseg,int unk_B,
          int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x)
{
   int i;
   double prefac_R,power_R,fac,resid;

   prefac_R = -1.0;    /*-x[unk_B][inode_box]*/;
   power_R = -(Nbond[pol_num][jseg]-2);

   fac=weight;
   for(i=0;i<nunk-1;i++) fac *=x[unk[i]][jnode_box];  /*Gs or Qs*/
   resid = fac*prefac_R*POW_DOUBLE_INT(x[unk[nunk-1]][jnode_box],power_R); /* Boltz Term */
   return(resid);
}
/****************************************************************************/
double CMS_Resid_Bulk_GCHAIN(int iunk,int pol_num,int jseg,int unk_B,
          int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x)
{
   int i;
   double prefac_R,power_R,fac,resid;

   prefac_R = -1.0;  /*x[unk_B][inode_box];*/
   power_R = -(Nbond[pol_num][jseg]-2);

   fac=weight;
   for(i=0;i<nunk-1;i++) fac *=constant_boundary(unk[i],jnode_box);  /*Gs or Qs*/
   resid = fac*prefac_R*POW_DOUBLE_INT(constant_boundary(unk[nunk-1],jnode_box),power_R); /* Boltz Term */
   return(resid);
}
/****************************************************************************/
