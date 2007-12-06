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
 *  FILE: dft_fill_wjdc.c
 *
 *  This file contains the matrix fill routine associated with the WJDC
 *  physics functionals.
 */

#include "dft_fill_wjdc.h"

/****************************************************************************/
double load_WJDC_field(int iunk, int loc_inode, int inode_box, int *ijk_box, 
                       int izone, double **x,struct  RB_Struct *dphi_drb,
                       int mesh_coarsen_flag_i,int resid_only_flag)
{
    double resid_EL;

    /* note that the call to load_euler_lagrange is found here because the
       euler-lagrange equation is packaged in the WJDC_FIELD variable rather
       than the DENSITY variable for this theory.  However, since the field variable
       we operate on is x[iunk][inode_box]=Xi=exp(D-Vext), it is quite convenient to
       map this to the usual EL equation where rho=exp(D-Vext+mu).  Just note that in this
       functional there is no chemical potential term - rather, that term is found
       in the expression for the density which is very similar to the CMS density expression. 
       see Jain, Dominik, and Chapman, 2007*/

    resid_EL=load_euler_lagrange(iunk,loc_inode,inode_box,ijk_box,izone,
                        x,dphi_drb,mesh_coarsen_flag_i,resid_only_flag);
    return(resid_EL);
}
/****************************************************************************/
double load_WJDC_density(int iunk, int loc_inode, int inode_box, double **x,int resid_only_flag)
{
   int itype_mer,unk_B,unkIndex[2],numEntries,iseg;
   double resid_R,resid,values[2];

   resid_R=0.0;

   iseg = iunk-Phys2Unk_first[DENSITY];
   itype_mer=Unk2Comp[iseg]; /* note that itype_mer is also known as icomp */

   if (Zero_density_TF[inode_box][itype_mer] || Vext[loc_inode][itype_mer] == VEXT_MAX){
         resid_R=fill_zero_value(iunk,loc_inode,inode_box,x,resid_only_flag);
   }
   else{
      unk_B=Phys2Unk_first[WJDC_FIELD]+iseg; /* Boltzmann factor for this segment */
      resid_R+=resid_and_Jac_ChainDensity (G_CHAIN,x,iunk,unk_B,loc_inode,inode_box,
                                           resid_only_flag, &prefactor_rho_wjdc);
   }
   return(resid_R);
}                                  
/****************************************************************************/
double prefactor_rho_wjdc(int iunk)
{
  int iseg,pol_number,jseg;
  double mu,fac;

  iseg = iunk-Phys2Unk_first[DENSITY];
  for (pol_number=0; pol_number<Npol_comp; ++pol_number){
     for (jseg=0;jseg<Nmer[pol_number];jseg++){
         if (SegChain2SegAll[pol_number][jseg]==iseg) mu=Betamu_chain[pol_number];
     }
  }
  fac=exp(mu);

  return (fac);
}
/****************************************************************************/
double load_WJDC_Geqns(int iunk, int loc_inode, int inode_box, int *ijk_box, int izone, double **x,int resid_only_flag)
{
    int Njacobian_types;
    int Njacobian_sums;
    double resid_G;
    void (*funcArray_Jac[3])(int,int,int,int,int,int,int,int,int *,double,double **);
    double (*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **);
    double (*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **);

    Njacobian_types=3;   /* derivatives with respect to field and G are needed */
    Njacobian_sums=2;    /* a sum(integral) must be computed for the field derivatives */

    funcArray_Jac[0]=&WJDC_Jacobian_GCHAIN_derivG;
    funcArray_Jac[1]=&WJDC_Jacobian_GCHAIN_derivFIELD;
    funcArray_Jac[2]=&WJDC_Jacobian_GCHAIN_derivCAVITY;
    fp_ResidG=&WJDC_Resid_GCHAIN;
    fp_ResidG_Bulk=&WJDC_Resid_Bulk_GCHAIN;
    
    resid_G=load_Chain_Geqns(WJDC_FIELD,Njacobian_types,Njacobian_sums,
                             funcArray_Jac,fp_ResidG,fp_ResidG_Bulk,
                             iunk,loc_inode,inode_box,
                             ijk_box,izone,x, resid_only_flag);

    return(resid_G);
}                                           
/****************************************************************************/
void WJDC_Jacobian_GCHAIN_derivG(int iunk,int loc_inode,int pol_num,int jseg,int unk_B,
                                int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x)
{
    int i,j,icomp,jcomp;
    double prefac,power,fac,mat_val;
    prefac= -1.0;
    power = -(Nbond[pol_num][jseg]-2);

    icomp=Unk2Comp[unk_B-Phys2Unk_first[WJDC_FIELD]];
    jcomp=Unk2Comp[jseg];

    for (i=0;i < nunk-1; i++){
       fac=weight*yterm_wjdc(icomp,jcomp,jnode_box,x);
       for (j=0; j<nunk-1; j++){
          if (j != i)  fac *= x[unk[j]][jnode_box];  /*Gs or Qs*/
       }
       mat_val = fac*prefac*POW_DOUBLE_INT(x[unk[nunk-1]][jnode_box],power);
       dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk[i],jnode_box,mat_val);
    }
}
/****************************************************************************/
void WJDC_Jacobian_GCHAIN_derivFIELD(int iunk,int loc_inode,int pol_num,int jseg,int unk_B,
                    int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x)
{
    int i,icomp,jcomp;
    double prefac,power,fac,mat_val;

    prefac = (Nbond[pol_num][jseg]-2);
    power = -(Nbond[pol_num][jseg]-1);

    icomp=Unk2Comp[unk_B-Phys2Unk_first[WJDC_FIELD]];
    jcomp=Unk2Comp[jseg];
    fac=weight*yterm_wjdc(icomp,jcomp,jnode_box,x);

    for(i=0;i<nunk-1;i++) fac *=x[unk[i]][jnode_box];  /*Gs or Qs*/
    mat_val = fac*prefac*POW_DOUBLE_INT(x[unk[nunk-1]][jnode_box],power); /* Boltz Field Term */
    dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk[nunk-1],jnode_box,mat_val);
    return;
}
/****************************************************************************/
void WJDC_Jacobian_GCHAIN_derivCAVITY(int iunk,int loc_inode,int pol_num,int jseg,int unk_B,
                    int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x)
{
    int i,icomp,jcomp,unk_xi;
    double prefac_R,power_R,fac,yterm,xi_2,xi_3,y,dydxi,mat_val;

    prefac_R = -1.0;
    power_R = -(Nbond[pol_num][jseg]-2);

    icomp=Unk2Comp[unk_B-Phys2Unk_first[WJDC_FIELD]];
    jcomp=Unk2Comp[jseg];

    yterm=0.5/yterm_wjdc(icomp,jcomp,jnode_box,x);

    for (unk_xi=Phys2Unk_first[CAVWTC];unk_xi<Phys2Unk_first[CAVWTC]+2;unk_xi++){
                       /* dydxi derivatives to be computed at jnode box */
        xi_2=x[Phys2Unk_first[CAVWTC]][jnode_box];
        xi_3=x[Phys2Unk_first[CAVWTC]+1][jnode_box];
        if (unk_xi==Phys2Unk_first[CAVWTC]){
             dydxi=dy_dxi2_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],xi_2,xi_3);
             /*dydxi=dy_dxi2_cav(Bond_ff[icomp][icomp],Bond_ff[jcomp][jcomp],xi_2,xi_3);*/
        }
        else {
             dydxi=dy_dxi3_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],xi_2,xi_3);
             /*dydxi=dy_dxi3_cav(Bond_ff[icomp][icomp],Bond_ff[jcomp][jcomp],xi_2,xi_3);*/
        }
        fac=weight;
        for(i=0;i<nunk-1;i++) fac *=x[unk[i]][jnode_box];  /*Gs or Qs*/
        mat_val=fac*prefac_R*POW_DOUBLE_INT(x[unk[nunk-1]][jnode_box],power_R)*(yterm*dydxi);
        dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_xi,jnode_box,mat_val);
    }           
    return;
}
/****************************************************************************/
double WJDC_Resid_GCHAIN(int iunk,int pol_num,int jseg,int unk_B,
          int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x)
{
   int i,icomp,jcomp;
   double prefac_R,power_R,fac,resid;

   prefac_R = -1.0;
   power_R = -(Nbond[pol_num][jseg]-2);

   icomp=Unk2Comp[unk_B-Phys2Unk_first[WJDC_FIELD]];
   jcomp=Unk2Comp[jseg];
   fac=weight*yterm_wjdc(icomp,jcomp,jnode_box,x);
   for(i=0;i<nunk-1;i++) fac *=x[unk[i]][jnode_box];  /*Gs or Qs*/
   resid = fac*prefac_R*POW_DOUBLE_INT(x[unk[nunk-1]][jnode_box],power_R); /* Boltz Term */
   return(resid);
}
/****************************************************************************/
double WJDC_Resid_Bulk_GCHAIN(int iunk,int pol_num,int jseg,int unk_B,
          int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x)
{
   int i,icomp,jcomp;
   double prefac_R,power_R,fac,resid;

   prefac_R = -1.0;
   power_R = -(Nbond[pol_num][jseg]-2);

   icomp=Unk2Comp[unk_B-Phys2Unk_first[WJDC_FIELD]];
   jcomp=Unk2Comp[jseg];
   fac=weight*yterm_wjdc(icomp,jcomp,jnode_box,x);
   for(i=0;i<nunk-1;i++){
       fac *=constant_boundary(unk[i],jnode_box);  /*Gs or Qs*/
   }
   resid = fac*prefac_R*POW_DOUBLE_INT(constant_boundary(unk[nunk-1],jnode_box),power_R); /* Boltz Term */
   return(resid);
}
/****************************************************************************/
double yterm_wjdc(int icomp, int jcomp,int jnode_box,double **x)
{
  int unk_xi2,unk_xi3;
  double xi_2,xi_3,y1,y2,y,term;

     unk_xi2=Phys2Unk_first[CAVWTC];
     unk_xi3=Phys2Unk_first[CAVWTC]+1;

     if (jnode_box <0){
        xi_2=constant_boundary(unk_xi2,jnode_box);
        xi_3=constant_boundary(unk_xi3,jnode_box);
     }
     else{
        xi_2=x[unk_xi2][jnode_box];
        xi_3=x[unk_xi3][jnode_box];
     }
     term=sqrt(y_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],xi_2,xi_3));
     /*term=sqrt(y_cav(Bond_ff[icomp][icomp],Bond_ff[jcomp][jcomp],xi_2j,xi_3j));*/

     return(term);
}
/****************************************************************************/
