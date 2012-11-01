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
   int itype_mer,unk_B,iseg;
   double resid_R;

   resid_R=0.0;

   if (Type_poly==WJDC || Type_poly==WJDC2){
      iseg = iunk-Phys2Unk_first[DENSITY];
      itype_mer=Unk2Comp[iseg]; /* note that itype_mer is also known as icomp */
   }
   else if (Type_poly==WJDC3){
      itype_mer=iunk-Phys2Unk_first[DENSITY];
   }

   if (Zero_density_TF[inode_box][itype_mer] || Vext[loc_inode][itype_mer] == VEXT_MAX){
         resid_R=fill_zero_value(iunk,loc_inode,inode_box,x,resid_only_flag);
   }
   else{
      if (Type_poly==WJDC2 || Type_poly==WJDC3) unk_B=Phys2Unk_first[WJDC_FIELD]+itype_mer; /* Boltzmann factor for this segment */
      else if (Type_poly==WJDC)                 unk_B=Phys2Unk_first[WJDC_FIELD]+iseg; 

/*      resid_R+=resid_and_Jac_ChainDensity_WJDC2 (G_CHAIN,x,iunk,unk_B,loc_inode,inode_box,
                                            resid_only_flag, &prefactor_rho_wjdc);*/

      resid_R+=resid_and_Jac_ChainDensity (G_CHAIN,x,iunk,unk_B,loc_inode,inode_box,
                                           resid_only_flag, &prefactor_rho_wjdc);
   }
   return(resid_R);
}                                  
/****************************************************************************/
double prefactor_rho_wjdc(int iseg,int inode_box,double **x)
{
  int pol_number,icomp;
  double mu,fac,scale_term;

  pol_number=SegAll_to_Poly[iseg];
  scale_term=0.0;

  if(Type_poly==WJDC3 && Grafted_Logical==TRUE && Grafted[pol_number]!=FALSE)  {
      if (iseg==Grafted_SegIDAll[pol_number]){
          if (Grafted[pol_number]==GRAFT_DENSITY) fac=Rho_g[pol_number];
          else fac=Rho_g[pol_number]/Total_area_graft[pol_number];
      }
      else{
          fac = Rho_g[pol_number]/Gsum_graft[pol_number]; 
      }
   }
   else{
      for (icomp=0;icomp<Ncomp;icomp++){
         scale_term-=Scale_fac_WJDC[pol_number][icomp]*Nseg_type_pol[pol_number][icomp];
      }

      if (Type_interface==DIFFUSIVE_INTERFACE){
         mu=x[Phys2Unk_first[DIFFUSION]+pol_number][inode_box];
      }
      else mu=Betamu_chain[pol_number];

      fac=exp(mu+scale_term);
   }

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

    if (Type_poly==WJDC2 || Type_poly==WJDC3) icomp=unk_B-Phys2Unk_first[WJDC_FIELD];
    else if (Type_poly==WJDC) icomp=Unk2Comp[unk_B-Phys2Unk_first[WJDC_FIELD]]; 
    jcomp=Type_mer[pol_num][jseg];
/*    jcomp=Unk2Comp[jseg];*/

    for (i=0;i < nunk-1; i++){
       fac=weight*yterm_wjdc(icomp,jcomp,jnode_box,x);
       for (j=0; j<nunk-1; j++){
          if (j != i)  fac *= x[unk[j]][jnode_box];  /*Gs or Qs*/
       }
       mat_val = fac*prefac*POW_DOUBLE_INT(x[unk[nunk-1]][jnode_box],power);
       if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[jnode_box]+Solver_Unk[unk[i]]*Nnodes]+=mat_val;
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

    if (Type_poly==WJDC2 || Type_poly==WJDC3)  icomp=unk_B-Phys2Unk_first[WJDC_FIELD];
    else if (Type_poly==WJDC) icomp=Unk2Comp[unk_B-Phys2Unk_first[WJDC_FIELD]]; 
    jcomp=Type_mer[pol_num][jseg];
/*    jcomp=Unk2Comp[jseg];*/
    fac=weight*yterm_wjdc(icomp,jcomp,jnode_box,x);

    for(i=0;i<nunk-1;i++) fac *=x[unk[i]][jnode_box];  /*Gs or Qs*/
    mat_val = fac*((double)prefac)*POW_DOUBLE_INT(x[unk[nunk-1]][jnode_box],power); /* Boltz Field Term */
    if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[jnode_box]+Solver_Unk[unk[nunk-1]]*Nnodes]+=mat_val;
    dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk[nunk-1],jnode_box,mat_val);
    return;
}
/****************************************************************************/
void WJDC_Jacobian_GCHAIN_derivCAVITY(int iunk,int loc_inode,int pol_num,int jseg,int unk_B,
                    int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x)
{
    int i,icomp,jcomp,unk_xi,power_R;
    double prefac_R,fac,yterm,xi_2,xi_3,dydxi,mat_val;

    prefac_R = -1.0;
    power_R = -(Nbond[pol_num][jseg]-2);

    if (Type_poly==WJDC2 || Type_poly==WJDC3)  icomp=unk_B-Phys2Unk_first[WJDC_FIELD];
    else if (Type_poly==WJDC) icomp=Unk2Comp[unk_B-Phys2Unk_first[WJDC_FIELD]]; 
    jcomp=Type_mer[pol_num][jseg];
/*    jcomp=Unk2Comp[jseg];*/

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
        if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[jnode_box]+Solver_Unk[unk_xi]*Nnodes]+=mat_val;
        dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_xi,jnode_box,mat_val);
    }           
    return;
}
/****************************************************************************/
double WJDC_Resid_GCHAIN(int iunk,int pol_num,int jseg,int unk_B,
          int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x)
{
   int i,icomp,jcomp,power_R;
   double prefac_R,fac,resid;

   prefac_R = -1.0;
   power_R = -(Nbond[pol_num][jseg]-2);

   if (Type_poly==WJDC2 || Type_poly==WJDC3) icomp=unk_B-Phys2Unk_first[WJDC_FIELD];
   else /*if (Type_poly==WJDC)*/ icomp=Unk2Comp[unk_B-Phys2Unk_first[WJDC_FIELD]]; 

   jcomp=Type_mer[pol_num][jseg];

   fac=weight*yterm_wjdc(icomp,jcomp,jnode_box,x);
   for(i=0;i<nunk-1;i++){
        fac *=x[unk[i]][jnode_box];  /*Gs or Qs*/
   }
   if (x[unk[nunk-1]][jnode_box] > 1.e-15)
        resid = fac*prefac_R*POW_DOUBLE_INT(x[unk[nunk-1]][jnode_box],power_R); /* Boltz Term */
   else resid=0.0;

   return(resid);
}
/****************************************************************************/
double WJDC_Resid_Bulk_GCHAIN(int iunk,int pol_num,int jseg,int unk_B,
          int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x)
{
   int i,icomp,jcomp,power_R;
   double prefac_R,fac,resid;

   prefac_R = -1.0;
   power_R = -(Nbond[pol_num][jseg]-2);

   if (Type_poly==WJDC2 || Type_poly==WJDC3) icomp=unk_B-Phys2Unk_first[WJDC_FIELD];
   else if (Type_poly==WJDC) icomp=Unk2Comp[unk_B-Phys2Unk_first[WJDC_FIELD]]; 

   jcomp=Type_mer[pol_num][jseg];
   fac=weight*yterm_wjdc(icomp,jcomp,jnode_box,x);
   for(i=0;i<nunk-1;i++){
       fac *=constant_boundary(unk[i],jnode_box);  /*Gs or Qs*/
   }
   if (constant_boundary(unk[nunk-1],jnode_box) > 1.e-15) 
        resid = fac*prefac_R*POW_DOUBLE_INT(constant_boundary(unk[nunk-1],jnode_box),power_R); /* Boltz Term */
   else resid=0.0;
   return(resid);
}
/****************************************************************************/
double yterm_wjdc(int icomp, int jcomp,int jnode_box,double **x)
{
  int unk_xi2,unk_xi3;
  double xi_2,xi_3,term;

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

     return(term);
}
/****************************************************************************/
/*load_polyWJDC_cavityEL:  Here add the bond contributions of the WTC association/bonding (cavity term) functionals                       to the Euler-Lagrange equation for the Wertheim-Tripathi-Chapman theory */
double load_polyWJDC_cavityEL(int iunk,int loc_inode,int inode_box,int icomp,int izone,int *ijk_box,double **x,int resid_only_flag)
{
  int jseg,kseg,kbond,jcomp,unk_xi2,unk_xi3,unk_rho,kcomp,unk_B,unk_GQ,jseg_tmp,i;
  double xi_2,xi_3,s1,s2,y,dy_dxi2,dy_dxi3,prefac2,prefac3;
  double d2y_dxi2_2,d2y_dxi3_2,d2y_dxi2_dxi3;
  int jzone,jnode_box,jlist,reflect_flag[3],jnode_boxJ;
  int   **sten_offset, *offset, isten;
  int   **sten_offsetJ, *offsetJ;
  int npol,jbond,iwall,jsurf_node,ipol;
  double *sten_weightJ,weightJ;
  double *sten_weight,  weight;
  struct Stencil_Struct *sten;
  struct Stencil_Struct *stenJ;
  double resid,mat_val,resid_sum,first_deriv,first_deriv_sum,dens,dens_Gderiv,resid_sum_tethered[NCOMP_MAX];

/* int iseg;
   if (Type_poly==WJDC) iseg = iunk-Phys2Unk_first[WJDC_FIELD];
  else                 iseg = iunk - Phys2Unk_first[DENSITY];*/

  resid_sum=resid=0.0;
  for (i=0;i<Npol_comp;i++) resid_sum_tethered[i]=0.0;
  unk_xi2 = Phys2Unk_first[CAVWTC];
  unk_xi3 = Phys2Unk_first[CAVWTC]+1;
  jzone=izone;

  sten = &(Stencil[THETA_FN_SIG][izone][icomp]);
  sten_offset = sten->Offset;
  sten_weight = sten->Weight;

  stenJ = &(Stencil[THETA_FN_SIG][jzone][icomp]);
  sten_offsetJ = stenJ->Offset;
  sten_weightJ = stenJ->Weight;

  npol = 0;
  while (Nmer_t[npol][icomp]==0){ npol++; }


  for (isten = 0; isten < sten->Length; isten++) {
    offset = sten_offset[isten];
    weight = sten_weight[isten];
    jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);

    if (jnode_box >= 0) {
         xi_2=x[unk_xi2][jnode_box]; xi_3=x[unk_xi3][jnode_box];
    }
    else{
         xi_2=constant_boundary(unk_xi2,jnode_box);
         xi_3=constant_boundary(unk_xi3,jnode_box);
    }



    resid=0.0;
    for (jseg=0;jseg<Nseg_tot;jseg++){
      
       jseg_tmp=jseg; 
       if (Pol_Sym_Seg[jseg] != -1) jseg_tmp=Pol_Sym_Seg[jseg];
       unk_rho = Unk2Comp[jseg_tmp]+Phys2Unk_first[DENSITY];
       jcomp=Unk2Comp[jseg];

       s1=Sigma_ff[jcomp][jcomp];
/*       s1=Bond_ff[jcomp][jcomp];*/
       if (Nlists_HW <= 2) jlist = 0;
       else                jlist = jcomp;

       if (jnode_box >=0 && !Zero_density_TF[jnode_box][jcomp]) {
             if (resid_only_flag==INIT_GUESS_FLAG)         dens=x[unk_rho][jnode_box]/Nseg_type[jcomp];
             else                                          dens=calc_dens_seg(jseg,jnode_box,x,FALSE);   
             if (Lhard_surf && Nodes_2_boundary_wall[jlist][jnode_box]!=-1){
                 weight = HW_boundary_weight (jcomp,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
             }
       }
       else if (jnode_box==-1 ||jnode_box==-3 ||jnode_box==-4)  {
            dens = constant_boundary(unk_rho,jnode_box)/Nseg_type[jcomp];
       }
       else                                                     dens=0.0;

       for (kbond=0; kbond<Nbonds_SegAll[jseg]; kbond++){
            if (Bonds_SegAll[jseg][kbond] !=-1){	
               kseg = Bonds_SegAll[jseg][kbond];
               kcomp=Unk2Comp[kseg];
               s2=Sigma_ff[kcomp][kcomp];
/*               s2=Bond_ff[kcomp][kcomp];*/

               y = y_cav(s1,s2,xi_2,xi_3);
               dy_dxi2=dy_dxi2_cav(s1,s2,xi_2,xi_3);
               dy_dxi3=dy_dxi3_cav(s1,s2,xi_2,xi_3);

               prefac2 = (PI/6.)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],2);
               prefac3 = (PI/6.)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],3);
               resid = -0.5*Fac_overlap[jcomp][kcomp]*weight*dens*(1./y)*(prefac2*dy_dxi2 + prefac3*dy_dxi3);
               if (resid_only_flag !=INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
               resid_sum += resid;
               if (Grafted[SegAll_to_Poly[jseg]]==TRUE) resid_sum_tethered[SegAll_to_Poly[jseg]]+=resid; 
            }

        } /* end of bond pair loop */
    }       /*  end of jseg loop */


    if (resid_only_flag==FALSE) {
       if (isten < stenJ->Length){
          offsetJ = sten_offsetJ[isten];
          weightJ = sten_weightJ[isten];
          jnode_boxJ = offset_to_node_box(ijk_box, offsetJ, reflect_flag);
          if (jnode_boxJ >= 0) {
               xi_2=x[unk_xi2][jnode_boxJ]; xi_3=x[unk_xi3][jnode_boxJ];
          }
          else {
               xi_2=constant_boundary(unk_xi2,jnode_boxJ);
               xi_3=constant_boundary(unk_xi3,jnode_boxJ);
          }

          for (jseg=0;jseg<Nseg_tot;jseg++){
              jseg_tmp = jseg;
              if (Pol_Sym_Seg[jseg] != -1) jseg_tmp = Pol_Sym_Seg[jseg];
              unk_rho = Unk2Comp[jseg_tmp]+Phys2Unk_first[DENSITY];
              jcomp=Unk2Comp[jseg];
              if (Type_poly==WJDC2 || Type_poly==WJDC3) unk_B=Phys2Unk_first[WJDC_FIELD]+jcomp;
              else if (Type_poly==WJDC) unk_B=Phys2Unk_first[WJDC_FIELD]+jseg; 

              s1=Sigma_ff[jcomp][jcomp];
/*              s1=Bond_ff[jcomp][jcomp];*/
              if (Nlists_HW <= 2) jlist = 0;
              else                jlist = jcomp;

              if (jnode_boxJ >=0 && !Zero_density_TF[jnode_boxJ][jcomp]){
                if (Lhard_surf) {
                    if (Nodes_2_boundary_wall[jlist][jnode_boxJ]!=-1)
                       weightJ = HW_boundary_weight (jcomp, jlist,
                         stenJ->HW_Weight[isten], jnode_boxJ, reflect_flag);
                }


                /* first compute density and first derivative terms 
                   that will be needed below */
/*                dens=x[unk_rho][jnode_boxJ]/Nseg_type[jcomp];*/
                dens = calc_dens_seg(jseg,jnode_boxJ,x,FALSE);
                first_deriv_sum=0.0;
                for (kbond=0; kbond<Nbonds_SegAll[jseg]; kbond++){
                  if (Bonds_SegAll[jseg][kbond] != -1){
                  kseg = Bonds_SegAll[jseg][kbond];
                  kcomp=Unk2Comp[kseg];
                  s2=Sigma_ff[kcomp][kcomp];
/*                  s2=Bond_ff[kcomp][kcomp];*/

                  y = y_cav(s1,s2,xi_2,xi_3);
                  dy_dxi2=dy_dxi2_cav(s1,s2,xi_2,xi_3);
                  dy_dxi3=dy_dxi3_cav(s1,s2,xi_2,xi_3);
   
                  /*d2y_dxi2_2=d2y_dxi2_sq(s1,s2,xi_2,xi_3);
                  d2y_dxi3_2=d2y_dxi3_sq(s1,s2,xi_2,xi_3);
                  d2y_dxi2_dxi3=d2y_dxi3_dxi2(s1,s2,xi_2,xi_3);*/

                  prefac2 = (PI/6.)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],2);
                  prefac3 = (PI/6.)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],3);

                  first_deriv_sum += Fac_overlap[jcomp][kcomp]*(prefac2*dy_dxi2 + prefac3*dy_dxi3)/y;
                  }
                }
 
                /* CAVITY variable derivatives */

                for (kbond=0; kbond<Nbonds_SegAll[jseg]; kbond++){
                  unk_GQ  = Phys2Unk_first[G_CHAIN] + Poly_to_Unk_SegAll[jseg][kbond];
                  if (Bonds_SegAll[jseg][kbond] != -1){
                  kseg = Bonds_SegAll[jseg][kbond];
                  kcomp=Unk2Comp[kseg];
                  s2=Sigma_ff[kcomp][kcomp];
/*                     s2=Bond_ff[kcomp][kcomp];*/

                  y = y_cav(s1,s2,xi_2,xi_3);
                  dy_dxi2=dy_dxi2_cav(s1,s2,xi_2,xi_3);
                  dy_dxi3=dy_dxi3_cav(s1,s2,xi_2,xi_3);

                  d2y_dxi2_2=d2y_dxi2_sq(s1,s2,xi_2,xi_3);
                  d2y_dxi3_2=d2y_dxi3_sq(s1,s2,xi_2,xi_3);
                  d2y_dxi2_dxi3=d2y_dxi3_dxi2(s1,s2,xi_2,xi_3);

                  prefac2 = (PI/6.)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],2);
                  prefac3 = (PI/6.)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],3);

                  first_deriv = (prefac2*dy_dxi2 + prefac3*dy_dxi3)/y;
 
                   /* Matrix entries for dR_Field/dcav_2 */
                   mat_val = -0.5*Fac_overlap[jcomp][kcomp]*weightJ*dens* (
                             (prefac2*d2y_dxi2_2 + prefac3*d2y_dxi2_dxi3)/y
                             - first_deriv*dy_dxi2/y );

                   if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[jnode_boxJ]+Solver_Unk[unk_xi2]*Nnodes]+=mat_val;
                   dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                             unk_xi2,jnode_boxJ,mat_val);

                   /* Matrix entries for dR_Field/dcav_3 */
                   mat_val = -0.5*Fac_overlap[jcomp][kcomp]*weightJ*dens* (
                              (prefac2*d2y_dxi2_dxi3 + prefac3*d2y_dxi3_2)/y
                              - first_deriv*dy_dxi3/y );

                   if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[jnode_boxJ]+Solver_Unk[unk_xi3]*Nnodes]+=mat_val;
                   dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                             unk_xi3,jnode_boxJ,mat_val);
                   }
               }

               /*Approximate matrix entries using a mean value of the segment density .....*/
                  if (Analyt_WJDC_Jac==FALSE){
                     mat_val = -0.5*weightJ*first_deriv_sum/Nseg_type[jcomp];

                     if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[jnode_boxJ]+Solver_Unk[unk_rho]*Nnodes]+=mat_val;
                     dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_rho,jnode_boxJ,mat_val);
                  }
                  else{
                
                   /* Analytic Matrix entries for dR_Field/dXi_alpha */   
                     mat_val = 0.5*weightJ*first_deriv_sum*(dens*((double)Nbonds_SegAll[jseg]-1.0)/x[unk_B][jnode_boxJ]);
                     if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[jnode_boxJ]+Solver_Unk[unk_B]*Nnodes]+=mat_val;
                     dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_B,jnode_boxJ,mat_val);

                   /* Analytic Matrix entries for dR_Field/dG_alpha */   
                 
                     for (kbond=0; kbond<Nbonds_SegAll[jseg]; kbond++){
                        unk_GQ  = Phys2Unk_first[G_CHAIN] + Poly_to_Unk_SegAll[jseg][kbond];
                         dens_Gderiv=calc_dens_seg_Gderiv(jseg,jnode_boxJ,kbond,x,FALSE);
                         mat_val = -0.5*weightJ*first_deriv_sum*dens_Gderiv;
                         if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[jnode_boxJ]+Solver_Unk[unk_GQ]*Nnodes]+=mat_val;
                         dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode, unk_GQ,jnode_boxJ,mat_val);
                     }


                 } /* end of analytical Jacobian section */

            }  /* check that the node has nonzero density */
          }  /* loop over all jseg */
       }  /* check on Jacobian fill */
     }  /* end Jacobian fill */

  } /* end of loop over stencil !! */

  /* Analytic Matrix entries for grafted chains */


  if(Analyt_WJDC_Jac!=FALSE && resid_only_flag==FALSE && Grafted_Logical==TRUE) { 
     for (iwall=0;iwall<Nwall;iwall++){
        for (ipol=0;ipol<Npol_comp;ipol++){
          if (WallType[iwall]==Graft_wall[ipol]){
            for (jsurf_node=0;jsurf_node<Nodes_Surf_Gsum[ipol][iwall];jsurf_node++){
                  jnode_box=Index_SurfNodes_Gsum[iwall][jsurf_node];

                  if (Nbonds_SegAll[Grafted_SegIDAll[ipol]]>2){
                         unk_B=Index_UnkB_Gsum[iwall][jsurf_node];
                         if (Grafted[ipol]==GRAFT_DENSITY) mat_val=resid_sum_tethered[ipol]*(GsumPrefac_XiDerivs[iwall][jsurf_node]/(Total_area_graft[ipol]*Gsum_graft[ipol]));
                         else                              mat_val=resid_sum_tethered[ipol]*(GsumPrefac_XiDerivs[iwall][jsurf_node]/Gsum_graft[ipol]);

                         if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node_extra[jnode_box]+Solver_Unk[unk_B]*Nnodes]+=mat_val;
                         dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_B,jnode_box,mat_val);
                   }

                   for (jbond=0;jbond<Nbonds_SegAll[Grafted_SegIDAll[ipol]];jbond++){
                         if (Bonds[ipol][Grafted_SegID[ipol]][jbond] != -1 ){
                           unk_GQ=Index_UnkGQ_Gsum[iwall][jsurf_node][jbond];
                           if (Grafted[ipol]==GRAFT_DENSITY) mat_val=-resid_sum_tethered[ipol]*(GsumPrefac_GDerivs[iwall][jsurf_node][jbond]/(Total_area_graft[ipol]*Gsum_graft[ipol]));
                           else                              mat_val=-resid_sum_tethered[ipol]*(GsumPrefac_GDerivs[iwall][jsurf_node][jbond]/Gsum_graft[ipol]);

                           if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node_extra[jnode_box]+Solver_Unk[unk_GQ]*Nnodes]+=mat_val;
                           dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_GQ,jnode_box,mat_val);
                         }
                    }
                 }

             }
          }
        }
  }

  return resid_sum;
}
/********************************************************************************************/
double calc_dens_seg(int iseg,int inode_box,double **x,int flag)
{
   int boltz_pow,unk_GQ,ibond,itype_mer,unk_B;
   double fac1,dens;
/* int unk_GQ_test; */

   boltz_pow = -(Nbonds_SegAll[iseg]-1);

   itype_mer=Unk2Comp[iseg]; /* note that itype_mer is also known as icomp */
   if (Type_poly==WJDC2 || Type_poly==WJDC3) unk_B=Phys2Unk_first[WJDC_FIELD]+itype_mer;
   else /*if (Type_poly==WJDC)*/ unk_B=Phys2Unk_first[WJDC_FIELD]+iseg; 

   fac1 = prefactor_rho_wjdc(iseg,inode_box,x);

   for (ibond=0; ibond<Nbonds_SegAll[iseg]; ibond++) {
              unk_GQ  = Phys2Unk_first[G_CHAIN] + Poly_to_Unk_SegAll[iseg][ibond];
/*              unk_GQ_test = unk_GQ-Phys2Unk_first[G_CHAIN];
              if (Pol_Sym[unk_GQ_test] != -1) unk_GQ=Pol_Sym[unk_GQ_test] + Phys2Unk_first[G_CHAIN];*/
              fac1 *= x[unk_GQ][inode_box];

    }
   dens = fac1*POW_DOUBLE_INT(x[unk_B][inode_box],boltz_pow);

   return(dens);
}
/********************************************************************************************/
double calc_dens_seg_Gderiv(int iseg,int inode_box,int kbond, double **x,int flag)
{
   int boltz_pow,unk_GQ,ibond,itype_mer,unk_B;
   double fac1,gderiv;
/* int unk_GQ_test; */

   boltz_pow = -(Nbonds_SegAll[iseg]-1);

   itype_mer=Unk2Comp[iseg]; /* note that itype_mer is also known as icomp */
   if (Type_poly==WJDC2 || Type_poly==WJDC3) unk_B=Phys2Unk_first[WJDC_FIELD]+itype_mer;
   else /*if (Type_poly==WJDC)*/ unk_B=Phys2Unk_first[WJDC_FIELD]+iseg; 

   fac1 = prefactor_rho_wjdc(iseg,inode_box,x);
   for (ibond=0; ibond<Nbonds_SegAll[iseg]; ibond++) {
        if (ibond != kbond){    /* taking a first derivative eliminates one G from the calcualtion */
              unk_GQ  = Phys2Unk_first[G_CHAIN] + Poly_to_Unk_SegAll[iseg][ibond];
/*              unk_GQ_test = unk_GQ-Phys2Unk_first[G_CHAIN];
              if (Pol_Sym[unk_GQ_test] != -1) unk_GQ=Pol_Sym[unk_GQ_test] + Phys2Unk_first[G_CHAIN];*/
              fac1 *= x[unk_GQ][inode_box];
         }
    }
   gderiv = fac1*POW_DOUBLE_INT(x[unk_B][inode_box],boltz_pow);
   return(gderiv);
}
/********************************************************************************************/
