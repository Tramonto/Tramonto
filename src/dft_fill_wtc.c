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

/*
 *  FILE: dft_fill_wtc.c
 *
 *  This file contains the fill for the Wertheim-Tripathi-Chapman functionals for 
 *  bonded systems
 */

/*#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"*/
#include "dft_fill_wtc.h"

/***********************************************************************************************/
/*load_polyTC_diagEL:  Here add the position diagonals of the WTC association/bonding contribution 
                       to the Euler-Lagrange equation for the Wertheim-Tripathi-Chapman theory */
double load_polyTC_diagEL(int iunk,int loc_inode,int inode_box, int icomp,
                          int izone,int *ijk_box,double **x,int resid_only_flag)
{
  double y,n,resid,resid_sum,xi_2,xi_3,dy_dxi2,dy_dxi3,mat_val;
  int iseg,jseg,ibond,jcomp,unk_xi2,unk_xi3,unk_bond;

  iseg = iunk-Phys2Unk_first[DENSITY];


  resid_sum=0.0;
  for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
     jseg = Bonds_SegAll[iseg][ibond];
     jcomp = Unk2Comp[jseg];
     unk_xi2=Phys2Unk_first[CAVWTC];
     unk_xi3=Phys2Unk_first[CAVWTC]+1;

     xi_2=x[unk_xi2][inode_box];
     xi_3=x[unk_xi3][inode_box];

     y=y_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],xi_2,xi_3);
     dy_dxi2=dy_dxi2_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],xi_2,xi_3);
     dy_dxi3=dy_dxi3_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],xi_2,xi_3);
/*     y=y_cav(Bond_ff[icomp][icomp],Bond_ff[jcomp][jcomp],xi_2,xi_3);
     dy_dxi2=dy_dxi2_cav(Bond_ff[icomp][icomp],Bond_ff[jcomp][jcomp],xi_2,xi_3);
     dy_dxi3=dy_dxi3_cav(Bond_ff[icomp][icomp],Bond_ff[jcomp][jcomp],xi_2,xi_3);*/
     unk_bond = Poly_to_Unk_SegAll[iseg][ibond];
     if (Pol_Sym[unk_bond] != -1) unk_bond=Pol_Sym[Poly_to_Unk_SegAll[iseg][ibond]];
     unk_bond += Phys2Unk_first[BONDWTC];
     n=x[unk_bond][inode_box];

                        /* FIRST ADD TERMS THAT ARE ON DIAGONAL WITH REPECT TO POSITION. 
                          (ALL JACOBIAN ENTRIES FOUND AT INODE_BOX...)*/
     /*resid = 0.5*Fac_overlap[icomp][jcomp]*(1.-log(y)-log(n));*/
     resid = 0.5*(1.-Fac_overlap[icomp][jcomp]*log(y)-log(n));
     resid_sum+=resid;
     dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

     mat_val = -0.5*Fac_overlap[icomp][jcomp]*(1./y)*dy_dxi2;
     dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_xi2,inode_box,mat_val);

     mat_val = -0.5*Fac_overlap[icomp][jcomp]*(1./y)*dy_dxi3;
     dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_xi3,inode_box,mat_val);

     /*mat_val = -0.5*Fac_overlap[icomp][jcomp]*(1./n);*/
     mat_val = -0.5*(1./n);
     dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_bond,inode_box,mat_val);
   }
   return resid_sum;
}
/***********************************************************************************************/
/*load_polyTC_bondEL:  Here add the bond contributions of the WTC association/bonding functionals
                       to the Euler-Lagrange equation for the Wertheim-Tripathi-Chapman theory */
double load_polyTC_bondEL(int iunk,int loc_inode,int inode_box,int icomp,int izone,int *ijk_box,double **x,int resid_only_flag)
{
  int iseg,jseg,ibond,jcomp,jzone,jnode_box,reflect_flag[3],junk_rho,unk_bond,jnode_boxJ,jlist;
  int   **sten_offset, *offset, isten;
  int   **sten_offsetJ, *offsetJ;
  double *sten_weightJ,weightJ;
  double *sten_weight,  weight,fac;
  struct Stencil_Struct *sten;
  struct Stencil_Struct *stenJ;
  double resid,mat_val,resid_sum=0.0;

  iseg = iunk-Phys2Unk_first[DENSITY];
  jzone=izone;

  for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
     jseg = Bonds_SegAll[iseg][ibond];
     junk_rho = jseg;
     if (Pol_Sym_Seg[jseg] != -1) junk_rho = Pol_Sym_Seg[jseg];
     junk_rho += Phys2Unk_first[DENSITY];
     unk_bond = Poly_to_Unk_SegAll[iseg][ibond];
     if (Pol_Sym[unk_bond] != -1) unk_bond=Pol_Sym[Poly_to_Unk_SegAll[iseg][ibond]];
     unk_bond += Phys2Unk_first[BONDWTC];
 
     jcomp = Unk2Comp[jseg];

     if (Nlists_HW <= 2) jlist = 0;
     else                jlist = jcomp;


                        /* NOW DO OFF DIAGONAL ENTRIES ASSOCIATED WITH THE BONDING NONLOCAL DENSITIES */

      sten = &(Stencil[DELTA_FN_BOND][izone][icomp+Ncomp*jcomp]);
      sten_offset = sten->Offset;
      sten_weight = sten->Weight;

      stenJ = &(Stencil[DELTA_FN_BOND][jzone][icomp+Ncomp*jcomp]);
      sten_offsetJ = stenJ->Offset;
      sten_weightJ = stenJ->Weight;

      for (isten = 0; isten < sten->Length; isten++) {
        offset = sten_offset[isten];
        weight = sten_weight[isten];
        jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);

        if (jnode_box >= 0 && !Zero_density_TF[jnode_box][jcomp]) {
            if (Lhard_surf) {
                if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1)
                   weight = HW_boundary_weight
                    (jcomp,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
            }
            resid = -0.5*weight*x[junk_rho][jnode_box]/x[unk_bond][jnode_box];
        }
        else if (jnode_box == -1 || jnode_box ==-3 || jnode_box == -4){
            resid = -0.5*weight*constant_boundary(junk_rho,jnode_box)/constant_boundary(unk_bond,jnode_box);
        }
        else resid = 0.0;
        dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
        resid_sum += resid;

        if (!resid_only_flag) {
           if (isten < stenJ->Length){
              offsetJ = sten_offsetJ[isten];
              weightJ = sten_weightJ[isten];
              jnode_boxJ = offset_to_node_box(ijk_box, offsetJ, reflect_flag);

              if (jnode_boxJ >=0 && !Zero_density_TF[jnode_boxJ][jcomp]){
                    if (Lhard_surf) {
                       if (Nodes_2_boundary_wall[jlist][jnode_boxJ]!=-1)
                          weightJ = HW_boundary_weight (jcomp, jlist,
                            stenJ->HW_Weight[isten], jnode_boxJ, reflect_flag);
                    }
                    mat_val = -0.5*weight/x[unk_bond][jnode_boxJ];
                    dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                                 junk_rho,jnode_boxJ,mat_val);
                    mat_val = 0.5*weight*x[junk_rho][jnode_box]/(x[unk_bond][jnode_boxJ]*x[unk_bond][jnode_boxJ]);
                    dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                                unk_bond,jnode_boxJ,mat_val);
              }
           }
       }

      }  /* end of loop over the bonding stencils */

  } /* end of loop over segments bonded to the segment of interest iseg */
  return(resid_sum);
}
/***********************************************************************************************/
/*load_polyTC_cavityEL:  Here add the bond contributions of the WTC association/bonding (cavity term) functionals
                       to the Euler-Lagrange equation for the Wertheim-Tripathi-Chapman theory */
double load_polyTC_cavityEL(int iunk,int loc_inode,int inode_box,int icomp,int izone,int *ijk_box,double **x,int resid_only_flag)
{
  int iseg,jseg,kseg,kbond,jcomp,unk_xi2,unk_xi3,unk_rho,junk_rho,kunk_rho,kcomp;
  double xi_2,xi_3,s1,s2,y,dy_dxi2,dy_dxi3,prefac2,prefac3;
  double d2y_dxi2_2,d2y_dxi3_2,d2y_dxi2_dxi3;
  int jzone,jnode_box,jlist,reflect_flag[3],jnode_boxJ;
  int   **sten_offset, *offset, isten;
  int   **sten_offsetJ, *offsetJ;
  double *sten_weightJ,weightJ;
  double *sten_weight,  weight,fac;
  struct Stencil_Struct *sten;
  struct Stencil_Struct *stenJ;
  double resid,mat_val,resid_sum,first_deriv,dens;

  iseg = iunk - Phys2Unk_first[DENSITY];

  resid_sum=0.0;
  unk_xi2 = Phys2Unk_first[CAVWTC];
  unk_xi3 = Phys2Unk_first[CAVWTC]+1;
  jzone=izone;

  sten = &(Stencil[THETA_FN_SIG][izone][icomp]);
  sten_offset = sten->Offset;
  sten_weight = sten->Weight;

  stenJ = &(Stencil[THETA_FN_SIG][jzone][icomp]);
  sten_offsetJ = stenJ->Offset;
  sten_weightJ = stenJ->Weight;

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
       unk_rho = jseg; 
       if (Pol_Sym_Seg[jseg] != -1) unk_rho=Pol_Sym_Seg[jseg];
       unk_rho += Phys2Unk_first[DENSITY]; 

       jcomp=Unk2Comp[jseg];
       s1=Sigma_ff[jcomp][jcomp];
/*       s1=Bond_ff[jcomp][jcomp];*/
       if (Nlists_HW <= 2) jlist = 0;
       else                jlist = jcomp;
       
       if (Lhard_surf && jnode_box>=0 && !Zero_density_TF[jnode_box][jcomp]) {
            if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1)
               weight = HW_boundary_weight
                (jcomp,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
       }
       for (kbond=0; kbond<Nbonds_SegAll[jseg]; kbond++){
            kseg = Bonds_SegAll[jseg][kbond];
            kunk_rho = kseg;
            if (Pol_Sym_Seg[kseg] != -1) kunk_rho=Pol_Sym_Seg[kseg];
            kunk_rho += Phys2Unk_first[DENSITY];
            kcomp=Unk2Comp[kseg];
            s2=Sigma_ff[kcomp][kcomp];
/*            s2=Bond_ff[kcomp][kcomp];*/

            y = y_cav(s1,s2,xi_2,xi_3);
            dy_dxi2=dy_dxi2_cav(s1,s2,xi_2,xi_3);
            dy_dxi3=dy_dxi3_cav(s1,s2,xi_2,xi_3);
             
            prefac2 = (PI/6.)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],2);
            prefac3 = (PI/6.)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],3);

            if (jnode_box >=0 && !Zero_density_TF[jnode_box][jcomp]) dens = x[unk_rho][jnode_box];
            else if (jnode_box==-1 ||jnode_box==-3 ||jnode_box==-4)  dens = constant_boundary(unk_rho,jnode_box);
            else                                                     dens=0.0;

            resid = -0.5*Fac_overlap[jcomp][kcomp]*weight*dens*(1./y)*(prefac2*dy_dxi2 + prefac3*dy_dxi3);
            dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
            resid_sum += resid;

        } /* end of bond pair loop */
    }        /* end of jseg loop */


    if (!resid_only_flag) {
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
              unk_rho = jseg; 
              if (Pol_Sym_Seg[jseg] != -1) unk_rho = Pol_Sym_Seg[jseg];
              unk_rho += Phys2Unk_first[DENSITY]; 
              jcomp=Unk2Comp[jseg];
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

                for (kbond=0; kbond<Nbonds_SegAll[jseg]; kbond++){
                  kseg = Bonds_SegAll[jseg][kbond];
                  kunk_rho = kseg;
                  if (Pol_Sym_Seg[kseg] != -1) kunk_rho = Pol_Sym_Seg[kseg];
                  kunk_rho += Phys2Unk_first[DENSITY];
                  kcomp=Unk2Comp[kseg];
                  s2=Sigma_ff[kcomp][kcomp];
/*                  s2=Bond_ff[kcomp][kcomp];*/
 
                  y = y_cav(s1,s2,xi_2,xi_3);
                  dy_dxi2=dy_dxi2_cav(s1,s2,xi_2,xi_3);
                  dy_dxi3=dy_dxi3_cav(s1,s2,xi_2,xi_3);
 
                  d2y_dxi2_2=d2y_dxi2_sq(s1,s2,xi_2,xi_3);
                  d2y_dxi3_2=d2y_dxi3_sq(s1,s2,xi_2,xi_3);
                  d2y_dxi2_dxi3=d2y_dxi3_dxi2(s1,s2,xi_2,xi_3);
   
                  prefac2 = (PI/6.)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],2);
                  prefac3 = (PI/6.)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],3);

                  if (jnode_boxJ >=0 && !Zero_density_TF[jnode_boxJ][jcomp]) {
                      dens = x[unk_rho][jnode_boxJ];
                      first_deriv = (prefac2*dy_dxi2 + prefac3*dy_dxi3)/y;
  
                      mat_val = -0.5*Fac_overlap[jcomp][kcomp]*weight*first_deriv;
                      dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                                unk_rho,jnode_boxJ,mat_val);
  
                      mat_val = -0.5*Fac_overlap[jcomp][kcomp]*weight*dens* ( 
                                (prefac2*d2y_dxi2_2 + prefac3*d2y_dxi2_dxi3)/y 
                                - first_deriv*dy_dxi2/y );
  
                      dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                                unk_xi2,jnode_boxJ,mat_val);
  
                      mat_val = -0.5*Fac_overlap[jcomp][kcomp]*weight*dens* (
                                 (prefac2*d2y_dxi2_dxi3 + prefac3*d2y_dxi3_2)/y
                                 - first_deriv*dy_dxi3/y );
  
                      dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                                unk_xi3,jnode_boxJ,mat_val);
                  }
                }  /*loop over kbonds */
            }  /* check that the node has nonzero density */
          }  /* loop over all jseg */
       }  /* check on Jacobian fill */
     }  /* end Jacobian fill */

  } /* end of loop over stencil !! */

  return resid_sum;
}
/********************************************************************************************/
   double d2y_dxi2_sq(double sigma_1,double sigma_2,double xi_2, double xi_3)
{
   double one_m_xi3_sq,one_m_xi3_cb,sig_sum,sig_m,d2y_dxi2sq;

   one_m_xi3_sq= (1.-xi_3)*(1.-xi_3); 
   one_m_xi3_cb= (1.-xi_3)*(1.-xi_3)*(1.-xi_3); 
   sig_sum=sigma_1+sigma_2;
   sig_m=sigma_1*sigma_2;

  d2y_dxi2sq = 4.*(sig_m/sig_sum)*(sig_m/sig_sum)*(1./one_m_xi3_cb);

   return d2y_dxi2sq;
}
/*******************************************************************************************/
   double d2y_dxi3_dxi2(double sigma_1,double sigma_2,double xi_2, double xi_3)
{
   double one_m_xi3_sq,one_m_xi3_cb,one_m_xi3_4th,sig_sum,sig_m,d2y_dxi3_dxi2;

   one_m_xi3_sq= (1.-xi_3)*(1.-xi_3); 
   one_m_xi3_cb= (1.-xi_3)*(1.-xi_3)*(1.-xi_3); 
   one_m_xi3_4th= (1.-xi_3)*one_m_xi3_cb; 
   sig_sum=sigma_1+sigma_2;
   sig_m=sigma_1*sigma_2;

  d2y_dxi3_dxi2 = 6.*(sig_m/sig_sum)*(1.0/one_m_xi3_cb)
                 + 12.*(sig_m/sig_sum)*(sig_m/sig_sum)*(xi_2/one_m_xi3_4th);

   return d2y_dxi3_dxi2;
}
/********************************************************************************************/
   double d2y_dxi3_sq(double sigma_1,double sigma_2,double xi_2, double xi_3)
{
   double one_m_xi3_sq,one_m_xi3_cb,one_m_xi3_4th,one_m_xi3_5th,sig_sum,sig_m,d2y_dxi3sq;

   one_m_xi3_sq= (1.-xi_3)*(1.-xi_3); 
   one_m_xi3_cb= (1.-xi_3)*(1.-xi_3)*(1.-xi_3); 
   one_m_xi3_4th= (1.-xi_3)*one_m_xi3_cb; 
   one_m_xi3_5th= (1.-xi_3)*one_m_xi3_4th; 
   sig_sum=sigma_1+sigma_2;
   sig_m=sigma_1*sigma_2;

  d2y_dxi3sq = 2./one_m_xi3_cb +(18.*sig_m/sig_sum)*(xi_2/one_m_xi3_4th)+
       24.*(sig_m/sig_sum)*(sig_m/sig_sum)*(xi_2*xi_2/one_m_xi3_5th);

   return d2y_dxi3sq;
}
/*********************************************************************************************/
/* load_bond_wtc: 
             In this routine we load the Jacobian entries for
             freely-jointed polymer bonds - WTC theory.                      */
double load_bond_wtc(int iunk, int loc_inode, int inode_box,int *ijk_box,
                    int izone, double **x,int resid_only_flag)
{
  int junk,unk_bond,pol_num,iseg,bond_num,jseg,jcomp,icomp,jzone_flag,ibond;
  double resid,mat_val,resid_bondwtc;

  resid=-x[iunk][inode_box];
  mat_val=-1.0;
  dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
  dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);
  resid_bondwtc=resid;

  unk_bond=iunk-Phys2Unk_first[BONDWTC];
  iseg=BondAll_to_isegAll[unk_bond];

  if (Pol_Sym[unk_bond]==-1){

     ibond=BondAll_to_ibond[unk_bond];
     jseg=Bonds_SegAll[iseg][ibond];

     jcomp = Unk2Comp[jseg];
     icomp = Unk2Comp[iseg];
     junk = jseg;
     if (Pol_Sym_Seg[jseg] != -1) junk=Pol_Sym_Seg[jseg];
     junk += Phys2Unk_first[DENSITY];

     jzone_flag=FALSE;

     resid_bondwtc+= resid_and_Jac_sten_fill(DELTA_FN_BOND,x,iunk,junk,
                   icomp,jcomp,loc_inode,inode_box,izone,
                   ijk_box,resid_only_flag,jzone_flag,
                    NULL, &resid_rho_bar,&jac_rho_bar);
  }
  else{
/*      if (Proc==0){
         printf("Polymer symmetries not working for WTC functionals: see code in load_bond_wtc\n");
         exit(-1);
      }*/
      junk=Pol_Sym[unk_bond]+Phys2Unk_first[BONDWTC];
      resid = x[junk][inode_box];
      resid_bondwtc+=resid;
      dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
      mat_val=1.0; 
      dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,inode_box,mat_val);
  }

  return(resid_bondwtc);
}
/*********************************************************************************************/
/* load_cavity_wtc: 
             In this routine we load the Jacobian entries for
             the cavity correlation function xi variables....WTC theory   */
double load_cavity_wtc(int iunk, int loc_inode, int inode_box, int *ijk_box,
                    int izone, double **x,int resid_only_flag)
{
  double resid,mat_val,resid_cavity;
  int jzone_flag;

  jzone_flag=FALSE;

  resid=-x[iunk][inode_box];
  mat_val=-1.0;
  dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);
  dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
  resid_cavity=resid;

  resid_cavity+=resid_and_Jac_sten_fill_sum_Ncomp(THETA_FN_SIG,x,iunk,loc_inode,inode_box,izone,
                   ijk_box,resid_only_flag,jzone_flag,
                    &prefactor_cavity_wtc, &resid_rho_bar,&jac_rho_bar);

  return(resid_cavity);
}
/*****************************************************************************/
double prefactor_cavity_wtc(int iunk,int icomp,int *offset)
{
  double fac;
  int ipow;

  if (iunk-Phys2Unk_first[CAVWTC]==0) ipow=2;
  else                                    ipow=3;

  fac=(PI/6.0)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],ipow);

  return (fac);
}
/*****************************************************************************/

