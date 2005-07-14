/*
 *  FILE: dft_fill_wtc.c
 *
 *  This file contains the fill for the Wertheim-Tripathi-Chapman functionals for 
 *  bonded systems
 */

#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"
double d2y_dxi3_sq(double,double,double,double);
double d2y_dxi2_sq(double,double,double,double);
double d2y_dxi3_dxi2(double,double,double,double);

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
     unk_xi2=Phys2Unk_first[CAVITY_WTC];
     unk_xi3=Phys2Unk_first[CAVITY_WTC]+1;

     xi_2=x[unk_xi2][inode_box];
     xi_3=x[unk_xi3][inode_box];

     y=y_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],xi_2,xi_3);
     dy_dxi2=dy_dxi2_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],xi_2,xi_3);
     dy_dxi3=dy_dxi3_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],xi_2,xi_3);
     unk_bond = Phys2Unk_first[BOND_WTC]+Poly_to_Unk_SegAll[iseg][ibond];
     n=x[unk_bond][inode_box];

                        /* FIRST ADD TERMS THAT ARE ON DIAGONAL WITH REPECT TO POSITION. 
                          (ALL JACOBIAN ENTRIES FOUND AT INODE_BOX...)*/
     resid = 0.5*(1.-log(y)-log(n));
     resid_sum+=resid;
     dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

     mat_val = -0.5*(1./y)*dy_dxi2;
     dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_xi2,inode_box,mat_val);

     mat_val = -0.5*(1./y)*dy_dxi3;
     dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_xi3,inode_box,mat_val);

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
     junk_rho = jseg+Phys2Unk_first[DENSITY];
     unk_bond = Poly_to_Unk_SegAll[iseg][ibond]+Phys2Unk_first[BOND_WTC];
     jcomp = Unk2Comp[junk_rho];

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
  double resid,mat_val,resid_sum;

  iseg = iunk - Phys2Unk_first[DENSITY];

  resid_sum=0.0;
  unk_xi2 = Phys2Unk_first[CAVITY_WTC];
  unk_xi3 = Phys2Unk_first[CAVITY_WTC]+1;
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

    resid=0.0;
    for (jseg=0;jseg<Nseg_tot;jseg++){
       unk_rho = Phys2Unk_first[DENSITY]+jseg; 
       jcomp=Unk2Comp[unk_rho];
       s1=Sigma_ff[jcomp][jcomp];
       if (Nlists_HW <= 2) jlist = 0;
       else                jlist = jcomp;

       if (jnode_box >= 0 && !Zero_density_TF[jnode_box][jcomp]) {
         xi_2=x[unk_xi2][jnode_box]; xi_3=x[unk_xi3][jnode_box];
       }
       else if (jnode_box == -1 || jnode_box ==-3 || jnode_box == -4){
          if (jnode_box == -1) {
             xi_2 = Xi_cav_b[2]; xi_3=Xi_cav_b[3];
          } 
          else if (jnode_box == -3) {xi_2 = Xi_cav_LBB[2]; xi_3=Xi_cav_LBB[3];} 
          else if (jnode_box == -4) {xi_2 = Xi_cav_RTF[2]; xi_3=Xi_cav_RTF[3];} 
          
       }
       
       if (Lhard_surf && jnode_box>=0 && !Zero_density_TF[jnode_box][jcomp]) {
            if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1)
               weight = HW_boundary_weight
                (jcomp,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
       }
       for (kbond=0; kbond<Nbonds_SegAll[jseg]; kbond++){
            kseg = Bonds_SegAll[jseg][kbond];
            kunk_rho = Phys2Unk_first[DENSITY]+kseg;
            kcomp=Unk2Comp[kunk_rho];
            s2=Sigma_ff[kcomp][kcomp];

            y = y_cav(s1,s2,xi_2,xi_3);
            dy_dxi2=dy_dxi2_cav(s1,s2,xi_2,xi_3);
            dy_dxi3=dy_dxi3_cav(s1,s2,xi_2,xi_3);
             
            prefac2 = (PI/6.)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],2);
            prefac3 = (PI/6.)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],3);

            if (jnode_box >=0 && !Zero_density_TF[jnode_box][jcomp]){
                resid = -0.5*(prefac2*dy_dxi2 + prefac3*dy_dxi3)*weight*(x[unk_rho][jnode_box]/y);
            }
            else if (jnode_box==-1 ||jnode_box==-3 ||jnode_box==-4){
                resid = -0.5*(prefac2*dy_dxi2 + prefac3*dy_dxi3)*weight*(constant_boundary(unk_rho,jnode_box)/y);
            }
            else resid=0.0;

            dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
            resid_sum += resid;

        } /* end of bond pair loop */
    }        /* end of jseg loop */


    if (!resid_only_flag) {
       if (isten < stenJ->Length){
          offsetJ = sten_offsetJ[isten];
          weightJ = sten_weightJ[isten];
          jnode_boxJ = offset_to_node_box(ijk_box, offsetJ, reflect_flag);
          xi_2=x[unk_xi2][jnode_boxJ];
          xi_3=x[unk_xi3][jnode_boxJ];

          for (jseg=0;jseg<Nseg_tot;jseg++){
              junk_rho = Phys2Unk_first[DENSITY]+jseg; 
              jcomp=Unk2Comp[junk_rho];
              s1=Sigma_ff[jcomp][jcomp];
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
                   kunk_rho = Phys2Unk_first[DENSITY]+kseg;
                   kcomp=Unk2Comp[kunk_rho];
                   s2=Sigma_ff[kcomp][kcomp];
  
                   y = y_cav(s1,s2,xi_2,xi_3);
                   dy_dxi2=dy_dxi2_cav(s1,s2,xi_2,xi_3);
                   dy_dxi3=dy_dxi3_cav(s1,s2,xi_2,xi_3);
  
                   d2y_dxi2_2=d2y_dxi2_sq(s1,s2,xi_2,xi_3);
                   d2y_dxi3_2=d2y_dxi3_sq(s1,s2,xi_2,xi_3);
                   d2y_dxi2_dxi3=d2y_dxi3_dxi2(s1,s2,xi_2,xi_3);
   
                   prefac2 = (PI/6.)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],2);
                   prefac3 = (PI/6.)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],3);
  
                   mat_val = -0.5*(prefac2*dy_dxi2 + prefac3*dy_dxi3)*weight/y;
                   dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                             junk_rho,jnode_boxJ,mat_val);
  
                   mat_val = -0.5*((prefac2*d2y_dxi2_2 + prefac3*d2y_dxi2_dxi3)*weight*(x[junk_rho][jnode_boxJ]/y)-
                             0.5*(prefac2*dy_dxi2 + prefac3*dy_dxi3)*weight*(x[junk_rho][jnode_boxJ]/(y*y))*dy_dxi2);
  
                   dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                             unk_xi2,jnode_boxJ,mat_val);
  
                   mat_val = -0.5*((prefac2*d2y_dxi2_dxi3 + prefac3*d2y_dxi3_2)*weight*(x[junk_rho][jnode_boxJ]/y)-
                             0.5*(prefac2*dy_dxi2 + prefac3*dy_dxi3)*weight*(x[junk_rho][jnode_boxJ]/(y*y))*dy_dxi3);
  
                   dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                             unk_xi3,jnode_boxJ,mat_val);
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
double load_bond_wtc(int iunk, int loc_inode, int inode_box,
                    int izone,int *ijk_box, double **x)
{
  int   **sten_offset, *offset, isten;
  double *sten_weight,  weight;
  struct Stencil_Struct *sten;

  int jlist,junk;
  int reflect_flag[NDIM_MAX];
  int i,j,jnode_box;
  int unk_bond,pol_num,iseg,bond_num,jseg,jcomp,icomp;
  double resid,resid_sum,mat_val;

  unk_bond=iunk-Phys2Unk_first[BOND_WTC];
  pol_num = Unk_to_Poly[unk_bond];
  iseg = Unk_to_Seg[unk_bond];
  bond_num = Unk_to_Bond[unk_bond];
  jseg = Bonds[pol_num][iseg][bond_num];
  jcomp = Type_mer[pol_num][jseg];
  icomp = Type_mer[pol_num][iseg];
  junk = SegChain2SegAll[pol_num][jseg]+Phys2Unk_first[DENSITY];

  if (Nlists_HW <= 2) jlist = 0;
  else                jlist = jcomp; 

  sten = &(Stencil[DELTA_FN_BOND][izone][icomp+Ncomp*jcomp]);
  sten_offset = sten->Offset;
  sten_weight = sten->Weight;

  resid_sum=0.0;
  for (isten = 0; isten < sten->Length; isten++) {
     offset = sten_offset[isten];
     weight = sten_weight[isten];

     /* Find the Stencil point */
     jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);
     if (jnode_box >= 0 && !Zero_density_TF[jnode_box][junk-Phys2Unk_first[DENSITY]]) {
        if (Lhard_surf) {
           if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1) 
           weight = HW_boundary_weight 
                    (icomp+Ncomp*jcomp,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
        }

        resid_sum += weight*x[junk][jnode_box]; 
        mat_val = weight;
        dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,jnode_box,mat_val);
           
     }
     else if ( jnode_box == -1 || jnode_box ==-3 || jnode_box == -4) 
             resid_sum += weight*constant_boundary(junk,jnode_box);
  }
  return(resid_sum);
}
/**********************************************************************************************/
/* load_cavity_wtc: 
             In this routine we load the Jacobian entries for
             the cavity correlation function xi variables....WTC theory   */
double load_cavity_wtc(int iunk, int loc_inode, int inode_box,
                    int izone,int *ijk_box, double **x)
{
  int ipow,iseg,unk_rho;
  int   **sten_offset, *offset, isten;
  double *sten_weight,  weight;
  struct Stencil_Struct *sten;

  int ilist,icomp;
  int reflect_flag[NDIM_MAX];
  int i,j,jnode_box;
  double resid,resid_sum,mat_val;

  if (iunk-Phys2Unk_first[CAVITY_WTC]==0) ipow=2;
  else                                    ipow=3;

  resid_sum=0.0;
  for (iseg=0; iseg<Nseg_tot; iseg++){
      unk_rho = iseg+Phys2Unk_first[DENSITY];

      icomp = Unk2Comp[iseg];
      if (Nlists_HW <= 2) ilist = 0;
      else                ilist = icomp;

     sten = &(Stencil[THETA_FN_SIG][izone][icomp]);
     sten_offset = sten->Offset;
     sten_weight = sten->Weight;

     for (isten = 0; isten < sten->Length; isten++) {
        offset = sten_offset[isten];
        weight = sten_weight[isten];

        /* Find the Stencil point */
        jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);
        if (jnode_box >= 0 && !Zero_density_TF[jnode_box][icomp]) {
           if (Lhard_surf) {
             if (Nodes_2_boundary_wall[ilist][jnode_box]!=-1) 
             weight = HW_boundary_weight(icomp,ilist,sten->HW_Weight[isten], jnode_box, reflect_flag);
            }

        resid_sum += (PI/6.0)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],ipow)*weight*x[unk_rho][jnode_box]; 
        mat_val = (PI/6.0)*weight*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],ipow);
        dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_rho,jnode_box,mat_val);
     }
     else if ( jnode_box == -1 || jnode_box ==-3 || jnode_box == -4) {
        resid_sum += (PI/6.0)*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],ipow)*weight*constant_boundary(unk_rho,jnode_box);
     }
  }
  }  
  return(resid_sum);
}
