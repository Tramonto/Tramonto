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

static struct RB_Struct d2phi_drb2_delta_rb(int, int, double **,double, 
					    int *,double *,double,double,
					    double);

static struct RB_Struct d2phi_drb2_theta_rb(int, int, double **,double,int *);

static struct RB_Struct d2phi_drb2_delta2_rb(int, int, double **,double, 
				      int *,double *,double, double,double);

static struct RB_Struct d2phi_drb2_theta2_rb(int, int, double **,double,int *);


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
  jzone = find_jzone(izone);

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
            junk=Phys2Unk_first[RHOBAR_ROSEN];

            for (idim = 0; idim<Ndim; idim++){
               if (reflect_flag[idim] == FALSE) sign[idim] = 1.0;
               else sign[idim] = -1.0;
            }

            if (sten_type == DELTA_FN) {
               if (Type_func == 0) tmp = 
                          d2phi_drb2_delta_rb(junk,jnode_boxJ,x,weightJ,offsetJ, 
			  sign,Inv_rad[icomp],Inv_4pir[icomp], 
			  Inv_4pirsq[icomp]);

               else       tmp = 
                          d2phi_drb2_delta2_rb(junk,jnode_boxJ,x,weightJ,offsetJ, 
			  sign,Inv_rad[icomp],Inv_4pir[icomp], 
			  Inv_4pirsq[icomp]);
            }
            else if (sten_type == THETA_FN) {
               if (Type_func == 0) tmp = 
                             d2phi_drb2_theta_rb(junk,jnode_boxJ,x,weightJ,offsetJ);
               else    tmp = d2phi_drb2_theta2_rb(junk,jnode_boxJ,x,weightJ,offsetJ);
            }
            numEntries=4;
            values[0]=Fac_overlap_hs[icomp]*tmp.S3; values[1]=Fac_overlap_hs[icomp]*tmp.S2; values[2]=Fac_overlap_hs[icomp]*tmp.S1; values[3]=Fac_overlap_hs[icomp]*tmp.S0;
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
/* load_rho_bar_s:  Load scalar rho_bar definition equations  

   WITH JACOBIAN COARSENING .....*/

double load_rho_bar_s(int sten_type,double **x, int iunk,
                     int loc_inode, int inode_box, int izone,int *ijk_box,
                     int resid_only_flag)
{
  double resid_sum=0.0,resid,mat_val;
  int jzone_flag,junk;

  jzone_flag=FALSE;

  resid =-x[iunk][inode_box];
  resid_sum+=resid;
  mat_val=-1.0;
  dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
  dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);
 
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

  resid =-x[iunk][inode_box]; 
  resid_sum+=resid;
  mat_val=-1.0;
  dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
  dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);

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
/* pre_calc_dphi_drb_rb1: rho_bars are calculated at each wall & fluid node
   		original rosenfeld functionals. */

void pre_calc_dphi_drb_rb1(struct RB_Struct *dphi_drb, double **x)
{
 double rb0,rb1,rb2,rb3,rb1v[3],rb2v[3],inv_one_m_rb3,inv_one_m_rb3_sq,
        inv_one_m_rb3_3rd,DOT_rho12,DOT_rho22;
 double rb0_l=0.0, rb1_l=0.0, rb2_l=0.0, rb3_l=0.0;
 double rb0_r=0.0, rb1_r=0.0, rb2_r=0.0, rb3_r=0.0;
 int idim, icomp, loc_inode;
 int inode_box;
 int i,imax,loc_i,loc_iv,junk;
 double rho;

 for (inode_box=0; inode_box < Nnodes_box; inode_box++) {

       junk = Phys2Unk_first[RHOBAR_ROSEN];
       rb3 = x[junk][inode_box];
       rb2 = x[junk+1][inode_box];
       rb1 = x[junk+2][inode_box];
       rb0 = x[junk+3][inode_box];

       for (idim = 0; idim<Ndim; idim++) {
          rb2v[idim] = x[junk+Nrho_bar_s+idim][inode_box];
          rb1v[idim] = x[junk+Nrho_bar_s+Ndim+idim][inode_box];
       }

       inv_one_m_rb3 = 1.0 / (1.0 - rb3);
       inv_one_m_rb3_sq = inv_one_m_rb3*inv_one_m_rb3;
       inv_one_m_rb3_3rd = inv_one_m_rb3_sq*inv_one_m_rb3;

       dphi_drb[inode_box].S0 =  log(inv_one_m_rb3);
       dphi_drb[inode_box].S1 =  rb2*inv_one_m_rb3;
       dphi_drb[inode_box].S2 =  rb1*inv_one_m_rb3 +
                             rb2*rb2*inv_one_m_rb3_sq / (8.0*PI);
       dphi_drb[inode_box].S3 = rb0*inv_one_m_rb3 +
                            rb1*rb2 * inv_one_m_rb3_sq +
                            rb2*rb2*rb2*inv_one_m_rb3_3rd / (12.0*PI);
       DOT_rho12 = 0.0;
       DOT_rho22 = 0.0;
       for (idim = 0; idim < Ndim; idim++) {
           DOT_rho12 += rb1v[idim] * rb2v[idim];
           DOT_rho22 += rb2v[idim] * rb2v[idim];
   
           dphi_drb[inode_box].V1[idim] = -rb2v[idim]*inv_one_m_rb3;
           dphi_drb[inode_box].V2[idim] = -rb1v[idim]*inv_one_m_rb3 -
                                    rb2*rb2v[idim]*inv_one_m_rb3_sq/(4.0*PI);
       }

       dphi_drb[inode_box].S2 += -DOT_rho22*inv_one_m_rb3_sq/(8.0*PI);
       dphi_drb[inode_box].S3 += -DOT_rho12*inv_one_m_rb3_sq -
                                 rb2*DOT_rho22*inv_one_m_rb3_3rd/(4.0*PI);

  }

  return;
}
/*****************************************************************************/
/* pre_calc_dphi_drb_rb2: rho_bars are calculated at each wall & fluid node
        more recent FM functional with better 0D crossover.*/

void pre_calc_dphi_drb_rb2(struct RB_Struct *dphi_drb, double **x)
{
 double rb0,rb1,rb2,rb3,rb1v[3],rb2v[3],inv_one_m_rb3,inv_one_m_rb3_sq,
        inv_one_m_rb3_3rd,DOT_rho12,DOT_rho22;
 double rb0_l=0.0,rb1_l=0.0,rb2_l=0.0,rb3_l=0.0;
 double rb0_r=0.0,rb1_r=0.0,rb2_r=0.0,rb3_r=0.0;
 int idim, icomp, loc_inode;
 int inode_box;
 int i,imax,loc_i,loc_iv,junk;
 double rho,alpha,alpha_sq,alpha_cb,beta,gamma[3];

 for (inode_box=0; inode_box < Nnodes_box; inode_box++) {

       junk = Phys2Unk_first[RHOBAR_ROSEN];
       rb3 = x[junk][inode_box];
       rb2 = x[junk+1][inode_box];
       rb1 = x[junk+2][inode_box];
       rb0 = x[junk+3][inode_box];

       for (idim = 0; idim<Ndim; idim++) {
         rb2v[idim] = x[junk+Nrho_bar_s+idim][inode_box];
         rb1v[idim] = x[junk+Nrho_bar_s+Ndim+idim][inode_box];
       }

       DOT_rho12 = 0.0;
       DOT_rho22 = 0.0;
       for (idim = 0; idim < Ndim; idim++) {
           DOT_rho12 += rb1v[idim] * rb2v[idim];
           DOT_rho22 += rb2v[idim] * rb2v[idim];
       }

       inv_one_m_rb3 = 1.0 / (1.0 - rb3);
       inv_one_m_rb3_sq = inv_one_m_rb3*inv_one_m_rb3;
       inv_one_m_rb3_3rd = inv_one_m_rb3_sq*inv_one_m_rb3;

       /* same as all old rosenfeld functional contributions */
       dphi_drb[inode_box].S0 =  log(inv_one_m_rb3);
       dphi_drb[inode_box].S1 =  rb2*inv_one_m_rb3;
       dphi_drb[inode_box].S2 =  rb1*inv_one_m_rb3;
       dphi_drb[inode_box].S3 = rb0*inv_one_m_rb3 +
                                (rb1*rb2-DOT_rho12) * inv_one_m_rb3_sq;

       for (idim = 0; idim < Ndim; idim++) {
           dphi_drb[inode_box].V1[idim] = -rb2v[idim]*inv_one_m_rb3;
           dphi_drb[inode_box].V2[idim] = -rb1v[idim]*inv_one_m_rb3;
       }

       /* new rosenfeld functional contributions */
       if (rb2 > 1.e-15){
             alpha=rb2-DOT_rho22/rb2;
             beta = 1.0+DOT_rho22/(rb2*rb2);
             for (idim = 0; idim < Ndim; idim++){
                 gamma[idim] = rb2v[idim]/rb2;
             }
       }
       else{
           alpha=rb2;
           beta = 1.0;
           for (idim = 0; idim < Ndim; idim++) gamma[idim] = 0.0;
       }
       alpha_sq=alpha*alpha;
       alpha_cb=alpha_sq*alpha;
  
       dphi_drb[inode_box].S2 += alpha_sq*beta*inv_one_m_rb3_sq/(8.0*PI);
       dphi_drb[inode_box].S3 += alpha_cb*inv_one_m_rb3_3rd/(12.0*PI);

       for (idim = 0; idim < Ndim; idim++) 
           dphi_drb[inode_box].V2[idim] -= 
                inv_one_m_rb3_sq*alpha_sq*gamma[idim]/(4.0*PI);

  }

  return;
}
/*****************************************************************************/
/* d2phi_drb2_delta_rb:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                 for the dphi_drb that use Delta_Fn Stencils (all but S3) */

static struct RB_Struct d2phi_drb2_delta_rb(int junk, int jnode_box,double **x, 
					    double weight, int *offset, double *sign,
					    double inv_rad, double inv_4pir, 
					    double inv_4pirsq)

{
  struct RB_Struct tmp;
  double rb0, rb1, rb2, rb3, rb1v[NDIM_MAX], rb2v[NDIM_MAX];
  double inv_one_m_rb3, inv_one_m_rb3_sq, inv_one_m_rb3_3rd;
  double vector[NDIM_MAX];
  int idim,loc_js,loc_jv;

  rb3=x[junk][jnode_box];
  rb2=x[junk+1][jnode_box];
  rb1=x[junk+2][jnode_box];
  rb0=x[junk+3][jnode_box];
   
  for (idim = 0; idim<Ndim; idim++) {
     rb2v[idim] = sign[idim]*x[junk+Nrho_bar_s+idim][jnode_box];
     rb1v[idim] = sign[idim]*x[junk+Nrho_bar_s+Ndim+idim][jnode_box];
  }

  inv_one_m_rb3 = 1.0 / (1.0 - rb3);
  inv_one_m_rb3_sq = inv_one_m_rb3*inv_one_m_rb3;
  inv_one_m_rb3_3rd = inv_one_m_rb3_sq*inv_one_m_rb3;
  for (idim = 0; idim<Ndim; idim++) 
    vector[idim] = offset[idim] * Esize_x[idim] * inv_rad;

  tmp.S2 = (inv_one_m_rb3*inv_4pir
         + rb2*Inv_4pi*inv_one_m_rb3_sq)*weight;

  for (idim = 0; idim<Ndim; idim++)
     tmp.S2 += weight * vector[idim]
              * rb2v[idim]*Inv_4pi * inv_one_m_rb3_sq;

  tmp.S3 = weight * (inv_4pirsq*inv_one_m_rb3 
         + rb2*inv_4pir*inv_one_m_rb3_sq
         + rb1*inv_one_m_rb3_sq+rb2*rb2*Inv_4pi*inv_one_m_rb3_3rd);

  for (idim = 0; idim<Ndim; idim++)
     tmp.S3 += weight * inv_one_m_rb3_sq *
           ( -rb2v[idim]*rb2v[idim]*Inv_4pi*inv_one_m_rb3
           + vector[idim]
           * ( rb2v[idim]* inv_4pir + rb1v[idim]
              + rb2*rb2v[idim]*inv_one_m_rb3/(2*PI) ) );

  tmp.S0 = 0.0;  
  tmp.S1 = inv_one_m_rb3*weight;  

  for (idim = 0; idim<Ndim; idim++)
    tmp.V1[idim] = sign[idim]*(weight * inv_one_m_rb3) * vector[idim];
  for (idim = 0; idim<Ndim; idim++)
    tmp.V2[idim] = sign[idim]*weight *(- rb2v[idim]*Inv_4pi*inv_one_m_rb3_sq
                    - (-inv_4pir*inv_one_m_rb3 
                    - rb2* Inv_4pi*inv_one_m_rb3_sq) * vector[idim]);
  return(tmp);
}
/****************************************************************************/
/* d2phi_drb2_delta2_rb:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                 for the dphi_drb that use Delta_Fn Stencils (all but S3) */

static struct RB_Struct d2phi_drb2_delta2_rb(int junk, int jnode_box,double **x, 
					     double weight, int *offset, double *sign,
					     double inv_rad, double inv_4pir, 
					     double inv_4pirsq)
{
  struct RB_Struct tmp;
  double rb0, rb1, rb2, rb3, rb1v[NDIM_MAX], rb2v[NDIM_MAX];
  double inv_one_m_rb3, inv_one_m_rb3_sq, inv_one_m_rb3_3rd, inv_one_m_rb3_4th;
  double vector[NDIM_MAX];
  int idim,loc_js,loc_jv;
  double DOT_rho22,alpha,alpha_sq,alpha_cb,beta,gamma[3],DOT_gamma,eps;
  
  rb3 = x[junk][jnode_box];
  rb2 = x[junk+1][jnode_box];
  rb1 = x[junk+2][jnode_box];
  rb0 = x[junk+3][jnode_box];

  for (idim = 0; idim<Ndim; idim++) {
    rb2v[idim] = sign[idim]*x[junk+Nrho_bar_s+idim][jnode_box];
    rb1v[idim] = sign[idim]*x[junk+Nrho_bar_s+Ndim+idim][jnode_box];
  }

  DOT_rho22 = 0.0;
  for (idim = 0; idim < Ndim; idim++) {
      DOT_rho22 += rb2v[idim] * rb2v[idim];
  }
  if (rb2 > 1.e-15){
     alpha=rb2-DOT_rho22/rb2;
     beta = 1.0+DOT_rho22/(rb2*rb2);
     DOT_gamma=0.0;
     for (idim = 0; idim < Ndim; idim++){
          gamma[idim] = rb2v[idim]/rb2;
          DOT_gamma+=gamma[idim]*gamma[idim];
     }
     eps=alpha/rb2;
  }
  else{
      alpha=rb2;
      beta = 1.0;
      for (idim = 0; idim < Ndim; idim++) gamma[idim] = 0.0;
      DOT_gamma=0.0;
      eps=1.0;
  }
  alpha_sq=alpha*alpha;
  alpha_cb=alpha_sq*alpha;

  inv_one_m_rb3 = 1.0 / (1.0 - rb3);
  inv_one_m_rb3_sq = inv_one_m_rb3*inv_one_m_rb3;
  inv_one_m_rb3_3rd = inv_one_m_rb3_sq*inv_one_m_rb3;
  inv_one_m_rb3_4th = inv_one_m_rb3_sq*inv_one_m_rb3_sq;
  for (idim = 0; idim<Ndim; idim++) 
    vector[idim] = offset[idim] * Esize_x[idim] * inv_rad;

  tmp.S2 = inv_one_m_rb3*inv_4pir*weight;      /*old rosen rb2-rb1*/

  tmp.S2 += weight*Inv_4pi*inv_one_m_rb3_sq*     /* new rosen rb2-rb2*/
            alpha*(beta*beta-eps*DOT_gamma);

  for (idim = 0; idim<Ndim; idim++)              /*new rosen rb2-rb2v*/
      tmp.S2 -= weight*Inv_4pi*inv_one_m_rb3_sq*
                alpha*(eps -2*beta)*gamma[idim]*vector[idim];

  tmp.S3 = weight * (inv_4pirsq*inv_one_m_rb3   /* old rosen*/
         + rb2*inv_4pir*inv_one_m_rb3_sq
         + rb1*inv_one_m_rb3_sq);

  for (idim = 0; idim<Ndim; idim++)             /* old rosen rb3-rbv1/v2*/
     tmp.S3 += weight * inv_one_m_rb3_sq *
               vector[idim] * (rb2v[idim]*inv_4pir + rb1v[idim]);
  
  tmp.S3 += weight*Inv_4pi*inv_one_m_rb3_3rd*alpha_sq*beta; /*new rosen rb3-rb2*/

  for (idim = 0; idim<Ndim; idim++)
     tmp.S3 += weight * 2.0*Inv_4pi*inv_one_m_rb3_3rd *    /*new rosen rb3-rb2v*/
               alpha_sq*gamma[idim]*vector[idim];

  tmp.S0 = 0.0;                     /*old rosen rb0-(rb0-rb2,rbv1,rbv2)*/
  tmp.S1 = inv_one_m_rb3*weight;    /*old rosen rb1-rb2*/

  for (idim = 0; idim<Ndim; idim++){                               /*old rosen*/
    tmp.V1[idim] = sign[idim]*(weight*inv_one_m_rb3)*vector[idim];     /*rb1v-rb2v*/
    tmp.V2[idim] = sign[idim]*(weight*inv_4pir*inv_one_m_rb3)*vector[idim]; /*rb2v-rb1v*/
  }

  for (idim = 0; idim<Ndim; idim++){                         
    tmp.V2[idim] += sign[idim]*weight*Inv_4pi*inv_one_m_rb3_sq*  /*new rb2v-rb2*/
                    alpha*(eps-2.0*beta)*gamma[idim];
                                                                
    tmp.V2[idim] -= sign[idim]*weight*Inv_4pi*inv_one_m_rb3_sq*  /*new rb2v-rb2v*/
                    alpha*(4.0*DOT_gamma-eps)*vector[idim];
  }
  return(tmp);
}
/****************************************************************************/
/* d2phi_drb2_theta_rb:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                    for the dphi_drb that use Theta_Fn Stencils (S3)      */

static struct RB_Struct d2phi_drb2_theta_rb(int junk, int jnode_box,double **x,double weight,
					    int *offset)
{
  struct RB_Struct tmp;
  double rb0, rb1, rb2, rb3, rb1v[NDIM_MAX], rb2v[NDIM_MAX];
  double inv_one_m_rb3, inv_one_m_rb3_sq, inv_one_m_rb3_3rd, inv_one_m_rb3_4th;
  int idim,loc_js,loc_jv;


  rb3 = x[junk][jnode_box];
  rb2 = x[junk+1][jnode_box];
  rb1 = x[junk+2][jnode_box];
  rb0 = x[junk+3][jnode_box];
  
  for (idim = 0; idim<Ndim; idim++) {
    rb2v[idim] = x[junk+Nrho_bar_s+idim][jnode_box];   
    rb1v[idim] = x[junk+Nrho_bar_s+Ndim+idim][jnode_box];   
  }

  inv_one_m_rb3 = 1.0 / (1.0 - rb3);
  inv_one_m_rb3_sq = inv_one_m_rb3*inv_one_m_rb3;
  inv_one_m_rb3_3rd = inv_one_m_rb3_sq*inv_one_m_rb3;
  inv_one_m_rb3_4th = inv_one_m_rb3_3rd*inv_one_m_rb3;

  tmp.S2 = weight * (rb1*inv_one_m_rb3_sq 
         + rb2*rb2*Inv_4pi*inv_one_m_rb3_3rd);
  for (idim = 0; idim<Ndim; idim++)
     tmp.S2 += - weight * rb2v[idim]*rb2v[idim]
                         *Inv_4pi*inv_one_m_rb3_3rd;
  tmp.S3 = weight * (rb0*inv_one_m_rb3_sq 
                     + 2*rb1*rb2*inv_one_m_rb3_3rd
                     + rb2*rb2*rb2*Inv_4pi*inv_one_m_rb3_4th);
  for (idim = 0; idim<Ndim; idim++)
     tmp.S3 += -weight * (2*rb1v[idim]*rb2v[idim]*inv_one_m_rb3_3rd
            +3*rb2*rb2v[idim]*rb2v[idim]*Inv_4pi*inv_one_m_rb3_4th);

  tmp.S0 = weight * inv_one_m_rb3;             
  tmp.S1 = weight * rb2 * inv_one_m_rb3_sq;    

  for (idim = 0; idim<Ndim; idim++)
      tmp.V1[idim] = - weight * rb2v[idim]*inv_one_m_rb3_sq;
  for (idim = 0; idim<Ndim; idim++)
      tmp.V2[idim] = - weight * (rb1v[idim]*inv_one_m_rb3_sq
                     + rb2*rb2v[idim]*inv_one_m_rb3_3rd/(2.0*PI) );
  return (tmp);
}
/****************************************************************************/
/* d2phi_drb2_theta2_rb:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                    for the dphi_drb that use Theta_Fn Stencils (S3)      */

static struct RB_Struct d2phi_drb2_theta2_rb(int junk, int jnode_box,double **x,double weight,
					     int *offset)
{
  struct RB_Struct tmp;
  double rb0, rb1, rb2, rb3, rb1v[NDIM_MAX], rb2v[NDIM_MAX];
  double inv_one_m_rb3, inv_one_m_rb3_sq, inv_one_m_rb3_3rd, inv_one_m_rb3_4th;
  int idim,loc_js,loc_jv;
  double DOT_rho22,DOT_rho12,alpha,alpha_sq,alpha_cb,beta,gamma[3];

  rb3 = x[junk][jnode_box];
  rb2 = x[junk+1][jnode_box];
  rb1 = x[junk+2][jnode_box];
  rb0 = x[junk+3][jnode_box];
  
  for (idim = 0; idim<Ndim; idim++) {
    rb2v[idim] = x[junk+Nrho_bar_s+idim][jnode_box];   
    rb1v[idim] = x[junk+Nrho_bar_s+Ndim+idim][jnode_box];   
  }

  inv_one_m_rb3 = 1.0 / (1.0 - rb3);
  inv_one_m_rb3_sq = inv_one_m_rb3*inv_one_m_rb3;
  inv_one_m_rb3_3rd = inv_one_m_rb3_sq*inv_one_m_rb3;
  inv_one_m_rb3_4th = inv_one_m_rb3_3rd*inv_one_m_rb3;

  DOT_rho22 = 0.0;
  DOT_rho12 = 0.0;
  for (idim = 0; idim < Ndim; idim++) {
      DOT_rho22 += rb2v[idim] * rb2v[idim];
      DOT_rho12 += rb1v[idim] * rb2v[idim];
  }

  if (rb2 > 1.e-15){
        alpha=rb2-DOT_rho22/rb2;
        beta = 1.0+DOT_rho22/(rb2*rb2);
        for (idim = 0; idim < Ndim; idim++){
            gamma[idim] = rb2v[idim]/rb2;
        }
  }
  else{
      alpha=rb2;
      beta = 1.0;
      for (idim = 0; idim < Ndim; idim++) gamma[idim] = 0.0;
  }
  alpha_sq=alpha*alpha;
  alpha_cb=alpha_sq*alpha;

  tmp.S2 = weight *rb1*inv_one_m_rb3_sq;                  /*old rb2-rb3*/
  tmp.S2 += weight * Inv_4pi*inv_one_m_rb3_3rd*alpha_sq*beta;      /*new rb2-rb3*/

  tmp.S3 = weight *(rb0*inv_one_m_rb3_sq                   /*old rb3-rb3*/
                   + 2.0*(rb1*rb2 - DOT_rho12)*inv_one_m_rb3_3rd);
  tmp.S3 += weight*Inv_4pi*inv_one_m_rb3_4th*alpha_cb;    /*new rb3-rb3*/

  tmp.S0 = weight * inv_one_m_rb3;            /*old rb0-rb3*/
  tmp.S1 = weight * rb2 * inv_one_m_rb3_sq;   /*old rb1-rb3*/

  for (idim = 0; idim<Ndim; idim++){         
      tmp.V1[idim] = - weight * rb2v[idim]*inv_one_m_rb3_sq;  /*old rbv1-rb3*/
      tmp.V2[idim] = - weight * rb1v[idim]*inv_one_m_rb3_sq;  /*old rbv2-rb3*/
  }

  for (idim = 0; idim<Ndim; idim++)                            /*new rbv2-rb3*/
      tmp.V2[idim] -= weight * 2.0*Inv_4pi*inv_one_m_rb3_3rd*alpha_sq*gamma[idim];

  return (tmp);
}
/****************************************************************************/
