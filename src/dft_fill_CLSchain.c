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

/* ---------------------------------------------------------
dft_fill_CLSchain.c:

Here are the generic fill routines for bonded systems.
These routines use function pointers to specify which functional
is being loaded -
------------------------------------------------------------*/

#include "dft_fill_CLSchain.h"

/****************************************************************************/
double resid_and_Jac_ChainDensity (int func_type, double **x, int iunk, int unk_B,
  int loc_inode, int inode_box, int resid_only_flag, double (*fp_prefactor)(int,int,double **))
{
  int i,loop_start,loop_end,itype_mer,npol,iseg,unk_GQ,unk_GQ_test,iref;
  int boltz_pow,boltz_pow_J,jbond,ibond,unkIndex[2],numEntries,unk_GQ_j,unk_GQ_j_test;
  int loc_jnode,jnode_box,ijk_box[3];
  int unk_xi2,unk_xi3;
  int loc_inode2,iwall,pwall,pwall_type,inode_w;
  int jsurf_node,icomp,icomp_graft,ilist;
  double nodepos[3],nodeposc[3],prefac;
  double fac1,fac2,mat_val,resid=0.0,resid_sum=0.0,resid_sum2=0.0,values[2];
  int inodel, inode_boxl,izone,graft_seg,graft_bond,gbond,jtmp,idim;
  double y,ysqrt,xi_2,xi_3,dummy=0.0;


  if (Lconstrain_interface && Type_interface==PHASE_INTERFACE && B2G_node[inode_box]==(int)(0.5*Size_x[Grad_dim]/Esize_x[Grad_dim]) && 
           iunk==Phys2Unk_first[DENSITY]){
         fill_constant_density_chain(iunk,0,0,x[unk_B][inode_box],loc_inode,inode_box,x,resid_only_flag);
  }
  else{

  if (resid_only_flag !=INIT_GUESS_FLAG){
     icomp=iunk-Phys2Unk_first[DENSITY];
     if (Grafted_Logical==FALSE || icomp != Grafted_TypeID[Icomp_to_polID[icomp]]) {
        resid = x[iunk][inode_box]*x[unk_B][inode_box];
        resid_sum=resid;
        if (resid_only_flag != CALC_RESID_ONLY)  dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
        if(resid_only_flag==FALSE){
            values[0]=x[iunk][inode_box]; values[1]=x[unk_B][inode_box];
            unkIndex[0]=unk_B; unkIndex[1]=iunk;
            numEntries=2;
         
            if (Iwrite_files==FILES_DEBUG_MATRIX) {
               for (jtmp=0;jtmp<numEntries;jtmp++) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[unkIndex[jtmp]]*Nnodes]+=values[jtmp];
            }
            dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                 unkIndex, inode_box, values, numEntries);
        }
      }
      else{
        resid = x[iunk][inode_box];
        resid_sum=resid;
        if (resid_only_flag != CALC_RESID_ONLY)  dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
        if(resid_only_flag==FALSE){
            values[0]=1.0;
            unkIndex[0]=iunk;
            numEntries=1;
         
            if (Iwrite_files==FILES_DEBUG_MATRIX) {
               for (jtmp=0;jtmp<numEntries;jtmp++) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[unkIndex[jtmp]]*Nnodes]+=values[jtmp];
            }
            dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                 unkIndex, inode_box, values, numEntries);
        }
      }
  }

  if (Type_poly==WJDC || Type_poly==WJDC2){
     itype_mer=Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
     loop_start=iunk;   /* don't sum over segments of a given type for WJDC or WJDC2 where segments rhos are explicit */
     loop_end=iunk+1;
  }
  else if (Type_poly==CMS || Type_poly==CMS_SCFT || Type_poly==WJDC3){
     itype_mer=iunk-Phys2Unk_first[DENSITY];
     npol = 0;
     while (Nmer_t[npol][itype_mer]==0){ npol++; }

                           /* for CMS or WJDC3 polymers where we are concerned with component densitie,  
                           we need to perform a sum over segment of a given type */
     loop_start=0;      
     loop_end=Nmer[npol];
  }

  for (i=loop_start; i<loop_end; i++){
     if (Lseg_densities || (Type_mer[npol][i] == itype_mer) ) {

        if (!Lseg_densities) iseg=SegChain2SegAll[npol][i];
        else                 iseg=iunk-Phys2Unk_first[DENSITY];

        boltz_pow = -(Nbonds_SegAll[iseg]-2);
        boltz_pow_J = -(Nbonds_SegAll[iseg]-1);

        if (fp_prefactor!=NULL) {
            if (Type_poly==CMS || Type_poly==CMS_SCFT) iref=itype_mer;
            if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3) iref=iseg;
            fac1 = (*fp_prefactor)(iref,inode_box,x);
	    if (Type_poly==CMS_SCFT) fac1 /= Gsum[npol];
        }
        else  fac1=1.0;

        if (Grafted_Logical==FALSE || Grafted[npol]==FALSE ||iseg !=Grafted_SegIDAll[npol]) {
        for (ibond=0; ibond<Nbonds_SegAll[iseg]; ibond++) {
            unk_GQ  = Phys2Unk_first[func_type] + Poly_to_Unk_SegAll[iseg][ibond];
            unk_GQ_test = unk_GQ-Phys2Unk_first[func_type];
            if (Pol_Sym[unk_GQ_test] != -1) unk_GQ=Pol_Sym[unk_GQ_test] + Phys2Unk_first[func_type];
            fac1 *= x[unk_GQ][inode_box];

            if(resid_only_flag==FALSE){   
               if (fp_prefactor!=NULL) {
                   if (Type_poly==CMS || Type_poly==CMS_SCFT) iref=itype_mer;
                   if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3) iref=iseg;
                   fac2 = (*fp_prefactor)(iref,inode_box,x);
                   if (Type_poly==CMS_SCFT) fac2 /=Gsum[npol];
               }
               else                     fac2=1.0;

               for (jbond=0; jbond<Nbonds_SegAll[iseg]; jbond++) {
                  if (jbond != ibond){
                     unk_GQ_j  = Phys2Unk_first[func_type]+Poly_to_Unk_SegAll[iseg][jbond];
                     unk_GQ_j_test=unk_GQ_j-Phys2Unk_first[func_type];
                     if (Pol_Sym[unk_GQ_j_test] != -1) unk_GQ_j=Pol_Sym[unk_GQ_j_test] + Phys2Unk_first[func_type];
                     fac2 *= x[unk_GQ_j][inode_box];
                   }
               }
               if (-boltz_pow > 0) mat_val = -fac2*POW_DOUBLE_INT(x[unk_B][inode_box],boltz_pow);
               else mat_val=-fac2;
               if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[unk_GQ]*Nnodes]+=mat_val;
               dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_GQ,inode_box,mat_val);
             }
        } /* end loop over ibond */			

        if(resid_only_flag==FALSE){
           if (-boltz_pow > 0){
              mat_val = -fac1*((double) boltz_pow)*POW_DOUBLE_INT(x[unk_B][inode_box],boltz_pow_J);
              if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[unk_B]*Nnodes]+=mat_val;
              dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_B,inode_box,mat_val);
           }
         }
         if (-boltz_pow >0){
               resid = -fac1*POW_DOUBLE_INT(x[unk_B][inode_box],boltz_pow); 
         }
         else{
              resid = -fac1;
         }

         if(resid_only_flag==FALSE && Grafted_Logical==TRUE && Grafted[npol]!=FALSE){    /* one more Jacobian entry to account for G_0^1 terms in grafted prefactor */
           icomp_graft=Grafted_TypeID[npol];
           for (iwall=0;iwall<Nwall;iwall++){
             if (WallType[iwall]==Graft_wall[npol]){
               for (jsurf_node=0;jsurf_node<Nodes_Surf_Gsum[npol][iwall];jsurf_node++){
                 jnode_box=Index_SurfNodes_Gsum[iwall][jsurf_node];

                 if (Nbonds_SegAll[Grafted_SegIDAll[npol]]>2){
                    unk_B=Index_UnkB_Gsum[iwall][jsurf_node];
                    if (Grafted[npol]==GRAFT_DENSITY) mat_val=resid*(GsumPrefac_XiDerivs[iwall][jsurf_node]/(Total_area_graft[npol]*Gsum_graft[npol]));
                    else                              mat_val=resid*(GsumPrefac_XiDerivs[iwall][jsurf_node]/Gsum_graft[npol]);
           
                    if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node_extra[jnode_box]+Solver_Unk[unk_B]*Nnodes]+=mat_val;
                    dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_B,jnode_box,mat_val);
                 }

                 for (jbond=0;jbond<Nbonds_SegAll[Grafted_SegIDAll[npol]];jbond++){
                    if (Bonds[npol][Grafted_SegID[npol]][jbond] != -1){
                      unk_GQ=Index_UnkGQ_Gsum[iwall][jsurf_node][jbond];
                      if (Grafted[npol]==GRAFT_DENSITY) mat_val=-resid*(GsumPrefac_GDerivs[iwall][jsurf_node][jbond]/(Total_area_graft[npol]*Gsum_graft[npol]));
                      else                              mat_val=-resid*(GsumPrefac_GDerivs[iwall][jsurf_node][jbond]/Gsum_graft[npol]);
           
                      if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node_extra[jnode_box]+Solver_Unk[unk_GQ]*Nnodes]+=mat_val;
                      dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_GQ,jnode_box,mat_val);
                    }
                 }
               }
              }
            }
         }
         resid_sum+=resid;
         resid_sum2+=resid;
         if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) 
                  dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
         }
         else{    /* end segments for grafted chains */
              resid = -fac1;
              resid_sum+=resid;
              if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) 
                  dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
         }
      } /* end if statement */
   } /* end loop over length of chain */

   if (resid_only_flag==INIT_GUESS_FLAG && (Grafted[npol]==FALSE ||iseg !=Grafted_SegIDAll[npol])) resid_sum /= x[unk_B][inode_box];

   /* extra term in Jacobian for SCF */
   if(Type_poly==CMS_SCFT && resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) {
      iseg = Nmer[npol]-1;
      unk_GQ = Geqn_start[npol] + Poly_to_Unk[npol][iseg][0];
      for(loc_jnode=0; loc_jnode<Nnodes_per_proc; loc_jnode++) {
         jnode_box = L2B_node[loc_jnode];
         mat_val = resid_sum*x[unk_GQ][jnode_box];
         if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[jnode_box]+Solver_Unk[unk_GQ]*Nnodes]+=mat_val;
         dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_GQ,jnode_box,mat_val);
      }
    }
	
    }
    return(resid_sum);
}
/****************************************************************************/
double resid_and_Jac_ChainDensity_WJDC2 (int func_type, double **x, int iunk, int unk_B,
  int loc_inode, int inode_box, int resid_only_flag, double (*fp_prefactor)(int))
{
  int itype_mer,iseg,unk_GQ,unk_GQ_test,jtmp;
  int jbond,ibond,unkIndex[2],numEntries,unk_GQ_j,unk_GQ_j_test;
  double fac1,fac2,mat_val,resid=0.0,resid_sum=0.0,values[2],fac1deriv;

  iseg=iunk-Phys2Unk_first[DENSITY];
  itype_mer=Unk2Comp[iunk-Phys2Unk_first[DENSITY]];

  /* compute bulk and local products of G functions associated with this segment */
  fac1=1.0;
  fac2=1.0;
  for (ibond=0; ibond<Nbonds_SegAll[iseg]; ibond++) {
       unk_GQ  = Phys2Unk_first[func_type] + Poly_to_Unk_SegAll[iseg][ibond];
       unk_GQ_test = unk_GQ-Phys2Unk_first[func_type];
       if (Pol_Sym[unk_GQ_test] != -1) unk_GQ=Pol_Sym[unk_GQ_test] + Phys2Unk_first[func_type];
       fac1 *= x[unk_GQ][inode_box];
       fac2 *= G_WJDC_b[unk_GQ_test];
  }

  resid = x[iunk][inode_box]*x[unk_B][inode_box]-Rho_seg_b[iseg]*Field_WJDC_b[itype_mer]*fac1/fac2;
  resid_sum=resid;
  dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
  if(!resid_only_flag){
      values[0]=x[iunk][inode_box]; 
      values[1]=x[unk_B][inode_box];
      unkIndex[0]=unk_B; 
      unkIndex[1]=iunk;
      numEntries=2;
      if (Iwrite_files==FILES_DEBUG_MATRIX){
          for (jtmp=0;jtmp<numEntries;jtmp++) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[unkIndex[jtmp]]*Nnodes]+=values[jtmp];
      }
      dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                              unkIndex, inode_box, values, numEntries);

     fac1deriv=1.0;
     for (ibond=0; ibond<Nbonds_SegAll[iseg]; ibond++) {
        unk_GQ  = Phys2Unk_first[func_type] + Poly_to_Unk_SegAll[iseg][ibond];
        for (jbond=0; jbond<Nbonds_SegAll[iseg]; jbond++) {
          if (jbond != ibond){
          unk_GQ_j  = Phys2Unk_first[func_type] + Poly_to_Unk_SegAll[iseg][ibond];
          unk_GQ_j_test = unk_GQ_j-Phys2Unk_first[func_type];
          if (Pol_Sym[unk_GQ_j_test] != -1) unk_GQ_j=Pol_Sym[unk_GQ_j_test] + Phys2Unk_first[func_type];
          fac1deriv *= x[unk_GQ_j][inode_box];
          }
        }
        mat_val=Rho_seg_b[iseg]*Field_WJDC_b[itype_mer]*fac1deriv/fac2;
        if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[unk_GQ]*Nnodes]+=mat_val;
        dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_GQ,inode_box,mat_val);
     }
  }
  return(resid_sum);
}
/****************************************************************************/
double load_Chain_Geqns(int func_type_field,int Njacobian_types, int Njacobian_sums,
       void (*funcArray_Jac[3])(int,int,int,int,int,int,int,int,int *,double,double **),
       double (*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **),
       double (*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **),
       int iunk, int loc_inode, int inode_box, 
       int *ijk_box, int izone, double **x,
       int resid_only_flag)
{
    int unk_GQ,unk_B,pol_num,seg_num,bond_num,junk,unkIndex[4],numEntries,sten;
    int itype_mer,jseg,jtype_mer,unk_xi2,unk_xi3,inode,jtmp;
	int inode_box2,inode2;
	int iwall,iwall_type;
    double nodepos[3],nodepos2[3],xbound,xbound_lim;
    double resid_G=0.0,resid,values[4],gint_tmp; 
 /* double mat_val;*/
    double y,ysqrt,dydxi,xi_2,xi_3,xi_2_0,xi_3_0;
	double y2,ysqrt2;
	   
       /* use arrays that are indexed with the geqns ordered starting at 0 */ 
    unk_GQ=iunk-Phys2Unk_first[G_CHAIN];
    pol_num = Unk_to_Poly[unk_GQ];
    seg_num = Unk_to_Seg[unk_GQ];
    bond_num = Unk_to_Bond[unk_GQ];
    itype_mer = Type_mer[pol_num][seg_num];                   /*itype_mer=icomp*/

     /* from iunk and the bond number find the jseg to which we are looking */

    jseg = Bonds[pol_num][seg_num][bond_num];
    if (jseg != -1 && jseg != -2) {jtype_mer = Type_mer[pol_num][jseg];}  /* ALF: fix*/
    else {jtype_mer=itype_mer;}

    resid_G=0.0;

    sten=DELTA_FN_BOND;

    if (Zero_density_TF[inode_box][itype_mer] || Vext[loc_inode][itype_mer] == VEXT_MAX) {
	resid_G=fill_zero_value(iunk,loc_inode,inode_box,x,resid_only_flag);
    }
    else {
		
        if (Pol_Sym[unk_GQ] != -1 && resid_only_flag != INIT_GUESS_FLAG){  /* FILL IN G's FOR SYMMETRIC BONDS */
            junk = Pol_Sym[unk_GQ] + Phys2Unk_first[G_CHAIN];
            resid = x[iunk][inode_box]-x[junk][inode_box];
            resid_G+=resid;
            if (resid_only_flag !=CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
            if (!resid_only_flag){
                unkIndex[0]=iunk; unkIndex[1]=junk;
                values[0]=1.0; values[1]=-1.0;
                numEntries=2;
                if (Iwrite_files==FILES_DEBUG_MATRIX){
                    for (jtmp=0;jtmp<numEntries;jtmp++) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[unkIndex[jtmp]]*Nnodes]+=values[jtmp];
                }
		dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode, unkIndex, inode_box, values, numEntries);
             }
        }
	else {          /* FILL IN G's FOR UNIQUE BONDS */
	   if (Bonds[pol_num][seg_num][bond_num] == -1){        /* fill end segment equation */
		/* Boltz unk for 1st seg */
		if (Type_poly==CMS || Type_poly==CMS_SCFT || Type_poly==WJDC2 || Type_poly==WJDC3) 
			unk_B = Phys2Unk_first[func_type_field] + itype_mer;
		else if (Type_poly==WJDC) unk_B = Phys2Unk_first[func_type_field] + SegChain2SegAll[pol_num][seg_num]; 
		
		if (Type_poly==CMS || Type_poly==CMS_SCFT){
			if (resid_only_flag==INIT_GUESS_FLAG) resid=-1.0;
			else  resid = x[iunk][inode_box]/x[unk_B][inode_box]-1.0;
			resid_G+=resid;
			if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) 
				dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
			if (resid_only_flag==FALSE){
				unkIndex[0]=iunk; unkIndex[1]=unk_B;
				values[0]=1.0/x[unk_B][inode_box]; 
				values[1]=-x[iunk][inode_box]/(x[unk_B][inode_box]*x[unk_B][inode_box]);
				numEntries=2;
                                if (Iwrite_files==FILES_DEBUG_MATRIX){
                                    for (jtmp=0;jtmp<numEntries;jtmp++) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[unkIndex[jtmp]]*Nnodes]+=values[jtmp];
                                }
				dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
									  unkIndex, inode_box, values, numEntries);
			}
		}
		else if (Type_poly==WJDC || Type_poly==WJDC2  || Type_poly==WJDC3){
			unk_xi2=Phys2Unk_first[CAVWTC]; unk_xi3=Phys2Unk_first[CAVWTC]+1;
			xi_2=x[unk_xi2][inode_box]; xi_3=x[unk_xi3][inode_box];
			y=y_cav(Sigma_ff[itype_mer][itype_mer],Sigma_ff[jtype_mer][jtype_mer],xi_2,xi_3);
			ysqrt=sqrt(y);
			if (resid_only_flag==INIT_GUESS_FLAG){
                              resid = -(1.0/ysqrt);
                        }
			else resid = x[iunk][inode_box]/(x[unk_B][inode_box]*ysqrt)-(1.0/ysqrt);
			resid_G+=resid;
			if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY)dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
					
			if (resid_only_flag==FALSE){
				unkIndex[0]=iunk; unkIndex[1]=unk_B; unkIndex[2]=unk_xi2; unkIndex[3]=unk_xi3;
				values[0]=1.0/(ysqrt*x[unk_B][inode_box]); 
				values[1]=-x[iunk][inode_box]/(ysqrt*x[unk_B][inode_box]*x[unk_B][inode_box]);
				
				dydxi=dy_dxi2_cav(Sigma_ff[itype_mer][itype_mer],Sigma_ff[jtype_mer][jtype_mer],xi_2,xi_3);
				values[2]=-(dydxi/(y*ysqrt))*((x[iunk][inode_box]/x[unk_B][inode_box]) -1.0);
				dydxi=dy_dxi3_cav(Sigma_ff[itype_mer][itype_mer],Sigma_ff[jtype_mer][jtype_mer],xi_2,xi_3);
				values[3]=-(dydxi/(y*ysqrt))*((x[iunk][inode_box]/x[unk_B][inode_box]) -1.0);
				numEntries=4;
                                if (Iwrite_files==FILES_DEBUG_MATRIX){
                                    for (jtmp=0;jtmp<numEntries;jtmp++) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[unkIndex[jtmp]]*Nnodes]+=values[jtmp];
                                }
				dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
									  unkIndex, inode_box, values, numEntries);
			}
		}
	   }
	   else{                                            /* fill G_seg eqns */

	   /* First calculate the residual contributions */

		if (Type_poly==CMS || Type_poly==CMS_SCFT || Type_poly==WJDC2 || Type_poly==WJDC3) 
			unk_B = Phys2Unk_first[func_type_field] + itype_mer;   /* Boltz unk for this seg */
		else if (Type_poly==WJDC) unk_B = Phys2Unk_first[func_type_field] + SegChain2SegAll[pol_num][seg_num];  

		if (Type_poly==CMS || Type_poly==CMS_SCFT){
		    if (resid_only_flag != INIT_GUESS_FLAG){
			resid = x[iunk][inode_box]/x[unk_B][inode_box];

			resid_G+=resid;
			if (resid_only_flag !=CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
		    }
			
		    if (resid_only_flag==FALSE){
		        unkIndex[0]=iunk; unkIndex[1]=unk_B;
                        values[0]=1.0/x[unk_B][inode_box]; 
		        values[1]=-x[iunk][inode_box]/(x[unk_B][inode_box]*x[unk_B][inode_box]);
		        numEntries=2;
                        if (Iwrite_files==FILES_DEBUG_MATRIX){
                            for (jtmp=0;jtmp<numEntries;jtmp++) 
                                Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[unkIndex[jtmp]]*Nnodes]+=values[jtmp];
                        }
		        dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
	   						  unkIndex, inode_box, values, numEntries);
		     }
                 }
                 else if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){
                     unk_xi2=Phys2Unk_first[CAVWTC]; 
                     unk_xi3=Phys2Unk_first[CAVWTC]+1;
                     xi_2=x[unk_xi2][inode_box]; xi_3=x[unk_xi3][inode_box];
                     y=y_cav(Sigma_ff[itype_mer][itype_mer],Sigma_ff[jtype_mer][jtype_mer],xi_2,xi_3);
                     ysqrt=sqrt(y);
                     if (resid_only_flag != INIT_GUESS_FLAG){
                        resid = x[iunk][inode_box]/(x[unk_B][inode_box]*ysqrt);
                        resid_G+=resid;
                        if (resid_only_flag !=CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                     }

                     if (resid_only_flag==FALSE){
                        unkIndex[0]=iunk; unkIndex[1]=unk_B; unkIndex[2]=unk_xi2; unkIndex[3]=unk_xi3;
                        values[0]=1.0/(ysqrt*x[unk_B][inode_box]); 
                        values[1]=-x[iunk][inode_box]/(ysqrt*x[unk_B][inode_box]*x[unk_B][inode_box]);
                        dydxi=dy_dxi2_cav(Sigma_ff[itype_mer][itype_mer],Sigma_ff[jtype_mer][jtype_mer],xi_2,xi_3);
                        values[2]=-0.5*x[iunk][inode_box]*dydxi/(y*ysqrt*x[unk_B][inode_box]);
                        dydxi=dy_dxi3_cav(Sigma_ff[itype_mer][itype_mer],Sigma_ff[jtype_mer][jtype_mer],xi_2,xi_3);
                        values[3]=-0.5*x[iunk][inode_box]*dydxi/(y*ysqrt*x[unk_B][inode_box]);
                        numEntries=4; 
                        if (Iwrite_files==FILES_DEBUG_MATRIX){
                            for (jtmp=0;jtmp<numEntries;jtmp++) 
                                Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[unkIndex[jtmp]]*Nnodes]+=values[jtmp];
                        }
                        dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                 unkIndex, inode_box, values, numEntries);
                     }
                 }

            /* Now Finish loading the Jacobian... */
                 gint_tmp = load_polymer_recursion(sten,func_type_field,
                            Njacobian_types,Njacobian_sums,
                            funcArray_Jac, fp_ResidG,fp_ResidG_Bulk,
                            iunk,loc_inode,inode_box,unk_B,
                            itype_mer,izone,ijk_box,x, resid_only_flag);
                 resid_G += gint_tmp;
         }
         if (resid_only_flag==INIT_GUESS_FLAG){
            if (Type_poly==CMS || Type_poly==CMS_SCFT)                        resid_G*=(-x[unk_B][inode_box]);
            else if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){
                   resid_G*=(-x[unk_B][inode_box]*ysqrt);
            }
         }

       } /*end of fill unique G bonds */
    }  /* end of fill something */
    return(resid_G);
}   
/****************************************************************************/
/* G eqtns for SCF theory--very similar to CMS */
/* currently implementing only for polymer brushes grafted to the x=0 plane */
double load_Chain_Geqns_SCF(int func_type_field,int Njacobian_types, int Njacobian_sums,
       void (*funcArray_Jac[3])(int,int,int,int,int,int,int,int,int *,double,double **),
       double (*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **),
       double (*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **),
       int iunk, int loc_inode, int inode_box, 
       int *ijk_box, int izone, double **x,
       int resid_only_flag)
{
    int unk_GQ,unk_B,pol_num,seg_num,bond_num,junk,unkIndex[4],numEntries,sten;
    int itype_mer,jseg,inode,jtmp;
    double nodepos[3],xbound;
    double resid_G=0.0,resid,mat_val,values[4],gint_tmp; 
   
       /* use arrays that are indexed with the geqns ordered starting at 0 */ 
    unk_GQ=iunk-Phys2Unk_first[G_CHAIN];
    pol_num = Unk_to_Poly[unk_GQ];
    seg_num = Unk_to_Seg[unk_GQ];
    bond_num = Unk_to_Bond[unk_GQ];
    itype_mer = Type_mer[pol_num][seg_num];                   /*itype_mer=icomp*/

     /* from iunk and the bond number find the jseg to which we are looking */

    jseg = Bonds[pol_num][seg_num][bond_num];

    resid_G=0.0;

    sten=DELTA_FN_BOND;
	
	inode = L2G_node[loc_inode];
	node_to_position(inode,nodepos);
	xbound = WallPos[0][0] + Poly_graft_dist[0];

    if (Zero_density_TF[inode_box][itype_mer] || Vext[loc_inode][itype_mer] == VEXT_MAX) {
		resid_G=fill_zero_value(iunk,loc_inode,inode_box,x,resid_only_flag);
    }
/*    else if (fabs(x[iunk][inode_box]) < 1.e-12  && resid_only_flag==FALSE){    make G stationary 
        resid = 0.0;
        if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
        if (!resid_only_flag){
          mat_val =1.0; 
          if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
          dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);
        }
    }*/
	else {
			/* need to modify for SCF */
		if (Pol_Sym[unk_GQ] != -1 && resid_only_flag != INIT_GUESS_FLAG){  /* FILL IN G's FOR SYMMETRIC BONDS */
			junk = Pol_Sym[unk_GQ] + Phys2Unk_first[G_CHAIN];
			resid = x[iunk][inode_box]-x[junk][inode_box];
			resid_G+=resid;
			if (resid_only_flag !=CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
			if (!resid_only_flag){
				unkIndex[0]=iunk; unkIndex[1]=junk;
				values[0]=1.0; values[1]=-1.0;
				numEntries=2;
                                if (Iwrite_files==FILES_DEBUG_MATRIX){
                                    for (jtmp=0;jtmp<numEntries;jtmp++) 
                                        Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[unkIndex[jtmp]]*Nnodes]+=values[jtmp];
                                }
				dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode, unkIndex, inode_box, values, numEntries);
			}
		}
		else {          /* FILL IN G's FOR UNIQUE BONDS */
			if (Bonds[pol_num][seg_num][bond_num] == -1){        /* fill end segment equation */
				/* Boltz unk for 1st seg */
				unk_B = Phys2Unk_first[func_type_field] + itype_mer;

				/* hard wall boundary condition, G = 0 on wall */
				if(nodepos[0] == xbound) {
					if (resid_only_flag==INIT_GUESS_FLAG) resid=0.0;
					else resid = x[iunk][inode_box]/x[unk_B][inode_box];
					resid_G+=resid;
					if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) 
						dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
					if (resid_only_flag==FALSE){
						mat_val=1.0/x[unk_B][inode_box]; 
                                                if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
						dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode, iunk, inode_box, mat_val);
					}										
				}
				else {
					if (resid_only_flag==INIT_GUESS_FLAG) resid=-1.0;
					else  resid = x[iunk][inode_box]/x[unk_B][inode_box]-1.0;
					resid_G+=resid;
					if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) 
						dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
					if (resid_only_flag==FALSE){
						unkIndex[0]=iunk; unkIndex[1]=unk_B;
						values[0]=1.0/x[unk_B][inode_box]; 
						values[1]=-x[iunk][inode_box]/(x[unk_B][inode_box]*x[unk_B][inode_box]);
						numEntries=2;
                                                if (Iwrite_files==FILES_DEBUG_MATRIX){
                                                    for (jtmp=0;jtmp<numEntries;jtmp++) 
                                                        Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[unkIndex[jtmp]]*Nnodes]+=values[jtmp];
                                                }
						dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode, unkIndex, inode_box, values, numEntries);
					}
				}
			}
			else if (Bonds[pol_num][seg_num][bond_num] == -2) {		/* for grafted ends */
				unk_B = Phys2Unk_first[func_type_field] + itype_mer;
				if (nodepos[0] == xbound) {
					if (resid_only_flag==INIT_GUESS_FLAG) resid=-1.0;
					else resid = x[iunk][inode_box]/x[unk_B][inode_box] - 1.0;
					resid_G+=resid;
					if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) 
						dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
					if (resid_only_flag==FALSE){
						unkIndex[0]=iunk; unkIndex[1]=unk_B;
						values[0]=1.0/x[unk_B][inode_box]; 
						values[1]=-(x[iunk][inode_box]-1.0)/(x[unk_B][inode_box]*x[unk_B][inode_box]);
						numEntries=2;
                                                if (Iwrite_files==FILES_DEBUG_MATRIX){
                                                    for (jtmp=0;jtmp<numEntries;jtmp++) 
                                                        Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[unkIndex[jtmp]]*Nnodes]+=values[jtmp];
                                                }
						dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode, unkIndex, inode_box, values, numEntries);
					}
				}
				else {			/* for grafted end, G = 0 not on wall */
					if (resid_only_flag==INIT_GUESS_FLAG) resid=0.0;
					else resid = x[iunk][inode_box]/x[unk_B][inode_box];	
					resid_G+=resid;
					if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) 
						dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
					if (resid_only_flag==FALSE){
						mat_val=1.0/x[unk_B][inode_box]; 
                                                if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
						dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode, iunk, inode_box, mat_val);
					}					
				}

			}
			else{                                            /* fill G_seg eqns */

			/* First calculate the residual contributions */
				unk_B = Phys2Unk_first[func_type_field] + itype_mer;   /* Boltz unk for this seg */

				if (resid_only_flag != INIT_GUESS_FLAG){
					resid = x[iunk][inode_box]/x[unk_B][inode_box];

					resid_G+=resid;
					if (resid_only_flag !=CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
				}

				if (resid_only_flag==FALSE){
					unkIndex[0]=iunk; unkIndex[1]=unk_B;
					values[0]=1.0/x[unk_B][inode_box]; 
					values[1]=-x[iunk][inode_box]/(x[unk_B][inode_box]*x[unk_B][inode_box]);
					numEntries=2;
                                        if (Iwrite_files==FILES_DEBUG_MATRIX){
                                           for (jtmp=0;jtmp<numEntries;jtmp++) 
                                               Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[unkIndex[jtmp]]*Nnodes]+=values[jtmp];
                                        }
					dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode, unkIndex, inode_box, values, numEntries);
				}

				if(nodepos[0] != xbound) {
					/* Now Finish loading the Jacobian... */
					gint_tmp = load_polymer_recursion(sten,func_type_field,
                            Njacobian_types,Njacobian_sums,
                            funcArray_Jac, fp_ResidG,fp_ResidG_Bulk,
                            iunk,loc_inode,inode_box,unk_B,
                            itype_mer,izone,ijk_box,x, resid_only_flag);
					resid_G += gint_tmp;
				}
			}
			if (resid_only_flag==INIT_GUESS_FLAG){
				resid_G*=(-x[unk_B][inode_box]);			/* put Boltzmann factor in guess */
			}
		}
    }  /* end of fill something */
    return(resid_G);
}   



/****************************************************************************/
/* load_polymer_recursion:
             In this routine we load the Jacobian for
             gaussian (or freely-jointed) polymer terms.                      */
double load_polymer_recursion(int sten_type,int func_type_field, int Njacobian_types, int Njacobian_sums,
       void (*funcArray_Jac[3])(int,int,int,int,int,int,int,int,int *,double,double **),
       double (*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **),
       double (*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **),
       int iunk,int loc_inode, int inode_box, int unk_B,
       int itype_mer, int izone,int *ijk_box, double **x, int resid_only_flag)
{
  int   **sten_offset, *offset, isten;
  double *sten_weight,  weight;
  struct Stencil_Struct *sten;

  int jlist,unk_GQ,nunk,unk[20],jseg,jtype_mer;
  int reflect_flag[NDIM_MAX];
  int i,jnode_box,jcomp,isum,iterm;
  int pol_num,seg_num,bond_num,ibond,unk_test;
  double resid,resid_sum,jac_sum[10],wt_sum;

  if (Nlists_HW <= 2) jlist = 0;
  else                jlist = itype_mer;
  for (i=0; i<20; i++) unk[i]=-1;
  wt_sum=0.0;

  /* As in precalc routine, we first need to assemble lists for the various G
     products that we care about ... note that there are two flavors of the
     products now.  one excludes k(beta)=i, the other also excludes k(beta=i2)*/

   /* (1) from iunk get the segment number i and bond number*/
  unk_GQ = iunk-Phys2Unk_first[G_CHAIN];
  pol_num  = Unk_to_Poly[unk_GQ];
  seg_num  = Unk_to_Seg[unk_GQ];
  bond_num = Unk_to_Bond[unk_GQ];

  /*(2) from iunk and the bond number find the jseg to which we are looking */

  jseg = Bonds[pol_num][seg_num][bond_num];
  jtype_mer = Type_mer[pol_num][jseg];

  /* (3) locate all of the unknown numbers corresponding to the
         G or Q functions that will be needed for the Jacobian contributions of
         derivatives with respect to the Boltzmann term. */

  nunk=0;
  for (ibond=0; ibond<Nbond[pol_num][jseg];ibond++){
      if (Bonds[pol_num][jseg][ibond] != seg_num){
        unk[nunk]     = Poly_to_Unk[pol_num][jseg][ibond]+Geqn_start[pol_num];
        unk_test=unk[nunk]-Phys2Unk_first[G_CHAIN];
        if (Pol_Sym[unk_test] != -1) unk[nunk]= Pol_Sym[unk_test]+Phys2Unk_first[G_CHAIN];
        nunk++;
      }
  }
  /* don't forget Boltzman factor */
  unk[nunk++]=Phys2Unk_first[func_type_field]+jtype_mer;

  sten = &(Stencil[sten_type][izone][itype_mer+Ncomp*jtype_mer]);
  sten_offset = sten->Offset;
  sten_weight = sten->Weight;

  resid_sum=0.0;
  for (isum=0;isum<Njacobian_sums;isum++) jac_sum[isum]=0.0;
  for (isten = 0; isten < sten->Length; isten++) {
     offset = sten_offset[isten];
     weight = sten_weight[isten];

     /* Find the Stencil point */
     jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);
     jcomp=jtype_mer;
     if (jnode_box >= 0 && !Zero_density_TF[jnode_box][jcomp]) {
        if (Lhard_surf) {
           if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1)
           weight = HW_boundary_weight
                    (jtype_mer,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
        }
        /* first compute the residual */
        resid=(*fp_ResidG)(iunk,pol_num,jseg,unk_B,inode_box,jnode_box,nunk,unk,weight,x);
        resid_sum+=resid;
        if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

        if (resid_only_flag==FALSE){

           /* now compute local entries to the Jacobian - at the jnode_box position */
           /* note that the loading of these terms to the linprobmgr is done in the function */
           for (iterm=0;iterm<Njacobian_types;iterm++){
               (*funcArray_Jac[iterm])(iunk,loc_inode,pol_num,jseg,unk_B,inode_box,jnode_box,nunk,unk,weight,x);
           }
        }
     }
     else if (jnode_box <0){
         resid = (*fp_ResidG_Bulk)(iunk,pol_num,jseg,unk_B,inode_box,jnode_box,nunk,unk,weight,x);
         resid_sum+=resid;
         if (resid_only_flag!=INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
     }
wt_sum+=weight;
  }

  return(resid_sum);
}
/******************************************************************************************/
double fill_constant_density_chain(int iunk, int icomp, int iseg, double fac_FIELD,int loc_inode, int inode_box, double **x,int resid_only_flag)
{
  double resid,mat_val;

  if (resid_only_flag != INIT_GUESS_FLAG){
     resid = x[iunk][inode_box] ;
     mat_val = 1.0;
     if (resid_only_flag !=CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
     if (resid_only_flag==FALSE){
       if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
        dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);
     }
   }

   if (Lseg_densities)   resid = -(0.5*(Rho_seg_LBB[iseg]+Rho_seg_RTF[iseg]));
   else {                 resid = -(0.5*(Rho_b_LBB[icomp]+Rho_b_RTF[icomp]));
   }

   if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) {
                      dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
   }

   return resid;
}
/******************************************************************************************/

