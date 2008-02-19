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
  int loc_inode, int inode_box, int resid_only_flag, double (*fp_prefactor)(int))
{
  int i,loop_start,loop_end,itype_mer,npol,iseg,unk_GQ,unk_GQ_test,iref;
  int boltz_pow,boltz_pow_J,jbond,ibond,unkIndex[2],numEntries,unk_GQ_j,unk_GQ_j_test;
  double fac1,fac2,mat_val,resid,resid_sum=0.0,values[2];

  resid = x[iunk][inode_box]*x[unk_B][inode_box];
  resid_sum=resid;
  dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
  if(!resid_only_flag){
      values[0]=x[iunk][inode_box]; values[1]=x[unk_B][inode_box];
      unkIndex[0]=unk_B; unkIndex[1]=iunk;
      numEntries=2;
      dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                              unkIndex, inode_box, values, numEntries);
  }

  if (Type_poly==WJDC){
     itype_mer=Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
     loop_start=iunk;  /* don't sum over segments of a given type for WJDC */
     loop_end=iunk+1;
  }
  else if (Type_poly==CMS){
     itype_mer=iunk-Phys2Unk_first[DENSITY];
     npol = 0;
     while (Nmer_t[npol][itype_mer]==0){
         npol++;
     }
     loop_start=0;      /* for CMS we _do_ sum over segment of a given type */
     loop_end=Nmer[npol];
  }

  for (i=loop_start; i<loop_end; i++){
     if (Lseg_densities || (Type_mer[npol][i] == itype_mer) ) {

        if (!Lseg_densities) iseg=SegChain2SegAll[npol][i];
        else                 iseg=iunk-Phys2Unk_first[DENSITY];


        boltz_pow = -(Nbonds_SegAll[iseg]-2);
        boltz_pow_J = -(Nbonds_SegAll[iseg]-1);

        if (fp_prefactor!=NULL) {
             if (Type_poly==CMS) iref=itype_mer;
             if (Type_poly==WJDC) iref=iseg;
             fac1 = (*fp_prefactor)(iref);
        }
        else                     fac1=1.0;

        for (ibond=0; ibond<Nbonds_SegAll[iseg]; ibond++) {
              unk_GQ  = Phys2Unk_first[func_type] + Poly_to_Unk_SegAll[iseg][ibond];
              unk_GQ_test = unk_GQ-Phys2Unk_first[func_type];
              if (Pol_Sym[unk_GQ_test] != -1) unk_GQ=Pol_Sym[unk_GQ_test] + Phys2Unk_first[func_type];
              fac1 *= x[unk_GQ][inode_box];

              if(!resid_only_flag){
                 if (fp_prefactor !=NULL) fac2= (*fp_prefactor)(iunk);
                 else                     fac2=1.0;
                 for (jbond=0; jbond<Nbonds_SegAll[iseg]; jbond++) {
                   if (jbond != ibond){
                     unk_GQ_j  = Phys2Unk_first[func_type]+Poly_to_Unk_SegAll[iseg][jbond];
                     unk_GQ_j_test=unk_GQ_j-Phys2Unk_first[func_type];
                     if (Pol_Sym[unk_GQ_j_test] != -1) unk_GQ_j=Pol_Sym[unk_GQ_j_test] + Phys2Unk_first[func_type];
                     fac2 *= x[unk_GQ_j][inode_box];
                   }
                 }
                 mat_val = -fac2*POW_DOUBLE_INT(x[unk_B][inode_box],boltz_pow);
                 dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_GQ,inode_box,mat_val);
               }
           }
           if(!resid_only_flag){
              mat_val = -fac1*((double) boltz_pow)*POW_DOUBLE_INT(x[unk_B][inode_box],boltz_pow_J);
              dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_B,inode_box,mat_val);
           }
           resid = -fac1*POW_DOUBLE_INT(x[unk_B][inode_box],boltz_pow);
           resid_sum+=resid;
           dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
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
    int itype_mer,jseg,jtype_mer,unk_xi2,unk_xi3;
    double nodepos[3];
    double resid_G,resid,mat_val,values[4],gint_tmp; 
    double y,ysqrt,dydxi,xi_2,xi_3;
   
       /* use arrays that are indexed with the geqns ordered starting at 0 */ 
    unk_GQ=iunk-Phys2Unk_first[G_CHAIN];
    pol_num = Unk_to_Poly[unk_GQ];
    seg_num = Unk_to_Seg[unk_GQ];
    bond_num = Unk_to_Bond[unk_GQ];
    itype_mer = Type_mer[pol_num][seg_num];                   /*itype_mer=icomp*/

     /* from iunk and the bond number find the jseg to which we are looking */

    jseg = Bonds[pol_num][seg_num][bond_num];
    jtype_mer = Type_mer[pol_num][jseg];

    resid_G=0.0;
    
    sten=DELTA_FN_BOND;
    if (Zero_density_TF[inode_box][itype_mer] || Vext[loc_inode][itype_mer] == VEXT_MAX){
         resid_G=fill_zero_value(iunk,loc_inode,inode_box,x,resid_only_flag);
    }
    else {
       
       if (Pol_Sym[unk_GQ] != -1){                            /* FILL IN G's FOR SYMMETRIC BONDS */
          junk = Pol_Sym[unk_GQ] + Phys2Unk_first[G_CHAIN];
          resid = x[iunk][inode_box]-x[junk][inode_box];
          resid_G+=resid;
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          if (!resid_only_flag){
             unkIndex[0]=iunk; unkIndex[1]=junk;
             values[0]=1.0; values[1]=-1.0;
             numEntries=1;
             dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                             unkIndex, inode_box, values, numEntries);
          }
       }
       else{                                                  /* FILL IN G's FOR UNIQUE BONDS */

         if (Bonds[pol_num][seg_num][bond_num] == -1){        /* fill end segment equation */
                                                /* Boltz unk for 1st seg */
            unk_B = Phys2Unk_first[func_type_field] + itype_mer;

            if (Type_poly==CMS){
               resid = x[iunk][inode_box]/x[unk_B][inode_box]-1.0;
               resid_G+=resid;
               dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
               if (!resid_only_flag){
                  unkIndex[0]=iunk; unkIndex[1]=unk_B;
                  values[0]=1.0/x[unk_B][inode_box]; 
                  values[1]=-x[iunk][inode_box]/(x[unk_B][inode_box]*x[unk_B][inode_box]);
                  numEntries=2;
                  dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                 unkIndex, inode_box, values, numEntries);
               }
            }
            else if (Type_poly==WJDC){
               unk_xi2=Phys2Unk_first[CAVWTC]; unk_xi3=Phys2Unk_first[CAVWTC]+1;
               xi_2=x[unk_xi2][inode_box]; xi_3=x[unk_xi3][inode_box];
               y=y_cav(Sigma_ff[itype_mer][itype_mer],Sigma_ff[jtype_mer][jtype_mer],xi_2,xi_3);
               ysqrt=sqrt(y);
               resid = x[iunk][inode_box]/(x[unk_B][inode_box]*ysqrt)-(1.0/ysqrt);
               resid_G+=resid;
               dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

               if (!resid_only_flag){
                  unkIndex[0]=iunk; unkIndex[1]=unk_B; unkIndex[2]=unk_xi2; unkIndex[3]=unk_xi3;
                  values[0]=1.0/(ysqrt*x[unk_B][inode_box]); 
                  values[1]=-x[iunk][inode_box]/(ysqrt*x[unk_B][inode_box]*x[unk_B][inode_box]);

                  dydxi=dy_dxi2_cav(Sigma_ff[itype_mer][itype_mer],Sigma_ff[jtype_mer][jtype_mer],xi_2,xi_3);
                  values[2]=-(dydxi/(y*ysqrt))*((x[iunk][inode_box]/x[unk_B][inode_box]) -1.0);
                  dydxi=dy_dxi3_cav(Sigma_ff[itype_mer][itype_mer],Sigma_ff[jtype_mer][jtype_mer],xi_2,xi_3);
                  values[3]=-(dydxi/(y*ysqrt))*((x[iunk][inode_box]/x[unk_B][inode_box]) -1.0);
                  numEntries=4;
                  dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                 unkIndex, inode_box, values, numEntries);
               }
            }
         }
         else{                                            /* fill G_seg eqns */

                        /* First calculate the residual contributions */
            unk_B = Phys2Unk_first[func_type_field] + itype_mer;   /* Boltz unk for this seg */

            if (Type_poly==CMS){
            resid = x[iunk][inode_box]/x[unk_B][inode_box];
            resid_G+=resid;
            dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

            if (!resid_only_flag){
               unkIndex[0]=iunk; unkIndex[1]=unk_B;
               values[0]=1.0/x[unk_B][inode_box]; 
               values[1]=-x[iunk][inode_box]/(x[unk_B][inode_box]*x[unk_B][inode_box]);
               numEntries=2;
               dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                 unkIndex, inode_box, values, numEntries);
            }
            }
            else if (Type_poly==WJDC){
               unk_xi2=Phys2Unk_first[CAVWTC]; 
               unk_xi3=Phys2Unk_first[CAVWTC]+1;
               xi_2=x[unk_xi2][inode_box]; xi_3=x[unk_xi3][inode_box];
               y=y_cav(Sigma_ff[itype_mer][itype_mer],Sigma_ff[jtype_mer][jtype_mer],xi_2,xi_3);
               ysqrt=sqrt(y);
               resid = x[iunk][inode_box]/(x[unk_B][inode_box]*ysqrt);
               resid_G+=resid;
               dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

               if (!resid_only_flag){
                  unkIndex[0]=iunk; unkIndex[1]=unk_B; unkIndex[2]=unk_xi2; unkIndex[3]=unk_xi3;
                  values[0]=1.0/(ysqrt*x[unk_B][inode_box]); 
                  values[1]=-x[iunk][inode_box]/(ysqrt*x[unk_B][inode_box]*x[unk_B][inode_box]);
                  dydxi=dy_dxi2_cav(Sigma_ff[itype_mer][itype_mer],Sigma_ff[jtype_mer][jtype_mer],xi_2,xi_3);
                  values[2]=-0.5*x[iunk][inode_box]*dydxi/(y*ysqrt*x[unk_B][inode_box]);
                  dydxi=dy_dxi3_cav(Sigma_ff[itype_mer][itype_mer],Sigma_ff[jtype_mer][jtype_mer],xi_2,xi_3);
                  values[3]=-0.5*x[iunk][inode_box]*dydxi/(y*ysqrt*x[unk_B][inode_box]);
                  numEntries=4;
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
  int i,j,jnode_box,jcomp,isum,iterm;
  int pol_num,seg_num,bond_num,ibond,unk_test;
  double resid,resid_sum,mat_val,jac_sum[10];

  if (Nlists_HW <= 2) jlist = 0;
  else                jlist = itype_mer;
  for (i=0; i<20; i++) unk[i]=-1;

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
        dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

        if (!resid_only_flag){

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
         dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
     }
  }

  return(resid_sum);
}
/****************************************************************************/
