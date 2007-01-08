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
 *  FILE: dft_fill_cms.c
 *
 *  This file contains the matrix fill routine associated with the CMS
 *  physics functionals.
 */

/*#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"*/
#include "dft_fill_CMS.h"

/****************************************************************************/
double load_CMS_field(int iunk, int loc_inode, int inode_box, int *ijk_box, int izone, double **x,int resid_only_flag)
{
   double resid_B,resid,mat_val;
   int itype_mer,junk;

   itype_mer = iunk - Phys2Unk_first[CMS_FIELD];

    if (Zero_density_TF[inode_box][itype_mer] || Vext[loc_inode][itype_mer] == VEXT_MAX){
         resid_B=fill_zero_value(iunk,loc_inode,inode_box,x);
    }
    else{

       resid_B = load_mean_field(POLYMER_CR,iunk,loc_inode,itype_mer,izone,ijk_box,x,resid_only_flag); 
       resid = Vext[loc_inode][itype_mer]+log(x[iunk][inode_box]);
       resid_B+=resid;
       dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
       mat_val = 1.0/x[iunk][inode_box];  
       dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);

       /* add a mean field electrostatic contribution to the CMS field if there are charges in the problem */
       if (Type_coul != NONE){
          junk = Phys2Unk_first[POISSON];
          resid = Charge_f[itype_mer]*x[junk][inode_box];
          resid_B+=resid;
          mat_val = Charge_f[itype_mer];
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,inode_box,mat_val);
       }          
    }
    return(resid_B);
}
/****************************************************************************/
double load_CMS_density(int iunk, int loc_inode, int inode_box, double **x,int resid_only_flag)
{
   int unk_B,boltz_pow,boltz_pow_J,unk_GQ_j,unk_GQ_j_test,jbond;
   double resid_R,resid,mat_val,values[2],k_poly,fac1,fac2;
   int unkIndex[2],numEntries,ibond,unk_GQ,unk_GQ_test,itype_mer,npol,i;

   resid_R=0.0;

   itype_mer = iunk-Phys2Unk_first[DENSITY];

   if (Zero_density_TF[inode_box][itype_mer] || Vext[loc_inode][itype_mer] == VEXT_MAX){
         resid_R=fill_zero_value(iunk,loc_inode,inode_box,x);
   }
   else{

      /* Boltzmann factor for this i */
      unk_B=Phys2Unk_first[CMS_FIELD]+itype_mer;

      resid = x[iunk][inode_box]*x[unk_B][inode_box];
      resid_R+=resid;
      values[0]=x[iunk][inode_box]; values[1]=x[unk_B][inode_box];
      unkIndex[0]=unk_B; unkIndex[1]=iunk;
      numEntries=2;
      dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
      dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                        unkIndex, inode_box, values, numEntries);

      npol = 0;
      while (Nmer_t[npol][itype_mer]==0) npol++;
      k_poly = Rho_b[itype_mer]/Nmer_t[npol][itype_mer];

      for (i=0; i<Nmer[npol]; i++){
        if (Type_mer[npol][i] == itype_mer) {

           boltz_pow = -(Nbond[npol][i]-2);
           boltz_pow_J = -(Nbond[npol][i]-1);
           fac1= k_poly;
           for (ibond=0; ibond<Nbond[npol][i]; ibond++) {
              unk_GQ  = Geqn_start[npol] + Poly_to_Unk[npol][i][ibond];
              unk_GQ_test = unk_GQ-Phys2Unk_first[CMS_G];
              if (Pol_Sym[unk_GQ_test] != -1) unk_GQ=Pol_Sym[unk_GQ_test] + Phys2Unk_first[CMS_G]; 
              fac1 *= x[unk_GQ][inode_box];

              fac2= k_poly;
              for (jbond=0; jbond<Nbond[npol][i]; jbond++) {
                if (jbond != ibond){       
                  unk_GQ_j  = Geqn_start[npol] + Poly_to_Unk[npol][i][jbond]; 
                  unk_GQ_j_test=unk_GQ_j-Phys2Unk_first[CMS_G];
                  if (Pol_Sym[unk_GQ_j_test] != -1) unk_GQ_j=Pol_Sym[unk_GQ_j_test] + Phys2Unk_first[CMS_G];
                  fac2 *= x[unk_GQ_j][inode_box];
                }        
              }        
              mat_val = -fac2*POW_DOUBLE_INT(x[unk_B][inode_box],boltz_pow);
              dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_GQ,inode_box,mat_val);
           }        
           mat_val = -fac1*((double) boltz_pow)*POW_DOUBLE_INT(x[unk_B][inode_box],boltz_pow_J);
           dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_B,inode_box,mat_val);
           resid = -fac1*POW_DOUBLE_INT(x[unk_B][inode_box],boltz_pow);
           resid_R+=resid;
           dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
        }        
     }        
  }
  return(resid_R);
}                                  
/****************************************************************************/
double load_CMS_Geqns(int iunk, int loc_inode, int inode_box, int *ijk_box, int izone, double **x,int resid_only_flag)
{
    int unk_GQ,unk_B,pol_num,seg_num,bond_num,junk,unkIndex[2],numEntries,sten;
    int itype_mer,npol,nseg;
    double nodepos[3];
    double resid_G,resid,mat_val,values[2],boltz,gint_tmp;

    npol = Unk_to_Poly[iunk-Phys2Unk_first[CMS_G]];
    nseg = Unk_to_Seg[iunk-Phys2Unk_first[CMS_G]];
    itype_mer = Type_mer[npol][nseg];

    resid_G=0.0;

    unk_GQ=iunk-Phys2Unk_first[CMS_G];
    pol_num = Unk_to_Poly[unk_GQ];
    seg_num = Unk_to_Seg[unk_GQ];
    bond_num = Unk_to_Bond[unk_GQ];

    sten=DELTA_FN;

    if (Zero_density_TF[inode_box][itype_mer] || Vext[loc_inode][itype_mer] == VEXT_MAX){
         resid_G=fill_zero_value(iunk,loc_inode,inode_box,x);
    }
    else {

       if (Pol_Sym[unk_GQ] != -1){                            /* FILL IN G's FOR SYMMETRIC BONDS */ 
          junk = Pol_Sym[unk_GQ] + Phys2Unk_first[CMS_G];
          resid = x[iunk][inode_box]-x[junk][inode_box];
          resid = x[iunk][inode_box];
          resid_G+=resid;
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          unkIndex[0]=iunk; unkIndex[1]=junk;
          values[0]=1.0; values[1]=-1.0;
          numEntries=1;
          dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                             unkIndex, inode_box, values, numEntries);
       }        
       else{                                                  /* FILL IN G's FOR UNIQUE BONDS */ 

         if (Bonds[pol_num][seg_num][bond_num] == -1){        /* fill end segment equation */
            node_to_position(B2G_node[inode_box],nodepos);
                                                /* Boltz unk for 1st seg */
            unk_B = Phys2Unk_first[CMS_FIELD] + Type_mer[npol][seg_num];     
            resid = x[iunk][inode_box]-x[unk_B][inode_box];
            resid_G+=resid;
            dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
            unkIndex[0]=iunk; unkIndex[1]=unk_B;
            values[0]=1.0; values[1]=-1.0;
            numEntries=2;
            dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                 unkIndex, inode_box, values, numEntries);
         }        
         else{                                            /* fill G_seg eqns */

                        /* First calculate the residual contributions */

            unk_B = Phys2Unk_first[CMS_FIELD] + itype_mer;   /* Boltz unk for this seg */
            boltz = x[unk_B][inode_box];
            resid = x[iunk][inode_box]; 
            resid_G+=resid;
            dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
            mat_val = 1.;
            dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);

            /* Now Finish loading the Jacobian... */
            gint_tmp = load_polymer_recursion(sten,iunk,loc_inode,inode_box,unk_B,itype_mer,izone,ijk_box,x);
            resid_G += gint_tmp;

            mat_val = gint_tmp / boltz; 
            dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_B,inode_box,mat_val);

         }        
       }        
    }
    return(resid_G);
}                                           
/****************************************************************************/
/* load_polymer_recursion: 
             In this routine we load the Jacobian for 
             gaussian (or freely-jointed) polymer terms.                      */
double load_polymer_recursion(int sten_type,int iunk,int loc_inode, int inode_box,
                    int unk_B,int itype_mer, int izone,int *ijk_box, double **x)
{
  int   **sten_offset, *offset, isten;
  double *sten_weight,  weight;
  struct Stencil_Struct *sten;

  int jlist,unk_GQ,nunk,unk[20],jseg,jtype_mer;
  int boltz_pow_1,boltz_pow_2,boltz_pow_R;
  double boltz_prefac_1,boltz_prefac_2,boltz_prefac_R,fac1,fac2;
  int reflect_flag[NDIM_MAX];
  int i,j,jnode_box;
  int pol_num,seg_num,bond_num,ibond,unk_test;
  double resid,resid_sum,mat_val;

  if (Nlists_HW <= 2) jlist = 0;
  else                jlist = itype_mer; 
  for (i=0; i<20; i++) unk[i]=-1;

  /* As in precalc routine, we first need to assemble lists for the various G 
     products that we care about ... note that there are two flavors of the
     products now.  one excludes k(beta)=i, the other also excludes k(beta=i2)*/

   /* (1) from iunk get the segment number i and bond number*/
  unk_GQ = iunk-Phys2Unk_first[CMS_G];
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
        unk_test=unk[nunk]-Phys2Unk_first[CMS_G];
        if (Pol_Sym[unk_test] != -1) unk[nunk]= Pol_Sym[unk_test]+Phys2Unk_first[CMS_G];
        nunk++;
      }
  }
  /* don't forget Boltzman factor */
  unk[nunk++]=jtype_mer+Phys2Unk_first[CMS_FIELD];

  boltz_pow_R = -(Nbond[pol_num][jseg]-2);
  boltz_prefac_R = -x[unk_B][inode_box];
  boltz_pow_1 = -(Nbond[pol_num][jseg]-1);
  boltz_prefac_1 = x[unk_B][inode_box]*(Nbond[pol_num][jseg]-2);
  boltz_pow_2 = -(Nbond[pol_num][jseg]-2);
  boltz_prefac_2 = -x[unk_B][inode_box];

  sten = &(Stencil[sten_type][izone][itype_mer+Ncomp*jtype_mer]);
  sten_offset = sten->Offset;
  sten_weight = sten->Weight;

  resid_sum=0.0;
  for (isten = 0; isten < sten->Length; isten++) {
     offset = sten_offset[isten];
     weight = sten_weight[isten];

     /* Find the Stencil point */
     jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);
     if (jnode_box >= 0 && !Zero_density_TF[jnode_box][unk[nunk-1]-Phys2Unk_first[CMS_FIELD]]) {
        if (Lhard_surf) {
           if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1) 
           weight = HW_boundary_weight 
                    (jtype_mer,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
        }
        /* first load the Boltzman factor derivatives */
        fac1=weight; 
        for(i=0;i<nunk-1;i++) fac1 *=x[unk[i]][jnode_box];  /*Gs or Qs*/
        resid = fac1*boltz_prefac_R*POW_DOUBLE_INT(x[unk[nunk-1]][jnode_box],boltz_pow_R); /* Boltz Term */
        resid_sum += resid;

        dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
        mat_val = fac1*boltz_prefac_1*POW_DOUBLE_INT(x[unk[nunk-1]][jnode_box],boltz_pow_1); /* Boltz Term */
        dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk[nunk-1],jnode_box,mat_val);
           
          
        /* now load the G/Q derivatives */

        for (i=0;i < nunk-1; i++){
           fac2=weight;
           for (j=0; j<nunk-1; j++){
              if (j != i)  fac2 *= x[unk[j]][jnode_box];  /*Gs or Qs*/
           }
           mat_val = fac2*boltz_prefac_2*POW_DOUBLE_INT(x[unk[nunk-1]][jnode_box],boltz_pow_2); /*Boltz Term*/
           dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk[i],jnode_box,mat_val);
        }
     }
  }
  return(resid_sum);
}
/****************************************************************************/
