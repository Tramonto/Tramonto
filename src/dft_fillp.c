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
 *  FILE: dft_fill.c
 *
 *  This file contains the fill of the residual equations and Jacobian
 *  matrix.
 */

#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"
double load_polymer_cr(int, int,int, int, int, int, int *, double **);
double load_polymer_G(int,int,int,int,int,int,int,int *,double **);

/****************************************************************************/
void fill_resid_and_matrix_P (double **x, int iter, int resid_only_flag, int unk_flag)
{
 /*
  * Local variable declarations
  */

  int     i, j, iunk,unk_GQ,unk_B,iunk_start,iunk_end,iloop;
  int     reflect_flag[3],junk;
  int     izone, mesh_coarsen_flag_i;
  int     loc_i, loc_inode, itype_mer;
  int     jtype_mer;
  int     nseg,ibond,pol_num,seg_num,bond_num;

  double resid_B=0.0,resid_R=0.0,resid_G=0.0,resid_P=0.0,k_poly,boltz=0.0;
  double nodepos[3];
 
  double gint_tmp;
  int ipol,iseg,numEntries,unkIndex[2];
  int boltz_pow,boltz_pow_J;
  double fac1,fac2;
  int jbond,unk_GQ_j,node_start;

  double  resid_tmp=0.0,resid,mat_val,values[2];

  /* the 6 offset patterns for nearest neighbors */
  int offset_idim_pm[18] = {1,0,0,  0,1,0,  0,0,1,  -1,0,0,  0,-1,0,  0,0,-1};
  int *offset_ptr; /* pointer into above */

  int i_box, inode_box,jnode_box, ijk_box[3];
  int npol=0,sten;

  int inode,ijk[3];

  /********************** BEGIN EXECUTION ************************************/

  if (Sten_Type[POLYMER_GAUSS]) sten=POLYMER_GAUSS;
  else                          sten=DELTA_FN;

  if (unk_flag == NODAL_FLAG){
      iunk_start = 0;
      iunk_end = Nunk_per_node;
  }
  else{
      iunk_start = unk_flag;
      iunk_end = unk_flag+1;
  }

  /* Load residuals and matrix */

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {    

    /* convert local node to global */

    inode_box = L2B_node[loc_inode];
    node_box_to_ijk_box(inode_box, ijk_box);

    if (Mesh_coarsening && Nwall_type >0) mesh_coarsen_flag_i = Mesh_coarsen_flag[inode_box];
    else                                  mesh_coarsen_flag_i = 0;

    for (iunk=iunk_start; iunk<iunk_end; iunk++) {               

    resid_B=0.0;resid_R=0.0;resid_G=0.0;resid_P=0.0;

      /* i_box is the equation number (matrix row) that we are filling now */

      if (Unk2Phys[iunk]==DENSITY) itype_mer = iunk-Phys2Unk_first[DENSITY]; 
      else if (Unk2Phys[iunk] == CMS_FIELD) itype_mer = iunk - Phys2Unk_first[CMS_FIELD];
      else if  (Unk2Phys[iunk] == CMS_G) {                   
                               npol = Unk_to_Poly[iunk-Phys2Unk_first[CMS_G]];
                               nseg = Unk_to_Seg[iunk-Phys2Unk_first[CMS_G]]; 
                               itype_mer = Type_mer[npol][nseg];
      }

                                                               /*DO ZERO DENSITY FILL*/
      if ( ((Zero_density_TF[inode_box][itype_mer] ||    
            Vext[loc_inode][itype_mer] == VEXT_MAX) && Unk2Phys[iunk]!=POISSON)
            /*|| (Unk2Phys[iunk]==DENSITY 
               && -log(x[loc_find(Phys2Unk_first[CMS_FIELD]+itype_mer,loc_inode,LOCAL)]) > VEXT_MAX )*/ ){
           resid = x[iunk][inode_box];
           mat_val=1.0;
           dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
           dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);
      }
      else if (mesh_coarsen_flag_i < 0) {                   /*DO COARSENED MESH FILL*/

         resid= x[iunk][inode_box];
         mat_val=1.0;
         dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);    

         for (iloop=0;iloop<2;iloop++){

          if (iloop==0) offset_ptr = &offset_idim_pm[3*(-mesh_coarsen_flag_i - 1)];  /* one node higher */
          else          offset_ptr += 9;                                             /* one node lower */

          jnode_box = offset_to_node_box(ijk_box, offset_ptr, reflect_flag);
  
          if (jnode_box >= 0) {
             resid= - 0.5*x[iunk][jnode_box];
             mat_val=-0.5;
             dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,jnode_box,mat_val);
          }
          else{  resid= - 0.5*constant_boundary(iunk,jnode_box); }
        }
        dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
      }
      else {                                           /* FILL REAL EQUATIONS */

         /* izone is set based on value of mesh_coarsen_flag       */
         /* izone is the zone number at the node of interest which */
         /* indicated which quadrature scheme is to be used.       */

         izone = mesh_coarsen_flag_i;
               /**************************************************************/
         if (Unk2Phys[iunk] == CMS_FIELD){             /*BOLTZMANN EQNS*/

             resid_B = load_polymer_cr(POLYMER_CR,iunk,loc_inode,inode_box,itype_mer,izone,ijk_box,x); 
             resid = Vext[loc_inode][itype_mer]+log(x[iunk][inode_box]);
             resid_B+=resid;
             dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
             mat_val = 1.0/x[iunk][inode_box];  
             dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);

             if (Ipot_ff_c == COULOMB){
                 junk = Phys2Unk_first[POISSON];
                 resid = Charge_f[itype_mer]*x[junk][inode_box];
                 resid_B+=resid;
                 mat_val = Charge_f[itype_mer];
                 dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                 dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,inode_box,mat_val);
              }
         }                                           /* END BOLTZMAN EQNS */
               /**************************************************************/
         else if (Unk2Phys[iunk]==DENSITY){   /* DENSITY EQNS FOR POLYMER SEGS*/

	   /* Boltzmann factor for this i */
            unk_B=Phys2Unk_first[CMS_FIELD]+itype_mer;

            if (Type_poly == 2 || Type_poly == 1) {
              resid = x[iunk][inode_box];
              resid_R+=resid;
              mat_val = 1.;
              dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
              dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);
            }
            else{
              resid = x[iunk][inode_box]*x[unk_B][inode_box];
              resid_R+=resid;
              values[0]=x[iunk][inode_box]; values[1]=x[unk_B][inode_box];
              unkIndex[0]=unk_B; unkIndex[1]=iunk;
              numEntries=2;
              dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
              dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                   unkIndex, inode_box, values, numEntries);
            }

            npol = 0;
            while (Nmer_t[npol][itype_mer]==0) npol++;
            k_poly = Rho_b[itype_mer]/Nmer_t[npol][itype_mer];

            for (i=0; i<Nmer[npol]; i++){
              if (Type_mer[npol][i] == itype_mer) {

                 if (Type_poly == 0 ||Type_poly==3){
                      boltz_pow = -(Nbond[npol][i]-2);
                      boltz_pow_J = -(Nbond[npol][i]-1);
                 }
                 else if (Type_poly == 1){
                      boltz_pow = -(Nbond[npol][i]-1);
                      boltz_pow_J = -(Nbond[npol][i]);
                 }
                 else {
                      boltz_pow = 1;
                      boltz_pow_J = 0;
                 }
                 fac1= k_poly;
                 for (ibond=0; ibond<Nbond[npol][i]; ibond++) {
                    unk_GQ  = Geqn_start[npol] + Poly_to_Unk[npol][i][ibond]; 
                    fac1 *= x[unk_GQ][inode_box];
 
                    fac2= k_poly;
                    for (jbond=0; jbond<Nbond[npol][i]; jbond++) {
                      if (jbond != ibond){       
                        unk_GQ_j  = Geqn_start[npol] + Poly_to_Unk[npol][i][jbond]; 
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
         }                                  /* END POLYMER DENSITY EQNS */
               /**************************************************************/
         else if (Unk2Phys[iunk] == CMS_G){          /* G EQNS */

             unk_GQ=iunk-Phys2Unk_first[CMS_G];
             pol_num = Unk_to_Poly[unk_GQ];
             seg_num = Unk_to_Seg[unk_GQ];
             bond_num = Unk_to_Bond[unk_GQ];
      
             if (Pol_Sym[unk_GQ] != -1){                            /* FILL IN G's FOR SYMMETRIC BONDS */
                junk = Pol_Sym[unk_GQ] + Phys2Unk_first[CMS_G];
                resid = x[iunk][inode_box]-x[junk][inode_box];
                resid_G+=resid;
                dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                unkIndex[0]=iunk; unkIndex[1]=junk;
                values[0]=1.0; values[1]=-1.0;
                numEntries=2;
                dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                   unkIndex, inode_box, values, numEntries);
             } 
             else{                                                  /* FILL IN G's FOR UNIQUE BONDS */

               if (Bonds[pol_num][seg_num][bond_num] == -1){        /* fill end segment equation */
                  if (Type_poly == 0 || Type_poly == 1 || Type_poly == 3){
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
                  else{
                        resid = x[iunk][inode_box]-1.0;
                        resid_G+=resid;
                        dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                        mat_val = 1.0;
                        dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);
                  }

               }
               else{                                            /* fill G_seg eqns */

                                 /* First calculate the residual contributions */

                  unk_B = Phys2Unk_first[CMS_FIELD] + itype_mer;   /* Boltz unk for this seg */
                  if (Type_poly == 0 || Type_poly == 3 ){
                     boltz = x[unk_B][inode_box];
                     resid = x[iunk][inode_box]; 
                     resid_G+=resid;
                     dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                     mat_val = 1.;
                     dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);
                  }
                  else if (Type_poly == 1) {
                     boltz = 1.0;
                     resid = x[iunk][inode_box]/x[unk_B][inode_box];
                     resid_G+=resid;
                     dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                     values[0] += 1./x[unk_B][inode_box];
                     values[1] -= x[iunk][inode_box]/(x[unk_B][inode_box]*x[unk_B][inode_box]);
                     unkIndex[0]=iunk; unkIndex[1]=unk_B;
                     numEntries=2;
                     dft_linprobmgr_insertmultiphysicsmatrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                             unkIndex, inode_box, values, numEntries);
                   }
                   else if (Type_poly == 2) {
                     boltz = 1.0;
                     resid = x[iunk][inode_box];
                     resid_G+=resid;
                     mat_val = 1.;
                     dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                     dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);
                   }

                  /* Now Finish loading the Jacobian... */
                  gint_tmp = load_polymer_G(sten,iunk,loc_inode,inode_box,unk_B,itype_mer,izone,ijk_box,x);
                  resid_G += gint_tmp;
                  if (Type_poly==0 || Type_poly==3){
                     mat_val = gint_tmp / boltz;
                     dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_B,inode_box,mat_val);
                  }

               }
             }
         }                                  /* END G EQNS */
               /**************************************************************/
         else if (Unk2Phys[iunk] == POISSON){ /* LOAD POISSON'S EQUATION*/
              resid_P=load_poissons_eqn(iunk,loc_inode,inode_box,ijk_box,x);
              resid_P+=load_poisson_bc(iunk,loc_inode,inode_box);
         }
               /**************************************************************/
         else {
            printf("Problem with unknowns in fill !!!\n");
            exit (-1);
         }
      }      /* end of else (not Zero_density and mesh_coarsen_flag_i >= 0) */

/*     if (fabs(resid_B)>1.e-3 ||fabs(resid_R)>1.e-3 ||fabs(resid_G)>1.e-3){
       printf("%d %d %d  %9.6f %9.6f %9.6f %9.6f \n",iunk+Nunk_per_node*loc_inode,inode_box,iunk, resid_B,resid_R,resid_G,resid_P);
     } */

    } /* end ofloop over number of unknowns per node */
  } /* end of loop over local nodes */
  return;
}
/*****************************************************************************/
/* load_polymer_G: 
             In this routine we load the Jacobian for 
             gaussian (or freely-jointed) polymer terms.                      */
double load_polymer_G(int sten_type,int iunk,int loc_inode, int inode_box,
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
  int pol_num,seg_num,bond_num,ibond;
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
        unk[nunk++]     = Poly_to_Unk[pol_num][jseg][ibond]+Phys2Unk_first[CMS_G];
      }
  }
  /* don't forget Boltzman factor */
  unk[nunk++]=jtype_mer;

  if (Type_poly == 0 || Type_poly == 3) {
      boltz_pow_R = -(Nbond[pol_num][jseg]-2);
      boltz_prefac_R = -x[unk_B][inode_box];
      boltz_pow_1 = -(Nbond[pol_num][jseg]-1);
      boltz_prefac_1 = x[unk_B][inode_box]*(Nbond[pol_num][jseg]-2);
      boltz_pow_2 = -(Nbond[pol_num][jseg]-2);
      boltz_prefac_2 = -x[unk_B][inode_box];
  }
  else if (Type_poly == 1){
      boltz_pow_R = -(Nbond[pol_num][jseg]-2);
      boltz_prefac_R = -1.0;
      boltz_pow_1 = -(Nbond[pol_num][jseg]-1);
      boltz_prefac_1 = (Nbond[pol_num][jseg]-2);
      boltz_pow_2 = -(Nbond[pol_num][jseg]-2);
      boltz_prefac_2 = -1.0;
  }
  else{
       boltz_pow_R= 1;
       boltz_prefac_R = -1.0;
       boltz_pow_1 = 0;
       boltz_prefac_1 = -1.0;
       boltz_pow_2 = 1;
       boltz_prefac_2 = -1.0;
  }

  sten = &(Stencil[sten_type][izone][itype_mer+Ncomp*jtype_mer]);
  sten_offset = sten->Offset;
  sten_weight = sten->Weight;

  resid_sum=0.0;
  for (isten = 0; isten < sten->Length; isten++) {
     offset = sten_offset[isten];
     weight = sten_weight[isten];

     /* Find the Stencil point */
     jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);
     if (jnode_box >= 0 && !Zero_density_TF[jnode_box][unk[nunk-1]]) {
        if (Lhard_surf) {
           if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1) 
           weight = HW_boundary_weight 
                    (itype_mer+Ncomp*jtype_mer,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
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
/*****************************************************************************/
/* load_polymer_cr: 
             In this routine we load the residual and Jacobian for 
             the polymer fields.                      */
double load_polymer_cr(int sten_type,int iunk,int loc_inode,int inode_box,int itype_mer, 
                int izone, int *ijk_box, double **x)
{
  int   **sten_offset, *offset, isten;
  double *sten_weight,  weight, weight_bulk;
  struct Stencil_Struct *sten;

  double sign,resid_sum,resid,mat_val;
  int jtype_mer, jlist,junk;
  int reflect_flag[NDIM_MAX];
  int jnode_box,node_start;

  sign = 1.0;
  resid_sum=0.0;

  for (jtype_mer=0; jtype_mer<Ncomp; jtype_mer++){
      junk=Phys2Unk_first[DENSITY]+jtype_mer;
      if (Nlists_HW <= 2) jlist = 0;
      else                jlist = jtype_mer;
           
      sten = &(Stencil[sten_type][izone][itype_mer+Ncomp*jtype_mer]);
      sten_offset = sten->Offset;
      sten_weight = sten->Weight;

      for (isten = 0; isten < sten->Length; isten++) {
        offset = sten_offset[isten];
        weight_bulk = sten_weight[isten];

         /* Find the Stencil point */
         jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);

         resid=0.0;
         if (jnode_box >= 0 ) {
            if (Lhard_surf) {
                if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1) 
                   weight = HW_boundary_weight 
                    (jtype_mer,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
                else weight = weight_bulk;
            }
            else weight = weight_bulk;

                                            /* density of this component type */
             resid =  -sign*(weight*x[junk][jnode_box]- weight_bulk*Rho_b[jtype_mer]);
             mat_val = -sign*weight;
             dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,jnode_box,mat_val);
         }
         else if ( jnode_box == -2){  /*in wall*/
              resid =  sign*weight_bulk*Rho_b[jtype_mer];
         }
         resid_sum+=resid;
         dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
         /* jnode_box == -1 = in bulk contributes nothing  */
      }
  }

  return(resid_sum);
}
/*****************************************************************************/
