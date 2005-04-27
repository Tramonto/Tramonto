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
#define HIT_FLAG 999
/****************************************************************************/
void fill_resid_and_matrix (double *x, double *resid,
                            int **bindx_2d, double *fill_time, int fill_flag,
                            int iter, int resid_only_flag,int unk_flag)
{
 /*
  * Local variable declarations
  */

  char   *yo = "fill_resid_and_matrix";
  int     i, j, icomp,idim,iunk,iunk_start,iunk_end;
  int     reflect_flag[3];
  int     izone, mesh_coarsen_flag_i;
  int    *bindx_tmp=NULL;/*Full storage of 1 row of bindx for MSR_PREPROCESS*/
  int     loc_i, loc_j, loc_inode, loc_i_charge, loc_mu,unk_mu,j_box;
  struct  RB_Struct *rho_bar = NULL;
  struct  RB_Struct *dphi_drb=NULL, dphi_drb_bulk,
                     dphi_drb_bulk_left, dphi_drb_bulk_right;
  double  t_lj=0.0, t_hs=0.0, t_rhobar=0.0, t_uatt=0.0, t_charge=0.0, t_psi=0.0, 
          t_put=0.0, t_all=0.0, t_precalc=0.0; /* time counters */
  double  t_lj_max,t_hs_max=0.0, t_rhobar_max=0.0, t_uatt_max=0.0, t_charge_max=0.0, 
          t_psi_max=0.0, t_put_max=0.0, t_all_max=0.0, t_precalc_max=0.0; /* time counters */
  double  t_lj_min,t_hs_min=0.0, t_rhobar_min=0.0, t_uatt_min=0.0, t_charge_min=0.0,
          t_psi_min=0.0, t_put_min=0.0, t_all_min=0.0,t_precalc_min=0.0; /* time counters */
  double *mat_row=NULL; /* full storage of 1 row of matrix */
   double resid_hs1,resid_hs2,resid_old,resid_rhobars,resid_rhobarv,resid_uatt;
   double resid_ig,resid_vext,resid_mu,resid_charge,resid_poisson,resid_transport,resid_el;

  int loc_i_charge_up,loc_i_charge_down,loc_i_charge_up2,loc_i_charge_down2,blocked;
  double fac_temp,gradphi;

  double  nodepos[3];

  /* the 6 offset patterns for nearest neighbors */
  int offset_idim_pm[18] = {1,0,0,  0,1,0,  0,0,1,  -1,0,0,  0,-1,0,  0,0,-1};
  int *offset_ptr; /* pointer into above */
  /* Jacobian weight precalc stuff */

  static double ***jac_weights_hs=NULL;
  static int    ***jac_columns_hs=NULL;
  int max_len, *jac_col;
  int node_count=0;
  int l_elec_RTF, l_elec_LBB;

  int i_box, inode_box,jnode_box, ijk_box[3], ijk[3],ijk_tmp[3],loc_jnode;


  /********************** BEGIN EXECUTION ************************************/

  if (Proc == 0 && !resid_only_flag && Iwrite != NO_SCREEN) printf("\n\t%s: Doing fill of residual and matrix\n",yo);

  /* Allocate and precalculate rho_bars */

  /*fflush(NULL);*/

  if (fill_flag != MSR_PREPROCESS) {

    if (Sten_Type[DELTA_FN] && Sten_Type[THETA_FN]) {
      rho_bar = (struct RB_Struct *) array_alloc
                          (1, Nnodes_1stencil, sizeof(struct RB_Struct));
      t_precalc -=MPI_Wtime();
      pre_calc_rho_bar(rho_bar, x, fill_flag, iter, NULL, NULL,fill_time);
/*      if (Mesh_coarsening != FALSE && Nwall_type >0 || L1D_bc) pre_calc_coarse_rho_bar(rho_bar);*/

      dphi_drb = (struct RB_Struct *) array_alloc
                      (1, Nnodes_box, sizeof(struct RB_Struct));
 
      if (Matrix_fill_flag >= 3 && Ipot_ff_n != IDEAL_GAS){
         if (Type_func==0)
            pre_calc_dphi_drb_rb1(dphi_drb, x, &dphi_drb_bulk, 
                        &dphi_drb_bulk_left,
                        &dphi_drb_bulk_right,rho_bar,fill_time);
         else
            pre_calc_dphi_drb_rb2(dphi_drb, x, &dphi_drb_bulk, 
                        &dphi_drb_bulk_left,
                        &dphi_drb_bulk_right,rho_bar,fill_time);
      }
      else{
         if (Type_func==0){
            pre_calc_dphi_drb(dphi_drb, rho_bar, &dphi_drb_bulk, 
                        &dphi_drb_bulk_left,
                        &dphi_drb_bulk_right);
         }
         else
            pre_calc_dphi_drb2(dphi_drb, rho_bar, &dphi_drb_bulk, 
                        &dphi_drb_bulk_left,
                        &dphi_drb_bulk_right);
      }

      t_precalc += MPI_Wtime();

      /* for debugging print out profiles on each iteration */
      if (Iwrite==VERBOSE) {
         print_rho_bar(rho_bar, "rb.out");
         print_rho_bar(dphi_drb, "dphi.out");
         if (Proc == 0){
             X_old = (double *) array_alloc (1, Nnodes*Nunk_per_node, 
                                                     sizeof(double));
             Vext_old = (double *) array_alloc (1, Nnodes*Ncomp, sizeof(double));
         }
         collect_x_old(x,0);
         collect_vext_old();
         if (Proc==0) {
              print_profile("dens_iter.dat");
              safe_free((void *) &X_old);
              safe_free((void *) &Vext_old);
         }
      }
      
    }

    mat_row = (double *) array_alloc(1, Nunk_int_and_ext, sizeof(double));

    for (j=0; j<Nunk_int_and_ext; j++) mat_row[j] =0.0;

    /* rename saved columns from global to local numbering scheme */

    if ((fill_flag==JAC_SAVE_FILL_1 || fill_flag==JAC_SAVE_FILL_2) && 
          (iter ==1 || (iter==2 && Load_Bal_Flag == LB_TIMINGS)) ) {
       for (j=0; j < Nnodes_1stencil; j++) {
         jac_col = jac_columns_hs[j][0];
         for (i=1; i<jac_col[0]; i++) 
             jac_col[i] = B2L_unknowns[jac_col[i]];
         jac_col = jac_columns_hs[j][1];
         for (i=1; i<jac_col[0]; i++) 
             jac_col[i] = B2L_unknowns[jac_col[i]];
       }
    }

  }

  else {
    bindx_tmp = (int *) array_alloc(1, Nunknowns_box, sizeof(int));
    for (j=0; j<Nunknowns_box; j++) bindx_tmp[j] = FALSE;

    if (Sten_Type[DELTA_FN] && Sten_Type[THETA_FN]) {
       /* set up B2L_1stencil for list of nodes where rho_bars are calculated */

       B2L_1stencil = (int *) array_alloc(1, Nnodes_box, sizeof(int));
       for (i=0; i<Nnodes_box; i++) B2L_1stencil[i] = -1;

       pre_calc_rho_bar(NULL, x, fill_flag, iter, NULL, NULL,NULL); 

       Nnodes_1stencil = 0;
       for (i=0; i<Nnodes_box; i++) {
         if (B2L_1stencil[i] == HIT_FLAG) {
            B2L_1stencil[i] = Nnodes_1stencil;
            Nnodes_1stencil++;
         }
       }

       if (Matrix_fill_flag==JAC_SAVE_FILL_1 || Matrix_fill_flag==JAC_SAVE_FILL_2) {
          if (jac_weights_hs != NULL) {
            safe_free((void *) &jac_weights_hs);
            safe_free((void *) &jac_columns_hs);
          }

          /* calculate the maximum length of a rhobar row of the matrix */
          max_len = 0;
          i = 0;
          izone = 0;
          for (icomp=0; icomp<Ncomp; icomp++) {
               max_len +=  Stencil[THETA_FN][izone][icomp].Length;
               i +=  Stencil[DELTA_FN][izone][icomp].Length;
          }
          if (i > max_len) max_len = i;

          /* allocate 1 longer, so real length can be put in first entry */

          jac_weights_hs = (double ***)
             array_alloc(3, Nnodes_1stencil, 2, max_len+1, sizeof(double));
          jac_columns_hs = (int ***)
             array_alloc(3, Nnodes_1stencil, 2, max_len+1, sizeof(int));

          /* load up the jac_save arrays */
 
          pre_calc_rho_bar(NULL, x, RHOBAR_JAC_SAVE, iter,
                           jac_weights_hs, jac_columns_hs,NULL); 

       }
    }
  }

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

    /* start timer for this node */
    if (fill_time != NULL) fill_time[loc_inode] -= MPI_Wtime();

    /* convert local node to global */

    inode_box = L2B_node[loc_inode];
    node_box_to_ijk_box(inode_box, ijk_box);
    ijk_box_to_ijk(ijk_box,ijk);

    if ( ((Mesh_coarsening != FALSE) && (Nwall_type >0)) || L1D_bc)
      mesh_coarsen_flag_i = Mesh_coarsen_flag[inode_box];
    else
      mesh_coarsen_flag_i = 0;
   

    for (iunk=iunk_start; iunk<iunk_end; iunk++) {

/*if (mesh_coarsen_flag_i != FLAG_BULK && !Zero_density_TF[inode_box][iunk]){
        printf("Proc: %d loc_inode: %d of %d : inode_box=%d and mesh_coarsen_flag: %d\n",
                 Proc,loc_inode,Nnodes_per_proc,inode_box,mesh_coarsen_flag_i);
}*/
       /*resid_old=resid_old2=resid_hs=resid_hs1=resid_hs2=
                   resid_uatt=0.0;*/
       resid_ig=resid_vext=resid_mu=resid_charge=resid_poisson=
                   resid_transport=resid_rhobars=resid_rhobarv=0.0;

      /* i_box is the equation number (matrix row) that we are filling now */
      i_box = loc_find(iunk,inode_box,BOX);
      loc_i = loc_find(iunk,loc_inode,LOCAL);
      if (fill_flag != MSR_PREPROCESS)
            loc_i = Aztec.update_index[loc_i];

      if (mesh_coarsen_flag_i == FLAG_1DBC){
         node_box_to_ijk_box(inode_box,ijk_box);
         for (idim=0; idim<Ndim; idim++) ijk_tmp[idim]=0;
         ijk_tmp[Grad_dim]=ijk_box[Grad_dim];
         jnode_box = ijk_box_to_node_box(ijk_tmp);
         j_box=loc_find(iunk,jnode_box,BOX);
         loc_jnode = B2L_node[jnode_box];

         if (jnode_box <0 ){
             printf("Proc: %d ijk_box: %d %d %d PROBLEMS: jnode_box: %d  ijk_tmp: %d %d %d\n",
                Proc,ijk_box[0],ijk_box[1],ijk_box[2],jnode_box,ijk_tmp[0],ijk_tmp[1],ijk_tmp[2]);
                exit(-1);
         }

         if (fill_flag != MSR_PREPROCESS){
            loc_j = B2L_unknowns[j_box];
            resid[loc_i] = x[loc_i]-x[loc_j];
            mat_row[loc_i] =  1.0;
            mat_row[loc_j] = -1.0;
            put_row_in_msr(i_box, loc_i, Aztec.bindx, bindx_tmp, bindx_2d,
                                 Aztec.val, mat_row, fill_flag,x);
         }
         else{
            bindx_tmp[loc_find(iunk,jnode_box,BOX)] = TRUE;
            put_1Dsolution_in_msr(mesh_coarsen_flag_i,i_box,loc_i,jnode_box,bindx_2d);
         }
      } 

      /* do mesh coarsening if indicated .... for all unknowns ! */
      else if (mesh_coarsen_flag_i < 0 && mesh_coarsen_flag_i != FLAG_BULK) {
           if (fill_flag != MSR_PREPROCESS){

             /* Go to node 1 higher in the appropriate ijk direction */

             offset_ptr = &offset_idim_pm[3*(-mesh_coarsen_flag_i - 1)];
             jnode_box = offset_to_node_box(ijk_box, offset_ptr, reflect_flag);

             if (jnode_box >= 0) {
                loc_j = B2L_unknowns[loc_find(iunk,jnode_box,BOX)];
                resid[loc_i] = x[loc_i] - 0.5*x[loc_j];
                mat_row[loc_i] =  1.0;
                mat_row[loc_j] = -0.5;
             }
             else{
                resid[loc_i] = x[loc_i];
                mat_row[loc_i] =  1.0;
                resid[loc_i] -= 0.5*constant_boundary(iunk,jnode_box);

             }

             /* Go to node 1 lower in the appropriate ijk direction */

             offset_ptr += 9;  /* gets negative of offset used above */
             jnode_box = offset_to_node_box(ijk_box, offset_ptr, reflect_flag);

             if (jnode_box >= 0){
                loc_j = B2L_unknowns[loc_find(iunk,jnode_box,BOX)];
                resid[loc_i] -= 0.5*x[loc_j];
                mat_row[loc_j] = -0.5;
             }
             else{
                 resid[loc_i] -= 0.5*constant_boundary(iunk,jnode_box);
             }
             t_put -= MPI_Wtime();

             if (!resid_only_flag) 
                put_row_in_msr(i_box, loc_i, Aztec.bindx, bindx_tmp, bindx_2d,
                                 Aztec.val, mat_row, fill_flag,x);
             t_put += MPI_Wtime();
           }
           else {

             /* Go to node 1 higher in the appropriate ijk direction */

             offset_ptr = &offset_idim_pm[3*(-mesh_coarsen_flag_i - 1)];
             jnode_box = offset_to_node_box(ijk_box, offset_ptr, reflect_flag);

             if (jnode_box >= 0)
                bindx_tmp[loc_find(iunk,jnode_box,BOX)] = TRUE;
     
             /* Go to node 1 lower in the appropriate ijk direction */

             offset_ptr += 9;  /* gets negative of offset used above */
             jnode_box = offset_to_node_box(ijk_box, offset_ptr, reflect_flag);
             
             if (jnode_box >= 0)
                bindx_tmp[loc_find(iunk,jnode_box,BOX)] = TRUE;
     
             t_put -= MPI_Wtime();
                put_coarse_in_msr(mesh_coarsen_flag_i,i_box,loc_i,bindx_2d);
             t_put += MPI_Wtime();
           }
      }
      else{

      /* SET ZONE */
      /* izone is set based on value of mesh_coarsen_flag       */
      /* izone is the zone number at the node of interest which */
      /* indicated which quadrature scheme is to be used.       */

      izone = mesh_coarsen_flag_i;

      /**** LOAD EULER-LAGRANGE EQUATIONS *****/
      if (Unk2Phys[iunk]==DENSITY) {      
        icomp = iunk-Phys2Unk_first[DENSITY];

        /* Do trivial fill if it is a zero-density node */
        if (Zero_density_TF[inode_box][icomp] || Vext[loc_inode][icomp] == VEXT_MAX) {
           if (fill_flag != MSR_PREPROCESS){
             resid[loc_i] = x[loc_i];
             mat_row[loc_i] = 1.0;
           }
        }
        else {

 
           if (fill_flag != MSR_PREPROCESS){
             /* First load diagonal and constant terms to the residual */
 
             resid[loc_i] = log(x[loc_i]) ; 
             mat_row[loc_i] +=1.0/x[loc_i];
resid_ig=resid[loc_i];

             if (mesh_coarsen_flag_i == FLAG_BULK) resid[loc_i] -= log(Rho_b[icomp]);
             else {

/*             if(Lsteady_state)
                resid[loc_i] += (- 3.0*log(Sigma_ff[icomp][icomp]) 
                                  -1.5*log(Mass[icomp]*Temp));*/

             if (Nwall > 0) resid[loc_i] += Vext[loc_inode][icomp];
resid_vext=resid[loc_i]-resid_ig;
 
             if (Lsteady_state == FALSE) { 
               if (Iliq_vap < 10 ){
                  resid[loc_i] -= log(Rho_b[icomp]);
  
                  if (Ipot_ff_n != IDEAL_GAS) resid[loc_i] -= Betamu_hs_ex[icomp];
                  if (Ipot_ff_n == LJ12_6)    resid[loc_i] -= Betamu_att[icomp];
                  if (Ipot_ff_c == COULOMB && Sten_Type[THETA_CHARGE]) 
                                                 resid[loc_i] += Deltac_b[icomp];
/*                  if (fill_flag !=MSR_PREPROCESS) resid_old=resid[loc_i];*/
               }
               else resid[loc_i] -= Betamu[icomp];
             }
             else{ /* Lsteady_state == TRUE */
               unk_mu = loc_find(Phys2Unk_first[DIFFUSION] + icomp,loc_inode,LOCAL);
               loc_mu = Aztec.update_index[unk_mu];
               resid[loc_i] -= x[loc_mu];
               mat_row[loc_mu] -= 1.0;
             }
             }
resid_mu=resid[loc_i]-resid_ig-resid_vext;

             if (Ipot_ff_c == COULOMB){
                loc_i_charge = Aztec.update_index[loc_find(Phys2Unk_first[POISSON],loc_inode,LOCAL)];
                resid[loc_i]  += Charge_f[icomp]*x[loc_i_charge];
                mat_row[loc_i_charge] += Charge_f[icomp];
resid_charge=resid[loc_i]-resid_ig-resid_vext-resid_mu;

                if (Lpolarize[icomp]){

                 fac_temp=Temp_elec/(4.0*PI*KAPPA_H2O);
                 if (Nodes_2_boundary_wall[0][inode_box] != -1){
                    if (Surf_normal[0][loc_inode][0] < 0){
                           loc_i_charge_down = Aztec.update_index[loc_find(Phys2Unk_first[POISSON],loc_inode-1,LOCAL)];
                           loc_i_charge_down2 = Aztec.update_index[loc_find(Phys2Unk_first[POISSON],loc_inode-2,LOCAL)];
                           gradphi = 0.5*(3.0*x[loc_i_charge]-4.0*x[loc_i_charge_down]+x[loc_i_charge_down2])/Esize_x[0];

                           mat_row[loc_i_charge]       += 1.5*Pol[icomp]*gradphi*fac_temp/Esize_x[0];
                           mat_row[loc_i_charge_down]  -= 2.0*Pol[icomp]*gradphi*fac_temp/Esize_x[0];
                           mat_row[loc_i_charge_down2] += 0.5*Pol[icomp]*gradphi*fac_temp/Esize_x[0];
                    }
                    else{
                           loc_i_charge_up = Aztec.update_index[loc_find(Phys2Unk_first[POISSON],loc_inode+1,LOCAL)];
                           loc_i_charge_up2 = Aztec.update_index[loc_find(Phys2Unk_first[POISSON],loc_inode+2,LOCAL)];
                           gradphi = 0.5*(-3.0*x[loc_i_charge]+4.0*x[loc_i_charge_up]-x[loc_i_charge_up2])/Esize_x[0];

                           mat_row[loc_i_charge]     -= 1.5*Pol[icomp]*gradphi*fac_temp/Esize_x[0];
                           mat_row[loc_i_charge_up]  += 2.0*Pol[icomp]*gradphi*fac_temp/Esize_x[0];
                           mat_row[loc_i_charge_up2] -= 0.5*Pol[icomp]*gradphi*fac_temp/Esize_x[0];
                    }
                  }
                  else{
                      loc_i_charge_up = Aztec.update_index[loc_find(Phys2Unk_first[POISSON],loc_inode+1,LOCAL)];
                      loc_i_charge_down = Aztec.update_index[loc_find(Phys2Unk_first[POISSON],loc_inode-1,LOCAL)];
                      gradphi = 0.5*(x[loc_i_charge_up]-x[loc_i_charge_down])/Esize_x[0];

                      mat_row[loc_i_charge_up]   += 0.5*Pol[icomp]*gradphi*fac_temp/Esize_x[0];
                      mat_row[loc_i_charge_down] -= 0.5*Pol[icomp]*gradphi*fac_temp/Esize_x[0];
                  }
                  resid[loc_i] += 0.5*Pol[icomp]*gradphi*gradphi*fac_temp;

                } 
             
             }
           }
           else { /* MSR_PREPROCESS == TRUE */

             /* If transport, put in off-diagonal entry for mu in E-L equation */
             if (Lsteady_state == TRUE) {
               unk_mu = loc_find(Phys2Unk_first[DIFFUSION]+icomp,inode_box,BOX);
               bindx_tmp[unk_mu] = TRUE;
             }

             /* If charged, put in off-diagonal entry for Psi in E-L equation */
             if (Ipot_ff_c == COULOMB){
                  bindx_tmp[loc_find(Phys2Unk_first[POISSON],inode_box,BOX)] = TRUE;

                if (Lpolarize[icomp]){
                    if (Nodes_2_boundary_wall[0][inode_box] != -1){
                       if (Surf_normal[0][loc_inode][0]<0){
                          bindx_tmp[loc_find(Phys2Unk_first[POISSON],inode_box-1,BOX)] = TRUE;
                          bindx_tmp[loc_find(Phys2Unk_first[POISSON],inode_box-2,BOX)] = TRUE;
                       }
                       else{
                          bindx_tmp[loc_find(Phys2Unk_first[POISSON],inode_box+1,BOX)] = TRUE;
                          bindx_tmp[loc_find(Phys2Unk_first[POISSON],inode_box+2,BOX)] = TRUE;
                       }
                    }
                    else{
                       bindx_tmp[loc_find(Phys2Unk_first[POISSON],inode_box-1,BOX)] = TRUE;
                       bindx_tmp[loc_find(Phys2Unk_first[POISSON],inode_box+1,BOX)] = TRUE;
                    }
                }
             }
           }
        
        /* Now loop over all the stencils to load off diagonal contributions
           to the Jacobian. */

           if (Ipot_ff_n != IDEAL_GAS && mesh_coarsen_flag_i != FLAG_BULK) {
              t_hs -= MPI_Wtime();

              if (Matrix_fill_flag >= 3) {
                 load_nonlocal_hs_rosen_rb(DELTA_FN,loc_i,icomp,izone,
                               ijk_box,fill_flag,x,dphi_drb,
                               &dphi_drb_bulk, &dphi_drb_bulk_left,
                               &dphi_drb_bulk_right,resid,mat_row,bindx_tmp,
                               rho_bar,loc_inode, resid_only_flag);
/*                 if (fill_flag !=MSR_PREPROCESS) resid_hs1 = resid[loc_i]-resid_old;*/
                 load_nonlocal_hs_rosen_rb(THETA_FN,loc_i,icomp,izone,
                               ijk_box,fill_flag,x,dphi_drb,
                               &dphi_drb_bulk, &dphi_drb_bulk_left,
                               &dphi_drb_bulk_right,resid,mat_row,bindx_tmp,
                               rho_bar,loc_inode, resid_only_flag);
/*                 if (fill_flag !=MSR_PREPROCESS) resid_hs2 = resid[loc_i]-resid_old-resid_hs1;*/
              }
              else{
                 if (fill_flag !=MSR_PREPROCESS) resid_old=resid[loc_i];
                 t_lj += load_nonlocal_hs_rosen(DELTA_FN,loc_i,icomp,izone,ijk_box,
                                       fill_flag,rho_bar,dphi_drb,
                                       &dphi_drb_bulk, &dphi_drb_bulk_left,
                                       &dphi_drb_bulk_right,resid,mat_row,bindx_tmp,
                                       jac_weights_hs, jac_columns_hs,
                                       resid_only_flag);
                 if (fill_flag !=MSR_PREPROCESS) resid_hs1 = resid[loc_i]-resid_old;

                 if (fill_flag !=MSR_PREPROCESS) resid_old=resid[loc_i];
                 t_lj += load_nonlocal_hs_rosen(THETA_FN,loc_i,icomp,izone,ijk_box,
                                       fill_flag,rho_bar,dphi_drb,
                                       &dphi_drb_bulk, &dphi_drb_bulk_left,
                                       &dphi_drb_bulk_right,resid,mat_row,bindx_tmp,
                                       jac_weights_hs, jac_columns_hs,
                                       resid_only_flag);
                 if (fill_flag !=MSR_PREPROCESS) resid_hs2 = resid[loc_i]-resid_old;

              }

              t_hs += MPI_Wtime();

              if (Sten_Type[U_ATTRACT]) {
                 t_uatt -= MPI_Wtime();

                 if (fill_flag !=MSR_PREPROCESS) resid_old=resid[loc_i];
                 load_mean_field(U_ATTRACT,loc_i,icomp,izone,ijk_box,fill_flag,
                                 resid,mat_row,bindx_tmp,x, resid_only_flag);
                 if (fill_flag !=MSR_PREPROCESS) resid_uatt = resid[loc_i] - resid_old;
                 t_uatt += MPI_Wtime();

              }    /* end load attractions */
            
              if (Sten_Type[THETA_CHARGE]) {
                 t_charge -= MPI_Wtime();
                 load_mean_field(THETA_CHARGE,loc_i,icomp,izone,ijk_box,fill_flag,
                                 resid,mat_row,bindx_tmp,x, resid_only_flag);
                 t_charge += MPI_Wtime();
              }    /* end load electrostatic c(r) corrections */

              /* if total effective field is greater than VEXT_MAX *
               * ... solve the equation x=0 for this node          */
/*COMMENT OUT DYNAMIC ZEROING
              if (fill_flag != MSR_PREPROCESS && Sten_Type[U_ATTRACT]){
                 if (resid[loc_i]-log(x[loc_i]) > VEXT_MAX){ 
                    for (j=Aztec.bindx[loc_i]; j<Aztec.bindx[loc_i+1]; j++)
                        mat_row[Aztec.bindx[j]] = 0.0;
                        resid[loc_i] = x[loc_i]-Rho_b[icomp]*exp(-VEXT_MAX);
                        resid[loc_i] = x[loc_i];
                        mat_row[loc_i] = 1.0;
                  }
              }
*/
 
           }
        }

        /* now that fill of row i is finished, put in MSR format */

         if (fill_flag==MSR_PREPROCESS){
            if (Ipot_ff_n == IDEAL_GAS){
              t_put -= MPI_Wtime();
              put_euler_lag_in_msr(i_box,loc_i,bindx_2d);
              t_put += MPI_Wtime();
            }
            else{
              t_put -= MPI_Wtime();
              if ( (Unk2Phys[iunk] != DENSITY && icomp !=0 || loc_inode !=0) &&
                   ( Zero_density_TF[inode_box][icomp] || 
                    Vext[loc_inode][icomp] == VEXT_MAX ) )
                 put_zero_in_msr(loc_i,bindx_2d);
              else{

                 put_row_in_msr(i_box, loc_i, Aztec.bindx, bindx_tmp, bindx_2d,
                                     Aztec.val, mat_row, fill_flag,NULL);
              }
              t_put += MPI_Wtime();
            }
         }
         else{  /* Fill matrix */
              t_put -= MPI_Wtime();

              if (!resid_only_flag){
                 put_row_in_msr(i_box, loc_i, Aztec.bindx, bindx_tmp, bindx_2d,
                                     Aztec.val, mat_row, fill_flag,x);
              }

              t_put += MPI_Wtime();
         }
      }
      /***** END EULER-LAGRANGE *****/

      /*********************************/
      /**** LOAD RHO_BAR EQUATIONS *****/
      /*********************************/
      else if (Unk2Phys[iunk]==RHOBAR_ROSEN){
          t_rhobar -= MPI_Wtime();

          if (iunk == Phys2Unk_first[RHOBAR_ROSEN])
             load_rho_bar_s(THETA_FN,x,loc_i,iunk,loc_inode,izone,ijk_box, 
                            fill_flag,resid,mat_row,bindx_tmp, resid_only_flag);
          else if (iunk < Phys2Unk_first[RHOBAR_ROSEN]+Nrho_bar_s)
             load_rho_bar_s(DELTA_FN,x,loc_i,iunk,loc_inode,izone,ijk_box, 
                            fill_flag,resid,mat_row,bindx_tmp, resid_only_flag);
          if (fill_flag !=MSR_PREPROCESS) resid_rhobars = resid[loc_i] - resid_ig-resid_vext-resid_mu-resid_charge;

          t_rhobar += MPI_Wtime();
          t_put -= MPI_Wtime();
          if (!resid_only_flag)
             put_row_in_msr(i_box, loc_i, Aztec.bindx, bindx_tmp, bindx_2d,
                                  Aztec.val, mat_row, fill_flag,x);
          t_put += MPI_Wtime();

          if (iunk >= Phys2Unk_first[RHOBAR_ROSEN]+Nrho_bar_s){
              if (Matrix_fill_flag==3){
              t_rhobar -= MPI_Wtime();
              load_rho_bar_v(x,loc_i,iunk,loc_inode,izone,ijk_box, fill_flag,
                             resid,mat_row,bindx_tmp, resid_only_flag);
              t_rhobar += MPI_Wtime();
              t_put -= MPI_Wtime();
              if (!resid_only_flag)
                put_row_in_msr(i_box, loc_i, Aztec.bindx, bindx_tmp, bindx_2d,
                                     Aztec.val, mat_row, fill_flag,x);
              t_put += MPI_Wtime();
              }
              if (fill_flag !=MSR_PREPROCESS) 
                 resid_rhobarv = resid[loc_i] - resid_ig-resid_vext-resid_mu-resid_charge-resid_rhobars;
          }
      }
      /****** END RHOBAR *****/

      /***************************************************/
      /**** LOAD POISSON'S EQUATION FOR ELECTROSTATICS ***/
      /***************************************************/
      else if (Unk2Phys[iunk]==POISSON){

         if (fill_flag == MSR_PREPROCESS){  /* we know where all the columns are
                                             for poisson's equation so use a faster
                                             routine for the put_row_in_msr */
            t_put -= MPI_Wtime();
             put_poisson_in_msr(i_box,loc_i,bindx_2d);
            t_put += MPI_Wtime();
         }
         else {                             /* Fill Poisson's equation */

            l_elec_RTF=FALSE;
            l_elec_LBB=FALSE;
            if (Lsteady_state){
               node_to_position(node_box_to_node(inode_box),nodepos);
               idim = Grad_dim;
               if ( nodepos[idim] +0.5*Size_x[idim] - X_const_mu <= 0.00000001) l_elec_LBB=TRUE;
               else if (nodepos[idim] - 0.5*Size_x[idim] + X_const_mu >=-0.00000001) l_elec_RTF=TRUE;
            }

            if (l_elec_LBB){
               resid[loc_i] = x[loc_i]-Elec_pot_LBB;
               mat_row[loc_i] = 1.0;
            }
            if (l_elec_RTF){
               resid[loc_i] = x[loc_i]-Elec_pot_RTF;
               mat_row[loc_i] = 1.0;
            }
            else if (!l_elec_LBB && !l_elec_RTF) {
               t_psi -= MPI_Wtime();
               if(Type_coul==POLARIZE){
                 load_polarize_poissons_eqn(i_box, inode_box, loc_i, ijk_box, mat_row,
                                                       resid, x, bindx_tmp, fill_flag);
               }
               else load_poissons_eqn(i_box, inode_box, loc_i, ijk_box, mat_row,
                                    resid, x, bindx_tmp, fill_flag);
               load_poisson_bc(resid,inode_box,loc_inode,loc_i);

               resid_poisson=resid[loc_i]-resid_ig-resid_vext-resid_mu-
                                resid_charge-resid_rhobars-resid_rhobarv;
               t_psi += MPI_Wtime();
            }
            t_put -= MPI_Wtime();
            if (!resid_only_flag)
                put_row_in_msr(i_box, loc_i,Aztec.bindx, bindx_tmp, bindx_2d,
                               Aztec.val, mat_row, fill_flag,x);
            t_put += MPI_Wtime();
         }
      }
      /******** END POISSON ********/

      /*******************************/
      /*** LOAD TRANSPORT EQUATION ***/
      /*******************************/
      else if (Unk2Phys[iunk]==DIFFUSION){
 
         if (fill_flag == MSR_PREPROCESS){
            put_transport_in_msr(i_box,loc_i,bindx_2d);
         }
         else {

            if (Linear_transport)
                load_linear_transport_eqn(i_box, inode_box, loc_i, ijk_box, 
                           mat_row, resid, x, bindx_tmp, fill_flag, iunk);
            else   load_nonlinear_transport_eqn(i_box, inode_box, loc_i, ijk_box, 
                              mat_row, resid, x, bindx_tmp, fill_flag, iunk);
 
            if (!resid_only_flag)
                 put_row_in_msr(i_box, loc_i, Aztec.bindx, bindx_tmp, bindx_2d,
                                   Aztec.val, mat_row, fill_flag,x);
            resid_transport=resid[loc_i]-resid_ig-resid_vext-resid_mu-resid_charge-resid_poisson;
         }
      }
      /******** END TRANSPORT ********/
      else if (Unk2Phys[iunk]==DENSITY_SEG){
      }
      else if (Unk2Phys[iunk]==CAVITY_WTC){
      }
      else if (Unk2Phys[iunk]==BOND_WTC){
      }

      else {
         printf("Problem with unknowns in fill !!!\n");
         exit (-1);
      }
      }
     
/*if (fill_flag !=MSR_PREPROCESS) {
    if (Unk2Phys[iunk]==DENSITY){
       resid_el = resid_ig + resid_vext + resid_mu + resid_charge;
       printf("inode_box=%d  iunk=%d  loc_i=%d  resid_el=%9.6f  resid_hs=%9.6f resid=%9.6f",
                                 inode_box,iunk,loc_i,resid_el,resid[loc_i]-resid_el,resid[loc_i]);
    }
    else if (Unk2Phys[iunk]==RHOBAR_ROSEN)
         printf("inode_box=%d : iunk_rbar=%d loc_i=%d [  %g  %g] : %9.6f",
                 inode_box,iunk,loc_i,resid_rhobars,resid_rhobarv,resid[loc_i]);
    else if (Unk2Phys[iunk]==POISSON)   
         printf(" inode_box=%d  iunk_poisson=%d   resid=%9.6f ",
                   inode_box,iunk,resid_poisson);
    else if (Unk2Phys[iunk]==DIFFUSION) 
         printf(" inode_box=%d  iunk_diffusion=%d  resid=%9.6f",inode_box,iunk,resid_transport);
    printf("  \n");
}*/

    } /* end of loop over # of unknowns per node */

    if (fill_time != NULL) {
       fill_time[loc_inode] += MPI_Wtime();
    }
  } /* end of loop over local nodes */


  if (Iwrite==VERBOSE && fill_flag != MSR_PREPROCESS){
     if (Proc == 0){
         X_old = (double *) array_alloc (1, Nnodes*Nunk_per_node, 
                                                 sizeof(double));
         Vext_old = (double *) array_alloc (1, Nnodes*Ncomp, sizeof(double));
     }
     collect_x_old(resid,0);
     collect_vext_old();
     if (Proc==0) {
         print_profile("resid_iter.dat");
         safe_free((void *) &X_old);
         safe_free((void *) &Vext_old);
     }
  }


  /* Now add in constraint equation(s)  */
/*  if (Proc==0 && Lstoichiometry) {
     for (icomp = 0; icomp<Ncomp; icomp++){

         
         load_stoichimetric_constraint(i_box, inode_box, loc_i, ijk_box, mat_row,
                                           resid, x, bindx_tmp, fill_flag, icomp); 
         t_put -= MPI_Wtime();
         put_row_in_msr(i_box, loc_i, Aztec.bindx, bindx_tmp, bindx_2d,
                                        Aztec.val, mat_row, fill_flag,x);
         t_put += MPI_Wtime();
     }
  }*/


/*  if (fill_flag != MSR_PREPROCESS && Ipot_wf_c) {
     t_psi -= MPI_Wtime();
     load_poisson_bc(resid);
     t_psi += MPI_Wtime();
  }*/

  
  if (Mesh_coarsening == FALSE){ 
  t_all = t_precalc + t_hs + t_uatt + t_charge + t_psi + t_put;
  if (Matrix_fill_flag >=3 && Ipot_ff_n !=IDEAL_GAS) t_all += t_rhobar;

   /* print out load balancing info */

  if (!resid_only_flag){
    t_precalc_max     = AZ_gmax_double(t_precalc,Aztec.proc_config);
    t_hs_max          = AZ_gmax_double(t_hs,Aztec.proc_config);
    if (Matrix_fill_flag>=3 && Ipot_ff_n !=IDEAL_GAS)  t_rhobar_max      = AZ_gmax_double(t_rhobar,Aztec.proc_config);
    else                      t_lj_max          = AZ_gmax_double(t_lj,Aztec.proc_config);
    t_uatt_max        = AZ_gmax_double(t_uatt,Aztec.proc_config);
    t_charge_max      = AZ_gmax_double(t_charge,Aztec.proc_config);
    t_psi_max         = AZ_gmax_double(t_psi,Aztec.proc_config);
    t_put_max         = AZ_gmax_double(t_put,Aztec.proc_config);
    t_all_max         = AZ_gmax_double(t_all,Aztec.proc_config);

    t_precalc_min     = AZ_gmin_double(t_precalc,Aztec.proc_config);
    t_hs_min          = AZ_gmin_double(t_hs,Aztec.proc_config);
    if (Matrix_fill_flag>=3 && Ipot_ff_n !=IDEAL_GAS)  t_rhobar_min      = AZ_gmin_double(t_rhobar,Aztec.proc_config);
    else                      t_lj_min          = AZ_gmin_double(t_lj,Aztec.proc_config);
    t_uatt_min        = AZ_gmin_double(t_uatt,Aztec.proc_config);
    t_charge_min      = AZ_gmin_double(t_charge,Aztec.proc_config);
    t_psi_min         = AZ_gmin_double(t_psi,Aztec.proc_config);
    t_put_min         = AZ_gmin_double(t_put,Aztec.proc_config);
    t_all_min         = AZ_gmin_double(t_all,Aztec.proc_config);
  }

  if (fill_flag != MSR_PREPROCESS) { 
      T_av_precalc_max    += t_precalc_max;
      T_av_fill_max += t_all_max;

      T_av_precalc_min    += t_precalc_min;
      T_av_fill_min  += t_all_min;
  }

/*  if (iter == 1){
 *    Bin_size = (t_all_max-t_all_min)/100.;
 *    for (i=0; i<200; i++) Hist_time[i] = 0;
 *    Time_min_avg = 0.0;
 * }

 * else if (iter > 2) {
 *   if (Bin_size != 0.0) {
 *      ibin = (int)((t_all-t_all_min)/Bin_size);
 *      Hist_time[ibin]++;
 *      Time_min_avg += t_all_min;
 *   }
 * }
 */

  if (Proc == 0 && !resid_only_flag && Iwrite ==VERBOSE) {
    if (Ipot_ff_n !=IDEAL_GAS) printf( "\t\tTime for Precalcs (max): %g   (min) %g secs\n",t_precalc_max,t_precalc_min);

    if (Ipot_ff_n !=IDEAL_GAS) printf("\t\tTime loading euler-lagrange equations (max): %g   (min): %g secs\n",t_hs_max,t_hs_min);

    if (Matrix_fill_flag >=3 && Ipot_ff_n !=IDEAL_GAS)
         printf("\t\tTime loading rhobar equations (max): %g   (min): %g secs\n",t_rhobar_max,t_rhobar_min);

    if (Sten_Type[U_ATTRACT])
      printf("\t\tTime loading U_Attract term (max): %g    (min): %g secs\n",t_uatt_max,t_uatt_min);

    if (Sten_Type[THETA_CHARGE])
      printf("\t\tTime loading Delta_C term    (max): %g    (min): %g secs\n",t_charge_max,t_charge_min);

    if (Ipot_wf_c)
      printf("\t\tTime loading Poisson's Equation (max): %g    (min): %g secs\n", t_psi_max,t_psi_min);

    printf("\t\tTime in put_row_in_msr = (max): %g    (min): %g secs\n", t_put_max,t_put_min);

    printf("\n\t\tTotal load time (max): %g   (min): %g secs\n",t_all_max,t_all_min);
  }
  }
  if (fill_flag != MSR_PREPROCESS) {
     if (Ipot_ff_n != IDEAL_GAS) {
       safe_free((void *) &rho_bar);
       safe_free((void *) &dphi_drb);
     }
     safe_free((void *) &mat_row);
  }
  else {
       safe_free((void *) &bindx_tmp);
  }
}
/****************************************************************************/
