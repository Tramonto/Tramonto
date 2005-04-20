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
double load_polymer_cr(int, int, int, int, int *, int,
                     double *, double *, int *, double *);
double load_polymer_G(int,int,int,int,int,int *,int,double *,int *,double *);

#define HIT_FLAG 999
#define FIELD_MAX 20.
#define BOLTZ_MIN 0.000001
/****************************************************************************/
void fill_resid_and_matrix_P (double *x, double *resid,
                            int **bindx_2d, double *fill_time, int fill_flag,
                            int iter, int resid_only_flag, int unk_flag)
{
 /*
  * Local variable declarations
  */

  int     i, j, iunk,unk_GQ,unk_B,loc_B,iunk_start,iunk_end;
  int     reflect_flag[3],junk;
  int     izone, mesh_coarsen_flag_i;
  int    *bindx_tmp=NULL;/*Full storage of 1 row of bindx for MSR_PREPROCESS*/
  int     loc_i, loc_j=0, loc_l=0,loc_inode, itype_mer,loc_i_charge;
  int     jtype_mer;
  int     nseg,ibond,pol_num,seg_num,bond_num;

  double   t_put=0.0, t_all=0.0 ; /* time counters */
  double   t_put_max=0.0, t_all_max=0.0; /* time counters */
  double   t_put_min=0.0, t_all_min=0.0; /* time counters */
  double *mat_row=NULL; /* full storage of 1 row of matrix */
  double resid_B=0.0,resid_R=0.0,resid_G=0.0,resid_P=0.0,row_fact,k_poly,boltz=0.0;
  double nodepos[3];
 
  double gint_tmp;
  int ipol,iseg;
  int boltz_pow,boltz_pow_J;
  double fac1,fac2;
  int jbond,unk_GQ_j,loc_GQ_j,loc_GQ,node_start;

  double  resid_tmp=0.0;

  /* the 6 offset patterns for nearest neighbors */
  int offset_idim_pm[18] = {1,0,0,  0,1,0,  0,0,1,  -1,0,0,  0,-1,0,  0,0,-1};
  int *offset_ptr; /* pointer into above */

  int i_box, inode_box,jnode_box, ijk_box[3], max_field[NCOMP_MAX];
  int npol=0,sten;

  int inode,ijk[3];

  /********************** BEGIN EXECUTION ************************************/

  if (Sten_Type[POLYMER_GAUSS]) sten=POLYMER_GAUSS;
  else                          sten=DELTA_FN;

  if (fill_flag==MSR_PREPROCESS) {
    bindx_tmp = (int *) array_alloc(1, Nunknowns_box, sizeof(int));
    for (j=0; j<Nunknowns_box; j++) bindx_tmp[j] = FALSE;
  }
  else{
    mat_row = (double *) array_alloc(1, Nunk_int_and_ext, sizeof(double));
    for (j=0; j<Nunk_int_and_ext; j++) mat_row[j] =0.0;
    /* if (Proc==0 && Iwrite == VERBOSE)  fp31 = fopen ("resid.out","w");*/
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

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {    /*************/

    /* start timer for this node */
    if (fill_time != NULL) fill_time[loc_inode] -= MPI_Wtime();
    /* convert local node to global */

    inode_box = L2B_node[loc_inode];
    node_box_to_ijk_box(inode_box, ijk_box);

    if (Mesh_coarsening && Nwall_type >0) mesh_coarsen_flag_i = Mesh_coarsen_flag[inode_box];
    else                                  mesh_coarsen_flag_i = 0;

    if (Sten_Type[POLYMER_CR])
        for (iunk=0; iunk<Ncomp; iunk++) max_field[iunk] = FALSE;

    for (iunk=iunk_start; iunk<iunk_end; iunk++) {                /*************/

    resid_B=0.0;resid_R=0.0;resid_G=0.0;resid_P=0.0;

      /* i_box is the equation number (matrix row) that we are filling now */

      i_box = loc_find(iunk,inode_box,BOX);
      loc_i = loc_find(iunk,loc_inode,LOCAL);


      if (fill_flag != MSR_PREPROCESS){
             loc_i = Aztec.update_index[loc_i];
             resid[loc_i]=0.0;                      /*initialize resid to zero*/
      }

      if (Unk_to_eq_type[iunk]==DENSITY) itype_mer = iunk-Unk_start_eq[DENSITY]; 
      else if (Unk_to_eq_type[iunk] == CMS_FIELD) itype_mer = iunk - Unk_start_eq[CMS_FIELD];
      else if  (Unk_to_eq_type[iunk] == CMS_G) {                   
                               npol = Unk_to_Poly[iunk-Unk_start_eq[CMS_G]];
                               nseg = Unk_to_Seg[iunk-Unk_start_eq[CMS_G]]; 
                               itype_mer = Type_mer[npol][nseg];
      }

                                                               /*DO ZERO DENSITY FILL*/
      if ( fill_flag != MSR_PREPROCESS && ( ((Zero_density_TF[inode_box][itype_mer] ||    
            Vext[loc_inode][itype_mer] == VEXT_MAX) && Unk_to_eq_type[iunk]!=POISSON)
            /*|| (Unk_to_eq_type[iunk]==DENSITY 
               && -log(x[loc_find(Unk_start_eq[CMS_FIELD]+itype_mer,loc_inode,LOCAL)]) > VEXT_MAX )*/) ){
           resid[loc_i] = x[loc_i];
           mat_row[loc_i] = 1.0;
      }
/*      else if (fill_flag != MSR_PREPROCESS && Unk_to_eq_type[iunk]==CMS_FIELD 
               && x[loc_i] < exp(-VEXT_MAX) ){
               resid[loc_i] = x[loc_i]-exp(-VEXT_MAX);
               mat_row[loc_i] = 1.0;
      }*/
      /* ALF:  check this code sometime! */
      else if (mesh_coarsen_flag_i < 0) {                   /*DO COARSENED MESH FILL*/
         if (fill_flag != MSR_PREPROCESS){

           /* Go to node 1 higher in the appropriate ijk direction */

           offset_ptr = &offset_idim_pm[3*(-mesh_coarsen_flag_i - 1)];
           jnode_box = offset_to_node_box(ijk_box, offset_ptr, reflect_flag);

           resid[loc_i] = x[loc_i];
           mat_row[loc_i] =  1.0;
           if (jnode_box >= 0) {
              loc_j = B2L_unknowns[loc_find(iunk,jnode_box,BOX)];
              resid[loc_i] -= 0.5*x[loc_j];
              mat_row[loc_j] = -0.5;
           }
           else if(jnode_box == -1 ||jnode_box==-3 ||jnode_box==-4){  /* BOUNDARY */
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
           else if(jnode_box == -1 ||jnode_box==-3 ||jnode_box==-4){
              resid[loc_i] -= 0.5*constant_boundary(iunk,jnode_box);
           }
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

         }
      }
      else {                                           /* FILL REAL EQUATIONS */

         /* izone is set based on value of mesh_coarsen_flag       */
         /* izone is the zone number at the node of interest which */
         /* indicated which quadrature scheme is to be used.       */

         izone = mesh_coarsen_flag_i;

	 /* ALF: modify field eqns for SCF option, Type_poly==3  */
         if (Unk_to_eq_type[iunk] == CMS_FIELD){             /*BOLTZMANN EQNS*/

	    if(Type_poly != 3) {
	      row_fact = load_polymer_cr(POLYMER_CR,loc_i,itype_mer,
		 izone,ijk_box,fill_flag,resid,mat_row,bindx_tmp,x); 
	      /* row_fact = 1.0;*/
	   }
	   else {    /* ALF: new code */
	     for (jtype_mer=0; jtype_mer<Ncomp; jtype_mer++){
	       if (fill_flag != MSR_PREPROCESS) {
		 loc_l = Aztec.update_index[loc_find(jtype_mer + Unk_start_eq[DENSITY],loc_inode,LOCAL)];
		 resid[loc_i] -= Rism_cr[itype_mer][jtype_mer][0]*(x[loc_l]-Rho_b[jtype_mer]);
		 mat_row[loc_l] -= Rism_cr[itype_mer][jtype_mer][0];
	       }
	       else{ /* MSR_PREPROCESS  */
		 bindx_tmp[loc_find(jtype_mer + Unk_start_eq[DENSITY],inode_box,BOX)] = TRUE;
               }
	     }
	     row_fact = 1.;
	   }

            if (fill_flag != MSR_PREPROCESS){   /* diagonal terms */
              resid[loc_i] += Vext[loc_inode][itype_mer];
              resid[loc_i] += log(x[loc_i]);      
              mat_row[loc_i] += row_fact/x[loc_i];  
              if (Ipot_ff_c == COULOMB){
                    loc_i_charge = Aztec.update_index[loc_find(Unk_start_eq[POISSON],loc_inode,LOCAL)];
                    resid[loc_i]  += Charge_f[itype_mer]*x[loc_i_charge];
                    mat_row[loc_i_charge] += Charge_f[itype_mer];
              }

/*DEBUG .... NO max_field */
/*             if (resid[loc_i] < FIELD_MAX){
                  resid[loc_i] += log(x[loc_i]);      
                  mat_row[loc_i] += row_fact/x[loc_i];  * won't work in parallel *
               }
               else {
                  max_field[itype_mer] = TRUE;
                  if (fill_flag != MSR_PREPROCESS){
                     resid[loc_i] = x[loc_i] - BOLTZ_MIN;
                     mat_row[loc_i] = 1.0;
                  }
               }
*/
                resid_B=resid[loc_i]; 
            }
            else{
               if (Ipot_ff_c==COULOMB)  bindx_tmp[loc_find(Unk_start_eq[POISSON],inode_box,BOX)] = TRUE;
            }
         }                                  /* END BOLTZMAN EQNS */

         else if (Unk_to_eq_type[iunk]==DENSITY){   /* DENSITY EQNS FOR POLYMER SEGS*/

            if (fill_flag != MSR_PREPROCESS && max_field[itype_mer]){
               resid[loc_i] = x[loc_i];
               mat_row[loc_i] = 1.0;
            }

            /* reset one density in the domain to a constant value */
            if (fill_flag !=MSR_PREPROCESS && iunk == Ncomp 
                         && L2G_node[loc_inode]==Bupdate_iters){
                resid[loc_i] = x[loc_i]-Bupdate_fact;
                mat_row[loc_i] = 1.;
            }
            else {

            if (fill_flag != MSR_PREPROCESS){
	      /* Boltzmann factor for this i */
               loc_B = Aztec.update_index[loc_find(itype_mer + Unk_start_eq[CMS_FIELD],loc_inode,LOCAL)];
               if (Type_poly == 2 || Type_poly == 1) {
                 resid[loc_i] += x[loc_i];
                 mat_row[loc_i] += 1.;
               }
               else{
                 resid[loc_i] += x[loc_i]*x[loc_B];
                 mat_row[loc_i] += x[loc_B]; 
                 mat_row[loc_B] += x[loc_i]; 
               }
            }
            else{
               bindx_tmp[loc_find(itype_mer+Unk_start_eq[CMS_FIELD],inode_box,BOX)] = TRUE;
            }

            npol = 0;
            while (Nmer_t[npol][itype_mer]==0) npol++;
            k_poly = Rho_b[itype_mer]/Nmer_t[npol][itype_mer];

            for (i=0; i<Nmer[npol]; i++){
              if (Type_mer[npol][i] == itype_mer) {

                 if (fill_flag != MSR_PREPROCESS) {
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
                      loc_GQ = Aztec.update_index[loc_find(unk_GQ,loc_inode,LOCAL)]; 
                      fac1 *= x[loc_GQ];


                      fac2= k_poly;
                      for (jbond=0; jbond<Nbond[npol][i]; jbond++) {
                        if (jbond != ibond){       
                          unk_GQ_j  = Geqn_start[npol] + Poly_to_Unk[npol][i][jbond]; 
                          loc_GQ_j = Aztec.update_index[loc_find(unk_GQ_j,loc_inode,LOCAL)]; 
                          fac2 *= x[loc_GQ_j];
                        }
                      }
                      mat_row[loc_GQ] -= fac2*POW_DOUBLE_INT(x[loc_B],boltz_pow);
                   }
                   mat_row[loc_B] -= fac1*((double) boltz_pow)*POW_DOUBLE_INT(x[loc_B],boltz_pow_J);
                   resid[loc_i] -= fac1*POW_DOUBLE_INT(x[loc_B],boltz_pow);
                 }
                 else{
                    for (ibond=0; ibond<Nbond[npol][i]; ibond++){
                         unk_GQ  = Geqn_start[npol] + Poly_to_Unk[npol][i][ibond]; 
                         bindx_tmp[loc_find(unk_GQ,inode_box,BOX)] = TRUE;
                    }
                 }
              }
           }

           }  /* not filling constant density node */
            if (fill_flag != MSR_PREPROCESS){ resid_R=resid[loc_i]-resid_B;
            }
         }                                  /* END POLYMER DENSITY EQNS */


/* I don't think we use this code any more... will verify before deleting code 
         else if (iunk < 2*Ncomp){          * DENSITY EQNS FOR MONO-SOLVENT 
            if (fill_flag!=MSR_PREPROCESS){
               if (max_field[itype_mer] && fill_flag != MSR_PREPROCESS){
                  resid[loc_i] = x[loc_i];
                  mat_row[loc_i] = 1.0;
               }
               loc_B = Aztec.update_index[loc_find(itype_mer,loc_inode,LOCAL)];
               resid[loc_i] = x[loc_i] - Rho_b[itype_mer]*x[loc_B];
               mat_row[loc_i] += 1.;
               mat_row[loc_B] -= Rho_b[itype_mer];
            }
            else{  
                bindx_tmp[loc_find(itype_mer,inode_box,BOX)] = TRUE;
            }      
         }                                   * END SOLVENT DENSITY EQNS *
**********************************************************************/

         else if (Unk_to_eq_type[iunk] == CMS_G){          /* G EQNS */

            if (max_field[itype_mer] && fill_flag != MSR_PREPROCESS){  /* ACCOUNT FOR MAX FIELD */
               resid[loc_i] = x[loc_i];
               mat_row[loc_i] = 1.0;
            }
            else {

             unk_GQ=iunk-Unk_start_eq[CMS_G];
             pol_num = Unk_to_Poly[unk_GQ];
             seg_num = Unk_to_Seg[unk_GQ];
             bond_num = Unk_to_Bond[unk_GQ];
      
             if (Pol_Sym[unk_GQ] != -1){                            /* FILL IN G's FOR SYMMETRIC BONDS */
                junk = Pol_Sym[unk_GQ] + Unk_start_eq[CMS_G];
                if (fill_flag == MSR_PREPROCESS){
                  bindx_tmp[loc_find(junk,inode_box,BOX)] = TRUE;
                }
                else{
                  loc_j = Aztec.update_index[loc_find(junk,loc_inode,LOCAL)];
                  resid[loc_i] = x[loc_i]-x[loc_j];
                  mat_row[loc_i] +=1.;
                  mat_row[loc_j] -=1.;
                }
             } 
             else{                                                  /* FILL IN G's FOR UNIQUE BONDS */

               if (Bonds[pol_num][seg_num][bond_num] == -1){        /* fill end segment equation */
                  if (Type_poly == 0 || Type_poly == 1 || Type_poly == 3){
                     node_to_position(B2G_node[inode_box],nodepos);
                   /*  if (fabs(nodepos[0])>100.0 && pol_num==0){ exclude end segments from the outer regions of the box 
                           if (fill_flag != MSR_PREPROCESS){
                             resid[loc_i] = x[loc_i];
                             mat_row[loc_i] += 1.;
                           }
                     }
                     else{*/
                                                         /* Boltz unk for 1st seg */
                        unk_B = Unk_start_eq[CMS_FIELD] + Type_mer[npol][seg_num];     
                        if (fill_flag == MSR_PREPROCESS){
                           bindx_tmp[loc_find(unk_B,inode_box,BOX)] = TRUE;
                        }
                        else {
                           loc_j = Aztec.update_index[loc_find(unk_B,loc_inode,LOCAL)];
                           resid[loc_i] = x[loc_i] - x[loc_j];
                           mat_row[loc_i] += 1.;
                           mat_row[loc_j] -= 1.;
                        }
                     /*}*/
                  }
                  else{
                     if (fill_flag != MSR_PREPROCESS){
                        resid[loc_i] = x[loc_i] - 1.0;
                        mat_row[loc_i] += 1.;
                     }
                  }

               }
               else{                                            /* fill G_seg eqns */

                                 /* First calculate the residual contributions */

                  unk_B = Unk_start_eq[CMS_FIELD] + itype_mer;   /* Boltz unk for this seg */
                  if (fill_flag == MSR_PREPROCESS){
                     if (Type_poly !=2) bindx_tmp[loc_find(unk_B,inode_box,BOX)] = TRUE;
                  }
                  else {
                     loc_B = Aztec.update_index[loc_find(unk_B,loc_inode,LOCAL)];
                     if (Type_poly == 0 || Type_poly == 3 ){
                        boltz = x[loc_B];
                        resid[loc_i] = x[loc_i]; 
                        mat_row[loc_i] += 1.;
                     }
                     else if (Type_poly == 1) {
                        boltz = 1.0;
                        resid[loc_i] = x[loc_i]/x[loc_B];
                        mat_row[loc_i] += 1./x[loc_B];
                        mat_row[loc_B] -= x[loc_i]/(x[loc_B]*x[loc_B]);
                      }
                      else if (Type_poly == 2) {
                        boltz = 1.0;
                        resid[loc_i] = x[loc_i];
                        mat_row[loc_i] += 1.;
                      }
                  }
                  /* Now Finish loading the Jacobian... */
                  gint_tmp = load_polymer_G(sten,iunk,loc_B,itype_mer, izone,ijk_box,
                                                        fill_flag,mat_row,bindx_tmp,x);
                  if (fill_flag != MSR_PREPROCESS){
                    resid[loc_i] += gint_tmp;
                    if (Type_poly==0 || Type_poly==3) mat_row[loc_B] += gint_tmp / boltz;
                  }

               }
             }
            }
            if (fill_flag != MSR_PREPROCESS){ resid_G=resid[loc_i]-resid_B-resid_R;
            }
	    /*	    if(Proc==0) {
	      if(resid[loc_i] > 10000)
		printf("large resid in G eqns\n");
		}*/
         }                                  /* END G EQNS */
         else if (Unk_to_eq_type[iunk] == POISSON){ /* LOAD POISSON'S EQUATION*/
              load_poissons_eqn(i_box, inode_box, loc_i, ijk_box, mat_row,
                                              resid, x, bindx_tmp, fill_flag);
              if (fill_flag !=MSR_PREPROCESS){
                   load_poisson_bc(resid,inode_box,loc_inode,loc_i);
                   resid_P = resid[loc_i] -resid_B-resid_R-resid_G;
              }
         }
         else {
            printf("Problem with unknowns in fill !!!\n");
            exit (-1);
         }
      }     /* end of else (not Zero_density and mesh_coarsen_flag_i >= 0) */

     if (fill_flag != MSR_PREPROCESS && fabs(resid[loc_i])>.000001 && resid_P != 0. ){
      printf("loc_i=%d inode %d : iunk %d : resid_B %g : resid_R %g : resid_G %g : resid_P %g resid_tot=%g\n",
               loc_i,L2G_node[loc_inode],iunk,resid_B,resid_R,resid_G,resid_P,resid[loc_i]);

/*       printf("loc_inode=%d  iunk=%d loc_i=%d ",loc_inode,iunk,loc_i); 
       for (i=0;i<Nnodes_box;i++){
          for (j=0;j<Nunk_per_node;j++){
             loc_j = Aztec.update_index[loc_find(j,i,BOX)];
             if (fabs (mat_row[loc_j]) >= 1.e-12) printf("j=%d mr=%9.6f  ",loc_j,mat_row[loc_j]); 
          }
       }*/
       printf("\n");
     }
  
      /* now that fill of row i is finished, put in MSR format */
      if (!resid_only_flag){
      t_put -= MPI_Wtime();
          put_row_in_msr(i_box, loc_i, Aztec.bindx, bindx_tmp, bindx_2d,
                             Aztec.val, mat_row, fill_flag,x);
      t_put += MPI_Wtime();
      }

/*      if (fill_flag != MSR_PREPROCESS ) 
           printf(" %d   %d  %d  %g  \n",loc_i,inode_box,iunk,resid[loc_i]);*/
    } /* end of loop over # of unknowns per node */

    if (fill_time != NULL)  fill_time[loc_inode] += MPI_Wtime();
  } /* end of loop over local nodes */
/*  if (fill_flag != MSR_PREPROCESS) exit(-1);*/

  t_all = t_put;

  /* print out load balancing info */

  if (!resid_only_flag){
  t_put_max         = AZ_gmax_double(t_put,Aztec.proc_config);
  t_all_max         = AZ_gmax_double(t_all,Aztec.proc_config);
  t_put_min         = AZ_gmin_double(t_put,Aztec.proc_config);
  t_all_min         = AZ_gmin_double(t_all,Aztec.proc_config);
  }

  if (fill_flag != MSR_PREPROCESS) { 
      T_av_fill_max += t_all_max;
      T_av_fill_min  += t_all_min;
  }

  if (fill_flag != MSR_PREPROCESS) {
                                   safe_free((void *) &mat_row);
  }
  else                             safe_free((void *) &bindx_tmp);
  
  if (Proc == 0) {
    if (Ipot_wf_c)
    printf("\t\tTime in put_row_in_msr = (max): %g    (min): %g secs\n", t_put_max,t_put_min);
    if (!resid_only_flag && Iwrite==VERBOSE && MATRIX_FILL_NODAL) 
         printf("\n\t\tTotal load time (max): %g   (min): %g secs\n",t_all_max,t_all_min);
  }
}
/*****************************************************************************/
/* load_polymer_G: 
             In this routine we load the Jacobian for 
             gaussian (or freely-jointed) polymer terms.                      */
double load_polymer_G(int sten_type,int iunk,int loc_B,int itype_mer,
                    int izone,int *ijk_box, int fill_flag,
                    double *mat_row,int *bindx_tmp,double *x)
{
  int   **sten_offset, *offset, isten;
  double *sten_weight,  weight;
  struct Stencil_Struct *sten;

  int jlist,unk_GQ,nunk,unk[20],jseg,jtype_mer;
  int boltz_pow_1,boltz_pow_2,boltz_pow_R;
  double boltz_prefac_1,boltz_prefac_2,boltz_prefac_R,jac_entry,resid_entry,resid_sum,fac1,fac2;
  int reflect_flag[NDIM_MAX];
  int i,j,jnode_box;
  int pol_num,seg_num,bond_num,ibond,loc_i2_j2,box_i2_j2;

  if (Nlists_HW <= 2) jlist = 0;
  else                jlist = itype_mer; 
  for (i=0; i<20; i++) unk[i]=-1;

  /* As in precalc routine, we first need to assemble lists for the various G 
     products that we care about ... note that there are two flavors of the
     products now.  one excludes k(beta)=i, the other also excludes k(beta=i2)*/

   /* (1) from iunk get the segment number i and bond number*/
  unk_GQ = iunk-Unk_start_eq[CMS_G];
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
        unk[nunk++]     = Poly_to_Unk[pol_num][jseg][ibond]+Unk_start_eq[CMS_G];
      }
  }
  /* don't forget Boltzman factor */
  unk[nunk++]=jtype_mer;

  if (fill_flag != MSR_PREPROCESS){
     if (Type_poly == 0 || Type_poly == 3) {
         boltz_pow_R = -(Nbond[pol_num][jseg]-2);
         boltz_prefac_R = -x[loc_B];
         boltz_pow_1 = -(Nbond[pol_num][jseg]-1);
         boltz_prefac_1 = x[loc_B]*(Nbond[pol_num][jseg]-2);
         boltz_pow_2 = -(Nbond[pol_num][jseg]-2);
         boltz_prefac_2 = -x[loc_B];
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
        if (Lhard_surf && fill_flag != MSR_PREPROCESS) {
           if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1) 
           weight = HW_boundary_weight 
                    (itype_mer+Ncomp*jtype_mer,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
        }
        if (fill_flag != MSR_PREPROCESS){
           /* first load the Boltzman factor derivatives */
           fac1=weight; 
           for(i=0;i<nunk-1;i++) fac1 *=x[B2L_unknowns[loc_find(unk[i],jnode_box,BOX)]];  /*Gs or Qs*/

           resid_entry = fac1*boltz_prefac_R*POW_DOUBLE_INT(x[B2L_unknowns[loc_find(unk[nunk-1],jnode_box,BOX)]],boltz_pow_R); /* Boltz Term */
           resid_sum += resid_entry;
           jac_entry = fac1*boltz_prefac_1*POW_DOUBLE_INT(x[B2L_unknowns[loc_find(unk[nunk-1],jnode_box,BOX)]],boltz_pow_1); /* Boltz Term */
           mat_row[B2L_unknowns[loc_find(unk[nunk-1],jnode_box,BOX)]] += jac_entry;
           
          
           /* now load the G/Q derivatives */

           for (i=0;i < nunk-1; i++){
              loc_i2_j2 = B2L_unknowns[loc_find(unk[i],jnode_box,BOX)];
              fac2=weight;
              for (j=0; j<nunk-1; j++){
                 if (j != i)  fac2 *= x[B2L_unknowns[loc_find(unk[j],jnode_box,BOX)]];  /*Gs or Qs*/
              }
              jac_entry = fac2*boltz_prefac_2*POW_DOUBLE_INT(x[B2L_unknowns[loc_find(unk[nunk-1],jnode_box,BOX)]],boltz_pow_2); /*Boltz Term*/
              mat_row[loc_i2_j2] += jac_entry;
           }
        }
        else{
           bindx_tmp[loc_find(unk[nunk-1],jnode_box,BOX)]=TRUE;
           for (i=0;i < nunk-1; i++){
              box_i2_j2 = loc_find(unk[i],jnode_box,BOX);
              bindx_tmp[box_i2_j2]=TRUE;
           }
        }
     }
  }

  return(resid_sum);
}
/*****************************************************************************/
/* load_polymer_cr: 
             In this routine we load the residual and Jacobian for 
             the polymer fields.                      */
double load_polymer_cr(int sten_type, int loc_i, int itype_mer, 
                int izone, int *ijk_box, int fill_flag,
                double *resid, double *mat_row, int *bindx_tmp, double *x)
{
  int   **sten_offset, *offset, isten;
  double *sten_weight,  weight, weight_bulk;
  struct Stencil_Struct *sten;

  double sign,r_tmp,row_fact;
  int jtype_mer, jlist;
  int reflect_flag[NDIM_MAX];
  int j_box, jnode_box,node_start;

  sign = 1.0;
  r_tmp = 0.;
  for (jtype_mer=0; jtype_mer<Ncomp; jtype_mer++){
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

         if (jnode_box >= 0 ) {
            if (Lhard_surf && fill_flag != MSR_PREPROCESS) {
                if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1) 
                   weight = HW_boundary_weight 
                    (jtype_mer,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
                else weight = weight_bulk;
            }
            else weight = weight_bulk;

             j_box=loc_find(Unk_start_eq[DENSITY]+jtype_mer,jnode_box,BOX);  
                                            /* density of this component type */
             if (fill_flag != MSR_PREPROCESS){
                    r_tmp -=  sign*(weight*x[B2L_unknowns[j_box]]-
                                           weight_bulk*Rho_b[jtype_mer]);
                    mat_row[B2L_unknowns[j_box]] -= sign*weight;
             }
             else /* MSR_PREPROCESS */
                 bindx_tmp[j_box] = TRUE;
         }
         else if ( jnode_box == -2 && fill_flag != MSR_PREPROCESS){  /*in wall*/
              r_tmp +=  sign*weight_bulk*Rho_b[jtype_mer];

         }
         /* jnode_box == -1 = in bulk contributes nothing  */
      }
  }

  if (fill_flag != MSR_PREPROCESS){
     if ((Sten_Type[POLYMER_CR] == 1) || (r_tmp > 0.5)) {
       resid[loc_i] += r_tmp;
       row_fact = 1.;
     }
     else{
       row_fact = sqrt(1.-2.*r_tmp);
       resid[loc_i] += 1. - row_fact;
       resid[loc_i] *= row_fact;
     }
  }
  else row_fact = 1.;

  return row_fact;
}
/*****************************************************************************/
