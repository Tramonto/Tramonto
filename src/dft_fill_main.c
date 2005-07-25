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
void fill_resid_and_matrix (double **x, int iter, int resid_only_flag,int unk_flag)
{
 /*
  * Local variable declarations
  */

  char   *yo = "fill_resid_and_matrix";
  int     i, j, icomp,idim,iunk,junk,iunk_start,iunk_end,iseg,ibond;
  int     zero_density_bond_check,unk_bond;
  double  n;
  int     reflect_flag[3];
  int     izone, mesh_coarsen_flag_i;
  int     loc_i, loc_inode;
  struct  RB_Struct *dphi_drb=NULL, dphi_drb_bulk,
                     dphi_drb_bulk_left, dphi_drb_bulk_right;

   double resid,mat_value,values[3];
   int numEntries,nodeIndices[3],iloop;
   double resid_hs1,resid_hs2,resid_rhobars,resid_rhobarv,resid_uatt;
   double resid_ig,resid_vext,resid_mu,resid_charge,resid_deltac;
   double resid_poisson,resid_transport,resid_el,resid_cavity,resid_bondwtc,resid_WTC1;

  double fac_temp,gradphi,fac;

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

  if (Ipot_ff_n != IDEAL_GAS){
     dphi_drb = (struct RB_Struct *) array_alloc
                    (1, Nnodes_box, sizeof(struct RB_Struct));
     if (Type_func==0)
        pre_calc_dphi_drb_rb1(dphi_drb, x, &dphi_drb_bulk, 
                    &dphi_drb_bulk_left,
                    &dphi_drb_bulk_right);
     else
        pre_calc_dphi_drb_rb2(dphi_drb, x, &dphi_drb_bulk, 
                    &dphi_drb_bulk_left,
                    &dphi_drb_bulk_right);
  }

  /* for debugging print out profiles on each iteration */
  if (Iwrite==VERBOSE) print_profile_box(x, "dens_iter.dat");

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
       resid_ig = resid_vext = resid_hs1 = resid_hs2 = resid_uatt = resid_mu = resid_charge
                = resid_poisson = resid_deltac = resid_transport = resid_rhobars = resid_rhobarv 
                = resid_cavity= resid_bondwtc=resid_WTC1=0.0;

      if (mesh_coarsen_flag_i == FLAG_1DBC){
         node_box_to_ijk_box(inode_box,ijk_box);
         for (idim=0; idim<Ndim; idim++) ijk_tmp[idim]=0;
         ijk_tmp[Grad_dim]=ijk_box[Grad_dim];
         jnode_box = ijk_box_to_node_box(ijk_tmp);
         loc_jnode = B2L_node[jnode_box];

         if (jnode_box <0 ){
             printf("Proc: %d ijk_box: %d %d %d PROBLEMS: jnode_box: %d  ijk_tmp: %d %d %d\n",
                Proc,ijk_box[0],ijk_box[1],ijk_box[2],jnode_box,ijk_tmp[0],ijk_tmp[1],ijk_tmp[2]);
                exit(-1);
         }

         resid = x[iunk][inode_box]-x[iunk][jnode_box];
         dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
         numEntries=2;
         values[0]=1.0; values[1]=-1.0;
         nodeIndices[0]=inode_box; nodeIndices[1]=jnode_box;
         dft_linprobmgr_insertmultinodematrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                               iunk, nodeIndices,values,numEntries);
      } 

      /* do mesh coarsening if indicated .... for all unknowns ! */
      else if (mesh_coarsen_flag_i < 0 && mesh_coarsen_flag_i != FLAG_BULK) {

         resid= x[iunk][inode_box];
         mat_value=1.0;
         dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_value);         

         for (iloop=0;iloop<2;iloop++){

          if (iloop==0) offset_ptr = &offset_idim_pm[3*(-mesh_coarsen_flag_i - 1)];  /* one node higher */
          else          offset_ptr += 9;                                             /* one node lower */

          jnode_box = offset_to_node_box(ijk_box, offset_ptr, reflect_flag);

          if (jnode_box >= 0) {
             resid-= 0.5*x[iunk][jnode_box];
             mat_value=-0.5;
             dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,jnode_box,mat_value);         
          }
          else{  resid-= 0.5*constant_boundary(iunk,jnode_box); }
        }
        dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
     }
     else{

      /* SET ZONE */
      /* izone is set based on value of mesh_coarsen_flag       */
      /* izone is the zone number at the node of interest which */
      /* indicated which quadrature scheme is to be used.       */

      izone = mesh_coarsen_flag_i;

      /**** LOAD EULER-LAGRANGE EQUATIONS *****/
      if (Unk2Phys[iunk]==DENSITY) {      
        i = iunk-Phys2Unk_first[DENSITY];
        if (Type_poly_TC){
                          iseg=i;
                          icomp=Unk2Comp[iseg];
        }
        else              icomp=i;

        /* Do trivial fill if it is a zero-density node */
        zero_density_bond_check=FALSE;
        if (Type_poly_TC){
            for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
                unk_bond = Phys2Unk_first[BOND_WTC]+Poly_to_Unk_SegAll[iseg][ibond];
                n=x[unk_bond][inode_box];
                if (fabs(n)<1.e-8) zero_density_bond_check=TRUE;
            }
        }
        if (Zero_density_TF[inode_box][icomp] || Vext[loc_inode][icomp] == VEXT_MAX ||zero_density_bond_check) {
             resid= x[iunk][inode_box];
             mat_value = 1.0;
             dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
             dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_value);         
        }
        else {

             /* First load diagonal and constant terms to the residual */
             resid = log(x[iunk][inode_box]) ; 
             mat_value = 1.0/x[iunk][inode_box];
             resid_ig=resid;
             dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
             dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_value);         

             if (mesh_coarsen_flag_i == FLAG_BULK){
                  if (Type_poly_TC)   resid = -log(Rho_seg_b[iseg]);
                  else                resid = -log(Rho_b[icomp]);
                  dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
             }
             else {
/*             if(Lsteady_state){
                  resid = (- 3.0*log(Sigma_ff[icomp][icomp]) -1.5*log(Mass[icomp]*Temp));
                  dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
               }
*/
               if (Nwall > 0){
                 resid = Vext[loc_inode][icomp];
                 dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                 resid_vext=resid;
               }
               
               resid_mu = 0.0; 
               if (Lsteady_state == FALSE) { 
                  if (Iliq_vap < 10 ){
                     if (Type_poly_TC)  resid_mu -= log(Rho_seg_b[iseg]);
                     else               resid_mu -= log(Rho_b[icomp]);
                  
                   if (Ipot_ff_n != IDEAL_GAS) resid_mu -= Betamu_hs_ex[icomp];
                   if (Ipot_ff_n == LJ12_6)    resid_mu -= Betamu_att[icomp];
                   if (Type_poly_TC)           resid_mu -= Betamu_wtc[iseg];
                   if (Ipot_ff_c == COULOMB && Sten_Type[THETA_CHARGE]) resid_mu += Deltac_b[icomp];
                  }
                  else{
                      if (Type_poly_TC) resid_mu = -Betamu_seg[iseg];
                      else              resid_mu = -Betamu[icomp];
                  }
               }
               else{                             /* Lsteady_state == TRUE */
                  if(Type_poly_TC)  junk=Phys2Unk_first[DIFFUSION] + iseg;
                  else              junk=Phys2Unk_first[DIFFUSION] + icomp;
                  resid_mu = -x[junk][inode_box];
                  mat_value=-1.0;
                  dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,inode_box,mat_value);         
               }
               dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid_mu);
             }

             if (Ipot_ff_c == COULOMB){
                junk = Phys2Unk_first[POISSON];
                resid_charge = Charge_f[icomp]*x[junk][inode_box];
                dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid_charge);
                mat_value = Charge_f[icomp];
                dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,inode_box,mat_value);         

                if (Lpolarize[icomp] && Ndim==1){

                 fac_temp=Temp_elec/(4.0*PI*KAPPA_H2O);
                 if (Nodes_2_boundary_wall[0][inode_box] != -1){
                    if (Surf_normal[0][loc_inode][0] < 0){
                        gradphi = 0.5*(3.0*x[junk][inode_box]-4.0*x[junk][inode_box-1]+x[junk][inode_box-2])/Esize_x[0];
                        fac=Pol[icomp]*gradphi*fac_temp/Esize_x[0];

                        numEntries=3;
                        nodeIndices[0]=inode_box; nodeIndices[1]=inode_box-1; nodeIndices[2]=inode_box-2;
                        values[0]=1.5*fac; values[1]=-2.0*fac; values[2]=0.5*fac;
                    }
                    else{
                        gradphi = 0.5*(-3.0*x[junk][inode_box]+4.0*x[junk][inode_box+1]-x[junk][inode_box+2])/Esize_x[0];
                        fac=Pol[icomp]*gradphi*fac_temp/Esize_x[0];

                        numEntries=3;
                        nodeIndices[0]=inode_box; nodeIndices[1]=inode_box+1; nodeIndices[2]=inode_box+2;
                        values[0]=-1.5*fac; values[1]=2.0*fac; values[2]=-0.5*fac;
                    }
                  }
                  else{
                      gradphi = 0.5*(x[junk][inode_box+1]-x[junk][inode_box-1])/Esize_x[0];
                      fac=Pol[icomp]*gradphi*fac_temp/Esize_x[0];
                      numEntries=2;
                      nodeIndices[0]=inode_box+1; nodeIndices[1]=inode_box-1; 
                      values[0]=0.5*fac; values[1]=-0.5*fac; 
                  }
                  resid = 0.5*Pol[icomp]*gradphi*gradphi*fac_temp;
                  resid_charge += resid;
                  dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                  dft_linprobmgr_insertmultinodematrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                       junk,nodeIndices,values,numEntries);
                } 
             }
        
        /* Now loop over all the stencils to load hard sphere and mean field 
           (off diagonal)  contributions to the Euler-Lagrange equations in the Jacobian. */

           if (Ipot_ff_n != IDEAL_GAS && mesh_coarsen_flag_i != FLAG_BULK) {

              resid_hs1 =load_nonlocal_hs_rosen_rb(DELTA_FN,iunk,loc_inode,inode_box,
                            icomp,izone,
                            ijk_box,x,dphi_drb,
                            &dphi_drb_bulk, &dphi_drb_bulk_left,
                            &dphi_drb_bulk_right,
                            resid_only_flag);

              resid_hs2=load_nonlocal_hs_rosen_rb(THETA_FN,iunk,loc_inode,inode_box,
                            icomp,izone,
                            ijk_box,x,dphi_drb,
                            &dphi_drb_bulk, &dphi_drb_bulk_left,
                            &dphi_drb_bulk_right,
                            resid_only_flag);


              if (Sten_Type[U_ATTRACT]) {  /* load attractions */
                 resid_uatt=load_mean_field(U_ATTRACT,iunk,loc_inode,
                                 icomp,izone,ijk_box, x, resid_only_flag);
              }   
            
              if (Sten_Type[THETA_CHARGE]) {   /* load electrostatics deltac correlations */
                 resid_deltac=load_mean_field(THETA_CHARGE,iunk,loc_inode,
                                 icomp,izone,ijk_box,x, resid_only_flag);
              }   

              if (Type_poly_TC){
                   resid_WTC1+=load_polyTC_diagEL(iunk,loc_inode,inode_box,icomp,
                                                 izone,ijk_box,x,resid_only_flag);
                   resid_WTC1+=load_polyTC_bondEL(iunk,loc_inode,inode_box,icomp,
                                                 izone,ijk_box,x,resid_only_flag);
                   resid_WTC1+=load_polyTC_cavityEL(iunk,loc_inode,inode_box,icomp,
                                                 izone,ijk_box,x,resid_only_flag);
              }

           }
        }
      }
      /***** END EULER-LAGRANGE *****/

      /*********************************/
      /**** LOAD RHO_BAR EQUATIONS *****/
      /*********************************/
      else if (Unk2Phys[iunk]==RHOBAR_ROSEN){

          if (iunk == Phys2Unk_first[RHOBAR_ROSEN])
             resid_rhobars+=load_rho_bar_s(THETA_FN,x,iunk,loc_inode,inode_box,izone,ijk_box, 
                            resid_only_flag);
          else if (iunk < Phys2Unk_first[RHOBAR_ROSEN]+Nrho_bar_s)
             resid_rhobars+=load_rho_bar_s(DELTA_FN,x,iunk,loc_inode,inode_box,izone,ijk_box, 
                            resid_only_flag);

          if (iunk >= Phys2Unk_first[RHOBAR_ROSEN]+Nrho_bar_s){
              resid_rhobarv+=load_rho_bar_v(x,iunk,loc_inode,inode_box,izone,ijk_box, resid_only_flag);
          }
      }
      /****** END RHOBAR *****/

      /***************************************************/
      /**** LOAD POISSON'S EQUATION FOR ELECTROSTATICS ***/
      /***************************************************/
      else if (Unk2Phys[iunk]==POISSON){

            l_elec_RTF=FALSE;
            l_elec_LBB=FALSE;
            if (Lsteady_state){
               node_to_position(node_box_to_node(inode_box),nodepos);
               idim = Grad_dim;
               if ( nodepos[idim] +0.5*Size_x[idim] - X_const_mu <= 0.00000001) l_elec_LBB=TRUE;
               else if (nodepos[idim] - 0.5*Size_x[idim] + X_const_mu >=-0.00000001) l_elec_RTF=TRUE;
            }

            if (l_elec_LBB){
               resid = x[iunk][inode_box]-Elec_pot_LBB;
               mat_value = 1.0;
               dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
               dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_value);         
               resid_poisson = resid;
            }
            else if (l_elec_RTF){
               resid = x[iunk][inode_box]-Elec_pot_RTF;
               mat_value = 1.0;
               dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
               dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_value);         
               resid_poisson = resid;
            }
            else if (!l_elec_LBB && !l_elec_RTF) {
               if(Type_coul==POLARIZE){
                  resid_poisson=load_polarize_poissons_eqn(iunk,loc_inode, inode_box,ijk_box,x);
               }
               else{
                  resid_poisson=load_poissons_eqn(iunk,loc_inode,inode_box,ijk_box,x);
               }
               resid_poisson+=load_poisson_bc(iunk,loc_inode,inode_box);
            }
      }
      /******** END POISSON ********/

      /*******************************/
      /*** LOAD TRANSPORT EQUATION ***/
      /*******************************/
      else if (Unk2Phys[iunk]==DIFFUSION){
 
            if (Linear_transport)
                resid_transport=load_linear_transport_eqn(iunk,loc_inode,inode_box,ijk_box,x);
            else   
                resid_transport=load_nonlinear_transport_eqn(iunk,loc_inode,inode_box,ijk_box,x);
      }
      /******** END TRANSPORT ********/
      else if (Unk2Phys[iunk]==CAVITY_WTC){
      /***********************************************/
      /*** LOAD CAVITY EQUATIONS OF WTC FUNCTIONALS***/
      /***********************************************/
         resid=x[iunk][loc_inode];
         mat_value=1.0;
         dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_value);
         resid += load_cavity_wtc(iunk,loc_inode,inode_box,izone,ijk_box,x);
         resid_cavity=resid;
         dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
      }
      else if (Unk2Phys[iunk]==BOND_WTC){

      /********************************************/
      /*** LOAD BOND EQUATIONS OF WTC FUNCTIONS ***/
      /********************************************/
         resid=x[iunk][loc_inode];
         mat_value=1.0;
         dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_value);
         resid += load_bond_wtc(iunk,loc_inode,inode_box,izone,ijk_box,x);
         resid_bondwtc=resid;
         dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
      }

      else {
         printf("Problem with unknowns in fill !!!\n");
         exit (-1);
      }
      }
    

      /* PRINT STATEMENTS FOR PHYSICS DEBUGGING .... CHECK RESIDUALS INDEPENDENTLY  */
/*    if (Unk2Phys[iunk]==DENSITY){
       resid_el = resid_ig + resid_vext + resid_mu + resid_charge;
          printf("loc_inode=%d  iunk=%d  resid_el=%9.6f (%9.6f %9.6f %9.6f %9.6f) resid_hs1=%9.6f resid_hs2=%9.6f  resid_WTC1=%9.6f",
                                 loc_inode,iunk,resid_el,resid_ig,resid_vext,resid_mu,resid_charge,resid_hs1,resid_hs2,resid_WTC1);
    }
    else if (Unk2Phys[iunk]==RHOBAR_ROSEN){
       printf("loc_inode=%d : iunk_rbar=%d resid_rhobars=%9.6f  resid_rhobarv=%9.6f ",
                 loc_inode,iunk,resid_rhobars,resid_rhobarv);
    }
    else if (Unk2Phys[iunk]==POISSON)   
         printf(" loc_inode=%d  iunk_poisson=%d   resid=%9.6f", loc_inode,iunk,resid_poisson);
    else if (Unk2Phys[iunk]==DIFFUSION) 
         printf(" loc_inode=%d  iunk_diffusion=%d  resid=%9.6f",loc_inode,iunk,resid_transport);
    else if (Unk2Phys[iunk]==CAVITY_WTC) 
         printf(" loc_inode=%d  iunk_cavity=%d  resid=%9.6f",loc_inode,iunk,resid_cavity);
    else if (Unk2Phys[iunk]==BOND_WTC) 
         printf(" loc_inode=%d  iunk_bondwtc=%d  resid=%9.6f",loc_inode,iunk,resid_bondwtc);
    printf("  \n");*/

    } /* end of loop over # of unknowns per node */

  } /* end of loop over local nodes */



/*** resid is no longer a vector so this won't work.
 * Perhaps first do a getRhs and then a ImportR2C to get assembled
 * resid vector in box coordinates
  if (Iwrite==VERBOSE){ print_profile_box(resid,"dens_iter.dat"); }
 ***/


  /* Now add in constraint equation(s)  */
/*  if (Proc==0 && Lstoichiometry) {
     for (icomp = 0; icomp<Ncomp; icomp++){
         load_stoichimetric_constraint(i_box, inode_box, loc_i, ijk_box, mat_row,
                                           resid, x, bindx_tmp, fill_flag, icomp); 
     }
  }*/

/*  if (Ipot_wf_c) {
     load_poisson_bc(resid);
  }*/

  if (Ipot_ff_n != IDEAL_GAS) {
     safe_free((void *) &dphi_drb);
  }
}
/****************************************************************************/
