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

/*
 *  FILE: dft_picard.c
 *
 *  This file solves the dft problem using a Picard successive substitution method.
 */
#include "dft_picard.h"

/*******************************************************************************/
int solve_problem_picard(double **x, double **x2)
/*
 * x and x2 in dft_main are allocated [Nunk_per_node][Nnodes_box]
 * x2 only relevent for Lbinodal calculations, so can mostly be ignored
 */
{
  int iter,iunk,i;
  int loc_inode;
  double **xOwned, **x2Owned;

  /* this routine looks just like the initial guess routine where we used the density field
     to generate all other fields.  In this case, we use the fields at the ith iteration to
     generate the solution at the i+1 st iteration */

  /* Set initial guess on owned nodes and reconcile ghost nodes manually here -- note this is a 
      rather silly approach since we have local nodes being converted to box coordinates both
      here and in the initial guess routine - I'll try to take this out, but don't want to make
      changes to the Newton's method routine while implementing something here.  */

  xOwned = (double **) array_alloc(2, Nunk_per_node, Nnodes_per_proc, sizeof(double));
  for (iunk=0;iunk<Nunk_per_node;iunk++)
        for (loc_inode=0;loc_inode<Nnodes_per_proc;loc_inode++) xOwned[iunk][loc_inode]=0.0;
  set_initial_guess(Iguess, xOwned);
  (void) dft_linprobmgr_importr2c(LinProbMgr_manager, xOwned, x);

  /* If requested, write out initial guess */
   if (Iwrite_files == FILES_DEBUG) print_profile_box(x,"rho_init.dat");

  /* Do same for second solution vector when Lbinodal is true */
  if (Lbinodal) {
    x2Owned = (double **) array_alloc(2, Nunk_per_node, Nnodes_per_proc, sizeof(double));
    for (iunk=0;iunk<Nunk_per_node;iunk++)
        for (loc_inode=0;loc_inode<Nnodes_per_proc;loc_inode++) x2Owned[iunk][loc_inode]=0.0;
    set_initial_guess(BINODAL_FLAG, x2Owned);
    (void) dft_linprobmgr_importr2c(LinProbMgr_manager, x2Owned, x2);

    if (Iwrite_files == FILES_DEBUG) print_profile_box(x2,"rho_init2.dat");
  }

  if (NL_Solver==PICARD_NOX || NL_Solver==PICNEWTON_NOX) iter=NOXLOCA_Solver(x, xOwned, x2Owned, TRUE);
  else{
      iter=picard_solver(x,xOwned,-1); 
      if(Lbinodal) iter=picard_solver(x2,x2Owned,-1); 
  }

  safe_free((void **) &xOwned);
  if (Lbinodal)  safe_free((void **) &x2Owned);

  return(iter);
}
/****************************************************************************/
int picard_solver(double **x, double **xOwned, int subIters){
  /* subIters  of -1 means Picard is being called as the nonlinear solver alg.
     Otherwise, this is an inner iteration to another solver (NOX) and
     will just proceed for subIters iterations and return. */
  int iter=0,i, max_iters;
  int iunk, ibox;
    char tempfilename[FILENAME_LENGTH]; //LMH
  int converged=FALSE;
  int skip_convergence_test=FALSE,Lprint_screen;
  double **x_old, **delta_x;

  if (Type_poly==WJDC3 && Grafted_Logical){
     x_old = (double **) array_alloc(2, Nunk_per_node, Nnodes_box_extra, sizeof(double));
     delta_x = (double **) array_alloc(2, Nunk_per_node, Nnodes_box_extra, sizeof(double));
  }
  else{
     x_old = (double **) array_alloc(2, Nunk_per_node, Nnodes_box, sizeof(double));
     delta_x = (double **) array_alloc(2, Nunk_per_node, Nnodes_box, sizeof(double));
  }

  if (subIters == -1) max_iters = Max_NL_iter;
  else {
    max_iters = subIters;
    skip_convergence_test = TRUE;
  }

  do {
    iter++;

    Lprint_screen=FALSE;
    if (max_iters<30){
       if ((Iwrite_screen== SCREEN_BASIC && (iter%(max_iters/1)==0 || iter==1)) || Iwrite_screen==SCREEN_VERBOSE){
          print_resid_norm_picard(x,iter);
          Lprint_screen=TRUE;
          sprintf(tempfilename,"dft_dens%d.vtk",iter); //LMH
          print_profile_box_vtk(x,tempfilename); //LMH print density for visualizing
          sprintf(tempfilename,"dft_dens%i.0",iter); //LMH
          print_profile_box(x,tempfilename); //LMH also save for restarting
       }
    }
    else if (max_iters<5000){
       if ((Iwrite_screen== SCREEN_BASIC && (iter%(max_iters/20)==0 || iter==1)) || Iwrite_screen==SCREEN_VERBOSE){
          print_resid_norm_picard(x,iter);
          Lprint_screen=TRUE;
           sprintf(tempfilename,"dft_dens%i.vtk",iter); //LMH
           print_profile_box_vtk(x,tempfilename); //LMH print density for visualizing
           sprintf(tempfilename,"dft_dens%i.0",iter); //LMH
           print_profile_box(x,tempfilename); //LMH also save for restarting
       }
    }
    else{
       if ((Iwrite_screen== SCREEN_BASIC && (iter%(max_iters/50)==0 || iter==1)) || Iwrite_screen==SCREEN_VERBOSE){
          print_resid_norm_picard(x,iter);
          Lprint_screen=TRUE;
           sprintf(tempfilename,"dft_dens%i.vtk",iter); //LMH
           print_profile_box_vtk(x,tempfilename); //LMH print density for visualizing
           sprintf(tempfilename,"dft_dens%i.0",iter); //LMH
           print_profile_box(x,tempfilename); //LMH also save for restarting
       }
    }

     /* copy current fields to the x_old array */
     for (iunk=0; iunk<Nunk_per_node;iunk++){
        if (Type_poly==WJDC3 && Grafted_Logical){ for (ibox=0; ibox<Nnodes_box_extra;ibox++) x_old[iunk][ibox]=x[iunk][ibox];}
        else{ for (ibox=0; ibox<Nnodes_box;ibox++) x_old[iunk][ibox]=x[iunk][ibox];}
     }
	  
	  /* for grafted chains */
     if(Type_poly==CMS || Type_poly==WJDC3) calc_Gsum_new(x);

     /* use successive substitution to update density field, then compute all other fields */ 
     if ((L_HSperturbation || Type_coul != NONE) && Type_poly != WJDC && Type_poly !=WJDC2 && Type_poly!=WJDC3) 
                                                                                    calc_density_next_iter_HSperturb(x,xOwned);
     else if(Type_poly==CMS)                                                        calc_density_next_iter_CMS(x,xOwned);
     else if(Type_poly==CMS_SCFT)                                                   calc_density_next_iter_SCF(x,xOwned);
     else if(Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3)               calc_density_next_iter_WJDC(x,xOwned);

     for (i=0; i<Phys2Nunk[DENSITY]; i++){
        iunk=Phys2Unk_first[DENSITY]+i;
        (void) dft_linprobmgr_importsingleunknownr2c(LinProbMgr_manager, xOwned[iunk], x[iunk]);
     }

     for (iunk=0; iunk<Nunk_per_node;iunk++){
         if (Type_poly==WJDC3 && Grafted_Logical){ for (ibox=0; ibox<Nnodes_box_extra;ibox++) delta_x[iunk][ibox]=x[iunk][ibox]-x_old[iunk][ibox]; }
         else{ for (ibox=0; ibox<Nnodes_box;ibox++) delta_x[iunk][ibox]=x[iunk][ibox]-x_old[iunk][ibox]; }
     }


   /* if (con_ptr != NULL) converged2 = continuation_hook_conwrap(x, delta_x, con_ptr, NL_rel_tol, NL_abs_tol); */

    /* Do: x += delta_x, and check for convergence .... */
    converged = update_solution_picard(x_old, xOwned, delta_x, iter,Lprint_screen);

    if (converged==TRUE){ 
        if (Iwrite_screen != SCREEN_NONE && Iwrite_screen != SCREEN_ERRORS_ONLY) print_resid_norm_picard(x,iter);
        Lprint_screen=TRUE;
        converged = update_solution_picard(x_old, xOwned, delta_x, iter,Lprint_screen);
    }

    (void) dft_linprobmgr_importr2c(LinProbMgr_manager, xOwned, x_old);
    if (skip_convergence_test) converged=FALSE;

     for (iunk=0; iunk<Nunk_per_node;iunk++){
       if (Type_poly==WJDC3 && Grafted_Logical){ for (ibox=0; ibox<Nnodes_box_extra;ibox++) x[iunk][ibox]=x_old[iunk][ibox];}
       else { for (ibox=0; ibox<Nnodes_box;ibox++) x[iunk][ibox]=x_old[iunk][ibox];}
     }
     fix_symmetries(x);
 
     /* now update all other fields in the solution vector */
     /* update all other fields (other than density - based on new density profile) */
     for (i=0;i<NEQ_TYPE;i++){
     switch(i){
         case DENSITY: /* don't do anything for density here */
           break;

         case MF_EQ:
           if (Phys2Nunk[MF_EQ]>0){
              calc_init_mf_attract(x,xOwned); 
           }
           break;

         case HSRHOBAR:
           if (Phys2Nunk[HSRHOBAR]>0){
              calc_init_rho_bar(x,xOwned); 
           }
           break;

         case POISSON:
           if (Phys2Nunk[POISSON]>0){
              calc_init_elec_pot(x,xOwned);
           }
           break;

         case DIFFUSION:
           if (Phys2Nunk[DIFFUSION]>0){
              calc_init_chem_pot(x,xOwned);
           }
           break;

         case CAVWTC:
           if (Phys2Nunk[CAVWTC]>0){
              calc_init_Xi_cavWTC(x,xOwned);
            }
            break;

         case BONDWTC:
            if (Phys2Nunk[BONDWTC]>0){
               calc_init_BondWTC(x,xOwned);
             }
             break;

         case WJDC_FIELD:
            if (Phys2Nunk[WJDC_FIELD]>0){
               calc_init_WJDC_field(x,xOwned);
             }
             break;

         case CMS_FIELD:
            if (Phys2Nunk[CMS_FIELD]>0){
              calc_init_CMSfield(x,xOwned);
            }
            break;
         case SCF_FIELD:
	    if (Phys2Nunk[SCF_FIELD]>0) calc_init_SCFfield(x,xOwned);
            break;
	 case SCF_CONSTR:
             if (Phys2Nunk[SCF_CONSTR]>0) calc_init_lambda(x,xOwned);
              break;
         case G_CHAIN:
            if (Phys2Nunk[G_CHAIN]>0){
               if (Type_poly==CMS || Type_poly==CMS_SCFT) calc_init_polymer_G_CMS(x,xOwned);
               else if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3) calc_init_polymer_G_wjdc(x,xOwned);
            }
            break;
         default:
           if (Iwrite_screen != SCREEN_NONE) printf("problem with switch in initial guess\n");
           exit(-1);
           break;
     }
     }

  (void) dft_linprobmgr_importr2c(LinProbMgr_manager, xOwned, x);

  } while (iter < max_iters && !converged);

    (void) dft_linprobmgr_importr2c(LinProbMgr_manager, xOwned, x_old);
    if (skip_convergence_test) converged=FALSE;

  /* Skip printing if NOX is controlling convergence */
  if (!converged && !skip_convergence_test) {
    if (Proc==0 && Iwrite_screen !=SCREEN_NONE) printf("\tPicard Solver: Failed to converge in %d iterations\n",iter);
    iter = -iter;
  }
  else if (converged && !skip_convergence_test){
    if (Proc==0 && Iwrite_screen !=SCREEN_NONE && Iwrite_screen != SCREEN_ERRORS_ONLY) printf("\tPicard Solver: Successful convergence in %d iterations\n",iter);
  }
  else{ if (Proc==0 && Iwrite_screen==SCREEN_VERBOSE) printf("\treturn control to NOX after %d iterations\n",iter);}

  safe_free((void **) &x_old);
  safe_free((void **) &delta_x);

  return iter;
}
/****************************************************************************/
/*calc_density_next_iter_HSperturb(): compute a new density profile for cases where
  we are doing HSperturbation - but not WJDC - calculations */
void calc_density_next_iter_HSperturb(double **xInBox, double **xOwned)
{
  int loc_inode,inode_box,ijk_box[3],iloop,iunk,izone,mesh_coarsen_flag_i,nloop;
  double resid_EL, **x_InBoxOld;
  struct  RB_Struct *dphi_drb=NULL;
  izone=0;
  mesh_coarsen_flag_i=0;

  
  if (Type_poly==WJDC3 && Grafted_Logical) x_InBoxOld = (double **) array_alloc(2, Nunk_per_node, Nnodes_box_extra, sizeof(double));
  else x_InBoxOld = (double **) array_alloc(2, Nunk_per_node, Nnodes_box, sizeof(double));
  
  if (Type_func !=NONE){
     dphi_drb = (struct RB_Struct *) array_alloc (1, Nnodes_box, sizeof(struct RB_Struct));
     FMT1stDeriv_switch(xInBox,dphi_drb);
  }

  if (Lseg_densities) nloop=Nseg_tot;
  else                nloop=Ncomp;

  for (inode_box=0; inode_box<Nnodes_box; inode_box++) {
     for (iloop=0; iloop<Nunk_per_node; iloop++) x_InBoxOld[iloop][inode_box]=xInBox[iloop][inode_box];
  }
 
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box=L2B_node[loc_inode]; 
     node_box_to_ijk_box(inode_box, ijk_box);
     
     for (iloop=0; iloop<nloop; iloop++){
        iunk=Phys2Unk_first[DENSITY]+iloop;
        resid_EL=load_euler_lagrange(iunk,loc_inode,inode_box,ijk_box,izone,
                          x_InBoxOld,dphi_drb,mesh_coarsen_flag_i,INIT_GUESS_FLAG);
        xInBox[iunk][inode_box]=resid_EL;
        xOwned[iunk][loc_inode]=xInBox[iunk][inode_box];
     }
  }
  if (Type_func != NONE) safe_free((void *) &dphi_drb);
  if (Type_func != NONE) safe_free((void *) &x_InBoxOld);
  return;
}
/****************************************************************************/
/* calc_density_next_iter_CMS(x); compute a new density profile for cases where
   we are doing CMS-DFT calculations */
void calc_density_next_iter_CMS(double **xInBox,double **xOwned)
{
  int loc_inode,inode_box,ijk_box[3],iloop,iunk,izone,mesh_coarsen_flag_i,i;
  double resid;
  struct  RB_Struct *dphi_drb=NULL;
  izone=0;
  mesh_coarsen_flag_i=0;

  calc_init_CMSfield(xInBox,xOwned);

  for (i=0; i<Phys2Nunk[CMS_FIELD]; i++){
     iunk=Phys2Unk_first[CMS_FIELD]+i;
     (void) dft_linprobmgr_importsingleunknownr2c(LinProbMgr_manager, xOwned[iunk], xInBox[iunk]);
  }

  calc_init_polymer_G_CMS(xInBox,xOwned);

  for (i=0; i<Phys2Nunk[G_CHAIN]; i++){
     iunk=Phys2Unk_first[G_CHAIN]+i;
     (void) dft_linprobmgr_importsingleunknownr2c(LinProbMgr_manager, xOwned[iunk], xInBox[iunk]);
  }

/*  (void) dft_linprobmgr_importr2c(LinProbMgr_manager, xOwned, xInBox);*/
  
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box=L2B_node[loc_inode]; 
     node_box_to_ijk_box(inode_box, ijk_box);
     
     for (iloop=0; iloop<Ncomp; iloop++){
        iunk=Phys2Unk_first[DENSITY]+iloop;
        resid=load_CMS_density(iunk,loc_inode,inode_box,xInBox,INIT_GUESS_FLAG);
        xInBox[iunk][inode_box]=-resid;
        xOwned[iunk][loc_inode]=xInBox[iunk][inode_box];
     }
  }
  return;
}
/****************************************************************************/
/* calc_density_next_iter_WJDC(x); compute a new density profile for cases where
   we are doing WJDC-DFT calculations */
void calc_density_next_iter_WJDC(double **xInBox,double **xOwned)
{
  int loc_inode,inode_box,ijk_box[3],iloop,iunk,izone,mesh_coarsen_flag_i,iloop_max,i;
  double resid;
  struct  RB_Struct *dphi_drb=NULL;
  izone=0;
  mesh_coarsen_flag_i=0;
  if (Type_poly==WJDC3) iloop_max=Ncomp;
  else iloop_max=Nseg_tot;

  calc_init_WJDC_field(xInBox,xOwned);
  for (i=0; i<Phys2Nunk[WJDC_FIELD]; i++){
     iunk=Phys2Unk_first[WJDC_FIELD]+i;
     (void) dft_linprobmgr_importsingleunknownr2c(LinProbMgr_manager, xOwned[iunk], xInBox[iunk]);
  }

  calc_init_polymer_G_wjdc(xInBox,xOwned);
  for (i=0; i<Phys2Nunk[G_CHAIN]; i++){
     iunk=Phys2Unk_first[G_CHAIN]+i;
     (void) dft_linprobmgr_importsingleunknownr2c(LinProbMgr_manager, xOwned[iunk], xInBox[iunk]);
  }
/*  (void) dft_linprobmgr_importr2c(LinProbMgr_manager, xOwned, xInBox);*/

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box=L2B_node[loc_inode]; 
     node_box_to_ijk_box(inode_box, ijk_box);
     
     for (iloop=0; iloop<iloop_max; iloop++){
        iunk=Phys2Unk_first[DENSITY]+iloop;
        resid=load_WJDC_density(iunk,loc_inode,inode_box,xInBox,INIT_GUESS_FLAG);
        xInBox[iunk][inode_box]=-resid;
        xOwned[iunk][loc_inode]=xInBox[iunk][inode_box];
     }
  }
  return;
}
/****************************************************************************/
/* calc_density_next_iter_SCF(x); compute a new density profile for cases where
we are doing CMS-SCFT calculations */
void calc_density_next_iter_SCF(double **xInBox,double **xOwned)
{
	int loc_inode,inode_box,ijk_box[3],iloop,iunk,izone,mesh_coarsen_flag_i;
	double resid;
	struct  RB_Struct *dphi_drb=NULL;
	izone=0;
	mesh_coarsen_flag_i=0;
	
	for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
		inode_box=L2B_node[loc_inode]; 
		node_box_to_ijk_box(inode_box, ijk_box);
		
		for (iloop=0; iloop<Ncomp; iloop++){
			iunk=Phys2Unk_first[DENSITY]+iloop;
			resid=load_SCF_density(iunk,loc_inode,inode_box,xInBox,INIT_GUESS_FLAG);
			xInBox[iunk][inode_box]=-resid;
                        xOwned[iunk][loc_inode]=xInBox[iunk][inode_box];
		}
	}
	return;
}

/*****************************************************************************************************/
void print_resid_norm_picard(double **x, int iter)
{
  double sum_local,norm;

  sum_local=fill_resid_and_matrix_control(x,iter,CALC_RESID_ONLY);
  norm = gsum_double(sum_local);

  if (Proc==0) printf("\tIter=%d\t::::\tResid norm=%g",iter, sqrt(norm)); 

  return;
}
/*****************************************************************************************************/
int update_solution_picard(double** x, double **xOwned, double **delta_x, int iter,int Lprint_screen) {
/* Routine to update solution vector x using delta_x and
 * to test for convergence of the nonlinear solution (Picard) method.
 * Note that Picard updates require different heuristics than the Newton updates, and only
 * the density field is updated here.
 *
 * Note: Modifications to Newton's method, including damping, have not yet
 *       been translated form previous update_solution.
 *       iter  value may be used in some damping methods
 */
  
  int i,iunk, ibox, inode,inodeG,ijk[3],go_update,idim;
  int nloop;
  double updateNorm=0.0, temp,xold;
  char *yo = "newupdate solution";

  if (Lseg_densities) nloop=Nseg_tot;
  else                nloop=Ncomp;

  if (Proc==0 && Iwrite_screen != SCREEN_NONE && Iwrite_screen != SCREEN_ERRORS_ONLY && Lprint_screen==TRUE) {
        printf("\t::::\tUpdate percent = %g ",NL_update_scalingParam*100);
  }
  
  for (ibox=0; ibox<Nnodes_box; ibox++) {
    
    /* Increment updateNorm only for owned nodes (inode=-1 for ghosts) */
    inode = B2L_node[ibox];
    if (inode != -1) {
      for (iunk=0; iunk<Nunk_per_node; iunk++) {
        temp =(NL_update_scalingParam*delta_x[iunk][ibox])/(NL_rel_tol_picard*x[iunk][ibox] + NL_abs_tol_picard);
        updateNorm +=  temp*temp;
      }
    }
  /* For some cases, we need to be able to keep the solution values at the boundaries constant
     and set equal to the values that are read in from a file.  Do not update the solution
     vector at these points */
    inodeG=L2G_node[inode]; 
    node_to_ijk(inodeG,ijk);
    go_update=TRUE;
    for (idim=0; idim<Ndim;idim++){
       if (  (ijk[idim]==0 && Type_bc[idim][0] == LAST_NODE_RESTART) ||
          (ijk[idim]==Nodes_x[idim]-1 && Type_bc[idim][1] == LAST_NODE_RESTART)) go_update=FALSE;
    }
    
    /* Update all solution componenets */
    if (go_update){
       for (i=0; i<nloop; i++) {
          iunk=i+Phys2Unk_first[DENSITY];
          x[iunk][ibox] += NL_update_scalingParam*delta_x[iunk][ibox];
       }
       /*if (Type_coul != NONE) x[POISSON][ibox] += NL_update_scalingParam*delta_x[POISSON][ibox];*/
    }
    if (inode>=0){ for (iunk=0; iunk<Nunk_per_node; iunk++) xOwned[iunk][inode]=x[iunk][ibox];}
  }
  
  updateNorm = sqrt(gsum_double(updateNorm));
  
  if (Proc==0 && Iwrite_screen!= SCREEN_NONE && Iwrite_screen != SCREEN_ERRORS_ONLY && Lprint_screen==TRUE) {
       printf("\t::::\t Weighted norm of update vector =%g\n", updateNorm);
  }

  if (updateNorm > 1.0 || iter==1 ) return(FALSE);
  else                  return(TRUE);

}
/*****************************************************************************************************/

