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
 *      All Global Definitions occur in the include files here.
 */
#include <mpi.h>
#include "include_global/dft_globals.h"
#include "rf_allo.h"
/*#include "GUI/dft_GUI.h"*/

#include "dft_main.h"
/*****************************************************************************/

void dftmain(double * engptr)

     /*
      *  Initial main driver for 3D? Finite Element DFT code.
      *
      *    Date:                    12/11/96
      *
      *    Command Line Arguments:
      *   ---------------------------
      *
      *     The code takes one command line argument: the name of the
      *       ascii input file.  The default name for the command
      *       file is "dft_input.dat".
      */

{

  double    start_t,t_mesh;
  char     *input_file;
  char     *output_file1, *output_file2, *output_file4;
  char     *output_TF, *yo = "main";
  int       iend, match, idim, icomp, i, niters;
  double    Esize_x_final[3];
  double    **x,**x2=NULL;  /* The solution vector of densities */
  double    time_save;
  double    t_preprocess=0.0,t_solve=0.0,t_postprocess=0.0,t_total=0.0;
  double    t_pre_max,t_solve_max,t_post_max,t_total_max;
  double    t_pre_min,t_solve_min,t_post_min,t_total_min;
  double    t_linsolv_first_min,t_linsolv_av_min,t_linsolv_first_max,t_linsolv_av_max;
  double    t_manager_first_min,t_manager_av_min,t_manager_first_max,t_manager_av_max,t_NLSolve_max,t_NLSolve_min;
  double    t_fill_first_min,t_fill_av_min,t_fill_first_max,t_fill_av_max,t_InitGuess_min,t_InitGuess_max;
  double    *t_linsolv_first_array,*t_linsolv_av_array;
  double    *t_manager_first_array,*t_manager_av_array;
  double    *t_fill_first_array,*t_fill_av_array;
  int       min_nnodes_per_proc,max_nnodes_per_proc,min_nnodes_box,max_nnodes_box;
  double    min_nodesLoc_over_nodesBox,max_nodesLoc_over_nodesBox;
  FILE      *fp;
  int izone,isten,jcomp,jmax;
  struct Stencil_Struct *sten;
/*  char line[100],linecwd[100];*/
  int proper_bc;
  int argc=1;

  Time_linsolver_first=0.0;
  Time_linsolver_av=0.0;
  Time_manager_first=0.0;
  Time_manager_av=0.0;
  Time_fill_first=0.0;
  Time_fill_av=0.0;
  Time_NLSolve=0.0;

/*  gethostname(line,100);*/
/*  getcwd(linecwd,100);*/

  (void) MPI_Comm_size(MPI_COMM_WORLD, &Num_Proc);
  (void) MPI_Comm_rank(MPI_COMM_WORLD, &Proc);

  start_t = MPI_Wtime();
  t_preprocess = -MPI_Wtime();
 /*
  * Interpret command line arguments, get input file name or use default
  */

  switch (argc) {

    case 2:
      /*  input_file = argv[1];*/
       break;

    case 1:
       input_file = "dft_input.dat";
       break;

    default:
       printf("%s: ERROR in command line usuage:\n", yo);
       printf("\tExecutable can have 0 or 1 argument\n");
       exit(-1);
       break;
  }
  output_file1 = "dft_out.lis";
  output_file2 = "dft_vext.dat";
  output_file4 = "dft_vext_c.dat";
  output_TF = "dft_zeroTF.dat";

 /*
  * Read in the ascii input file from the front end 
  * This routine also sets up some flags
  * (see file dft_input.c)
  */

  read_input_file(input_file,output_file1);
/*  dft_OptikaGUI();*/

/*  printf("Ndim - SET IN GUI is %d\n",Ndim);*/
  setup_stencil_logicals();
  if (Type_attr != NONE) setup_stencil_uattr_core_properties();
  setup_nunk_per_node(output_file1);
  setup_pairPotentials(output_file1);

  if (Type_poly == CMS){
      Rism_cr = (double ***) array_alloc (3, Ncomp,Ncomp,N_NZCR_MAX,sizeof(double));
      setup_polymer_cr();
  }

 /* count_zero = 0;
  count_nonzero = 0;
  */

/*  if (Loca.method != -1) {
    Nruns = 1;
    if (Proc==0 && Iwrite!=NO_SCREEN) printf("\nWARNING: Loca library requested, disabling"
                        "Tramonto continuation loops\n");
  }*/

  for (Imain_loop=0; Imain_loop<Nruns; Imain_loop++){ /*only does mesh continuation now */
    if (Proc==0){  
         t_fill_first_array = (double *) array_alloc (1, Num_Proc,sizeof(double));
         t_fill_av_array = (double *) array_alloc (1, Num_Proc,sizeof(double));
         t_manager_first_array = (double *) array_alloc (1, Num_Proc,sizeof(double));
         t_manager_av_array = (double *) array_alloc (1, Num_Proc,sizeof(double));
         t_linsolv_first_array = (double *) array_alloc (1, Num_Proc,sizeof(double));
         t_linsolv_av_array = (double *) array_alloc (1, Num_Proc,sizeof(double));
     }

   if (Imain_loop > 0)  continuation_shift(); /* only mesh continuation is now possible outside of loca */ 

     for (idim=0; idim<Ndim; idim++){
        Esize_x_final[idim] = Esize_x[idim];
        if (Lmesh_refine==TRUE) Esize_x[idim] = 0.2;
     }
     iend=0;

     while (iend==0) {

     t_mesh = MPI_Wtime();

     if (Loca.method == -1 || Imain_loop==0) {
    /* we need the HS diameters before we compute the integration stencils. */
     if (Type_func != NONE) calc_HS_diams();

    /*
     * Stencils for non-local interactions are assembled next, and are
     * dependent on choice for integration schemes
     * (see file dft_stencil.c)
     */
     if (Imain_loop==0) calc_stencils();
    /*
     * do all the thermodynamics for the bulk fluid mixture
     */
     thermodynamics(output_file1);
     }

    /*
     * Set up the mesh based on input info, and allocate the solution vector
     * (see file dft_mesh.c)
     */


     if (Imain_loop > 0){
           free_mesh_arrays();
           boundary_free();
           if (Type_interface==DIFFUSIVE_INTERFACE && Ndim==1) safe_free((void *) &Area_IC);
           safe_free((void *) &Comm_node_proc);
           safe_free((void *) &Comm_unk_proc);
           safe_free((void *) &Comm_offset_node);
           safe_free((void *) &Comm_offset_unk);
           safe_free((void *)&Nel_hit);
           safe_free((void *)&Nel_hit2);
     }
     set_up_mesh(output_file1,output_file2);

    /*
     * Set up boundary and surface charge information -- this must come
     * after load-balancing because quantities are local on a Proc
     */
     boundary_setup(output_file1);

     /*
      * bvbw moved setup_integral after boundary_setup to eliminate seg fault
      */

     if (Imain_loop==0) setup_domain_multipliers();
     setup_integrals();

     /* check boudary conditions to see Coulomb external field can be computed
        for this run ..... we can't do it with periodic, multiply reflecting, or
        in wall boundary conditions due to need for infinite ewald sums */
     Vext_coul=NULL;
     if (Vol_charge_flag && Ndim==3){
         proper_bc = TRUE;
         for (idim=0; idim<Ndim; idim++){
            if (Type_bc[idim][0] == PERIODIC || 
                 (Type_bc[idim][0] == REFLECT && Type_bc[idim][1] != IN_BULK) ||
                 (Type_bc[idim][1] == REFLECT && Type_bc[idim][0] != IN_BULK) ||
                 (Type_bc[idim][0] == IN_WALL && Type_bc[idim][0] != IN_WALL)) 
             proper_bc=FALSE; 
         }
         if (proper_bc) setup_vext_coulomb_vol();
         else{
           if(Iwrite !=NO_SCREEN) printf("Not computing vext_coulomb due to boundary conditions\n");
         }
     }
     if (Iwrite==VERBOSE) {
        print_vext(Vext,output_file2);
        if (Restart_Vext == READ_VEXT_STATIC) print_vext(Vext_static,"dft_vext_static.dat");
        print_zeroTF(Zero_density_TF,output_TF);
        if (Vol_charge_flag && Ndim==3) print_vext(Vext_coul,output_file4);
     }


     /*
      * Use Newton's method to solve the problem
      * (see file dft_newton.c, which calls dft_fill.c and Aztec library)
      * The variable "niters" is positive if Newton's method converrged.
      */
      x = (double **) array_alloc (2, Nunk_per_node, Nnodes_box, sizeof(double));
      if (Lbinodal) x2 = (double **) array_alloc (2, Nunk_per_node, Nnodes_box, sizeof(double));
		 
      t_preprocess += MPI_Wtime();
      t_solve = -MPI_Wtime();
      if (NL_Solver == PICARD_BUILT_IN || NL_Solver==PICARD_NOX ||
           NL_Solver==PICNEWTON_NOX || NL_Solver==PICNEWTON_BUILT_IN){
          if (NL_Solver==PICNEWTON_NOX || NL_Solver==PICNEWTON_BUILT_IN) NL_update_scalingParam/=100.;
          if (Proc==0) printf("Calling Picard Solver!\n");
          niters = solve_problem_picard(x, x2);
          if (NL_Solver==PICNEWTON_NOX || NL_Solver==PICNEWTON_BUILT_IN){
               NL_update_scalingParam*=100.;
              if (Proc==0) printf("Calling Newton Solver!\n");
               niters = solve_problem(x, x2);
          }
      }
      else if (NL_Solver == NEWTON_BUILT_IN || NL_Solver==NEWTON_NOX){
            niters = solve_problem(x, x2);
      }
      else {
         printf("Problem with solver type...set to %d\n",NL_Solver);
         exit(-1);
      }   
      t_solve += MPI_Wtime();
     /*
      * Post-Process the results
      * (see file dft_output.c)
      */
      time_save =  MPI_Wtime() - t_mesh;
      t_postprocess = -MPI_Wtime();

      if (NL_Solver !=NEWTON_NOX && NL_Solver != PICARD_NOX){
      if (Loca.method == -1 || Nruns>1) {
        if (Lbinodal) post_process(x2, &niters, &time_save,Imain_loop, TRUE,FROM_MAIN);
        post_process(x, &niters, &time_save, Imain_loop, FALSE,FROM_MAIN);
      }
      }
      Nodes_old = Nnodes;
      for (idim=0; idim<Ndim; idim++) Nodes_x_old[idim] = Nodes_x[idim];


      t_postprocess += MPI_Wtime();
/*      printf("\n*********************************************\n");
      printf("The number of nonzeros counted is %d\n",count_nonzero);
      printf("The number of zeros counted is %d\n",count_zero);
      printf("*********************************************\n");
*/

      iend=0;
      match=0;
      for (idim = 0; idim<Ndim; idim++)
           if (Esize_x[idim]==Esize_x_final[idim]) match++;
      if(match==Ndim)  iend=1; /* we've done the final mesh */ 
 
      for (idim = 0; idim<Ndim; idim++){
        if (Esize_x[idim] > Esize_x_final[idim]) {
            if (Esize_x[idim] == 0.2) Esize_x[idim] = 0.1;
            else if (Esize_x[idim] == 0.1) Esize_x[idim] = 0.05;
            else if (Esize_x[idim] == 0.05) Esize_x[idim] = 0.02;
            else if (Esize_x[idim] == 0.02) Esize_x[idim] = 0.01;
            else Esize_x[idim] = Esize_x_final[idim];
        }
        if (Esize_x[idim] < Esize_x_final[idim]) 
            Esize_x[idim] = Esize_x_final[idim];
      }

     /*
      * Get total CPU time for run 
      * (see files md_timer_*.c)
      */

      t_total = MPI_Wtime() - start_t;
      t_pre_max= gmax_double(t_preprocess);
      t_solve_max=gmax_double(t_solve);
      t_post_max=gmax_double(t_postprocess);
      t_total_max=gmax_double(t_total);

      t_pre_min= gmin_double(t_preprocess);
      t_solve_min=gmin_double(t_solve);
      t_post_min=gmin_double(t_postprocess);
      t_total_min=gmin_double(t_total);

      if (L_Schur != 2){
      MPI_Gather(&Time_linsolver_first,1,MPI_DOUBLE,t_linsolv_first_array,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Gather(&Time_linsolver_av,1,MPI_DOUBLE,t_linsolv_av_array,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Gather(&Time_manager_first,1,MPI_DOUBLE,t_manager_first_array,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Gather(&Time_manager_av,1,MPI_DOUBLE,t_manager_av_array,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Gather(&Time_fill_first,1,MPI_DOUBLE,t_fill_first_array,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Gather(&Time_fill_av,1,MPI_DOUBLE,t_fill_av_array,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      if (Proc==0 && Iwrite==VERBOSE){
         fp  = fopen("dft_time.out","w");
            fprintf(fp,"\n Time histogram for the fill on the first iteration\n\n");
            for (i=0;i<Num_Proc;i++) fprintf(fp,"\t %d  %9.6f\n",i,t_fill_first_array[i]);
            fprintf(fp,"\n Time histogram for the fill averaging the 2-%d iterations\n\n",niters);
            for (i=0;i<Num_Proc;i++) fprintf(fp,"\t %d  %9.6f\n",i,t_fill_av_array[i]/((double)niters-1.0));

            fprintf(fp,"\n Time histogram for the linear solver manager on the first iteration\n\n");
            for (i=0;i<Num_Proc;i++) fprintf(fp,"\t %d  %9.6f\n",i,t_manager_first_array[i]);
            fprintf(fp,"\n Time histogram for the linear solver manager averaging the 2-%d iterations\n\n",niters);
            for (i=0;i<Num_Proc;i++) fprintf(fp,"\t %d  %9.6f\n",i,t_manager_av_array[i]/((double)niters-1.0));

            fprintf(fp,"\n Time histogram for the linear solve on the first iteration\n\n");
            for (i=0;i<Num_Proc;i++) fprintf(fp,"\t %d  %9.6f\n",i,t_linsolv_first_array[i]);
            fprintf(fp,"\n Time histogram for the linear solve averaging the 2-%d iterations\n\n",niters);
            for (i=0;i<Num_Proc;i++) fprintf(fp,"\t %d  %9.6f\n",i,t_linsolv_av_array[i]/((double)niters-1.0));
         fclose(fp); 
      }

      t_linsolv_first_min=gmin_double(Time_linsolver_first);
      t_linsolv_first_max=gmax_double(Time_linsolver_first);
      t_manager_first_min=gmin_double(Time_manager_first);
      t_manager_first_max=gmax_double(Time_manager_first);
      t_fill_first_min=gmin_double(Time_fill_first);
      t_fill_first_max=gmax_double(Time_fill_first);

      min_nnodes_per_proc=gmin_int(Nnodes_per_proc);
      max_nnodes_per_proc=gmax_int(Nnodes_per_proc);
      min_nnodes_box=gmin_int(Nnodes_box);
      max_nnodes_box=gmax_int(Nnodes_box);
      min_nodesLoc_over_nodesBox=gmin_double((double)Nnodes_box/(double)Nnodes_per_proc);
      max_nodesLoc_over_nodesBox=gmax_double((double)Nnodes_box/(double)Nnodes_per_proc);

      if(niters > 1){
        t_linsolv_av_min=gmin_double(Time_linsolver_av/((double)niters-1.0));
        t_linsolv_av_max=gmax_double(Time_linsolver_av/((double)niters-1.0));
        t_manager_av_min=gmin_double(Time_manager_av/((double)niters-1.0));
        t_manager_av_max=gmax_double(Time_manager_av/((double)niters-1.0));
        t_fill_av_min=gmin_double(Time_fill_av/((double)niters-1.0));
        t_fill_av_max=gmax_double(Time_fill_av/((double)niters-1.0));
      }
      else{
        t_linsolv_av_min = t_linsolv_first_min;
        t_linsolv_av_max = t_linsolv_first_max;
        t_manager_av_min = t_manager_first_min;
        t_manager_av_max = t_manager_first_max;
        t_fill_av_min = t_fill_first_min;
        t_fill_av_max = t_fill_first_max;
      }
      }
      t_NLSolve_min=gmin_double(Time_NLSolve);
      t_NLSolve_max=gmax_double(Time_NLSolve);
      t_InitGuess_min=gmin_double(Time_InitGuess);
      t_InitGuess_max=gmax_double(Time_InitGuess);


      if (Proc == 0 &&Iwrite !=NO_SCREEN) {
        printf ("\n\n\n\n");
        printf ("===================================================\n");
        printf ("MESH SUMMARY ....");
        printf ("                    min             max       \n");
        printf ("---------------------------------------------------\n");
        printf ("---------------------------------------------------\n");
        printf ("Nodes_per_proc                %d         %d       \n", min_nnodes_per_proc,max_nnodes_per_proc);
        printf ("Nodes_box                     %d         %d       \n", min_nnodes_box,max_nnodes_box);
        printf ("Nodes_local/Nnodes_box        %g         %g       \n", 1./min_nodesLoc_over_nodesBox,1./max_nodesLoc_over_nodesBox);
        printf ("===================================================\n");
        printf ("===================================================\n");
        printf ("\n\n");
        printf ("===================================================\n");
        if (NL_Solver !=NEWTON_NOX && NL_Solver !=PICARD_NOX) 
        printf ("TIMING SUMMARY .... averaged over all %d iterations\n",niters);
        printf ("===================================================\n");
        printf ("                         min(sec)        max(sec)  \n");
        printf ("---------------------------------------------------\n");
        printf ("---------------------------------------------------\n");
        printf ("MESH SETUP               %g         %g       \n",
                                                    t_pre_max,t_pre_min); 
        printf ("---------------------------------------------------\n");
        printf ("---------------------------------------------------\n");
        printf ("TOTAL SOLVE              %g         %g       \n",
                                                t_solve_min,t_solve_max);
        printf ("---------------------------------------------------\n");
        printf ("---------------------------------------------------\n");
        printf ("INNER NONLINEAR SOLVER       %g         %g       \n",
                                                t_NLSolve_min,t_NLSolve_max);
        printf ("---------------------------------------------------\n");
        printf ("---------------------------------------------------\n");
        printf ("TIME FOR INITIAL GUESSS      %g         %g       \n",
                                                t_InitGuess_min,t_InitGuess_max);
        printf ("---------------------------------------------------\n");
        printf ("---------------------------------------------------\n");
        printf ("First fill time          %g         %g       \n",
                                                t_fill_first_min,t_fill_first_max);
        printf ("2-niter avg fill time    %g         %g       \n",
                                                t_fill_av_min,t_fill_av_max);
        if (L_Schur != 2){
        printf ("First manager time       %g         %g       \n",
                                                t_manager_first_min,t_manager_first_max);
        printf ("2-niter avg manager time %g         %g       \n",
                                                t_manager_av_min,t_manager_av_max);
        printf ("First linsolv time       %g         %g       \n",
                                                t_linsolv_first_min,t_linsolv_first_max);
        printf ("2-niter avg linsolv time %g         %g       \n",
                                                t_linsolv_av_min,t_linsolv_av_max);
        }
        printf ("---------------------------------------------------\n");
        printf ("---------------------------------------------------\n");
        printf ("POST-PROCESSING          %g         %g       \n",
                                                  t_post_min,t_post_max);
        printf ("===================================================\n");
        printf("TOTAL TIME               %g         %g       \n",
                                                t_total_min,t_total_max);
        printf ("===================================================\n");

      }

      }   /* end of loop over different mesh densities (for automated
             mesh refinement */

      /*
       * release lots of memory to start a new mesh
       */


      /* thermo arrays */
      if (Ipot_ff_c == COULOMB) 
           safe_free((void *) &Deltac_b);

      /* solution arrays*/
      safe_free((void *) &B2L_unknowns);
      safe_free((void *) &x); 
      if (Lbinodal) safe_free((void *) &x2); 

    if (Proc==0){  
         safe_free((void *)&t_fill_first_array);
         safe_free((void *)&t_fill_av_array);
         safe_free((void *)&t_manager_first_array);
         safe_free((void *)&t_manager_av_array);
         safe_free((void *)&t_linsolv_first_array);
         safe_free((void *)&t_linsolv_av_array);
     }

     safe_free((void *)&Nel_hit);
     safe_free((void *)&Nel_hit2);
  } /* end of loop over continuation field #1 */

  if (Nwall_type !=0){
    safe_free((void *) &Link_list);
    safe_free((void *) &Nwall_this_link);
  }
  free_mesh_arrays();
  boundary_free();
  if (Type_interface==DIFFUSIVE_INTERFACE && Ndim==1) safe_free((void *) &Area_IC);
  safe_free((void *) &Comm_node_proc);
  safe_free((void *) &Comm_unk_proc);
  safe_free((void *) &Comm_offset_node);
  safe_free((void *) &Comm_offset_unk);

  safe_free((void *)  &Xtest_reflect_TF);
  safe_free((void *)  &Lsemiperm);
  safe_free((void *)  &Vext_membrane);
/*
  safe_free((void *)  &Unk_to_Poly);
  safe_free((void *)  &Unk_to_Seg);
  safe_free((void *)  &Unk_to_Bond);
  safe_free((void *)  &Bonds);
  safe_free((void *)  &Nbond);
*/


  /* free stencils */
  for (izone=0; izone<Nzone; izone++){
  for (isten=0; isten<NSTEN; isten++)
   if (Sten_Type[isten]==TRUE) {
     jmax=stencil_Njcomp_switch(isten);

     for (icomp=0; icomp<Ncomp; icomp++) {
       for (jcomp=0; jcomp<jmax; jcomp++) {
        sten = &(Stencil[isten][izone][icomp+Ncomp*jcomp]);
        safe_free((void **) &sten->Offset);
        safe_free((void **) &sten->Weight);
        if (Lhard_surf) safe_free((void **) &sten->HW_Weight);
       }
     }
   }
  }
  safe_free((void **) &Stencil);

  if (Type_poly == CMS) safe_free((void **) &Rism_cr);


  if (Proc == 0){
    safe_free((void *) &X_old); /* This is also done for each continuation run
              in set_initial_guess (dft_guess.c) called from solve_problem.
              The array is allocated in collect_x_old (dft_output.c) called
              from post_process.    */
    if (Lbinodal) safe_free((void *) &X2_old); 
  }
  *engptr = Energy;
  return;
} 
/************END of main() *************************************************/


/****************************************************************************/
/*continutation_shift:  In this routine we set up parameters for a continuation
                        where the walls have been shifted apart.*/
void continuation_shift()
{
  int idim,iwall,icharge,nadd;
  double size_y_tmp;

  /* special case of continuation for a hexagonal system */
  /* commenting out of main source for now */
  size_y_tmp= sqrt(3*(Size_x[0]+Del_1[0])*(Size_x[0]+Del_1[0]));
  nadd = round_to_int((size_y_tmp-Size_x[1])/Esize_x[1]);
  Del_1[1]=nadd*Esize_x[1];
  /* end special case */

  for (idim=0; idim<Ndim; idim++) {
      Size_x[idim] += Del_1[idim];
                                   /*comment out special case for hexagonal mesh change */
      if (Plane_new_nodes == idim || round_to_int(Del_1[idim]/Esize_x[idim])>0){
         if (Pos_new_nodes == 0) {    /*adding nodes to center*/
            for (iwall=0; iwall<Nwall; iwall++) {
               if (WallPos[idim][iwall]<0.0)
                 WallPos[idim][iwall] -= 0.5*Del_1[idim];
               else if (WallPos[idim][iwall] > 0.0)
                 WallPos[idim][iwall] += 0.5*Del_1[idim];
            }
            for (icharge=0;icharge<Nlocal_charge;icharge++){
               if (Charge_x[icharge][idim] <0.0)
                 Charge_x[icharge][idim] -=0.5*Del_1[idim];
               else if (Charge_x[icharge][idim] > 0.0)
                 Charge_x[icharge][idim] +=0.5*Del_1[idim];
            }
         }
         else if (Pos_new_nodes == -1) {    /*adding nodes to lbb*/
            for (iwall=0; iwall<Nwall; iwall++)
                 WallPos[idim][iwall] += 0.5*Del_1[idim];
            for (icharge=0;icharge<Nlocal_charge;icharge++)
                 Charge_x[icharge][idim] +=0.5*Del_1[idim];
         }
         else if (Pos_new_nodes == 1) {    /*adding nodes to rft*/
            for (iwall=0; iwall<Nwall; iwall++)
                  WallPos[idim][iwall] -= 0.5*Del_1[idim];
            for (icharge=0;icharge<Nlocal_charge;icharge++)
                 Charge_x[icharge][idim] -=0.5*Del_1[idim];
         }
      } /* end of check for idim = Plane_new_nodes */
  }  /* end of loop over dimentions */
  return;
}
/******************************************************************/
