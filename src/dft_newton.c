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
 *  FILE: dft_newton.c
 *
 *  This file solves the dft problem using Newton's method. The fill
 *  of the residual and matrix are in file dft_fill.c, and the linear
 *  solve is done using calls to the Aztec library.
 */
#include "mpi.h" 


#ifdef HAVE_DFT_SCHUR_SOLVER
#include "dft_schur_solver.h"
#endif

#include "dft_globals_const.h" 
#include "rf_allo.h"

int update_solution(int,double *, double *);
void full_to_msr(double **, int *, double *, int, int);
static void set_multiple_B2L(int *ijk_box, int iunk, int loc_i);
#ifdef NUMERICAL_JACOBIAN
static void do_numerical_jacobian(double *x, double *resid, double *resid_tmp);
#endif
double dotp(double *, double *);


int solve_problem(double **x_internal_ptr, double **x2_internal_ptr,
                  char *output_file1,int *niters)
/*
 * This routine drives Newton's method for finding the converged solution.
 * It returns TRUE if Newton's method converged successfully, else FALSE
 *
 * (*x_internal) is the array of unknowns. The address is passed so that
 * it can be re-allocated in case of load balancing
 */

{
 /*
  * Local variable declarations
  */
#ifdef HAVE_DFT_SCHUR_SOLVER
  DFT_SCHUR_SOLVER *schur_solver;
#endif
  void * aux_info = NULL;
  char *yo = "solve_problem";
  int   iter=0, i,j;  /* iteration counter for Newton's method */
  double   t1;
  int    **bindx_2d=NULL; /* Temp matrix indexing array before size is known*/
  double *x, *x2=NULL;  /* Solution vector, with space for externals */
  double *fill_time = NULL; /* Array to hold fill time for each node */
  double *fill_non0 = NULL; /* Array to hold fill unks for each node */
  int outer, outer_start, max_iter;
  int max_nonzeros=1, min_nonzeros, tot_nonzeros, max_ints,max_exts,min_ints,min_exts; 
  double time_preproc=0.0, t_solve_max=0.0;

  int ijk[3], ijk_box[3]= {0,0,0}, inode_box, iunk;
  int inode,loc_global;
  FILE *fp1=NULL;

#define REGULAR_LOOP 1
#define ONE_ITER_LOOP REGULAR_LOOP-1

 
  /********************** BEGIN EXECUTION ************************************/

  *niters = iter;
  if (Proc==0 && Iwrite != NO_SCREEN) printf("\n%s: Doing Newton's method to solve the problem\n",yo);

  Rhobar3_old = (double *) array_alloc (1, Nnodes_box, sizeof(double));

  /* put in logic for load balancing after 1 completed Newton iter */

  if (Load_Bal_Flag >= LB_TIMINGS) {
       outer_start = ONE_ITER_LOOP;
  }
  else outer_start = REGULAR_LOOP;

  for (outer=outer_start; outer <= REGULAR_LOOP; outer ++) {

  if (outer == ONE_ITER_LOOP) {
    max_iter = 1;
    fill_time = (double *) array_alloc(1, Nnodes_per_proc, sizeof(double));
    fill_non0 = (double *) array_alloc(1, Nnodes_per_proc, sizeof(double));
  }
  else {
     max_iter = Max_Newton_iter;

     if (Load_Bal_Flag >= LB_TIMINGS && !Lbinodal) { /* need Xold2 for Lbinodal */

       /* Save current solution after 1 Newton iter in global vec */
       if (Proc == 0)
          X_old = (double *) array_alloc (1, Nnodes*Nunk_per_node, sizeof(double));

       collect_x_old(*x_internal_ptr,1);

       safe_free((void *) &*x_internal_ptr);
       
       /* load balancing based on fill_time per node */
       if      (Load_Bal_Flag == LB_TIMINGS) load_balance(1, fill_time);
       else if (Load_Bal_Flag == LB_NON0)    load_balance(1, fill_non0);
       else if (Load_Bal_Flag == LB_MIXED) {
          /* weight fill time and # of unknowns by relative times of
             matrix fill and solve. Convert # unknowns to units of
             time using max_nonzeros t_solve_max */
          /* printf("Proc %d: max nonzeros %d and t_solve_max %g\n",Proc,max_nonzeros,t_solve_max);*/
           for (i=0; i<Nnodes_per_proc; i++)
              fill_time[i] += fill_non0[i] * t_solve_max / (double) max_nonzeros;
           load_balance(1, fill_time);
       }

       MPI_Gather(&Nnodes_per_proc,1,MPI_INT,
             Comm_node_proc,1,MPI_INT,0,MPI_COMM_WORLD);

       if (Proc == 0){
          for (i=0; i<Num_Proc; i++)
            Comm_unk_proc[i] = Comm_node_proc[i]*Nunk_per_node;
 
          Comm_offset_node[0] = 0; Comm_offset_unk[0] = 0;
          for (i=1; i<Num_Proc; i++){
             Comm_offset_node[i] = Comm_offset_node[i-1] + Comm_node_proc[i-1];
             Comm_offset_unk[i]  = Comm_offset_unk[i-1]  + Comm_unk_proc[i-1];
          }
       }

       /* And redo the mesh variables based on the new load balancing. */
       if (Proc==0) {
	 if(Iwrite != NO_SCREEN) printf("\n%s: Setting up the mesh again ... ",yo);
	 fp1 = fopen("dft_out.lis","a+");
       }
       free_mesh_arrays();
       control_mesh(fp1,"dft_vext2.dat",FALSE);

       /* *x_internal_ptr local array must be redistributed global->local */
      
       *x_internal_ptr = (double *)
                         array_alloc(1, Aztec.N_update, sizeof(double));

      /*
       * With new load balancing, need to free local arrays 
       * and redo the boundaries
       */

       safe_free((void *) &fill_time); fill_time = NULL;
       safe_free((void *) &fill_non0); fill_non0 = NULL;
       if (Ipot_ff_n != IDEAL_GAS)safe_free((void *) &B2L_1stencil);
       safe_free((void *) &B2L_unknowns);
       safe_free((void *) &Aztec.external);
       safe_free((void *) &Aztec.extern_index);
       safe_free((void *) &Aztec.update_index);
       safe_free((void *) &Aztec.data_org);

       boundary_free();
       boundary_setup(output_file1);


     }
  }

 /*
  *  If filling directly into MSR matrix, do pre-processing loop to get
  *  size of arrays. Allocate MSR arrays.
  */

   if (Proc==0 && Iwrite !=NO_SCREEN){
     printf("\n=====================================");
     printf("=======================================\n");
     printf("%s: Preprocessing MSR matrix\n",yo);
   }

   t1 = MPI_Wtime();
   /* Now continue with pre-processing step to allocate teh matrix */

   /* Form temporary array to keep track of nonzeros per row */

   bindx_2d = (int **) array_alloc(1, Aztec.N_update+1, sizeof(int *));
   bindx_2d[Aztec.N_update] = 
                   (int *) array_alloc(1, Aztec.N_update, sizeof(int));

   /* alert fill of first time through fill after a preprocessing  */

/*   if ((Load_Bal_Flag < LB_TIMINGS && iter==1) ||
       (Load_Bal_Flag >= LB_TIMINGS && iter==2)) first = TRUE;
   else                                          first = FALSE;

   if (Iwrite != NO_SCREEN)printf("FIRST=%d  iter: %d \n",first,iter);
*/

   /* Call fill, number of nonzeros returned in Aztec.nonzeros */

   T_msr_setup = -MPI_Wtime();
   fill_resid_and_matrix_control(NULL, NULL, bindx_2d, NULL, MSR_PREPROCESS, iter, FALSE);
   T_msr_setup += MPI_Wtime();


   /* Allocate MSR matrix with known size */

   min_nonzeros = AZ_gmin_int(Aztec.nonzeros,Aztec.proc_config);
   max_nonzeros = AZ_gmax_int(Aztec.nonzeros,Aztec.proc_config);
   tot_nonzeros = AZ_gsum_int(Aztec.nonzeros,Aztec.proc_config);

   time_preproc = MPI_Wtime()-t1;

   if (Proc ==0 && Iwrite != NO_SCREEN) {
      printf("\n\tMSR Preproc found %d to %d (total %d) nonzeros\n",
                                min_nonzeros,max_nonzeros, tot_nonzeros);
      printf("\tMSR Preproc took %g seconds\n", time_preproc);
   }

   Aztec.bindx  = (int *) array_alloc(1, Aztec.nonzeros+1, sizeof(int));
   Aztec.val = (double *) array_alloc(1, Aztec.nonzeros+1, sizeof(double));
   if (Aztec.val==NULL) {
      printf("%s ERROR: Not enough memory to allocate MSR matrix\n",yo);
      exit(-1);
   }


   /* load diagnal entries of bindx */

   Aztec.bindx[0] = Aztec.N_update + 1;
   for (i=0; i< Aztec.N_update; i++){
     Aztec.bindx[i+1] = Aztec.bindx[i] + bindx_2d[Aztec.N_update][i];
     /*printf("i: %d  Aztec.bindx[i]: %d  bindx_2d[Aztec.N_update][i]: %d\n",
              i,Aztec.bindx[i],bindx_2d[Aztec.N_update][i]);*/
   }

   /* load off-diagnal entries of bindx */

   for (i=0; i<Aztec.N_update; i++){
     for (j=0; j<bindx_2d[Aztec.N_update][i]; j++){
/*       printf("i=%d  j=%d trial1=%d \n", i,j,bindx_2d[i][j]);*/
       Aztec.bindx[Aztec.bindx[i] + j] = bindx_2d[i][j];
     } 
   }

   for (i=0; i<Aztec.N_update+1; i++) safe_free((void *) &bindx_2d[i]);
   safe_free((void *) &bindx_2d);

   /* FIll number of unknowns for each node into array for load balancing */

   if (Load_Bal_Flag >= LB_NON0 && iter==0) {
     for (i=0; i<Nnodes_per_proc; i++)
        fill_non0[i] = Aztec.bindx[(i+1)*Nunk_per_node] - Aztec.bindx[i*Nunk_per_node]
                       + Nunk_per_node; 
   }

   /*
    * Transform from global to local (on this Processor) numbering scheme,
    * matrix then resid. Change sign on resid for use in Newton's method.
    * Set the matrix "name" to 1 every Newton iter so Aztec doesn't  think
    * it's a different matrix and waste memory.
    */

/*   if (Num_Proc>1) MPI_Barrier(MPI_COMM_WORLD);
   printf("Proc %d: Aztec.N_update = %d, nonzeros %d\n",Proc,Aztec.N_update,
           Aztec.bindx[Aztec.bindx[0]-1]);
   for (i=0; i<Aztec.bindx[Aztec.bindx[0]-1]; i++)
     printf("Proc %d: %d. bindx = %d  \n", Proc,i, Aztec.bindx[i]);
*/
  
   AZ_transform(Aztec.proc_config, &(Aztec.external), Aztec.bindx, 
                Aztec.val, Aztec.update, &(Aztec.update_index),
                &(Aztec.extern_index), &(Aztec.data_org), Aztec.N_update, 
                NULL, NULL, NULL, NULL, AZ_MSR_MATRIX);

   Aztec.data_org[AZ_name] = 1;

#ifdef HAVE_DFT_SCHUR_SOLVER

   dft_create_schur_solver(Aztec.proc_config, Aztec.external, Aztec.bindx, 
			   Aztec.val, Aztec.update, Aztec.update_index,
			   Aztec.extern_index, Aztec.data_org, Aztec.N_update,
			   Nunk_per_node, &schur_solver);
   aux_info = (void *) schur_solver;
#endif

   /* Allocate solution vector x with space for externals */
   /* and communicate initial guess between procs         */

   Nunk_int_and_ext = Aztec.N_update + Aztec.data_org[AZ_N_external];
   x = (double *) array_alloc (1, Nunk_int_and_ext, sizeof(double));
   if (Lbinodal) x2 = (double *) array_alloc (1, Nunk_int_and_ext, sizeof(double));

   if (iter == 0 || (iter==1 && Load_Bal_Flag >= LB_TIMINGS )) {
      set_initial_guess(Iguess1,*x_internal_ptr,&iter);
      if (Lbinodal) set_initial_guess(BINODAL_FLAG,*x2_internal_ptr,&iter);
      if (Proc==0 && iter==1)  safe_free((void *) &X_old);

/*      if (Imain_loop == 0) 
 *        for (i=0; i<200; i++) Hist_time[i] = 0;
 */
   }

   for (i=0; i < Aztec.N_update; i++) {
      x[Aztec.update_index[i]] = (*x_internal_ptr)[Aztec.update_index[i]];
   }
   AZ_exchange_bdry(x, Aztec.data_org,Aztec.proc_config);

   if (Lbinodal){
     for (i=0; i < Aztec.N_update; i++) 
        x2[Aztec.update_index[i]] = (*x2_internal_ptr)[Aztec.update_index[i]];
     AZ_exchange_bdry(x2, Aztec.data_org,Aztec.proc_config);
   }

   /* Set up Box to local array */

   B2L_unknowns = (int *) array_alloc(1, Nunknowns_box, sizeof(int));
   Az2G_unknowns = (int *) array_alloc(1, Aztec.N_update+Aztec.data_org[AZ_N_external], sizeof(int));
   Az2G_unknowns_by_type = (int *) array_alloc(1, Aztec.N_update+Aztec.data_org[AZ_N_external], sizeof(int));
   Az2eq_type = (int *) array_alloc(1, Aztec.N_update+Aztec.data_org[AZ_N_external], sizeof(int));
   Nunk_eq_type = (int *) array_alloc(1, Ntype_blocks, sizeof(int));

   for (i=0; i < Nunknowns_box; i++) 
      B2L_unknowns[i] = -1;

   for (i=0;i<Ntype_blocks;i++){
      Nunk_eq_type[i]=0;
   }

   for (i=0; i < Aztec.N_update; i++) {
     if (MATRIX_FILL_NODAL){
        node_to_ijk(Aztec.update[i]/Nunk_per_node,ijk);
        iunk = Aztec.update[i] - Nunk_per_node * ijk_to_node(ijk);
     }
     else{
        iunk = Aztec.update[i]/Nnodes_per_proc;
        node_to_ijk(Aztec.update[i]-Nnodes_per_proc*iunk,ijk);
     }
     inode = ijk_to_node(ijk);
     loc_global = loc_find(iunk,inode,GLOBAL);
     Az2G_unknowns[i]=loc_global;
     Az2G_unknowns_by_type[Block_type[Unk_to_eq_type[iunk]][i]=loc_global;
     Az2eq_type[i]=Block_type[Unk_to_eq_type[iunk]];
     Nunk_eq_type[Block_type[Unk_to_eq_type[iunk]]]++;

     ijk_to_ijk_box(ijk, ijk_box);
     inode_box = ijk_box_to_node_box(ijk_box);
     B2L_unknowns[loc_find(iunk,inode_box,BOX)] = Aztec.update_index[i];

     /*
      * With periodic BC, the same global node can be in multiple
      * box locations -- find all periodic images.  Non_unique_G2B[i]
      * holds the maximum # of wraps in each dimension.
      */

     if (Non_unique_G2B[3])
       set_multiple_B2L(ijk_box, iunk, Aztec.update_index[i]);
   }
   /* loop over block types in this problem ... call the row map for this processor */
   for (i=0;i<Ntype_blocks;i++){
     dft_solvermanager_setrowmap(solver_manager, i, Nunk_eq_type[i], Az2G_unknowns_by_type[i]); 
   }

   for (i=0; i < Aztec.data_org[AZ_N_external]; i++){
     if (MATRIX_FILL_NODAL){
        node_to_ijk(Aztec.external[i]/Nunk_per_node,ijk);
        iunk = Aztec.external[i] - Nunk_per_node * ijk_to_node(ijk);
     }
     else{
        iunk = Aztec.external[i]/(Aztec.data_org[AZ_N_external]/Nunk_per_node);
        node_to_ijk(Aztec.external[i]-(Aztec.data_org[AZ_N_external]/Nunk_per_node)*iunk,ijk);
     }
     inode = ijk_to_node(ijk);
     loc_global = loc_find(iunk,inode,GLOBAL);
     Az2G_unknowns[i+Aztec.N_update]=loc_global;
     Az2G_unknowns_by_type[Block_type[Unk_to_eq_type[iunk]][i+Aztec.N_update]=loc_global;
     Az2eq_type[i+Aztec.N_update]=Unk_to_eq_type[iunk];
     Nunks_eq_type[Block_type[Unk_to_eq_type[iunk]]]++;

     ijk_to_ijk_box(ijk, ijk_box);
     inode_box = ijk_box_to_node_box(ijk_box);
     B2L_unknowns[loc_find(iunk,inode_box,BOX)] = Aztec.extern_index[i];

     /*
      * With periodic BC, the same global node can be in multiple
      * box locations -- find all periodic images.  Non_unique_G2B[i]
      * holds the maximum # of wraps in each dimension.
      */

     if (Non_unique_G2B[3]) 
       set_multiple_B2L(ijk_box, iunk, Aztec.extern_index[i]);
   }
   /* loop over block types in this problem ... call the row map for this processor */
   for (i=0;i<Ntype_blocks;i++){
     dft_solvermanager_setcolmap(solver_manager, i, Nunk_eq_type[i], Az2G_unknowns_by_type[i]); 
   }

   max_ints = AZ_gmax_int(Aztec.N_update,Aztec.proc_config);
   max_exts = AZ_gmax_int(Aztec.data_org[AZ_N_external],Aztec.proc_config);
   min_ints = AZ_gmin_int(Aztec.N_update,Aztec.proc_config);
   min_exts = AZ_gmin_int(Aztec.data_org[AZ_N_external],Aztec.proc_config);

   if (Proc ==0 && Iwrite != NO_SCREEN) { 
       printf("\n\tmax # externals %d : min # externals: %d\n",
                                        max_exts,min_exts);
       printf("\tmax # internals %d : min # internals: %d\n",
                                        max_ints,min_ints);
   }

   /* print out load balancing info */

   if (Iwrite == VERBOSE) {
     if (Num_Proc>1) MPI_Barrier(MPI_COMM_WORLD);
     printf("Proc %2d: nonzeros %d, internals %d, externals %d\n", Proc,
             Aztec.nonzeros, Aztec.N_update, Aztec.data_org[AZ_N_external]);
   }

/* Call LOCA continuation library if it is requested and compiled in */

  if (outer != ONE_ITER_LOOP && Loca.method != -1)
#ifdef LOCA
    iter += solve_continuation(x, x2, Sten_Type[POLYMER_CR], aux_info);
#endif
  else iter += newton_solver(x, x2, fill_time, NULL, max_iter, &t_solve_max, aux_info);

  /* return solution vector with only internals, unsorted */

  for (i=0; i < Aztec.N_update; i++)
       (*x_internal_ptr)[i] = x[Aztec.update_index[i]];
  if (Lbinodal) 
    for (i=0; i < Aztec.N_update; i++)
       (*x2_internal_ptr)[i] = x2[Aztec.update_index[i]];

  /* free memory from local arrays */

#ifdef HAVE_DFT_SCHUR_SOLVER
  dft_destroy_schur_solver(&schur_solver);
#endif

  safe_free((void *) &Aztec.bindx);
  safe_free((void *) &Aztec.val);
  safe_free((void *) &x);
  safe_free((void *) &Rhobar3_old);
  if (Lbinodal) safe_free((void *) &x2);

  } /* end of outer loop for possible re-load balancing */

  *niters = iter;

  return(iter);
}


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
int newton_solver(double *x, double *x2, double *fill_time, void *con_ptr,
		  int max_iter, double *t_solve_max_ptr, void * aux_info)
{
 /*
  * Local variable declarations
  */
#ifdef HAVE_DFT_SCHUR_SOLVER
  DFT_SCHUR_SOLVER * schur_solver;
#endif

  char *yo = "newton_solver";
  int   iter=0, i;  /* iteration counter for Newton's method */
  int   converged = FALSE, converged2=TRUE;  /* Convergence flag ste to True or False */
  double   t1, t2=0., l2_resid;
  double  *delta_x=NULL;/* Update of solution vector returned by matrix solve*/
  double  *resid, *resid_tmp;    /* Residual vector */
  double t_solve=0.0, t_solve_min, t_solve_max;
  int nodesave;
  int ijk[3], iunk,inode;

 /*
  *  Do Newton's method with a do ... while loop.
  */

  delta_x = (double *) array_alloc (1, Nunk_int_and_ext, sizeof(double));
  resid     = (double *) array_alloc(1, Aztec.N_update, sizeof(double));
  resid_tmp = (double *) array_alloc(1, Aztec.N_update, sizeof(double));

  do {

    iter++;

    if (Proc==0 && Iwrite != NO_SCREEN) {
      printf("\n=====================================");
      printf("=======================================\n");
      printf("%s: Starting Newton iteration #%d\n",yo,iter);
    }

    for (i=0; i<Aztec.N_update; i++) resid_tmp[i] = 0.0;

    /* 
     * Fill residual and matrix (in bindx and val), and print CPU time 
     */

    for (i=0; i<Aztec.nonzeros+1; i++) Aztec.val[i] = 0.0;

    t1 = MPI_Wtime();
    fill_resid_and_matrix_control(x,resid_tmp,NULL,fill_time,Matrix_fill_flag,iter,FALSE);
    t2 = MPI_Wtime() - t1;

    if (Proc == 0 && Iwrite != NO_SCREEN) printf("\n\t%s: Matrix fill took %g seconds\n", yo, t2);


 /*     if (iter < 3 && Iwrite) { 
  * 
  *    int *index, loc_inode;
  *    double *proc_store,*proc_global=NULL,*time_global=NULL;
  *    *
  *    * Write timings to a file for a graphical (AVS) output !!
  *    * 

  *    index = (int *) array_alloc (1, Nnodes_per_proc, sizeof(int));
  *    proc_store = (double *) array_alloc(1,Nnodes_per_proc,sizeof(double));
  *    if (Proc == 0){
  *      time_global = (double *) array_alloc (1, Nnodes, sizeof(double));
  *      proc_global = (double *) array_alloc (1, Nnodes, sizeof(double));
  *    }

  *    for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
  *      inode = Aztec.update[Nunk_per_node * loc_inode] / Nunk_per_node;
  *      index[loc_inode] = inode;
  *      proc_store[loc_inode] = (double)Proc;
  *    }

  *    if (Proc == 0){
  *      for (inode=0; inode<Nnodes; inode++) {
  *          time_global[inode] = 0.0;
  *          proc_global[inode] = 0.0;
  *      }
  *    }

  *    gather_global_vec(fill_time , index, Nnodes_per_proc, time_global);
  *    gather_global_vec(proc_store, index, Nnodes_per_proc, proc_global);
 *
  *    safe_free((void *) &index);
  *    safe_free((void *) &proc_store);
  *    if (Proc == 0) {
  *        safe_free((void *) &time_global);
  *        safe_free((void *) &proc_global);
  *    }
  * }
  */

    /* Change sign on resid vector for Newton;s method, calc resid norm */

    l2_resid = 0.0;
    for (i=0; i<Aztec.N_update; i++) {

       resid[i] = -resid_tmp[i];
       l2_resid += resid_tmp[i]*resid_tmp[i];
    }

    l2_resid = AZ_gsum_double(l2_resid, Aztec.proc_config);
    if (Proc==0 && Iwrite != NO_SCREEN) printf("\t\t Norm of resid vector = %g\n", sqrt(l2_resid));

/*#define NUMERICAL_JACOBIAN*/
#ifdef NUMERICAL_JACOBIAN
    do_numerical_jacobian(x, resid, resid_tmp);
#endif

    /*
     * Call Aztec to solve matrix, after allocating the solution vector
     */

    for (i=0; i<Nunk_int_and_ext; i++) delta_x[i] = 0.0;

    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();

    AZ_free_memory(Aztec.data_org[AZ_name]);
#ifdef HAVE_DFT_SCHUR_SOLVER

    schur_solver = (DFT_SCHUR_SOLVER *) aux_info;
    
    dft_update_schur_solver(Aztec.bindx, Aztec.val, Aztec.update_index, schur_solver);
    dft_apply_schur_solver((DFT_SCHUR_SOLVER *) schur_solver, Aztec.update_index, delta_x, resid);
#else

    AZ_solve(delta_x, resid, Aztec.options, Aztec.params, NULL, Aztec.bindx, 
             NULL, NULL, NULL, Aztec.val, Aztec.data_org, Aztec.status, 
             Aztec.proc_config);
#endif
    t_solve = MPI_Wtime() - t1;
    t_solve_max = AZ_gmax_double(t_solve,Aztec.proc_config);
    t_solve_min = AZ_gmin_double(t_solve,Aztec.proc_config);
    if (Proc==0 && Iwrite != NO_SCREEN)
      printf("\tAZTEC solve time = %9.6f sec\n", t_solve);

    T_av_solve_max += t_solve_max;
    T_av_solve_min += t_solve_min;


    if (con_ptr != NULL) {
      for (i=0; i<Nunk_int_and_ext; i++) delta_x[i] *= -1.0;
      converged2 = continuation_hook(x, delta_x, con_ptr, 
		                     Newton_rel_tol, Newton_abs_tol);
      for (i=0; i<Nunk_int_and_ext; i++) delta_x[i] *= -1.0;
    }

    /*** start binodal calc */
    else if (Lbinodal) 
    {
    double *a, *b, *c, *d, *f1, *f2, g, ads, del_rb0, dmu_drb;
    a = (double *) array_alloc(1, Nunk_int_and_ext, sizeof(double));
    b = (double *) array_alloc(1, Nunk_int_and_ext, sizeof(double));
    c = (double *) array_alloc(1, Nunk_int_and_ext, sizeof(double));
    d = (double *) array_alloc(1, Nunk_int_and_ext, sizeof(double));
    f1= (double *) array_alloc(1, Aztec.N_update, sizeof(double));
    f2= (double *) array_alloc(1, Aztec.N_update, sizeof(double));

    setup_integrals();  /*to get the Nel_hit factors*/

    for (i=0; i<Nunk_int_and_ext; i++) {
       a[i] = delta_x[i];
       b[i] = 0.0;
       c[i] = 0.0;
       d[i] = 0.0;
    }
    for (i=0; i<Aztec.N_update; i++){
         f1[i] = resid_tmp[i]*Nel_hit2[0][Aztec.update_index[i]]
                             *Vol_el/((double)Nnodes_per_el_V);
    }

 /*  for (i=0; i < Aztec.N_update; i++) 
      x[Aztec.update_index[i]] = (*x_internal_ptr)[Aztec.update_index[i]];*/

    /* set up d_f/d_mu (later change to numerical deriv) */
    dmu_drb= (dmu_drho_hs(Rho_b) + dmu_drho_att(Rho_b));
    for (i=0; i<Nnodes_per_proc; i++) {
      for (iunk=0; iunk<Nunk_per_node; iunk++) {
        if (iunk==Unk_start_eq[DENSITY])  resid_tmp[B2L_unknowns[loc_find(iunk,L2B_node[i],BOX)]] = 
                              dmu_drb*Nel_hit2[0][Aztec.update_index[loc_find(iunk,i,LOCAL)]]
                                      *Vol_el/((double)Nnodes_per_el_V);
        else             resid_tmp[B2L_unknowns[loc_find(iunk,L2B_node[i],BOX)]]= 0.0;
      }
    }
    
    AZ_solve(b, resid_tmp, Aztec.options, Aztec.params, NULL, Aztec.bindx, 
             NULL, NULL, NULL, Aztec.val, Aztec.data_org, Aztec.status, 
             Aztec.proc_config);

    for (i=0; i<Aztec.N_update; i++) resid_tmp[i] = 0.0;
    for (i=0; i<Aztec.nonzeros+1; i++) Aztec.val[i] = 0.0;

    fill_resid_and_matrix_control(x2,resid_tmp,NULL,fill_time,Matrix_fill_flag,iter,FALSE);
    l2_resid = 0.0;
    for (i=0; i<Aztec.N_update; i++) {
       resid[i] = -resid_tmp[i];
       l2_resid += resid_tmp[i]*resid_tmp[i];
    }
    l2_resid = AZ_gsum_double(l2_resid, Aztec.proc_config);
    if (Proc==0 && Iwrite != NO_SCREEN) printf("\t\tNorm of resid vector = %g\n", sqrt(l2_resid));

    for (i=0; i<Aztec.N_update; i++)  
      f2[i] = resid_tmp[i]*Nel_hit2[0][Aztec.update_index[i]]*Vol_el/((double)Nnodes_per_el_V);
    for (i=0; i<Aztec.N_update; i++)  resid_tmp[i] *= -1.0;
    
    AZ_solve(c, resid_tmp, Aztec.options, Aztec.params, NULL, Aztec.bindx, 
             NULL, NULL, NULL, Aztec.val, Aztec.data_org, Aztec.status, 
             Aztec.proc_config);

    /* set up d_f/d_mu (later change to numerical deriv) */
    for (i=0; i<Nnodes_per_proc; i++) {
      for (iunk=0; iunk<Nunk_per_node; iunk++) {
        if (iunk ==Unk_start_eq[DENSITY])  resid_tmp[B2L_unknowns[loc_find(iunk,L2B_node[i],BOX)]] = 
                           dmu_drb*Nel_hit2[0][Aztec.update_index[loc_find(iunk,i,LOCAL)]]
                                                           *Vol_el/((double)Nnodes_per_el_V);
        else            resid_tmp[B2L_unknowns[loc_find(iunk,L2B_node[i],BOX)]] = 0.0;
      }
    }
    
    AZ_solve(d, resid_tmp, Aztec.options, Aztec.params, NULL, Aztec.bindx, 
             NULL, NULL, NULL, Aztec.val, Aztec.data_org, Aztec.status, 
             Aztec.proc_config);

    g = calc_free_energy(NULL,x,1.0,1.0,FALSE)-calc_free_energy(NULL,x2,1.0,1.0,FALSE);
    ads = -dmu_drb*(calc_adsorption(NULL,x,1.0,1.0)-calc_adsorption(NULL,x2,1.0,1.0));

/*if (iter<=2){    SAVE THIS FOR NUMERICAL CHECKING OF RESIDUALS (FREE ENERGY DERIVS) FOR
                   NEW FUNCTIONALS ADDED TO CODE !!!
double g1, g2, ee=1.0e-8, bt, bh, ba;
    printf("Starting numerical loop for last matrix row: dmu_drb=%g\n",dmu_drb);
    for (i=0; i<Nunk_int_and_ext; i++)  a[i]=0;
    for (i=0; i<Nunk_int_and_ext; i++)  b[i]=0;
    g1 = calc_free_energy(NULL,x,1.0,1.0,FALSE);
    g2 = calc_free_energy(NULL,x2,1.0,1.0,FALSE);

    for (i=0; i<Nunk_int_and_ext; i++)  {
      if (x[i] > ee) {
      x[i] += ee;  
      a[i] = calc_free_energy(NULL,x,1.0,1.0,FALSE);  
      if (i==33) printf("NODE 33: x:%g  a: %g \t",x[i],a[i]);
      x[i] -= 2*ee; 
      a[i] = (a[i] - calc_free_energy(NULL,x,1.0,1.0,FALSE))/(2.0*ee); 
      if (i==33) printf("NODE 33  x2:%g  a2:%g \n",x[i],calc_free_energy(NULL,x,1.0,1.0,FALSE));
      x[i] += ee;
      }
      else a[i] = 0.0;

      if (x2[i] > ee) {
      x2[i] += ee; 
      b[i] = calc_free_energy(NULL,x2,1.0,1.0,FALSE);  
      x2[i] -= 2*ee; 
      b[i] = (b[i] - calc_free_energy(NULL,x2,1.0,1.0,FALSE))/(2.0*ee); 
      x2[i] += ee;
      }
      else b[i] = 0.0;
    }

    Rho_b[0] += ee; thermodynamics("dft_out.lis", FALSE);
    g1 = calc_free_energy(NULL,x,1.0,1.0,FALSE);
    g2 = calc_free_energy(NULL,x2,1.0,1.0,FALSE);
    bt = Betamu[0];
    bh = Betamu_hs_ex[0] + Betamu_id[0];
    ba = Betamu_att[0];
    Rho_b[0] -= 2*ee; thermodynamics("dft_out.lis", FALSE);
    bt = (bt-Betamu[0])/ (2.0*ee);
    bh = (bh-Betamu_hs_ex[0]-Betamu_id[0])/ (2.0*ee);
    ba = (ba-Betamu_att[0])/ (2.0*ee);
    Rho_b[0] += ee; thermodynamics("dft_out.lis", FALSE);
    printf("Numerical dmu_drb: total %g  hs %g  att %g\n\n",bt, bh,ba);

    printf("Compare ads:   %g  %g  diff \n\n",ads, (g1-g2-g)/ee);
    printf("Compare:  f1a  f1n f1r    f2a  f2n  f2r  Vext\n");
    for (i=0; i<Nunk_int_and_ext; i++)  {
      printf("%d  %g  %g   %g     %g  %g   %g  %g  %g  %g  %g\n",i, f1[i],a[i],f1[i]/(a[i]+1.e-15),
                                         f2[i],b[i],f2[i]/(b[i]+1.e-15), Vext[i][0],
                           Nel_hit[0][Aztec.update_index[i]]*Vol_el/((double)Nnodes_per_el_V),
                           g,ads);
    }

    if (iter==2) exit(-1);
}*/
    safe_free((void *)&Nel_hit);
    safe_free((void *)&Nel_hit2);

    del_rb0 = -(g + dotp(f1,a) - dotp(f2,c)) / (ads + dotp(f1,b) - dotp(f2,d));

    if ((Imain_loop==0 && iter<3) /*|| fabs(g)<1.e-5*/) del_rb0=0.0;
    if ((iter<3) /*|| fabs(g)<1.e-5*/) del_rb0=0.0;

    if (Rho_b[0] + del_rb0 < 0.75*Rho_b[0]) del_rb0 = - 0.25*Rho_b[0];

    Rho_b[0] += del_rb0;

    thermodynamics("dft_out.lis", FALSE);

    for (i=0; i<Nunk_int_and_ext; i++){  
      delta_x[i] = c[i] + del_rb0*d[i];
    }

    converged2 = FALSE;
    converged2 = update_solution(iter,x2, delta_x);

    if (Imain_loop==0 && iter<1) converged2=FALSE;

    /*printf("\nBinodal norms: a=%g, b=%g, c=%g, d=%g\n",dotp(a,a), dotp(b,b), dotp(c,c), dotp(d,d));
    printf("\t Proc: %d  f1=%g, f2=%g, g=%g, ads=%g\n",Proc,dotp(f1,f1), dotp(f2,f2), g, ads);
    printf("\t del_rb0=%g, delta_x2=%g",del_rb0, dotp(delta_x,delta_x));
    printf("\t dmu_drb=%g, ",dmu_drb);*/

    for (i=0; i<Nunk_int_and_ext; i++){
      delta_x[i] = a[i] + del_rb0*b[i];
    }

    /*printf("delta_x1=%g\n\n",dotp(delta_x,delta_x));*/
    safe_free((void *) &a);
    safe_free((void *) &b);
    safe_free((void *) &c);
    safe_free((void *) &d);
    safe_free((void *) &f1);
    safe_free((void *) &f2);
    }

    /* Update solution vector and return convergence flag */

    converged = update_solution(iter,x, delta_x);
    MPI_Barrier(MPI_COMM_WORLD);

  /* END of Newton iteration loop */
 
  } while (iter < max_iter && (!converged || !converged2));

  safe_free((void *) &delta_x);
  safe_free((void *) &resid);
  safe_free((void *) &resid_tmp);

  *t_solve_max_ptr = t_solve_max;
  /* print out message on success or failure of Newton's nmethod */

  if (Proc==0 && Iwrite!=NO_SCREEN) {
    printf("\n=====================================");
    printf("=======================================\n");
    if (converged && converged2) 
      printf("\n### %s: Newton's method SUCCEEDED in %d iterations\n",yo,iter);
    else
      printf("\n### %s: Newton's method FAILED after %d iterations\n",yo,iter);
  }

  if (converged && converged2) return iter;
  else return -iter;
}

/****************************************************************************/
int update_solution(int iter,double *x, double *delta_x)
/*
 * This function updates the solution vector, x, with the delta_x,
 * the vector returned from the matrix solve. Any tampering with
 * Newton's method, such as damping or backtracking, can be
 * added to this routine.
 *
 * This function also decides whether or not Newton's method is 
 * converged and returns TRUE if it is, FALSE if it is not.
 */
{
  int loc_i,loc_inode,iunk,itmp,eq_type;
  double l2_norm = 0.0, temp;
  char *yo = "update solution";
  double frac,frac_min;


  frac_min = 1.0; 
  for (iunk=0; iunk<Nunk_per_node; iunk++){
     if ( (Type_poly>-1 && Unk_to_eq_type[iunk]!=POISSON) || 
          (Type_poly==-1 && Unk_to_eq_type[iunk]==DENSITY && Ipot_ff_n !=IDEAL_GAS)){

        for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {
           loc_i=Aztec.update_index[loc_find(iunk,loc_inode,LOCAL)];
           if (x[loc_i]+delta_x[loc_i]<0.0){
                frac = min(1.0,x[loc_i]/(-delta_x[loc_i]));
                frac = max(frac,Min_update_frac);
           }
           else frac = 1.0;
           if (Type_poly==-1 && Sten_Type[U_ATTRACT] && Unk_to_eq_type[iunk]==DENSITY
               && x[loc_i] <= Rho_b[iunk-Unk_start_eq[DENSITY]]*exp(-VEXT_MAX)) frac = 1.0;

           if (frac < frac_min) frac_min = frac;
         }
     }
  }
  frac_min = AZ_gmin_double(frac_min, Aztec.proc_config);
  if (Proc==0 && Iwrite != NO_SCREEN) 
      printf("\tUPDATE FRAC = %g percent\n",frac_min*100);

  if (frac_min < 1.0){
    for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {
      for (iunk=0; iunk<Nunk_per_node; iunk++){
          loc_i=Aztec.update_index[loc_find(iunk,loc_inode,LOCAL)];
          delta_x[loc_i] *= frac_min;
      }
    }
  }
  
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {
     for (iunk=0; iunk<Nunk_per_node; iunk++){
     eq_type=Unk_to_eq_type[iunk];

      loc_i=Aztec.update_index[loc_find(iunk,loc_inode,LOCAL)];


      if (Type_poly != -1){
          if (eq_type==DENSITY) itmp=iunk-Unk_start_eq[DENSITY];
          else if (eq_type==CMS_FIELD) itmp=iunk-Unk_start_eq[CMS_FIELD];
          else if (eq_type==CMS_G){
              itmp = Type_mer[Unk_to_Poly[iunk-Unk_start_eq[CMS_G]]][Unk_to_Seg[iunk-Unk_start_eq[CMS_G]]];
          }
      }
                                    /* TAKE CARE OF ZERO DENSITY NODES */
      if (Type_poly==-1 && eq_type==DENSITY && 
                   Zero_density_TF[L2B_node[loc_inode]][iunk-Unk_start_eq[DENSITY]]) x[loc_i] = 0.0; 

      else if (Type_poly>=0 && eq_type != POISSON && 
                  Zero_density_TF[L2B_node[loc_inode]][itmp]) x[loc_i] = 0.0; 

      else{

        temp     = delta_x[loc_i]/(Newton_rel_tol*x[loc_i] + Newton_abs_tol);
        l2_norm += temp*temp;
    
        /* Don't allow solutions to go negative */
            /* polymer B,rho, & G's always >= 0, but need to change this line
                when charges are added */

        if (( (Type_poly==-1 && eq_type==DENSITY) || (Type_poly>=0 && eq_type != POISSON) ) && 
               x[loc_i]+delta_x[loc_i] <= 1.e-15){ 
              x[loc_i] = x[loc_i]*0.1;  
        }
        else if(Matrix_fill_flag >=3 && Type_coul==-1 && iunk==Unk_start_eq[RHOBAR_ROSEN] 
                          && x[loc_i] >=1.0 && Ipot_ff_n !=IDEAL_GAS){
                                        x[loc_i] += 0.5*(1.0-x[loc_i]);
        }
        else {
               x[loc_i] += delta_x[loc_i]; 
        }
        /* check for density out of bounds */
        if (Type_poly > -1 && eq_type==DENSITY && x[loc_i]>Rho_max)
                        x[loc_i]=Rho_max;
      }
    }
  }

  AZ_exchange_bdry(x,Aztec.data_org,Aztec.proc_config);

  l2_norm = sqrt(AZ_gsum_double(l2_norm, Aztec.proc_config));

  if (Proc==0 && Iwrite != NO_SCREEN) printf("\t\t%s: Weighted norm of update vector =  %g\n", yo, l2_norm);

  if (l2_norm > 1.0) return(FALSE);
  else               return(TRUE);
}

/****************************************************************************/
static void set_multiple_B2L(int *ijk_box, int iunk, int loc_i)
{
  int i2, j2, k2, ijk_box_wrap[3], inode_box;

  for (i2 = -Non_unique_G2B[0]; i2 <= Non_unique_G2B[0]; i2++) {

    ijk_box_wrap[0] = ijk_box[0]+i2*Nodes_x[0];

    if (ijk_box_wrap[0] >=0 && ijk_box_wrap[0] < Nodes_x_box[0]) {

      for (j2 = -Non_unique_G2B[1]; j2 <= Non_unique_G2B[1]; j2++) {

        ijk_box_wrap[1] = ijk_box[1]+j2*Nodes_x[1];

        if (ijk_box_wrap[1] >=0 && ijk_box_wrap[1] < Nodes_x_box[1]) {

          for (k2 = -Non_unique_G2B[2]; k2 <= Non_unique_G2B[2]; k2++) {

            ijk_box_wrap[2] = ijk_box[2]+k2*Nodes_x[2];

            if (ijk_box_wrap[2] >=0 && ijk_box_wrap[2] < Nodes_x_box[2]) {

              inode_box = ijk_box_to_node_box(ijk_box_wrap);
              B2L_unknowns[loc_find(iunk,inode_box,BOX)] = loc_i;
            }
          }
        }
      }
    }
  }
}
/****************************************************************************/

void full_to_msr(double **a, int *bindx, double *val, int maxrow, int maxcol)
/* This routine translates from a full matrix to MSR format */

{
  int i,j,k;

  bindx[0] = maxrow+1;

  /* put diagonal entries in val array */
  for (i=0; i<maxrow; i++)
    val[i] = a[i][i];
      
  /* put off-diagonal entries in val and bindx arrays */
  k = maxrow+1;
    
  for (i=0; i<maxrow; i++) {
    for (j=0; j<maxcol; j++) {

       if ((a[i][j] != 0.0) && (i != j)) {
         bindx[k] = j;
         val[k] = a[i][j];
         k++;
       }
    }
    /* first maxrow+1 bindx entries give starting position for off-diags */
    bindx[i+1] = k;
  }
  
}
/****************************************************************************/
double dotp(double *x, double *y) 
{

  double sum=0.0;
  int i;

  for (i=0; i<Aztec.N_update; i++)   sum += x[i]*y[i];

  sum = AZ_gsum_double(sum, Aztec.proc_config);

  return(sum);

}
/****************************************************************************/
#ifdef NUMERICAL_JACOBIAN
static void do_numerical_jacobian(double *x, double *resid, double *resid_tmp)
/* This routine compares the analytic and numerical jacobians for the     */
/* purpose of checking the accuracy of an analytic Jacobian. It is only   */
/* written for serial runs and will be slow, so only for small problems.  */
/* Three files are created: ja, jn, jd  containing the analytic jacobian, */
/* numerical jacobian, and the difference between them, printed as the    */
/* row, column, and value.                                                */
{
  double **full=NULL, del;
  double **diff=NULL;
  double error;
  int i, j, N=Aztec.N_update,count=0;
  char filename[20];
  FILE *ifp;

  if (Iwrite != NO_SCREEN) printf("Proc: %d N=%d\n",Proc,N);

  /* print out analytic matrix in MSR format */
  sprintf(filename, "ja%0d",Proc);
  ifp = fopen(filename,"w");
  for (i=0; i< N; i++) {
    fprintf(ifp,"%d  %d   %g\n",i,i,Aztec.val[i]);
    for (j= Aztec.bindx[i]; j<Aztec.bindx[i+1]; j++) {
      fprintf(ifp,"%d  %d  %9.6f   %g\n",i,Aztec.bindx[j],
                Esize_x[0]*(Aztec.bindx[j]-i),Aztec.val[j]);
    }
  }
  fclose(ifp);

  /* compute numerical jacobian by repeated residual fill calls */

  full = (double **) array_alloc(2,N,N,sizeof(double));
  diff = (double **) array_alloc(2,N,N,sizeof(double));
  if (full==NULL){printf("Not enough memory for full numerical jacobian\n");
                  exit(-1);}
  for (i=0; i<N; i++) {
     del=1.e-6*fabs(x[i])+1.e-8;
     x[i] += del;
     for (j=0; j<N; j++) resid_tmp[j] = 0.0;
     fill_resid_and_matrix_control(x, resid_tmp, NULL, NULL, Matrix_fill_flag, 1, TRUE);

     for (j=0; j<N; j++){ 
        full[j][i] = (resid[j] + resid_tmp[j])/del;
     }
     x[i] -= del;
     if (count==100 || i==N-1){
       if (Iwrite != NO_SCREEN)printf("Proc: %d :: varying ith unknown=%d, out of %d unknowns.\n",Proc,i,N-1);
       count=0;
     } count++;
  }
  if (Iwrite != NO_SCREEN)printf("Proc: %d finished computing numerical jacobian !!\n",Proc);

  /* print out nonzero entries of numerical jacobian */
  sprintf(filename, "jn%0d",Proc);
  ifp = fopen(filename,"w");
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      if (fabs(full[i][j]) > 1.0e-10)
       fprintf(ifp,"%d  %d   %g\n",i,j,full[i][j]);
    }
  }
  fclose(ifp);

  /* compute and print the difference between the two jacobans */
  sprintf(filename, "jd%0d",Proc);
  ifp = fopen(filename,"w");
  for (i=0; i< N; i++) {
     diff[i][i] = full[i][i]-Aztec.val[i];
     for (j= Aztec.bindx[i]; j<Aztec.bindx[i+1]; j++) {
       diff[i][Aztec.bindx[j]] = full[i][Aztec.bindx[j]]-Aztec.val[j];
     }
   }

  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      if (fabs(full[i][j])>1.e-5){
         error = fabs(diff[i][j]/full[i][j])*100;
         if (error>1.)
          fprintf(ifp,"%d  %d   error=%g  diff=%g  numeric=%g\n",i,j,error,diff[i][j],full[i][j]);
      }
      else if (fabs(diff[i][j])>1.e-6)
          fprintf(ifp,"%d  %d   diff=%g  numeric=%g\n",i,j,diff[i][j],full[i][j]);
    }
  }
  fclose(ifp);
  safe_free((void *) &full);
  if (Num_Proc>1) MPI_Barrier(MPI_COMM_WORLD);
    if (Iwrite != NO_SCREEN)printf("Proc: %d KILLING CODE AT END OF NUMERICAL JACOBIAN PRINT\n",Proc);
  exit(-1);
}
#endif
/****************************************************************************/
