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
 *      All Global Definitions occur in the include files here.
 */
#include <mpi.h>
#include <unistd.h>
#include "dft_globals.h"
#include "rf_allo.h"
/*****************************************************************************/
void continuation_shift();
void setup_nunk_per_node(char *);
void dftmain(double *);

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
  char     *output_file1, *output_file2, *output_file3, *output_file4;
  char     *output_TF, *yo = "main";
  int       iend, match, idim, icomp, i, niters;
  double    Esize_x_final[3];
  double    *x,*x2=NULL;  /* The solution vector of densities */
  double    time_save;
  double    t_preprocess=0.0,t_solve=0.0,t_postprocess=0.0,t_total=0.0;
  double    t_pre_max,t_solve_max,t_post_max,t_total_max,t_msr_max;
  double    t_pre_min,t_solve_min,t_post_min,t_total_min,t_msr_min;
  int izone,isten,jcomp,jmax;
  struct Stencil_Struct *sten;
  char line[100],linecwd[100];
  int *dummy,proper_bc;
  int argc=1;

  AZ_set_proc_config(Aztec.proc_config,MPI_COMM_WORLD);
  gethostname(line,100);
  getcwd(linecwd,100);
  /*  printf(" hello %s %s \n",line,linecwd); */
  MPI_Barrier(MPI_COMM_WORLD);

  /*  AZ_processor_info(Aztec.proc_config);
  Proc     = Aztec.proc_config[AZ_node];
  Num_Proc = Aztec.proc_config[AZ_N_procs];
  */
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
  output_file3 = "dft_output.dat";
  output_file4 = "dft_vext_c.dat";
  /* ALF */
  output_TF = "dft_zeroTF.dat";

  /*dummy = (int *) array_alloc (1, 500000000, sizeof(int));
  safe_free((void *) &dummy);*/

 /*
  * Read in the ascii input file from the front end 
  * This routine also sets up some flags
  * (see file dft_input.c)
  */

  read_input_file(input_file,output_file1);
  setup_nunk_per_node(output_file1);

  if (Mix_type == 0) pot_parameters(output_file1);
  if (Sten_Type[POLYMER_CR]){
      Rism_cr = (double ***) array_alloc (3, Ncomp,Ncomp,N_NZCR_MAX,sizeof(double));
      setup_polymer_cr();
  }

 /* count_zero = 0;
  count_nonzero = 0;
  */

  if (Loca.method != -1) {
    Nruns = 1;
    if (Proc==0 && Iwrite!=NO_SCREEN) printf("\nWARNING: Loca library requested, disabling"
                        "Tramonto continuation loops\n");
  }

  for (Imain_loop=0; Imain_loop<Nruns; Imain_loop++){

   if (Imain_loop > 0)  continuation_shift(); /* only mesh continuation is now possible outside of loca */ 

     for (idim=0; idim<Ndim; idim++){
        Esize_x_final[idim] = Esize_x[idim];
        if (Lmesh_refine==TRUE) Esize_x[idim] = 0.2;
     }
     iend=0;

     while (iend==0) {

     t_mesh = MPI_Wtime();
    /*
     * Stencils for non-local interactions are assembled next, and are
     * dependent on choice for integration schemes
     * (see file dft_stencil.c)
     */
     if (Imain_loop==0) calc_stencils();
    /*
     * do all the thermodynamics for the bulk fluid mixture
     */
     thermodynamics(output_file1, TRUE);

    /*
     * Set up the mesh based on input info, and allocate the solution vector
     * (see file dft_mesh.c)
     */


     if (Imain_loop > 0){
           free_mesh_arrays();
           boundary_free();
           if (Lsteady_state && Ndim==1) safe_free((void *) &Area_IC);
           safe_free((void *) &Comm_node_proc);
           safe_free((void *) &Comm_unk_proc);
           safe_free((void *) &Comm_offset_node);
           safe_free((void *) &Comm_offset_unk);
           safe_free((void *) &Aztec.update);

     }
     set_up_mesh(output_file1,output_file2);

     /* call load balancer if requested */

     /*if (Load_Bal_Flag == LB_TIMINGS) {
       load_balance(1, NULL);
     }
     */

    /*
     * Set up boundary and surface charge information -- this must come
     * after load-balancing because quantities are local on a Proc
     */
     boundary_setup(output_file1);

     /* check boudary conditions to see Coulomb external field can be computed
        for this run ..... we can't do it with periodic, multiply reflecting, or
        in wall boundary conditions due to need for infinite ewald sums */
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
        print_zeroTF(Zero_density_TF,output_TF);
        if (Vol_charge_flag && Ndim==3) print_vext(Vext_coul,output_file4);
     }


     /*
      * Use Newton's method to solve the problem
      * (see file dft_newton.c, which calls dft_fill.c and Aztec library)
      * The variable "niters" is positive if Newton's method converrged.
      */
      x = (double *) array_alloc (1, Nnodes_per_proc*Nunk_per_node, sizeof(double));
      if (Lbinodal) x2 = (double *) array_alloc (1, Nnodes_per_proc*Nunk_per_node, sizeof(double));


      t_preprocess += MPI_Wtime();
      t_solve = -MPI_Wtime();
      T_av_precalc_min = T_av_fill_min = T_av_solve_min = T_av_lj_min = 0.0;
      T_av_precalc_max = T_av_fill_max = T_av_solve_max = T_av_lj_max = 0.0;
      niters = solve_problem(&x, &x2);
      t_solve += MPI_Wtime();
     /*
      * Post-Process the results
      * (see file dft_output.c)
      */
      time_save =  MPI_Wtime() - t_mesh;
      t_postprocess = -MPI_Wtime();

      if (Loca.method == -1) {
        if (Lbinodal) post_process(x2, output_file3, &niters, &time_save,
                                   Imain_loop, TRUE);
        if (Proc==0 && Lbinodal){
          X2_old = (double *) array_alloc(1, Nodes_old*Nunk_per_node, sizeof(double));
          for (i=0; i<Nodes_old*Nunk_per_node; i++) X2_old[i] = X_old[i];
        }
        post_process(x, output_file3, &niters, &time_save, Imain_loop, FALSE);
      }

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
      t_msr_max= gmax_double(T_msr_setup);
      t_solve_max=gmax_double(t_solve);
      t_post_max=gmax_double(t_postprocess);
      t_total_max=gmax_double(t_total);

      t_pre_min= gmin_double(t_preprocess);
      t_msr_min= gmin_double(T_msr_setup);
      t_solve_min=gmin_double(t_solve);
      t_post_min=gmin_double(t_postprocess);
      t_total_min=gmin_double(t_total);

      if (Proc == 0 &&Iwrite !=NO_SCREEN) {
        printf ("\n\n\n\n");
        printf ("===================================================\n");
        printf ("TIMING SUMMARY .... averaged over all %d iterations\n",niters);
        printf ("===================================================\n");
        printf ("                         min(sec)        max(sec)  \n");
        printf ("---------------------------------------------------\n");
        printf ("MESH SETUP               %g         %g       \n",
                                                    t_pre_max,t_pre_min); 
        printf ("---------------------------------------------------\n");
        printf ("MSR_PREPROC              %g         %g       \n",
                                                    t_msr_min,t_msr_max);
        printf ("PRE_CALC ROUTINES        %g        %g       \n",
                        T_av_precalc_min/niters,T_av_precalc_max/niters);
        if (Matrix_fill_flag <3)
        printf ("JACOBIAN                 %g         %g       \n",
                                  T_av_lj_min/niters,T_av_lj_max/niters);
        printf ("FILL ROUTINE             %g         %g       \n",
                              T_av_fill_min/niters,T_av_fill_max/niters);
        printf ("---------------------------------------------------\n");
        printf ("AZTEC SOLVE              %g          %g       \n",
                            T_av_solve_min/niters,T_av_solve_max/niters);
        printf ("---------------------------------------------------\n");
        printf ("TOTAL SOLVE              %g         %g       \n",
                                                t_solve_min,t_solve_max);
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
      if ((Ipot_ff_n == LJ12_6) && !Sten_Type[POLYMER_CR]) {
        safe_free((void *) &Avdw);
        safe_free((void *) &Betamu_att);
      }
      if (Ipot_ff_c == COULOMB) 
           safe_free((void *) &Deltac_b);

      /* solution arrays*/
      safe_free((void *) &B2L_unknowns);
      safe_free((void *) &x); 
      if (Lbinodal) safe_free((void *) &x2); 

      /* Aztec arrays -- mostly allocated by AZ_transform */
      safe_free((void *) &Aztec.external);
      safe_free((void *) &Aztec.extern_index);
      safe_free((void *) &Aztec.update_index);
      /*AZ_free_memory(Aztec.data_org[AZ_name]);*/
      safe_free((void *) &Aztec.data_org);

  } /* end of loop over continuation field #1 */

  if (Nwall_type !=0){
    safe_free((void *) &Link_list);
    safe_free((void *) &Nwall_this_link);
  }
  free_mesh_arrays();
  boundary_free();
  if (Lsteady_state && Ndim==1) safe_free((void *) &Area_IC);
  safe_free((void *) &Comm_node_proc);
  safe_free((void *) &Comm_unk_proc);
  safe_free((void *) &Comm_offset_node);
  safe_free((void *) &Comm_offset_unk);
  safe_free((void *) &Aztec.update);

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
  /* ALF: fixed for changed stencils */
  for (izone=0; izone<Nzone; izone++){
  for (isten=0; isten<NSTEN; isten++)
   if (Sten_Type[isten]==1) {
     if (isten == U_ATTRACT || isten == THETA_CHARGE || isten == POLYMER_CR) jmax = Ncomp;
     else if (isten == DELTA_FN && Sten_Type[POLYMER_CR]) jmax = Ncomp;
     else  jmax = 1;

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

  if (Sten_Type[POLYMER_CR]) safe_free((void **) &Rism_cr);


  if (Proc == 0)
    safe_free((void *) &X_old); /* This is also done for each continuation run
              in set_initial_guess (dft_guess.c) called from solve_problem.
              The array is allocated in collect_x_old (dft_output.c) called
              from post_process.    */
  *engptr = Energy;
  return;
} 
/************END of main() *************************************************/


/****************************************************************************/
/*continutation_shift:  In this routine we set up parameters for a continuation
                        where the walls have been shifted apart.*/
void continuation_shift()
{
  int idim,iwall,icharge;

  for (idim=0; idim<Ndim; idim++) {
      Size_x[idim] += Del_1[idim];

      if (Plane_new_nodes == idim){
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
/*setup_polymer_cr: read in c(r) from file and add attractions
  ALF: modified 10/21 for Type_poly == 3 */
void setup_polymer_cr()
{
   FILE *fp7,*fp8,*fp9,*fp10;
   int i,j,lines,lines2,ir,ir2, nextChar;
   double r,u,cr_rad_max,cr_rad_max2,rsave,rsave2,dummy_read,crread;
   double crfac1,crfac2,crfac3,crfac4,xs;
   char c;

   /* not doing automated interpolation - just mixing two c(r) files */
/*   crfac1=Crfac;
   crfac2=(1.0-crfac1);
   crfac3=0.0; crfac4=0.0;*/

/*   if (fabs(Crfac+1.0) < 1.e-8) { *//* do automated interpolation */
      /* note that this particular interpolation is based on xs varying between 0 and 1....
         to interpolate other parameters with other limits will require a different implementation */
      xs=Rho_b[2]/(Rho_b[0]+Rho_b[1]+Rho_b[2]); 
      if (Ncr_files == 1){
             crfac1=Crfac;
      }
      else if (Ncr_files == 2){
             crfac1=xs;
             crfac2=1.0-xs;
      }
      else if (Ncr_files ==3){
         if (xs>0 && xs<Cr_break[0]){
            crfac1=xs/(Cr_break[0]);
            crfac2=1.0-crfac1;
            crfac3=0.0;
         }
         else{ 
           crfac1=0.0;
           crfac2=(xs-Cr_break[0])/(1.0-Cr_break[0]);
           crfac3=1.0-crfac2;
         }
      }
      else if (Ncr_files ==4){
         if (xs>0 && xs<Cr_break[0]){
            crfac1=1.0-xs/(Cr_break[0]);
            crfac2=1.0-crfac1;
            crfac3=0.0; crfac4=0.0;
         }
         else if (xs>=Cr_break[0] && xs<Cr_break[1]){
            crfac1=0.0; crfac4=0.0;
            crfac2 = 1.0-(xs-Cr_break[0])/(Cr_break[1]-Cr_break[0]);
            crfac3=1.0-crfac2;
         }
         else if (xs>=Cr_break[1]){
            crfac1=0.0; crfac2=0.0;
            crfac3 = 1.0 - (xs-Cr_break[1])/(1.0-Cr_break[1]);
            crfac4=1.0-crfac3;
         }
      }
  /* }*/
   if (Ncr_files == 1) printf("crfac1=%9.6f  ",crfac1);
   if (Ncr_files == 2) printf("crfac2=%9.6f  ",crfac2);
   if (Ncr_files == 3) printf("crfac3=%9.6f  ",crfac3);
   if (Ncr_files == 4) printf("crfac4=%9.6f  ",crfac4);
   printf("\n");

   /* reading in c(r) file */
   if(Proc==0) printf("reading in %d c(r) file(s)...\n",Ncr_files);
   if (Type_poly != 3) {

   if (Proc==0){
     if( (fp7  = fopen(Cr_file,"r")) == NULL) {
       printf("Can't open file %s\n", Cr_file);
       exit(1);
     }
     if (Ncr_files>=2){
       if( (fp8  = fopen(Cr_file2,"r")) == NULL) {
          printf("Can't open file %s\n", Cr_file2);
          exit(1);
        }
        fclose(fp8);
     }
     if (Ncr_files>=3){
       if( (fp9  = fopen(Cr_file3,"r")) == NULL) {
          printf("Can't open file %s\n", Cr_file3);
          exit(1);
        }
        fclose(fp9);
     }
     if (Ncr_files==4){
       if( (fp10  = fopen(Cr_file4,"r")) == NULL) {
          printf("Can't open file %s\n", Cr_file4);
          exit(1);
        }
        fclose(fp10);
     }
     for (ir=1; ir<=3; ir++){
       fscanf(fp7,"%lf",&r );
       fscanf(fp7,"%c",&c );
       for (i=0; i<Ncomp; i++)  /* for (i=0; i<Ntype_mer; i++)  */
	 fscanf(fp7,"%lf",&dummy_read);
       fscanf(fp7,"%c",&c );
       for (i=0; i<Ncomp; i++) {                /* "  */
	 for (j=i+1; j<Ncomp; j++) {                /* "  */
	   fscanf(fp7,"%lf",&dummy_read);
	   fscanf(fp7,"%c",&c );
         }
       }
       if (ir == 1) rsave = r;
       if (ir == 2) Deltar_cr = r-rsave;
     }
     fclose(fp7);
   }
   MPI_Bcast(&Deltar_cr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   lines = 0;
   cr_rad_max=0.0;
   for (i=0; i < Ncomp; ++i)
      for (j=0; j<Ncomp; ++j){
	 if (Cr_rad_hs[i][j] > cr_rad_max) cr_rad_max = Cr_rad_hs[i][j];
      }
   lines=(int)(cr_rad_max/Deltar_cr);

   if (Proc==0 && lines >= N_NZCR_MAX) {
      printf("Warning: Rism_cr array may be truncated :: Max Cr_rad_hs > allowed\n");
      lines = N_NZCR_MAX-1;
    }
    /*printf("Proc=%d lines=%d\n",Proc,lines); */

   Last_nz_cr = lines;
   if (Proc==0) {
     fp7  = fopen(Cr_file,"r");
     if (Ncr_files>=2) fp8  = fopen(Cr_file2,"r");
     if (Ncr_files>=3) fp9  = fopen(Cr_file3,"r");
     if (Ncr_files==4) fp10  = fopen(Cr_file4,"r");
    
     for (ir=1; ir<=lines; ir++){
       fscanf(fp7,"%lf",&r );
       fscanf(fp7,"%c",&c );
       if (ir == 1) rsave = r;
       if (ir == 2) Deltar_cr = r-rsave;

       if (Ncr_files>=2){
         fscanf(fp8,"%lf",&r );
         fscanf(fp8,"%c",&c );
       }
       if (Ncr_files>=3){
         fscanf(fp9,"%lf",&r );
         fscanf(fp9,"%c",&c );
       }
       if (Ncr_files==4){
         fscanf(fp10,"%lf",&r );
         fscanf(fp10,"%c",&c );
       }

       for (i=0; i<Ncomp; i++){  /* for (i=0; i<Ntype_mer; i++)  */
	 fscanf(fp7,"%lf",&crread);
         Rism_cr[i][i][ir]=Crfac*crfac1*crread;
	 if (Ncr_files>=2){
              fscanf(fp8,"%lf",&crread);
              Rism_cr[i][i][ir]+=Crfac*crfac2*crread;
         }
	 if (Ncr_files>=3){
              fscanf(fp9,"%lf",&crread);
              Rism_cr[i][i][ir]+=Crfac*crfac3*crread;
         }
	 if (Ncr_files==4){
              fscanf(fp10,"%lf",&crread);
              Rism_cr[i][i][ir]+=Crfac*crfac4*crread;
         }
       }
       fscanf(fp7,"%c",&c );
       if (Ncr_files>=2) fscanf(fp8,"%c",&c );
       if (Ncr_files>=3) fscanf(fp9,"%c",&c );
       if (Ncr_files==4) fscanf(fp10,"%c",&c );

       for (i=0; i<Ncomp; i++) {                
	 for (j=i+1; j<Ncomp; j++) {                
	   fscanf(fp7,"%lf",&crread);
           Rism_cr[i][j][ir]=Crfac*crfac1*crread;
	   fscanf(fp7,"%c",&c );
	   if (Ncr_files>=2){
                 fscanf(fp8,"%lf",&crread);
                 Rism_cr[i][j][ir]+=Crfac*crfac2*crread;
	         fscanf(fp8,"%c",&c );
           }
	   if (Ncr_files>=3){
                 fscanf(fp9,"%lf",&crread);
                 Rism_cr[i][j][ir]+=Crfac*crfac3*crread;
	         fscanf(fp9,"%c",&c );
           }
	   if (Ncr_files==4){
                 fscanf(fp10,"%lf",&crread);
                 Rism_cr[i][j][ir]+=Crfac*crfac4*crread;
	         fscanf(fp10,"%c",&c );
           }
	   Rism_cr[j][i][ir] = Rism_cr[i][j][ir];
	 }
       }
       while(c != '\n') c=getc(fp7);
       if(Ncr_files>=2) while(c != '\n') c=getc(fp8);
       if(Ncr_files>=3) while(c != '\n') c=getc(fp9);
       if(Ncr_files==4) while(c != '\n') c=getc(fp10);
     }
     fclose(fp7);
     if(Ncr_files>=2) fclose(fp8);
     if(Ncr_files>=3) fclose(fp9);
     if(Ncr_files==4) fclose(fp10);

     for (ir=lines+1; ir<N_NZCR_MAX; ir++)
       for (i=0; i<Ncomp; i++){
	 Rism_cr[i][i][ir] = 0.;
	 for (j=i+1; j<Ncomp; j++) {     
	   Rism_cr[i][j][ir] = 0.;
	   Rism_cr[j][i][ir] = 0.;
	 }
       }
     /* extrapolate c(r) to r = 0 */
     for (i=0; i<Ncomp; i++)                 /* "  */
       for (j=0; j<Ncomp; j++)                 /* "  */
	 Rism_cr[i][j][0] = 2.*Rism_cr[i][j][1] - Rism_cr[i][j][2];
   }
   MPI_Bcast(**Rism_cr,Ncomp*Ncomp*N_NZCR_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

  /* add attractions to polymer c(r)  */

   if (Iwrite == VERBOSE && Proc==0){
      fp7  = fopen("cr.out","w");
      for (ir=0; ir <= Last_nz_cr; ++ir){
         fprintf(fp7,"%lf   ",ir*Deltar_cr);
	 for (i=0; i < Ncomp; ++i)
	    for (j=0; j < Ncomp; ++j)
	       fprintf(fp7,"%lf   ",Rism_cr[i][j][ir]);
	 fprintf(fp7,"\n");
      }
      fclose(fp7);
   }

   if (Ipot_ff_n == LJ12_6){
     for (i=0; i < Ncomp; ++i)
       for (j=0; j<Ncomp; ++j){
	 if (Cut_ff[i][j] > Cr_rad_hs[i][j]) Cr_rad[i][j] = Cut_ff[i][j];
         else                                Cr_rad[i][j] = Cr_rad_hs[i][j];
         if (Cr_rad[i][j] > cr_rad_max) cr_rad_max=Cr_rad[i][j];
       }
     
   lines = cr_rad_max/Deltar_cr;
   if (lines >= N_NZCR_MAX) {
     if (Proc==0) printf("Need to increase N_NZCR_MAX\n");
     lines = N_NZCR_MAX-1;
   }
   if (lines > Last_nz_cr) Last_nz_cr = lines;

     for (i=0; i < Ncomp; ++i)
       for (j=i; j<Ncomp; j++) { 
	 for (ir=0; ir<=lines; ir++){
	   r = ir*Deltar_cr;
	   if (r-1.e-8 > Sigma_ff[i][j]*1.122462) { /* watch roundoffs*/
	     u = uLJatt_n(r,i,j);
	     Rism_cr[i][j][ir] -= u;
	     if (i != j) Rism_cr[j][i][ir] -= u;
	   }
	 }
       }
   }
   if (Iwrite == VERBOSE && Proc==0){
      fp7  = fopen("cr.lj.out","w");
      for (ir=0; ir <= Last_nz_cr; ++ir){
         fprintf(fp7,"%lf   ",ir*Deltar_cr);
	 for (i=0; i < Ncomp; ++i)
	    for (j=0; j < Ncomp; ++j)
	       fprintf(fp7,"%lf   ",Rism_cr[i][j][ir]);
	 fprintf(fp7,"\n");
      }
      fclose(fp7);
   }

   } /* end of if (Type_poly != 3) */

   /* for Type_poly == 3, just read in c(0) */
   else {
     if (Proc==0){
       if( (fp7  = fopen(Cr_file,"r")) == NULL) {
	 printf("Can't open file %s\n", Cr_file);
	 exit(1);
       }
 
       for (i=0; i<Ncomp; i++) {  /* for (i=0; i<Ntype_mer; i++)  */
	 fscanf(fp7,"%lf",&Rism_cr[i][i][0]);
	 fscanf(fp7,"%c",&c );
	 printf("cr[%d][%d] = %f\t", i, i, Rism_cr[i][i][0]);
       }
       for (i=0; i<Ncomp; i++) {                /* "  */
	 for (j=i+1; j<Ncomp; j++) {                /* "  */
	   fscanf(fp7,"%lf",&Rism_cr[i][j][0]);
	   fscanf(fp7,"%c",&c );
	   Rism_cr[j][i][0] = Rism_cr[i][j][0];
	   printf("cr[%d][%d] = %f\t", i,j,Rism_cr[i][j][0]);
	 }
       }
       fclose(fp7);
       for (ir=1; ir<N_NZCR_MAX; ir++){
	 for (i=0; i<Ncomp; i++){
	   Rism_cr[i][i][ir] = 0.;
	   for (j=i+1; j<Ncomp; j++) {     
	     Rism_cr[i][j][ir] = 0.;
	     Rism_cr[j][i][ir] = 0.;
	   }
	 }
       }
     } /* end of if(Proc==0) */

     MPI_Bcast(**Rism_cr,Ncomp*Ncomp*N_NZCR_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
   } /* end of Type_poly==3 */

   return;
}
/*************************************************************/
/*******************************************************************************/
/*setup_nunk_per_node:  here we just set up the basic parameters for the number
  of unknowns per node for the run of interest, and some arrays to move between
  unknown number and equation type easily. */
void setup_nunk_per_node(char *output_file1)
{
  int i,iunk,icomp;
  int NCMSField_unk;
  FILE *fp2=NULL;

  if( (fp2 = fopen(output_file1,"w+")) == NULL) {
      printf("Can't open file %s\n", output_file1);
      exit(1);
   }

/*in Makefile set switch for a polymer run an set unknowns accordingly*/
   
   for (i=0;i<NEQ_TYPE;i++){
     switch(i){
         case DENSITY:                  /* unknowns of Euler-Lagrange equation */
            Phys2Nunk[DENSITY]=Ncomp;
            break;

         case RHOBAR_ROSEN:     /* unknowns of Nonlocal Density Eqns for Rosenfeld Functionals */
            Nrho_bar = 0;
            if (Type_poly == -1 && Matrix_fill_flag >= 3 && Ipot_ff_n != IDEAL_GAS){
                 if (Matrix_fill_flag ==3) Nrho_bar = 4 + 2*Ndim;
                 else Nrho_bar = 4;
                 Nrho_bar_s = 4;
            }
            Phys2Nunk[RHOBAR_ROSEN]=Nrho_bar;
            break; 

         case POISSON:                             /* unknowns of Poisson Equation */
            Npoisson=0;
            if ( Type_coul !=NONE){
                 Npoisson=1;
            }
            Phys2Nunk[POISSON]=Npoisson;
            break;

         case DIFFUSION:                            /* unknowns of Diffusion Equation */
            Ndiffusion=0;
            if (Type_poly==-1 && Lsteady_state){
              if (!Type_poly_TC) Ndiffusion=Ncomp;
              if (Type_poly_TC) Ndiffusion=Nseg_tot;
            }
            Phys2Nunk[DIFFUSION]=Ndiffusion;
            break;

         case DENSITY_SEG:
            Ntype_unk=0;
            if (Type_poly_TC){
               Ntype_unk=Ncomp;
            }
            Phys2Nunk[DENSITY_SEG]=Ntype_unk;
            break;
             
         case CAVITY_WTC:
            Nrho_bar_cavity=0;
            if (Type_poly_TC){
                 Nrho_bar_cavity = 4;
                 Phys2Nunk[CAVITY_WTC]=Nrho_bar_cavity-2;  
                                   /* strange case because y function only uses two of the 
                                      defined xi's.  But, I'm leaving a placeholder for 4 */
            }
            else Phys2Nunk[CAVITY_WTC]=0;
            break;

         case BOND_WTC:
            Nrho_bar_bond=0;
            if (Type_poly_TC){
                 Nrho_bar_bond = Nbonds;
            }
            Phys2Nunk[BOND_WTC]=Nrho_bar_bond;
            break;

         case CMS_FIELD:
            NCMSField_unk = 0;
            if (Type_poly != NONE){
                 NCMSField_unk=Ncomp;
            }
            Phys2Nunk[CMS_FIELD]=NCMSField_unk;
            break;

         case CMS_G:
            if (Type_poly != NONE){
                 Phys2Nunk[CMS_G] = Ngeqn_tot;
            }
            break;

         default:
            printf("problems with defining equation type %d\n",i);
            exit(-1);
            break;
     }
     if (i==0) Phys2Unk_first[i]=0;
     else      Phys2Unk_first[i]=Phys2Unk_last[i-1];
     Phys2Unk_last[i]=Phys2Unk_first[i]+(Phys2Nunk[i]);
     for (iunk=Phys2Unk_first[i]; iunk< Phys2Unk_last[i]; iunk++) Unk2Phys[iunk]=i;
     Nunk_per_node += Phys2Nunk[i];
   }
   for (i=0;i<NEQ_TYPE;i++){
     if (Phys2Nunk[i]==0){
        Phys2Unk_first[i]=NO_UNK;
        Phys2Unk_last[i]=NO_UNK;
     }
     if (i==CMS_G && Type_poly != NONE){
        for (icomp=0;icomp<Npol_comp;icomp++) Geqn_start[icomp]+=Phys2Unk_first[i];
     }
   }
   if (Iwrite==VERBOSE){
        printf("\n******************************************************\n");
        printf("TOTAL Nunk_per_node=%d\n",Nunk_per_node);
        for (i=0;i<NEQ_TYPE;i++) printf("Phys2Nunk[%d]=%d  start_unk=%d  end_unk=%d\n",
                                   i,Phys2Nunk[i],Phys2Unk_first[i],Phys2Unk_last[i]);
        for (iunk=0;iunk<Nunk_per_node;iunk++) printf("iunk=%d equation_type=%d\n",iunk,Unk2Phys[iunk]);
        printf("******************************************************\n");
   }
   fprintf(fp2,"\n******************************************************\n");
   fprintf(fp2,"TOTAL Nunk_per_node=%d\n",Nunk_per_node);
   for (i=0;i<NEQ_TYPE;i++) fprintf(fp2,"Phys2Nunk[%d]=%d  start_unk=%d  end_unk=%d\n",
                                   i,Phys2Nunk[i],Phys2Unk_first[i],Phys2Unk_last[i]);
   for (iunk=0;iunk<Nunk_per_node;iunk++) fprintf(fp2,"iunk=%d equation_type=%d\n",iunk,Unk2Phys[iunk]);
   fprintf(fp2,"******************************************************\n");

   return;
}
/*******************************************************************************/
