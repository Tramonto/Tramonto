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
/* if dftnormal is 0 then compile for a library */
/* if dftnormal is 1 then compile as stand alone code */
#define DFTNORMAL 1
/*****************************************************************************/
void continuation_shift ();
void dftmain (double *);

#if DFTNORMAL
int
main (int argc, char *argv[])
{
  int i;
  double dumb;
  MPI_Init (&argc, &argv);
/*  for (i=0; i<1000; i++)*/ dftmain (&dumb);
  MPI_Finalize ();
  return (1);
}
#endif

void
dftmain (double *engptr)
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

  double start_t, t_mesh;
  char *input_file;
  char *output_file1, *output_file2, *output_file3, *output_file4;
  char *output_TF, *yo = "main";
  int successful;
  int iend, match, idim, icomp, i, niters;
  double Esize_x_final[3];
  double *x, *x2 = NULL;	/* The solution vector of densities */
  double time_save;
  double t_preprocess = 0.0, t_solve = 0.0, t_postprocess = 0.0, t_total =
    0.0;
  double t_pre_max, t_solve_max, t_post_max, t_total_max, t_msr_max;
  double t_pre_min, t_solve_min, t_post_min, t_total_min, t_msr_min;
  int izone, isten, jcomp, jmax;
  struct Stencil_Struct *sten;
  char line[100], linecwd[100];
  int *dummy, proper_bc;
  int argc = 1;

  AZ_set_proc_config (Aztec.proc_config, MPI_COMM_WORLD);
  gethostname (line, 100);
  getcwd (linecwd, 100);
  /*  printf(" hello %s %s \n",line,linecwd); */
  MPI_Barrier (MPI_COMM_WORLD);

  /*  AZ_processor_info(Aztec.proc_config); */
  Proc = Aztec.proc_config[AZ_node];
  Num_Proc = Aztec.proc_config[AZ_N_procs];

  start_t = MPI_Wtime ();
  t_preprocess = -MPI_Wtime ();
  /*
   * Interpret command line arguments, get input file name or use default
   */

  switch (argc)
    {

    case 2:
      /*  input_file = argv[1]; */
      break;

    case 1:
      input_file = "dft_input.dat";
      break;

    default:
      printf ("%s: ERROR in command line usuage:\n", yo);
      printf ("\tExecutable can have 0 or 1 argument\n");
      exit (-1);
      break;
    }
  output_file1 = "dft_out.lis";
  output_file2 = "dft_vext.dat";
  output_file3 = "dft_output.dat";
  output_file4 = "dft_vext_c.dat";
  /* ALF */
  output_TF = "dft_zeroTF.dat";

  /*dummy = (int *) array_alloc (1, 500000000, sizeof(int));
     safe_free((void *) &dummy); */

  /*
   * Read in the ascii input file from the front end 
   * This routine also sets up some flags
   * (see file dft_input.c)
   */

  read_input_file (input_file, output_file1);
  if (Mix_type == 0)
    pot_parameters (output_file1);
  if (Sten_Type[POLYMER_CR])
    {
      Rism_cr =
	(double ***) array_alloc (3, Ncomp, Ncomp, N_NZCR_MAX,
				  sizeof (double));
      setup_polymer_cr ();
    }

  /* count_zero = 0;
     count_nonzero = 0;
   */

  if (Loca.method != -1)
    {
      Nruns = 1;
      if (Proc == 0 && Iwrite != NO_SCREEN)
	printf ("\nWARNING: Loca library requested, disabling"
		"Tramonto continuation loops\n");
    }

  for (Imain_loop = 0; Imain_loop < Nruns; Imain_loop++)
    {

      if (Imain_loop > 0)
	continuation_shift ();	/* only mesh continuation is now possible outside of loca */

      for (idim = 0; idim < Ndim; idim++)
	{
	  Esize_x_final[idim] = Esize_x[idim];
	  if (Lmesh_refine == TRUE)
	    Esize_x[idim] = 0.2;
	}
      iend = 0;

      while (iend == 0)
	{

	  t_mesh = MPI_Wtime ();
	  /*
	   * Stencils for non-local interactions are assembled next, and are
	   * dependent on choice for integration schemes
	   * (see file dft_stencil.c)
	   */
	  if (Imain_loop == 0)
	    calc_stencils ();
	  /*
	   * do all the thermodynamics for the bulk fluid mixture
	   */

	  thermodynamics (output_file1, TRUE);

	  /*
	   * Set up the mesh based on input info, and allocate the solution vector
	   * (see file dft_mesh.c)
	   */


	  if (Imain_loop > 0)
	    {
	      free_mesh_arrays ();
	      boundary_free ();
	      if (Lsteady_state && Ndim == 1)
		safe_free ((void *) &Area_IC);
	      safe_free ((void *) &Comm_node_proc);
	      safe_free ((void *) &Comm_unk_proc);
	      safe_free ((void *) &Comm_offset_node);
	      safe_free ((void *) &Comm_offset_unk);
	      safe_free ((void *) &Aztec.update);

	    }
	  set_up_mesh (output_file1, output_file2);

	  /* call load balancer if requested */

	  /*if (Load_Bal_Flag == LB_TIMINGS) {
	     load_balance(1, NULL);
	     }
	   */

	  /*
	   * Set up boundary and surface charge information -- this must come
	   * after load-balancing because quantities are local on a Proc
	   */
	  printf ("calling boundary setup\n");
	  boundary_setup (output_file1);
	  printf ("after boundary setup\n");

	  /* check boudary conditions to see Coulomb external field can be computed
	     for this run ..... we can't do it with periodic, multiply reflecting, or
	     in wall boundary conditions due to need for infinite ewald sums */
	  if (Vol_charge_flag && Ndim == 3)
	    {
	      proper_bc = TRUE;
	      for (idim = 0; idim < Ndim; idim++)
		{
		  if (Type_bc[idim][0] == PERIODIC ||
		      (Type_bc[idim][0] == REFLECT
		       && Type_bc[idim][1] != IN_BULK)
		      || (Type_bc[idim][1] == REFLECT
			  && Type_bc[idim][0] != IN_BULK)
		      || (Type_bc[idim][0] == IN_WALL
			  && Type_bc[idim][0] != IN_WALL))
		    proper_bc = FALSE;
		}
	      if (proper_bc)
		setup_vext_coulomb_vol ();
	      else
		{
		  if (Iwrite != NO_SCREEN)
		    printf
		      ("Not computing vext_coulomb due to boundary conditions\n");
		}
	    }
	  if (Iwrite == VERBOSE)
	    {
	      print_vext (Vext, output_file2);
	      print_zeroTF (Zero_density_TF, output_TF);
	      if (Vol_charge_flag && Ndim == 3)
		print_vext (Vext_coul, output_file4);
	    }

	  printf ("after coputing vext_coulomb\n");

	  /*
	   * Use Newton's method to solve the problem
	   * (see file dft_newton.c, which calls dft_fill.c and Aztec library)
	   * The variable "successful" is TRUE if Newton's method converrged.
	   */
	  x =
	    (double *) array_alloc (1, Nnodes_per_proc * Nunk_per_node,
				    sizeof (double));
	  if (Lbinodal)
	    x2 =
	      (double *) array_alloc (1, Nnodes_per_proc * Nunk_per_node,
				      sizeof (double));


	  t_preprocess += MPI_Wtime ();
	  t_solve = -MPI_Wtime ();
	  T_av_precalc_min = T_av_fill_min = T_av_solve_min = T_av_lj_min =
	    0.0;
	  T_av_precalc_max = T_av_fill_max = T_av_solve_max = T_av_lj_max =
	    0.0;
	  successful = solve_problem (&x, &x2, output_file1, &niters);
	  t_solve += MPI_Wtime ();
	  /*
	   * Post-Process the results
	   * (see file dft_output.c)
	   */
	  time_save = MPI_Wtime () - t_mesh;
	  t_postprocess = -MPI_Wtime ();

	  if (Loca.method == -1)
	    {
	      if (Lbinodal)
		post_process (x2, output_file3, &niters, &time_save,
			      Imain_loop, TRUE);
	      if (Proc == 0 && Lbinodal)
		{
		  X2_old =
		    (double *) array_alloc (1, Nodes_old * Nunk_per_node,
					    sizeof (double));
		  for (i = 0; i < Nodes_old * Nunk_per_node; i++)
		    X2_old[i] = X_old[i];
		}
	      post_process (x, output_file3, &niters, &time_save, Imain_loop,
			    FALSE);
	    }

	  t_postprocess += MPI_Wtime ();
/*      printf("\n*********************************************\n");
      printf("The number of nonzeros counted is %d\n",count_nonzero);
      printf("The number of zeros counted is %d\n",count_zero);
      printf("*********************************************\n");
*/

	  iend = 0;
	  match = 0;
	  for (idim = 0; idim < Ndim; idim++)
	    if (Esize_x[idim] == Esize_x_final[idim])
	      match++;
	  if (match == Ndim)
	    iend = 1;		/* we've done the final mesh */

	  for (idim = 0; idim < Ndim; idim++)
	    {
	      if (Esize_x[idim] > Esize_x_final[idim])
		{
		  if (Esize_x[idim] == 0.2)
		    Esize_x[idim] = 0.1;
		  else if (Esize_x[idim] == 0.1)
		    Esize_x[idim] = 0.05;
		  else if (Esize_x[idim] == 0.05)
		    Esize_x[idim] = 0.02;
		  else if (Esize_x[idim] == 0.02)
		    Esize_x[idim] = 0.01;
		  else
		    Esize_x[idim] = Esize_x_final[idim];
		}
	      if (Esize_x[idim] < Esize_x_final[idim])
		Esize_x[idim] = Esize_x_final[idim];
	    }

	  /*
	   * Get total CPU time for run 
	   * (see files md_timer_*.c)
	   */

	  t_total = MPI_Wtime () - start_t;
	  t_pre_max = AZ_gmax_double (t_preprocess, Aztec.proc_config);
	  t_msr_max = AZ_gmax_double (T_msr_setup, Aztec.proc_config);
	  t_solve_max = AZ_gmax_double (t_solve, Aztec.proc_config);
	  t_post_max = AZ_gmax_double (t_postprocess, Aztec.proc_config);
	  t_total_max = AZ_gmax_double (t_total, Aztec.proc_config);

	  t_pre_min = AZ_gmin_double (t_preprocess, Aztec.proc_config);
	  t_msr_min = AZ_gmin_double (T_msr_setup, Aztec.proc_config);
	  t_solve_min = AZ_gmin_double (t_solve, Aztec.proc_config);
	  t_post_min = AZ_gmin_double (t_postprocess, Aztec.proc_config);
	  t_total_min = AZ_gmin_double (t_total, Aztec.proc_config);

	  if (Proc == 0 && Iwrite != NO_SCREEN)
	    {
	      printf ("\n\n\n\n");
	      printf
		("===================================================\n");
	      printf ("TIMING SUMMARY .... averaged over all %d iterations\n",
		      niters);
	      printf
		("===================================================\n");
	      printf
		("                         min(sec)        max(sec)  \n");
	      printf
		("---------------------------------------------------\n");
	      printf ("MESH SETUP               %g         %g       \n",
		      t_pre_max, t_pre_min);
	      printf
		("---------------------------------------------------\n");
	      printf ("MSR_PREPROC              %g         %g       \n",
		      t_msr_min, t_msr_max);
	      printf ("PRE_CALC ROUTINES        %g        %g       \n",
		      T_av_precalc_min / niters, T_av_precalc_max / niters);
	      if (Matrix_fill_flag < 3)
		printf ("JACOBIAN                 %g         %g       \n",
			T_av_lj_min / niters, T_av_lj_max / niters);
	      printf ("FILL ROUTINE             %g         %g       \n",
		      T_av_fill_min / niters, T_av_fill_max / niters);
	      printf
		("---------------------------------------------------\n");
	      printf ("AZTEC SOLVE              %g          %g       \n",
		      T_av_solve_min / niters, T_av_solve_max / niters);
	      printf
		("---------------------------------------------------\n");
	      printf ("TOTAL SOLVE              %g         %g       \n",
		      t_solve_min, t_solve_max);
	      printf
		("---------------------------------------------------\n");
	      printf ("POST-PROCESSING          %g         %g       \n",
		      t_post_min, t_post_max);
	      printf
		("===================================================\n");
	      printf ("TOTAL TIME               %g         %g       \n",
		      t_total_min, t_total_max);
	      printf
		("===================================================\n");

	    }

	}			/* end of loop over different mesh densities (for automated
				   mesh refinement */

      /*
       * release lots of memory to start a new mesh
       */


      /* thermo arrays */
      if ((Ipot_ff_n == LJ12_6) && !Sten_Type[POLYMER_CR])
	{
	  safe_free ((void *) &Avdw);
	  safe_free ((void *) &Betamu_att);
	}
      if (Ipot_ff_c == COULOMB)
	safe_free ((void *) &Deltac_b);

      /* solution arrays */
      if (Ipot_ff_n != IDEAL_GAS)
	safe_free ((void *) &B2L_1stencil);
      safe_free ((void *) &B2L_unknowns);
      safe_free ((void *) &x);
      if (Lbinodal)
	safe_free ((void *) &x2);

      /* Aztec arrays -- mostly allocated by AZ_transform */
      safe_free ((void *) &Aztec.external);
      safe_free ((void *) &Aztec.extern_index);
      safe_free ((void *) &Aztec.update_index);
      /*AZ_free_memory(Aztec.data_org[AZ_name]); */
      safe_free ((void *) &Aztec.data_org);

    }				/* end of loop over continuation field #1 */

  if (Nwall_type != 0)
    {
      safe_free ((void *) &Link_list);
      safe_free ((void *) &Nwall_this_link);
    }
  free_mesh_arrays ();
  boundary_free ();
  if (Lsteady_state && Ndim == 1)
    safe_free ((void *) &Area_IC);
  safe_free ((void *) &Comm_node_proc);
  safe_free ((void *) &Comm_unk_proc);
  safe_free ((void *) &Comm_offset_node);
  safe_free ((void *) &Comm_offset_unk);
  safe_free ((void *) &Aztec.update);

  safe_free ((void *) &Xtest_reflect_TF);
  safe_free ((void *) &Lsemiperm);
  safe_free ((void *) &Vext_membrane);
/*
  safe_free((void *)  &Unk_to_Poly);
  safe_free((void *)  &Unk_to_Seg);
  safe_free((void *)  &Unk_to_Bond);
  safe_free((void *)  &Bonds);
  safe_free((void *)  &Nbond);
*/


  /* free stencils */
  /* ALF: fixed for changed stencils */
  for (izone = 0; izone < Nzone; izone++)
    {
      for (isten = 0; isten < NSTEN; isten++)
	if (Sten_Type[isten] == 1)
	  {
	    if (isten == U_ATTRACT || isten == THETA_CHARGE
		|| isten == POLYMER_CR)
	      jmax = Ncomp;
	    else if (isten == DELTA_FN && Sten_Type[POLYMER_CR])
	      jmax = Ncomp;
	    else
	      jmax = 1;

	    for (icomp = 0; icomp < Ncomp; icomp++)
	      {
		for (jcomp = 0; jcomp < jmax; jcomp++)
		  {
		    sten = &(Stencil[isten][izone][icomp + Ncomp * jcomp]);
		    safe_free ((void **) &sten->Offset);
		    safe_free ((void **) &sten->Weight);
		    if (Lhard_surf)
		      safe_free ((void **) &sten->HW_Weight);
		  }
	      }
	  }
    }
  safe_free ((void **) &Stencil);

  if (Sten_Type[POLYMER_CR])
    safe_free ((void **) &Rism_cr);


  if (Proc == 0)
    safe_free ((void *) &X_old);	/* This is also done for each continuation run
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
void
continuation_shift ()
{
  int idim, iwall, icharge;

  for (idim = 0; idim < Ndim; idim++)
    {
      Size_x[idim] += Del_1[idim];

      if (Plane_new_nodes == idim)
	{
	  if (Pos_new_nodes == 0)
	    {			/*adding nodes to center */
	      for (iwall = 0; iwall < Nwall; iwall++)
		{
		  if (WallPos[idim][iwall] < 0.0)
		    WallPos[idim][iwall] -= 0.5 * Del_1[idim];
		  else if (WallPos[idim][iwall] > 0.0)
		    WallPos[idim][iwall] += 0.5 * Del_1[idim];
		}
	      for (icharge = 0; icharge < Nlocal_charge; icharge++)
		{
		  if (Charge_x[icharge][idim] < 0.0)
		    Charge_x[icharge][idim] -= 0.5 * Del_1[idim];
		  else if (Charge_x[icharge][idim] > 0.0)
		    Charge_x[icharge][idim] += 0.5 * Del_1[idim];
		}
	    }
	  else if (Pos_new_nodes == -1)
	    {			/*adding nodes to lbb */
	      for (iwall = 0; iwall < Nwall; iwall++)
		WallPos[idim][iwall] += 0.5 * Del_1[idim];
	      for (icharge = 0; icharge < Nlocal_charge; icharge++)
		Charge_x[icharge][idim] += 0.5 * Del_1[idim];
	    }
	  else if (Pos_new_nodes == 1)
	    {			/*adding nodes to rft */
	      for (iwall = 0; iwall < Nwall; iwall++)
		WallPos[idim][iwall] -= 0.5 * Del_1[idim];
	      for (icharge = 0; icharge < Nlocal_charge; icharge++)
		Charge_x[icharge][idim] -= 0.5 * Del_1[idim];
	    }
	}			/* end of check for idim = Plane_new_nodes */
    }				/* end of loop over dimentions */
  return;
}

/******************************************************************/
/*setup_polymer_cr: read in c(r) from file and add attractions
  ALF: modified 10/21 for Type_poly == 3 */
void
setup_polymer_cr ()
{
  FILE *fp7;
  int i, j, lines, ir, nextChar;
  double r, u, cr_rad_max, rsave, dummy_read;
  char c;

  /* reading in c(r) file */
  if (Proc == 0)
    printf ("reading in c(r) file...\n");
  if (Type_poly != 3)
    {

      if (Proc == 0)
	{
	  if ((fp7 = fopen (Cr_file, "r")) == NULL)
	    {
	      printf ("Can't open file %s\n", Cr_file);
	      exit (1);
	    }
	  for (ir = 1; ir <= 3; ir++)
	    {
	      fscanf (fp7, "%lf", &r);
	      fscanf (fp7, "%c", &c);
	      for (i = 0; i < Ncomp; i++)	/* for (i=0; i<Ntype_mer; i++)  */
		fscanf (fp7, "%lf", &dummy_read);
	      fscanf (fp7, "%c", &c);
	      for (i = 0; i < Ncomp; i++)
		{		/* "  */
		  for (j = i + 1; j < Ncomp; j++)
		    {		/* "  */
		      fscanf (fp7, "%lf", &dummy_read);
		      fscanf (fp7, "%c", &c);
		    }
		}
	      if (ir == 1)
		rsave = r;
	      if (ir == 2)
		Deltar_cr = r - rsave;
	    }
	  fclose (fp7);
	}
      MPI_Bcast (&Deltar_cr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


      lines = 0;
      cr_rad_max = 0.0;
      for (i = 0; i < Ncomp; ++i)
	for (j = 0; j < Ncomp; ++j)
	  {
	    if (Cr_rad_hs[i][j] > cr_rad_max)
	      cr_rad_max = Cr_rad_hs[i][j];
	  }

      lines = (int) (cr_rad_max / Deltar_cr);

      if (Proc == 0 && lines >= N_NZCR_MAX)
	{
	  printf
	    ("Warning: Rism_cr array may be truncated :: Max Cr_rad_hs > allowed\n");
	  lines = N_NZCR_MAX - 1;
	}
      /* printf("Proc=%d lines=%d\n",Proc,lines); */

      Last_nz_cr = lines;
      if (Proc == 0)
	{
	  fp7 = fopen (Cr_file, "r");
	  for (ir = 1; ir <= lines; ir++)
	    {
	      fscanf (fp7, "%lf", &r);
	      fscanf (fp7, "%c", &c);
	      if (ir == 1)
		rsave = r;
	      if (ir == 2)
		Deltar_cr = r - rsave;
	      for (i = 0; i < Ncomp; i++)
		{		/* for (i=0; i<Ntype_mer; i++)  */
		  fscanf (fp7, "%lf", &Rism_cr[i][i][ir]);
		  /*Rism_cr[i][i][ir] *=0.7; */
		  if (i == 0)
		    Rism_cr[i][i][ir] *= 1.0;	/*1.052758; */
		  /*else if (i==1) Rism_cr[i][i][ir] *=0.78;*//*0.750362672; */
		  /*else if (i==2) Rism_cr[i][i][ir] *=0.80; *//*0.7723047; */
		}
	      fscanf (fp7, "%c", &c);
	      for (i = 0; i < Ncomp; i++)
		{		/* "  */
		  for (j = i + 1; j < Ncomp; j++)
		    {		/* "  */
		      fscanf (fp7, "%lf", &Rism_cr[i][j][ir]);
		      if (i == 0 && j == 1)
			Rism_cr[i][j][ir] *= 1.0;	/*0.9100006; */
		      else if (i == 0 && j == 2)
			Rism_cr[i][j][ir] *= 1.0;	/*0.91902196; */
		      else if (i == 1 && j == 2)
			Rism_cr[i][j][ir] *= 1.0;	/*0.77002807; */
/*   Rism_cr[i][j][ir] *=0.7;*/
		      fscanf (fp7, "%c", &c);
		      Rism_cr[j][i][ir] = Rism_cr[i][j][ir];
		    }
		}
	      while (c != '\n')
		c = getc (fp7);
	    }
	  fclose (fp7);
	  for (ir = lines + 1; ir < N_NZCR_MAX; ir++)
	    for (i = 0; i < Ncomp; i++)
	      {
		Rism_cr[i][i][ir] = 0.;
		for (j = i + 1; j < Ncomp; j++)
		  {
		    Rism_cr[i][j][ir] = 0.;
		    Rism_cr[j][i][ir] = 0.;
		  }
	      }
	  /* extrapolate c(r) to r = 0 */
	  for (i = 0; i < Ncomp; i++)	/* "  */
	    for (j = 0; j < Ncomp; j++)	/* "  */
	      Rism_cr[i][j][0] = 2. * Rism_cr[i][j][1] - Rism_cr[i][j][2];
	}
      MPI_Bcast (**Rism_cr, Ncomp * Ncomp * N_NZCR_MAX, MPI_DOUBLE, 0,
		 MPI_COMM_WORLD);

      /* add attractions to polymer c(r)  */

      if (Iwrite == VERBOSE && Proc == 0)
	{
	  fp7 = fopen ("cr.out", "w");
	  for (ir = 0; ir <= Last_nz_cr; ++ir)
	    {
	      fprintf (fp7, "%lf   ", ir * Deltar_cr);
	      for (i = 0; i < Ncomp; ++i)
		for (j = 0; j < Ncomp; ++j)
		  fprintf (fp7, "%lf   ", Rism_cr[i][j][ir]);
	      fprintf (fp7, "\n");
	    }
	  fclose (fp7);
	}

      if (Ipot_ff_n == LJ12_6)
	{
	  for (i = 0; i < Ncomp; ++i)
	    for (j = 0; j < Ncomp; ++j)
	      {
		if (Cut_ff[i][j] > Cr_rad_hs[i][j])
		  Cr_rad[i][j] = Cut_ff[i][j];
		else
		  Cr_rad[i][j] = Cr_rad_hs[i][j];
		if (Cr_rad[i][j] > cr_rad_max)
		  cr_rad_max = Cr_rad[i][j];
	      }

	  lines = cr_rad_max / Deltar_cr;
	  if (lines >= N_NZCR_MAX)
	    {
	      if (Proc == 0)
		printf ("Need to increase N_NZCR_MAX\n");
	      lines = N_NZCR_MAX - 1;
	    }
	  if (lines > Last_nz_cr)
	    Last_nz_cr = lines;

	  for (i = 0; i < Ncomp; ++i)
	    for (j = i; j < Ncomp; j++)
	      {
		for (ir = 0; ir <= lines; ir++)
		  {
		    r = ir * Deltar_cr;
		    if (r - 1.e-8 > Sigma_ff[i][j] * 1.122462)
		      {		/* watch roundoffs */
			u = uLJatt_n (r, i, j);
			Rism_cr[i][j][ir] -= u;
			if (i != j)
			  Rism_cr[j][i][ir] -= u;
		      }
		  }
	      }
	}
      if (Iwrite == VERBOSE && Proc == 0)
	{
	  fp7 = fopen ("cr.lj.out", "w");
	  for (ir = 0; ir <= Last_nz_cr; ++ir)
	    {
	      fprintf (fp7, "%lf   ", ir * Deltar_cr);
	      for (i = 0; i < Ncomp; ++i)
		for (j = 0; j < Ncomp; ++j)
		  fprintf (fp7, "%lf   ", Rism_cr[i][j][ir]);
	      fprintf (fp7, "\n");
	    }
	  fclose (fp7);
	}

    }				/* end of if (Type_poly != 3) */

  /* for Type_poly == 3, just read in c(0) */
  else
    {
      if (Proc == 0)
	{
	  if ((fp7 = fopen (Cr_file, "r")) == NULL)
	    {
	      printf ("Can't open file %s\n", Cr_file);
	      exit (1);
	    }

	  for (i = 0; i < Ncomp; i++)
	    {			/* for (i=0; i<Ntype_mer; i++)  */
	      fscanf (fp7, "%lf", &Rism_cr[i][i][0]);
	      fscanf (fp7, "%c", &c);
	      printf ("cr[%d][%d] = %f\t", i, i, Rism_cr[i][i][0]);
	    }
	  for (i = 0; i < Ncomp; i++)
	    {			/* "  */
	      for (j = i + 1; j < Ncomp; j++)
		{		/* "  */
		  fscanf (fp7, "%lf", &Rism_cr[i][j][0]);
		  fscanf (fp7, "%c", &c);
		  Rism_cr[j][i][0] = Rism_cr[i][j][0];
		  printf ("cr[%d][%d] = %f\t", i, j, Rism_cr[i][j][0]);
		}
	    }
	  fclose (fp7);
	  for (ir = 1; ir < N_NZCR_MAX; ir++)
	    {
	      for (i = 0; i < Ncomp; i++)
		{
		  Rism_cr[i][i][ir] = 0.;
		  for (j = i + 1; j < Ncomp; j++)
		    {
		      Rism_cr[i][j][ir] = 0.;
		      Rism_cr[j][i][ir] = 0.;
		    }
		}
	    }
	}			/* end of if(Proc==0) */

      MPI_Bcast (**Rism_cr, Ncomp * Ncomp * N_NZCR_MAX, MPI_DOUBLE, 0,
		 MPI_COMM_WORLD);
    }				/* end of Type_poly==3 */

  return;
}

/*************************************************************/
