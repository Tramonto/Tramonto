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
 * $Name$
 *====================================================================*/

/* DOCUMENTATION FOR CONTINUATION LIBRARY:    10/26/1999
 *    Andrew Salinger  Org. 9221, MS-1111, 845-3523
 *
 * The continuation library currently consists of 4 files. 
 * (1) [rf_continuation.c]  The MPSalsa version of the interface
 *     to the continuation library. This file (and only this file)
 *     must be extensively modified for every new code that
 *     uses the library. It consists of solve_continuation, the
 *     top level routine called by the application code, and ~10
 *     wrapper routines, which must be filled for your specific
 *     code. The most noteworthy of these include
 *     nonlinear_solver_conwrap, linear_solver_conwrap,
 *     matrix_residual_fill_conwrap, and assign_parameter_conwrap.
 *     Each wrapper has its own comments.
 * (2) [rf_con_lib.c]  This is the stepping algorithm for zeroth-order,
 *     first-order, and arc-length continuation. The turning point
 *     tracking algorithm is only zeroth-order. This routine includes
 *     such things as step-size control and predictor steps. Ideally,
 *     this file does not need to be modified when linking to a new
 *     application code.
 * (3) [rf_con_bord.c] This file contains the bordering algorithms,
 *     currently for arc-length continuation and turning point tracking.
 *     These routines are accessed from within the Newton iteration,
 *     and make frequent calls to the wrapper routines (matrix fill
 *     and solve). Ideally, this file does not need to be modified when
 *     linking to a new application code.
 * (4) [rf_con_const.h] This header file includes definitions of the
 *     continuation structures, define statements for flags used
 *     by the continuation library, and prototypes for the wrapper
 *     routines. Ideally, this file does not need to be modified when
 *     linking to a new application code.
 *
 * How to interface to this library:
 * (0) Have a steady-state code that uses Newton's method
 * (1) Call solve_continuation from your code, in the place
 *     where you normally call the steady-state or transient
 *     drivers.
 * (2) In your current nonlinear solver, add "(void *) con_ptr"
 *     to the argument list. Pass NULL in that place for all
 *     calls to the nonlinear solver besides the one from
 *     nonlinear_solver_conwrap below. In your Newton loop,
 *     after the matrix equation has been solved for the
 *     update vector, delta_x, but before the solution vector, x,
 *     has been updated by delta_x, put in the following call:
 *      if (con_ptr != NULL) continuation_converged =
 *         continuation_hook(x, delta_x, con_ptr, Reltol, Abstol);
 *     (Reltol=1.0e-3, Abstol=1.0e-8 are good defaults.)
 *     When checking convergence of your Newton's method, also
 *     check that (continuation_converged == TRUE).
 * (3) Change the contents of all routines in this file. In
 *     solve_continuation, many fields of the con and set_con structures
 *     must be set. Follow the template and the comments.
 *     Also, the passdown structure can be defined and set with 
 *     any information that needs to be passed from the 
 *     solve_continuation routine to the wrapper routines below
 *     that isn't needed by the continuation library. All
 *     of the wrapper routines (all routines in this file besides
 *     solve_continuation, which all have names ending in _conwrap)
 *     must be filled in with the corresponding call to your
 *     application code. Each has its own comments.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Put include statements for your code here. */
#include "dft_globals_const.h"
#include "rf_allo.h"

/*****************************************************************************/

/* This include file is for the continuation structures, defines, prototypes */

#include "loca_const.h"
static void print_con_struct(const struct con_struct* con);

/* Define passdown structure: this structure is global to this file, and
 * provides a way to pass variables from the top solve_continuation routine
 * down to the various wrapper routines below, without passing through the
 * continuation library. For instance, the Jacobian matrix is never seen
 * in the continuation routines, but is needed by several of the wrapper
 * routines in this file. The passdown structure is a way make the
 * matrix global to this file.
 */

struct passdown_struct {
double epswf_previous; /* Save old value of Eps_wf for rescaling external field*/
double epsff_previous; /* Save old value of Eps_ff for rescaling stencils */
double temp_previous; /* Save old value of Temp for rescaling external field*/
double chg_scale_previous;
int    polymer_flag;  /* Flag for calling polymer vs.regular fill*/
} passdown;

static double get_init_param_value(int cont_type);
static void print_final(double param, int step_num);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int solve_continuation( double *x, double *x2, int polymer_flag, void * aux_info)

/* Interface routine to the continuation library.
 */
{

  /* Define 2 continuation structures */
  struct con_struct con;

  /* Local Variables */
  int nstep=0;

  /******************************* First Executable Statment *****************/

  /* 
   * Set passdown structure -- variables needed in the argument
   * lists to wrapped routines but not needed in the continuation
   * library.
   */


   passdown.polymer_flag = polymer_flag;

  /*
   * Fill all the necessary information into the 'con' structure
   * for running the continuation routines. (In this example, this
   * info was read into an identical structure, so the copying isn't
   * necessary. However, if a application has input the continuation
   * information in a different format, these lines below will be
   * needed to load the 'con' structure.)
   */

  /* First load the general info structure */

  switch (Loca.method) {
    case 0: con.general_info.method = ZERO_ORDER_CONTINUATION; break;
    case 1: con.general_info.method = FIRST_ORDER_CONTINUATION; break;
    case 2: con.general_info.method = ARC_LENGTH_CONTINUATION; break;
    case 3: con.general_info.method = TURNING_POINT_CONTINUATION; break;
    case 4: con.general_info.method = PHASE_TRANSITION_CONTINUATION; break;
    default:
       if (Proc==0) printf("\nERROR in solve_continuation: Unknown "
                           "continuation method: %d\n",Loca.method);
       exit(-1);
  }
 
  printf("\n###\nPROC = %d\n\n",Proc);

  con.general_info.param        = get_init_param_value(Loca.cont_type1);
  con.general_info.x            = x;
  con.general_info.numUnks      = Nunk_int_and_ext;
  con.general_info.numOwnedUnks = Aztec.N_update;
  if (Proc==0) {
     switch (Iwrite) {
	     case NO_SCREEN: con.general_info.printproc = 0; break;
	     case MINIMAL:   con.general_info.printproc = 1; break;
	     case DENSITIES: con.general_info.printproc = 5; break;
	     case VERBOSE:   con.general_info.printproc = 8; break;
     }
  }
  else         con.general_info.printproc = FALSE;
  con.general_info.nv_save   = FALSE;   /* Salsa saves null vector anyway */
  con.general_info.perturb  = 1.0e-7;


  /* Then load the stepping info structure */

  con.stepping_info.first_step     = Loca.step_size;
  con.stepping_info.base_step      = 0;
  con.stepping_info.max_steps      = Loca.num_steps;
  con.stepping_info.max_param      = 1.0e10;
  con.stepping_info.max_delta_p    = 1.0e10;
  con.stepping_info.min_delta_p    = 1.0e-10;
  con.stepping_info.step_ctrl      = Loca.aggr;
  con.stepping_info.max_newton_its = Max_Newton_iter;

  /* Then load one of the method dependent structures */

  switch (con.general_info.method) {
    case ZERO_ORDER_CONTINUATION:
    case FIRST_ORDER_CONTINUATION:
      break;
    case ARC_LENGTH_CONTINUATION:
      con.arclength_info.dp_ds2_goal = 0.1;
      con.arclength_info.dp_ds_max   = 0.6;
      con.arclength_info.tang_exp    = 2.0;
      con.arclength_info.tang_step_limit = 0.5;
      break;
    case TURNING_POINT_CONTINUATION:
      con.turning_point_info.bif_param = get_init_param_value(Loca.cont_type2);
      con.general_info.nv_restart       = FALSE;
      break;
    case PITCHFORK_CONTINUATION:
      con.pitchfork_info.bif_param = get_init_param_value(Loca.cont_type2);
      con.pitchfork_info.psi        = NULL;
      break;
    case HOPF_CONTINUATION:
      con.hopf_info.bif_param = get_init_param_value(Loca.cont_type2);
      con.hopf_info.omega       = 0.0;
      con.hopf_info.y_vec       = NULL;
      con.hopf_info.z_vec       = NULL;
      con.hopf_info.mass_flag   = 0;
      break;
    case PHASE_TRANSITION_CONTINUATION:
      con.phase_transition_info.bif_param
          = get_init_param_value(Loca.cont_type2);
      con.phase_transition_info.x2  = x2;
      break;
    default:
      printf("ERROR: Unknown LOCA input method: %d\n",con.general_info.method);
      exit(-1);
  }

  /* Finally, load the eigensolver structures */

  con.eigen_info.Num_Eigenvalues   = 0;
  con.eigen_info.Shift_Point[0]    = 5.0;
  con.eigen_info.Shift_Point[1]    = 50.0;
  con.eigen_info.Shift_Point[2]    = 1.0;
  con.eigen_info.Arnoldi           = 30;
  con.eigen_info.Residual_Tol[0]   = 1.0e-4;
  con.eigen_info.Residual_Tol[1]   = 1.0e-4;
  con.eigen_info.Max_Iter          = 1;
  con.eigen_info.Every_n_Steps     = 1;

  /* print out continuation structure */

  print_con_struct(&con);

  /* Now call continuation library and return */

  nstep = con_lib(&con, aux_info);

  if (con.general_info.printproc) print_final(con.general_info.param, nstep);

  return nstep;
} /**************** END of solve_continuation () *****************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int nonlinear_solver_conwrap(double *x, void *con_ptr, int step_num,
		             double lambda, double delta_s, void * aux_info)
/* Put the call to your nonlinear solver here.
 * Input:
 *    x         solution vector
 *    con_ptr   pointer to continuation structure, cast to (void *)
 *              must be passed to nonlinear solver and then passed
 *              to bordering algorithms.
 *    step_num  Continuation step number
 *
 * Output:
 *
 * Return Value:
 *    num_its  Number of Newton iterations needed for
 *             convergence, used to pick next step size.
 *             Negative value means nonlinear solver didn't converge.
 */
{
  int num_its;
  double t=0; /* dumm spot with solve time returned */

  num_its = newton_solver(x, NULL, NULL, con_ptr, Max_Newton_iter, &t, aux_info);

  return (num_its);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int linear_solver_conwrap(double *x, int jac_flag, double *tmp)
/* Put the call to your linear solver here. There are three options
 * about the reuse of the preconditioner. It is always safe to 
 * just solve the matrix from scratch.
 * Input:
 *    x          Right hand side
 *    jac_flag   Flag indicating the status of the Jacobian so that
 *               preconditioners can be used: 
 *               NEW_JACOBIAN:   recalculate preconditioner
 *               OLD_JACOBIAN:   reuse preconditioner
 *               SAME_BUT_UNSCALED_JACOBIAN: Must rescale the matrix and
 *                               can reuse preconditioner. This happens
 *                               when the matrix has been recalculated
 *                               at the same conditions as before.
 *    tmp        Work space array same length as x, only used for
 *               the SAME_BUT_UNSCALED_JACOBIAN option.
 *
 * Output:
 *    x          Solution vector
 *
 * Return Value:
 *    Negative value means linear solver didn't converge.
 */
{
  double *x_tmp;
  int i;

  x_tmp = (double *) array_alloc(1, Nunk_int_and_ext, sizeof(double));
  for (i=0; i< Nunk_int_and_ext; i++) x_tmp[i] = 0.0;

  if (jac_flag == OLD_JACOBIAN || jac_flag == CHECK_JACOBIAN) {

    /* reuse the preconditioner for this same matrix */
    Aztec.options[AZ_pre_calc] = AZ_reuse;

  }
  else if (jac_flag == SAME_BUT_UNSCALED_JACOBIAN) {

    /* This option is for when the matrix is the same as the last
     * solve but it has been recalculated. Therfore it must be
     * rescaled.
     */

     /* Currently, reuse preconditioner if there was no scaling,  */
     /* otherwise recomput it                                     */

     if (Aztec.options[AZ_scaling] == AZ_none) {
         Aztec.options[AZ_pre_calc] = AZ_reuse;

     /* In future, use AZ_MATRIX structure and switch to this coding */

/**
     scaling = AZ_scaling_create();
     AZ_scale_f(AZ_SCALE_MAT_RHS_SOL, m, passdown.options, tmp, NULL,
                passdown.proc_config, scaling);
     AZ_scaling_destroy(&scaling);

     Aztec.options[AZ_pre_calc] = AZ_reuse;
**/

    }
  }
  else if (jac_flag != NEW_JACOBIAN) {
    printf("ERROR: linear solve conwrap: unknown flag value %d\n",jac_flag);
    exit(-1);
  }


  AZ_solve(x_tmp, x, Aztec.options, Aztec.params, NULL, Aztec.bindx,
           NULL, NULL, NULL, Aztec.val, Aztec.data_org, Aztec.status,
           Aztec.proc_config);

  for (i=0; i< Nunk_int_and_ext; i++) x[i] = x_tmp[i];
  safe_free((void *) &x_tmp);
 
  Aztec.options[AZ_pre_calc] = AZ_calc;
 
  return 0;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int komplex_linear_solver_conwrap(double *c, double *d,
    int jac_flag, double *omega, double *tmp){return 0;}
void mass_matrix_fill_conwrap(double *x, double *rhs){}
void mass_matvec_mult_conwrap(double *x, double *y){}
void create_shifted_matrix_conwrap(){}
void shifted_matrix_fill_conwrap(double sigma){}
void shifted_linear_solver_conwrap(double *x, double *y,
                                   int jac_flag, double tol){}
void destroy_shifted_matrix_conwrap(){}

/*****************************************************************************/
void matrix_residual_fill_conwrap(double *x, double *rhs, int matflag)
/* Put the call to your matrix/residual fill routine here.
 * Input:
 *    x         Solution vector
 *    matflag   Flag indicating residual (RHS_ONLY), matrix (MATRIX_ONLY),
 *              or both (RHS_MATRIX) are requested.
 *
 * Output:
 *    rhs       Right hand side
 *
 * Return Value:
 */
{
  int i, resid_only_flag;
  double l2_resid;

  if (matflag == RHS_ONLY || matflag == RHS_MATRIX_SAVE) resid_only_flag = TRUE;
  else                     resid_only_flag = FALSE;

  /*fill_time not currently plugged in, iter hardwire above 2 */
  fill_resid_and_matrix_control(x, rhs, NULL, NULL, Matrix_fill_flag,
                                                 4, resid_only_flag);

  /* calculate and print resid norm */
  l2_resid=0.0;
  for (i=0; i<Aztec.N_update; i++) l2_resid += rhs[i]*rhs[i];

  l2_resid= sqrt(gsum_double_conwrap(l2_resid));
  if (Proc==0) printf("\t\tNorm of resid vector = %g\n", l2_resid);

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void matvec_mult_conwrap(double *x, double *y)
/* Put the call to your matrix-vector multiply here.
 * Input:
 *    x         Vector of length number of unknowns
 *
 * Output:
 *    y         Matrix times x.
 *
 * Return Value:
 */
{

 AZ_matvec_mult(Aztec.val, NULL, Aztec.bindx,NULL,
		NULL, NULL, x, y, 1, Aztec.data_org);

}
/*****************************************************************************/
/*************************************************************/
/* scale_elec_param: when changing eps wall only a
   simple scaling is required. */
static void scale_elec_param(double ratio)
{
   int iwall,iel,loc_inode,idim;

   if (Surf_charge_flag){ 
      for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++)
           for (idim=0; idim<Ndim; idim++) Charge_w_sum_els[loc_inode][idim] *=ratio;
   }
   if (Vol_charge_flag) for (iel=0; iel<Nelements_box; iel++) Charge_vol_els[iel]*=ratio;
   for (iwall=0; iwall<Nwall; iwall++) Elec_param_w[iwall]*=ratio;

   if (Iwrite == VERBOSE){
       if (Proc==0) printf ("PRINTING CHARGE DISTRIBUTIONS: Scale_fac=%9.6f\n",Scale_fac);
       if (Vol_charge_flag) print_charge_vol(Charge_vol_els,"dft_charge_vol.dat");
       if (Surf_charge_flag) print_charge_surf(Charge_w_sum_els,"dft_charge_surf.dat");
    }

   return;
}

/*****************************************************************************/
/* recalculate_stencils:  free the old Stencil structure and recompute stencil 
   with the new parameters */
static void recalculate_stencils()
{
   int izone, isten,jmax,icomp,jcomp;
   struct Stencil_Struct *sten;

   for (izone=0; izone<Nzone; izone++){
     for (isten=0; isten<NSTEN; isten++)
       if (Sten_Type[isten]) {
          if (isten == U_ATTRACT || isten == THETA_CHARGE 
             || isten == POLYMER_CR || (isten==DELTA_FN && Sten_Type[POLYMER_CR])) jmax = Ncomp;
          else                                                                     jmax = 1;
                         
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
  calc_stencils();
  return;
}
/*************************************************************/
/* scale_att_stencil_temp: when changing the temperature, a simple
   scaling of the attractive integration stencil is required.  */
static void scale_att_stencil_temp(double temp)
{
/* were going to rescale attractive stencil, but then realized that polymer stencil
   will require complete recalculation as formulated.  So, we opted to just perform
   a recalculation of the entire stencil for this option.  Note that the polymer
   stencils should really be reformulated so the hard chain part is separate from
   the attractive part.  In that case, we could easily just scale the attractive
   part of the stencil.*/

   double ratio;
   int icomp,jcomp,izone,isten;
   struct Stencil_Struct *sten;
   FILE *ifp=NULL;

   if (Iwrite == VERBOSE) ifp = fopen("stencil.out", "a");
   ratio = passdown.temp_previous/temp;
   for (izone=0; izone < Nzone; izone++){
      for (icomp=0; icomp<Ncomp; icomp++){
      for (jcomp=0; jcomp<Ncomp; jcomp++){
          sten = &(Stencil[U_ATTRACT][izone][icomp+Ncomp*jcomp]);
          for (isten = 0; isten < sten->Length; isten++) {
              sten->Weight[isten] *= ratio;
/*              if (Iwrite==VERBOSE) print_out_stencil(isten, izone,icomp, jcomp, ifp);*/
          }
      }
      }
   }
   if (Iwrite==VERBOSE) fclose(ifp);
   return;
}
/*****************************************************************************/
/*************************************************************/
/* scale_external_field: when changing eps wall only a
   simple scaling is required. */
static void scale_vext_temp(double ratio)
{
   int loc_inode,icomp,iwall,idim,iunk;

   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
      for (icomp=0; icomp<Ncomp; icomp++){
         Vext[loc_inode][icomp] /= ratio;
         if (!Lhard_surf && Restart !=4 ){
         for (iwall=0; iwall<Nwall; iwall++)
         for (idim=0; idim<Ndim; idim++){
            iunk = iwall*Ncomp + icomp;
            Vext_dash[loc_inode][iunk][idim] /= ratio;
         }
         }
      }
   }
   if (Iwrite==VERBOSE) {
      print_vext(Vext,"dft_vext.dat");
   }
   return;
}

/*****************************************************************************/
/*************************************************************/
/* scale_external_field: when changing eps wall-fluid for the
   first component, only a simple scaling is required. */

static void scale_vext_epswf(double ratio, int icomp,int iwall)
/* nc is the number of components whose parameter has changed,
 * and is either 1 or Ncomp */
{
   int loc_inode,idim,iunk;
   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
       Vext[loc_inode][icomp] *= ratio;
       if (!Lhard_surf && Restart !=4){
       iunk = iwall*Ncomp+icomp;
       for (idim=0; idim<Ndim; idim++) Vext_dash[loc_inode][iunk][idim] *= ratio;
       }
   }
   if (Iwrite==VERBOSE) {
      print_vext(Vext,"dft_vext.dat");
   }
   return;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void assign_parameter_tramonto(int cont_type, double param)
/* Note: Post_processing assumes the cont_type flags are the same as those
   used in Tramonto's own continuation */
{
  int i,j,icomp,jcomp,iw,inode;
  double ratio,temp_save,scale_save,eps_wf_save[NCOMP_MAX][NWALL_MAX_TYPE],param_save;
  char     *output_TF,*output_file1, *output_file2;
  
  output_file1 = "dft_out.lis";
  switch(cont_type){
     case CONT_MESH: 
       printf("ERROR: Continuation Library cannot do mesh changes\n");
       exit(-1); break;

      case CONT_TEMP: 
                      temp_save = Temp;
                      Temp      = param;
                      ratio = Temp/temp_save;
                      if (Ipot_ff_c == COULOMB ) Temp_elec *=ratio;
                      if (Mix_type==0){
                         if (Ipot_ff_n == LJ12_6){
                            for (icomp=0; icomp<Ncomp; icomp++) Eps_ff[icomp][icomp] /= ratio;
                         }
                         for (i=0; i<Nwall_type;i++){
                             if(Ipot_wf_n[i] != VEXT_HARD) Eps_w[i] /= ratio;
                         }
                         pot_parameters(NULL);
                      }
                      else if (Mix_type==1){
                           for (icomp=0; icomp<Ncomp; icomp++){
                             for(jcomp=0; jcomp<Ncomp; jcomp++) Eps_ff[icomp][jcomp] /= ratio;
                             for (i=0; i<Nwall_type;i++) Eps_wf[icomp][i] /= ratio;
                           }
                           for (i=0; i<Nwall_type;i++) {
                              for (j=0;j<Nwall_type;j++) Eps_ww[i][j] /= ratio;
                           }
                      }
                      if (Type_poly >=0) setup_polymer_cr();
                      recalculate_stencils();
                      scale_vext_temp(ratio);
                      break;

      case CONT_RHO_0:   Rho_b[2]= param;    
                         /*Rho_b[0] = (16./18.)*(Rho_b[2]);
                         Rho_b[1] = (2./18.)*(Rho_b[2]);*/

                         Rho_b[0] = (16./18.)*(0.59-Rho_b[2]);
                         Rho_b[1] = (2./18.)*(0.59-Rho_b[2]);

                         /*Rho_b[0] = (16./18.)*(0.825-Rho_b[2]);
                         Rho_b[1] = (2./18.)*(0.825-Rho_b[2]);*/
                         if (Type_poly >=0) setup_polymer_cr();
                         recalculate_stencils();
                         break;
      case CONT_RHO_ALL: for (i=0; i<Ncomp;i++)  Rho_b[i]= param;    break;
      case CONT_LOG_RHO_0: Rho_b[0]        = exp(param);    break;
      case CONT_LOG_RHO_ALL: for (i=0;i<Ncomp;i++) Rho_b[i]= exp(param);    break;
      case CONT_SCALE_RHO: scale_save = Scale_fac;
                           Scale_fac = param;
                           ratio = Scale_fac/scale_save;
                           for (i=0; i<Ncomp-1; i++) Rho_b[i] *= ratio;
                           Rho_b[Ncomp-1]=0.678;
                           for (i=0; i<Ncomp-1; i++) Rho_b[Ncomp-1] -= Rho_b[i];
                           recalculate_stencils();
                           break;
                           

      case CONT_EPSW_0:  
                      if (Mix_type ==0) {
                         Eps_w[0] = param;
                         for (i=0; i<Ncomp; i++){ 
                           for (iw=0; iw<Nwall_type; iw++) eps_wf_save[i][iw]=Eps_wf[i][iw];
                         }
                         pot_parameters(NULL);
                         for (i=0; i<Ncomp; i++){
                            ratio = Eps_wf[i][0]/eps_wf_save[i][WallType[0]];
                            scale_vext_epswf(ratio,i,0); 
                         }
                      }
                      else Eps_ww[0][0]=param;
                      break;

      case CONT_EPSW_ALL: 
                      if (Mix_type==0){
                          for (i=0;i<Nwall_type;i++) Eps_w[i] = param;
                          for (i=0; i<Ncomp; i++){ 
                            for (iw=0; iw<Nwall_type; iw++) eps_wf_save[i][iw]=Eps_wf[i][iw];
                          }
                          pot_parameters(NULL);
                          for (i=0; i<Ncomp; i++){
                            for (iw=0; iw<Nwall; iw++){
                         
                               ratio = Eps_wf[i][WallType[iw]]/eps_wf_save[i][WallType[iw]];
                               scale_vext_epswf(ratio,i,iw); 
                            }
                          }
                      }
                      else {
                            for (i=0;i<Nwall_type;i++) 
                                for (j=0;j<Nwall_type;j++) Eps_ww[i][j] = param;
                      }
                      break;

      case CONT_SCALE_EPSW: 
                      scale_save = Scale_fac;
                      Scale_fac = param;
                      ratio = Scale_fac/scale_save;
                      if (Mix_type==0){
                          for (iw=0; iw <Nwall_type;iw++) Eps_w[iw] *= ratio;
                          for (i=0; i<Ncomp; i++){ 
                            for (iw=0; iw<Nwall_type; iw++) eps_wf_save[i][iw]=Eps_wf[i][iw];
                          }
                          pot_parameters(NULL);
                          for (i=0; i<Ncomp; i++){
                            for (iw=0; iw<Nwall; iw++){
                               ratio = Eps_wf[i][WallType[iw]]/eps_wf_save[i][WallType[iw]];
                               scale_vext_epswf(ratio,i,iw); 
                            }
                          } 
                        }
                        else{
                            for (i=0;i<Nwall_type;i++) 
                                for (j=0;j<Nwall_type;j++) Eps_ww[i][j] *=ratio;
                        }
                        break;

      case CONT_EPSWF00: 
/*                         ratio = param/Eps_wf[0][0];
                         Eps_wf[0][0]    = param;
                         scale_vext_epswf(ratio,0,0); break;*/

                          /* component 2 wall 0 */
/*                         ratio = param/Eps_wf[2][0];
                         Eps_wf[2][0]    = param;
                         scale_vext_epswf(ratio,2,0); break;*/

                         ratio = param/Eps_wf[1][0];
                         Eps_wf[1][0]    = param;
                         Eps_wf[1][1]    = param;
                         Eps_wf[2][0]    = param;
                         Eps_wf[2][1]    = param;
                         scale_vext_epswf(ratio,1,0); break;
                         scale_vext_epswf(ratio,1,1); break;
                         scale_vext_epswf(ratio,2,0); break;
                         scale_vext_epswf(ratio,2,1); break;

      case CONT_EPSWF_ALL_0: 
                         for (i=0; i<Ncomp-1; i++){ 
                            for (iw=0; iw<Nwall_type; iw++) eps_wf_save[i][iw]=Eps_wf[i][iw];
                         }
                         for (i=0; i<Ncomp-1; i++){
                              Eps_wf[i][0] = param;
                              ratio = Eps_wf[i][0]/eps_wf_save[i][WallType[0]];
                              scale_vext_epswf(ratio,i,0);
                         }
                         break;
      
      case CONT_SCALE_EPSWF:
                          scale_save = Scale_fac;
                          Scale_fac = param;
                          ratio = Scale_fac/scale_save;
                          for (i=0; i<Ncomp; i++){ 
                            for (iw=0; iw<Nwall_type; iw++) eps_wf_save[i][iw]=Eps_wf[i][iw];
                          }
                          for (i=0; i<Ncomp; i++){
                             for (iw=0; iw<Nwall_type; iw++){
                                Eps_wf[i][iw] *= ratio;
                             }
                          }
                          for (i=0; i<Ncomp; i++){
                             for (iw=0; iw<Nwall; iw++){
                                ratio = Eps_wf[i][WallType[iw]]/eps_wf_save[i][WallType[iw]];
                                scale_vext_epswf(ratio,i,iw);
                             }
                          }
                          break;

      case CONT_EPSFF_00: Eps_ff[0][0]=param;
/*now do a special case where we change two of them at once */
                /*         Eps_ff[2][0]=param;
                         Eps_ff[0][2]=param;
                         Eps_ff[1][2]=param;
                         Eps_ff[2][1]=param;*/

                      /*   Eps_ff[2][2]=param;
                         Eps_ff[1][2]=param;
                         Eps_ff[2][1]=param;*/

                       /*  Eps_ff[2][0]=param;*/
                         if (Mix_type==0) {
                             for (iw=0; iw<Nwall_type; iw++) eps_wf_save[0][iw]=Eps_wf[0][iw];
                             pot_parameters("dft_out.lis"); 
                             for (iw=0; iw<Nwall; iw++){
                                 ratio = Eps_wf[0][WallType[iw]]/eps_wf_save[0][WallType[iw]];
                                 scale_vext_epswf(ratio,0,iw); 
                             }
                         }
                         if (Type_poly >=0) setup_polymer_cr();
                         recalculate_stencils();
                         break;

      case CONT_EPSFF_ALL: 
                         for (i=0; i<Ncomp; i++) 
                            for (j=0; j<Ncomp; j++) if (fabs(Eps_ff[i][j])>1.e-15) Eps_ff[i][j] = param;
                      
                         if (Mix_type==0) {
                             for (i=0; i<Ncomp; i++){ 
                              for (iw=0; iw<Nwall_type; iw++) eps_wf_save[i][iw]=Eps_wf[i][iw];
                             }
                             pot_parameters("dft_out.lis"); 
                             for (i=0; i<Ncomp; i++){
                               for (iw=0; iw<Nwall; iw++){
                             
                                ratio = Eps_wf[i][WallType[iw]]/eps_wf_save[i][WallType[iw]];
                                scale_vext_epswf(ratio,i,iw); 
                             }
                           }
                         }
                         if (Type_poly >=0) setup_polymer_cr();
                         recalculate_stencils();
                         break;

      case CONT_SCALE_EPSFF:
                         scale_save = Scale_fac;
                         Scale_fac = param;
                         ratio = Scale_fac/scale_save;
                         for (i=0; i<Ncomp; i++) Eps_ff[i][i] *= ratio;
                         if (Mix_type==0) {
                             for (i=0; i<Ncomp; i++){ 
                              for (iw=0; iw<Nwall_type; iw++) eps_wf_save[i][iw]=Eps_wf[i][iw];
                             }
                             pot_parameters("dft_out.lis"); 
                             for (i=0; i<Ncomp; i++){
                               for (iw=0; iw<Nwall; iw++){
                             
                                ratio = Eps_wf[i][WallType[iw]]/eps_wf_save[i][WallType[iw]];
                                scale_vext_epswf(ratio,i,iw); 
                             }
                           }
                         }
                         if (Type_poly >=0) setup_polymer_cr();
                         recalculate_stencils();
                         break;

      case CONT_SCALE_CHG: scale_save=Scale_fac;
                           Scale_fac=param;
                           ratio = Scale_fac/scale_save;
                           scale_elec_param(ratio); 
                           break;

      case CONT_SEMIPERM:
           param_save=Vext_membrane[0][0];
           for (inode=0;inode<Nnodes_per_proc;inode++){
                  if (fabs(Vext[inode][0]-param_save)<1.e-10) Vext[inode][0]=param;
           }
           Vext_membrane[0][0]=param;
           break;

      case CONT_WALLPARAM:
           WallParam[1]=param;
           WallPos[0][1] = 0.5*Size_x[0]-WallParam[WallType[1]]; 
           WallPos[0][0] = 0.5*Size_x[0]-2.*WallParam[WallType[1]]-WallParam[WallType[0]]; 
           free_mesh_arrays();
           boundary_free();
           if (Lsteady_state && Ndim==1) safe_free((void *) &Area_IC);
           safe_free((void *) &Comm_node_proc);
           safe_free((void *) &Comm_unk_proc);
           safe_free((void *) &Comm_offset_node);
           safe_free((void *) &Comm_offset_unk);
           safe_free((void *) &Aztec.update);
           output_file2 = "dft_vext.dat";
           output_TF = "dft_zeroTF.dat";
           set_up_mesh(output_file1,output_file2);
           boundary_setup(output_file1);
           if (Iwrite==VERBOSE) {
                 print_vext(Vext,output_file2);
                 print_zeroTF(Zero_density_TF,output_TF);
          }

           break;

     case CONT_CRFAC:
         Crfac=param;
         setup_polymer_cr();
         recalculate_stencils();
         break;

      default:
        printf("ERROR_apt: Unknown Continuation parameter %d\n",cont_type);
        exit(-1); break;
  }

  /* recalculate thermo based on new parameter */
  thermodynamics(output_file1, TRUE);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void assign_parameter_conwrap(double param)
/* Put the call to a routine to assign the continuation parameter here.
 * Input:
 *    param     New value of continuation parameter.
 *
 * Output:
 *
 * Return Value:
 */
{
  if (Proc==0) printf("\tContinuation parameter #%d set to %g\n",
                         Loca.cont_type1, param);
  assign_parameter_tramonto(Loca.cont_type1, param);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void assign_bif_parameter_conwrap(double tp_param)
/* Put the call to a routine to assign the continuation parameter here.
 * Input:
 *    tp_param     New value of continuation parameter.
 *
 * Output:
 *
 * Return Value:
 */
{
  if (Proc==0) printf("\tSecond (floating) parameter #%d set to %g\n",
                         Loca.cont_type2, tp_param);
  assign_parameter_tramonto(Loca.cont_type2, tp_param);

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static double get_init_param_value(int cont_type)
{
  int i,j; 
  switch(cont_type){
      case CONT_MESH: 
       printf("ERROR: Continuation Library cannot do mesh changes\n");
       exit(-1); break;

      case CONT_TEMP: return Temp;

      case CONT_RHO_0:   return Rho_b[2];
      case CONT_RHO_ALL:   
             for (i=0;i<Ncomp;i++) if (Rho_b[i] != Rho_b[0]){
                 printf("ERROR: need all Rho_b to be the same for CONT_RHO_ALL\n"); 
                 exit(-1);
             }
             return Rho_b[0];

      case CONT_LOG_RHO_0: return log(Rho_b[0]);

      case CONT_LOG_RHO_ALL:   
             for (i=0;i<Ncomp;i++) if (Rho_b[i] != Rho_b[0]){
                 printf("ERROR: need all Rho_b to be the same for CONT_LOG_RHO_ALL\n"); 
                 exit(-1);
             }
             return log(Rho_b[0]);

      case CONT_SCALE_RHO: return Scale_fac;

      case CONT_EPSW_0:   if (Mix_type==0) return Eps_w[0];
                          else             return Eps_ww[0][0];

      case CONT_EPSW_ALL: 
             if (Mix_type ==0 ){
                for (i=0;i<Nwall_type;i++) if (Eps_w[i] != Eps_w[0]){
                    printf("ERROR: need all Eps_w to be the same for CONT_EPS_W_ALLL\n"); 
                    exit(-1);
                }
                return Eps_w[0];
             }
             else{
                for (i=0;i<Nwall_type;i++){
                   for (j=0; j<Nwall_type;j++)  if (Eps_ww[i][j] != Eps_ww[0][0]){
                      printf("ERROR: need all Eps_ww to be the same for CONT_EPS_W_ALLL Mix_type=1\n"); 
                   }
                }
                return Eps_ww[0][0];
             }
      case CONT_SCALE_EPSW: return Scale_fac;

      case CONT_EPSWF00: /*return Eps_wf[0][0];*/ return Eps_wf[2][0];

      case CONT_EPSWF_ALL_0: 
           for (i=0;i<Ncomp-1;i++) if (Eps_wf[i][0] != Eps_wf[0][0]) {
                 printf("ERROR: all Eps_wf must be equal for CONT_EPSWF_ALL_0\n");
                 exit(-1);
           }
           return Eps_wf[0][0];
      case CONT_SCALE_EPSWF: return Scale_fac;

      case CONT_EPSFF_00:   return Eps_ff[0][0]; /*Eps_ff[0][2];*/  /*Eps_ff[2][2];*/
      case CONT_EPSFF_ALL:
           for (i=0; i<Ncomp; i++) if (Eps_ff[i][i] != Eps_ff[0][0]) {
                 printf("ERROR: need all Eps_ff[i][i] to be equal for  (CONT_EPSFF_ALL)\n");
                 exit(-1);
           }
           return Eps_ff[0][0];
      case CONT_SCALE_EPSFF: return Scale_fac;

      case CONT_SCALE_CHG:  return Scale_fac;

      case CONT_SEMIPERM: return Vext_membrane[0][0];

      case CONT_WALLPARAM: return WallParam[WallType[1]];

      case CONT_CRFAC:
              return Crfac;

      default:
        printf("ERROR: Unknown Continuation parameter %d\n",cont_type);
        exit(-1); break;
  }
  return 0.0;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void calc_scale_vec_conwrap(double *x, double *scale_vec, int numUnks)
/* Put the call to a routine to calculate a scaling vector here.
 * Input:
 *    x          New value of continuation parameter.
 *    numUnks    Number of unknowns on this proc, the length of x
 *               and scale_vec.
 *
 * Output:
 *    scale_vec  Vector of length number of unknowns used to scale
 *               variables so that one type of unknown (e.g. pressure)
 *               doesn't dominate over others. Used to balance the
 *               variables and the arc-length variable in arc-length
 *               continuation, and for scaling the null vector in
 *               turning point tracking. Using reciprocal of the average
 *               value of that variable type is a good choice. Vector
 *               of all ones should suffice for most problems.
 *
 *               For tramonto, I've seen better behavior with vector 
 *               of all ones instead of using Rho_b.
 *
 * Return Value:
 */
{
  int      i;

  /* This was the default -- changed for spinodal tracking */
  for (i=0; i<numUnks; i++) scale_vec[i] = 1.0;
  /*for (i=0; i<numUnks; i++) scale_vec[i] = 1.0/ (sqrt(fabs(x[i])) + 1.0e-5);*/

/* Scale vector using Rho_bulk seems to make things worse
  for (i=0; i<Aztec.N_update/Nunk_per_node; i++) 
    for (j=0; j<Ncomp; j++) 
      scale_vec[Aztec.update_index[i*Nunk_per_node + Phys2Unk_first[DENSITY]+j]] = 1.0/Rho_b[j];
*/

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
double gsum_double_conwrap(double sum)
/* Put the call to a routine to calculate a global sum.
 * Just return sum for single processor jobs.
 * Input:
 *    sum     Value of double on this processor to be summed on all procs.
 *
 * Output:
 *
 * Return Value:
 *    The global sum is returned on all processors.
 */
{
  return AZ_gsum_double(sum, Aztec.proc_config);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int gmax_int_conwrap(int max)
/* Put the call to a routine to calculate a global sum.
 * Just return sum for single processor jobs.
 * Input:
 *    max     Value of integer on this processor to be maxed on all procs.
 *
 * Output:
 *
 * Return Value:
 *    The global max is returned on all processors.
 *
 * Only used by Eigensolver
 */
{
  return AZ_gmax_int(max, Aztec.proc_config);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void random_vector_conwrap(double *x, int numOwnedUnks)
/* Put a routine to calculate a random vector.
 * Input:
 *    numOwnedUnks  Length of owned nodes part of x.
 *
 * Output:
 *    x             Random vector.
 *
 * Used by eigensolver only
 */
{
  AZ_random_vector(x, Aztec.data_org, Aztec.proc_config);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void perturb_solution_conwrap(double *x, double *x_old,
		              double *scale_vec, int numOwnedUnks)
/* Put a routine to perturb the solution vector a little bit here.
 * This is to move a converged solution at a singularity off
 * of the singularity before doing continuation. This ain't pretty
 * but has helped convergence on some turning point tracking problems.
 * Input:
 *    x_old         Current solution vector.
 *    scale_vec     Work space for a vector to scale x.
 *    numOwnedUnks  Length of owned nodes part of x, x_old, scale_vec
 *
 * Output:
 *    x             Solution vector perturbed a bit.
 *    
 * Return Value:
 */
{
  int i;

  random_vector_conwrap(x, numOwnedUnks);

  for (i=0; i<numOwnedUnks; i++) x[i] = x_old[i] * (1.0 + x[i]*1.0e-5);

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void solution_output_conwrap(int num_soln_flag, double *x, double param,
                             double *x2, double param2,
                             double *x3, double param3,
                             int step_num, int num_its, struct con_struct *con)

/* Put the call to your solution output (both file and screen) routines here.
 * Input:
 *    num_soln_flag  Flag for how many solution vectors are being passed for
 *                   output. For parameter continuation and turning point
 *                   tracking there is just 1, for pitchfork and phase
 *                   transitions there are 2, and for Hopfs there are 3.
 *                   The eigensolver sends -1 for a single eigenvector and
 *                   -2 fir a complex pair of eigenvectors.
 *    x            First solution vector for output (x or first eigen.
 *    param        Continuation parameter value (Real part ev when flag < 0)
 *    x2           Second solution vector for output (y_vec or x2)
 *    param2       Bifurcation parameter value (imaginary part ev if flag=-2)
 *    x3           Third solution vector for output (z_vec for Hopfs)
 *    param3       Third Parameter value (frequency Hopfs)
 *    step_num+1   Time index to output to (step_num is 0 based).
 *    num_its      Number of Newton iterations used for for convergence
 *    con          pointer to continuation structure, for passing to
 *                 the eigensolver.
 *
 * Output:
 *
 * Return Value:
 */
{
  char *output_file3 = "dft_output.dat";
  double time_save=0, *x_internal;
  int i;

  x_internal = (double *) array_alloc(1, Aztec.N_update, sizeof(double));

  for (i=0; i<Aztec.N_update; i++) x_internal[i] = x[Aztec.update_index[i]];

  post_process(x_internal, output_file3, &num_its, &time_save, 
               step_num, FALSE);

  /* Print second solution for Phase Trans, but not Spinodal tracking */

  if (num_soln_flag == 2 
     && con->general_info.method == PHASE_TRANSITION_CONTINUATION) {

    for (i=0; i<Aztec.N_update; i++) x_internal[i] = x2[Aztec.update_index[i]];

    post_process(x_internal, output_file3, &num_its, &time_save, 
                 step_num, TRUE);
  }

  safe_free((void *) &x_internal);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
double free_energy_diff_conwrap(double *x, double *x2)
/* Call to return the free energy difference betwen two solutions
 * Input:
 *    x    One solution vector
 *    x2   Second solution vector
 *
 * Output:
 *
 * Return Value:
 *    The difference in the free energy beween the two solutions
 */
{
  double energy1, energy2;

  setup_integrals();

  if (!Sten_Type[POLYMER_CR]){
    energy1 = calc_free_energy(NULL,x,1.0,1.0,FALSE);
    energy2 = calc_free_energy(NULL,x2,1.0,1.0,FALSE);
  }
  else {
    calc_adsorption(NULL,x,1.0,1.0);
    energy1=calc_free_energy_polymer(NULL,x,1.0,1.0);
    calc_adsorption(NULL,x2,1.0,1.0);
    energy2=calc_free_energy_polymer(NULL,x2,1.0,1.0);
  }

  safe_free((void *)&Nel_hit);
  safe_free((void *)&Nel_hit2);

  if (Proc==0)  printf("In energy calc: e1,e2,e1-e2 = %g %g %g\n",
        energy1, energy2, energy1-energy2);
  

  return energy1-energy2;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void print_con_struct(const struct con_struct* con)
/* Routine for printing out the con structure to the screen */

{
  printf("\n"); /* print_line("~", 80); */
  printf("\tcon->general_info.param=             %10.4g\n",
    con->general_info.param);
  printf("\tcon->general_info.numUnks=           %10d\n",
    con->general_info.numUnks);
  printf("\tcon->general_info.numOwnedUnks=      %10d\n",
    con->general_info.numOwnedUnks);
  printf("\tcon->general_info.printproc=         %10d\n",
    con->general_info.printproc);

  printf("\tcon->stepping_info.first_step=       %10.4g\n",
    con->stepping_info.first_step);
  printf("\tcon->stepping_info.max_steps=        %10d\n",
    con->stepping_info.max_steps);
  printf("\tcon->stepping_info.max_param=        %10.4g\n",
    con->stepping_info.max_param);
  printf("\tcon->stepping_info.max_delta_p=      %10.4g\n",
    con->stepping_info.max_delta_p);
  printf("\tcon->stepping_info.step_ctrl=        %10.4g\n",
    con->stepping_info.step_ctrl);
  printf("\tcon->stepping_info.max_newton_its=   %10d\n",
    con->stepping_info.max_newton_its);

  if(con->general_info.method==ARC_LENGTH_CONTINUATION) {
    printf("\tcon->arclength_info.dp_ds2_goal=     %10.4g\n",
      con->arclength_info.dp_ds2_goal);
    printf("\tcon->arclength_info.dp_ds_max=       %10.4g\n",
      con->arclength_info.dp_ds_max);
    printf("\tcon->arclength_info.tang_exp=        %10.4g\n",
      con->arclength_info.tang_exp);
    printf("\tcon->arclength_info.tang_step_limit= %10.4g\n",
      con->arclength_info.tang_step_limit);
  }

  if(con->general_info.method==TURNING_POINT_CONTINUATION)
    printf("\tcon->turning_point_info.bif_param=   %10.4g\n",
      con->turning_point_info.bif_param);

  if(con->general_info.method==PITCHFORK_CONTINUATION)
    printf("\tcon->pitchfork_info.bif_param=       %10.4g\n",
      con->pitchfork_info.bif_param);

  if(con->general_info.method==HOPF_CONTINUATION)
    printf("\tcon->hopf_info.bif_param=            %10.4g\n",
      con->hopf_info.bif_param);

  if(con->eigen_info.Num_Eigenvalues > 0) {
    printf("\tcon->eigen_info.Shift_Point[0]=      %10.4g\n",
      con->eigen_info.Shift_Point[0]);
    printf("\tcon->eigen_info.Shift_Point[1]=      %10.4g\n",
      con->eigen_info.Shift_Point[1]);
    printf("\tcon->eigen_info.Arnoldi=             %10d\n",
      con->eigen_info.Arnoldi);
    printf("\tcon->eigen_info.Residual_Tol[0]=     %10.4g\n",
      con->eigen_info.Residual_Tol[0]);
    printf("\tcon->eigen_info.Every_n_Steps=       %10d\n",
      con->eigen_info.Every_n_Steps);
  }
  /* print_line("~", 80); */ printf("\n");
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static void print_final(double param, int step_num)

/*
 * Print out the final results and counters
 */

{
  printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  printf("CONTINUATION ROUTINE HAS FINISHED: \n");
  printf("\tEnding Parameter value     = %g\n", param);
  printf("\tNnumber of steps           = %d\n", step_num);
  printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

} /************* END of print_final () ***************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
