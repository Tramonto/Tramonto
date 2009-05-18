/* This file was automatically generated.  Do not edit! */
#include <stdio.h>
int continuation_hook(double *x,double *delta_x,void *con_void,double reltol,double abstol);
int continuation_hook_conwrap(double **xx,double **delta_xx,void *con_ptr,double reltol,double abstol);
void box2owned(double **xBox,double **xOwned);
void box2owned(double **xBox,double **xOwned);
void setup_integrals();
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if defined(LOCA_MF) && !defined(_CON_CONST_H_)
#include <mf.h>
#endif
double free_energy_diff_conwrap(double *x,double *x2);
#if !defined(_CON_CONST_H_)
double free_energy_diff_conwrap(double *x,double *x2);
#endif
double calc_free_energy(FILE *fp,double **x);
double calc_free_energy_conwrap(double **xB);
#include <stdlib.h>
#include <math.h>
#if defined(HAS_VALUES_H)
#include <values.h>
#include <unistd.h>
#include <string.h>
#endif
#include "mpi.h"
#include "az_aztec.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define FROM_LOCA 0
void post_process(double **x,int *niters,double *time_save,int loop1,int binodal_flag,int call_from_flag);
#if !defined(_CON_CONST_H_)
typedef struct con_struct con_struct;
#endif
void solution_output_conwrap(int num_soln_flag,double *x,double param,double *x2,double param2,double *x3,double param3,int step_num,int num_its,struct con_struct *con);
#if !defined(_CON_CONST_H_)
void solution_output_conwrap(int type,double *x,double param,double *x2,double param2,double *x3,double param3,int stp,int nits,struct con_struct *con);
#endif
void perturb_solution_conwrap(double *x,double *x_old,double *scale_vec,int numOwnedUnks);
#if !defined(_CON_CONST_H_)
void perturb_solution_conwrap(double *x,double *p,double *s,int n);
#endif
void random_vector_conwrap(double *x,int numOwnedUnks);
#if !defined(_CON_CONST_H_)
void random_vector_conwrap(double *x,int n);
#endif
int gmax_int(int c);
int gmax_int_conwrap(int max);
#if !defined(_CON_CONST_H_)
int gmax_int_conwrap(int sum);
#endif
void *array_alloc_2d(size_t n1,size_t n2,size_t size);
void *array_alloc_2d(size_t n1,size_t n2,size_t size);
double **array_alloc_2d_conwrap(unsigned int ii,unsigned int jj,unsigned int kk);
void safe_free_conwrap(void **p);
void fill_resid_and_matrix_control_conwrap(double **xBox,int ii,int jj);
double gmax_double(double c);
double gmax_double_conwrap(double sum);
double gsum_double(double c);
void calc_scale_vec_conwrap(double *x,double *scale_vec,int numUnks);
#if !defined(_CON_CONST_H_)
void calc_scale_vec_conwrap(double *x,double *scale_vec,int numUnks);
#endif
void assign_bif_parameter_conwrap(double tp_param);
#if !defined(_CON_CONST_H_)
void assign_bif_parameter_conwrap(double bif_param);
#endif
void assign_parameter_tramonto(int cont_type,double param,int Loca_contID);
void assign_parameter_conwrap(double param);
#if !defined(_CON_CONST_H_)
void assign_parameter_conwrap(double param);
#endif
void matvec_mult_conwrap(double *x,double *y);
#if !defined(_CON_CONST_H_)
void matvec_mult_conwrap(double *x,double *y);
#endif
double gsum_double_conwrap(double sum);
double gsum_double_conwrap(double sum);
#if !defined(_CON_CONST_H_)
double gsum_double_conwrap(double sum);
#endif
double fill_resid_and_matrix_control(double **x,int iter,int resid_only_flag);
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
#if !defined(_CON_CONST_H_)
#define  RHS_ONLY        100
#endif
void matrix_residual_fill_conwrap(double *x,double *rhs,int matflag);
#if !defined(_CON_CONST_H_)
void matrix_residual_fill_conwrap(double *x,double *rhs,int matflag);
#endif
void destroy_shifted_matrix_conwrap();
#if !defined(_CON_CONST_H_)
void destroy_shifted_matrix_conwrap();
#endif
void shifted_linear_solver_conwrap(double *x,double *y,int jac_flag,double tol);
#if !defined(_CON_CONST_H_)
void shifted_linear_solver_conwrap(double *x,double *y,int jac_flag,double tol);
#endif
void shifted_matrix_fill_conwrap(double sigma);
#if !defined(_CON_CONST_H_)
void shifted_matrix_fill_conwrap(double sigma);
#endif
void create_shifted_matrix_conwrap();
#if !defined(_CON_CONST_H_)
void create_shifted_matrix_conwrap();
#endif
void mass_matvec_mult_conwrap(double *x,double *y);
#if !defined(_CON_CONST_H_)
void mass_matvec_mult_conwrap(double *x,double *y);
#endif
void mass_matrix_fill_conwrap(double *x,double *rhs);
#if !defined(_CON_CONST_H_)
void mass_matrix_fill_conwrap(double *x,double *rhs);
#endif
int komplex_linear_solver_conwrap(double *c,double *d,int jac_flag,double *omega,double *tmp);
#if !defined(_CON_CONST_H_)
int komplex_linear_solver_conwrap(double *x,double *y,int jac_flag,double *omega,double *tmp);
#define  NEW_JACOBIAN                200
#define  SAME_BUT_UNSCALED_JACOBIAN  202
#define  CHECK_JACOBIAN              204
#define  OLD_JACOBIAN                201
#endif
extern void *LinProbMgr_manager;
int linear_solver_conwrap(double *x,int jac_flag,double *tmp);
#if !defined(_CON_CONST_H_)
int linear_solver_conwrap(double *x,int jac_flag,double *tmp);
#endif
int newton_solver(double **x,void *con_ptr);
int newton_solver(double **x,void *con_ptr);
int nonlinear_solver_conwrap(double *x,void *con_ptr,int step_num,double lambda,double delta_s,void *aux_info);
#if !defined(_CON_CONST_H_)
int nonlinear_solver_conwrap(double *x,void *con,int step_num,double lambda,double delta_s,void *aux_info);
#endif
void safe_free(void **ptr);
void safe_free(void **ptr);
#if !defined(_CON_CONST_H_)
int con_lib(struct con_struct *con,void *aux_info);
#endif
int con_lib(struct con_struct *con,void *aux_info);
#if !defined(_CON_CONST_H_)
#define HOPF_CONTINUATION              5
#define PITCHFORK_CONTINUATION         4
#endif
extern int Max_NL_iter;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define VERBOSE      3 
#define DENSITIES    1 
#define MINIMAL      0
#define NO_SCREEN    2 
extern int Iwrite;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
#if !defined(_CON_CONST_H_)
#define PHASE_TRANSITION_CONTINUATION  6
#define TURNING_POINT_CONTINUATION     3
#define ARC_LENGTH_CONTINUATION        2
#define FIRST_ORDER_CONTINUATION       1
#define ZERO_ORDER_CONTINUATION        0
#endif
typedef struct Loca_Struct Loca_Struct;
struct Loca_Struct {
  int    method;      /* Continuation method                          */
  int    cont_type1;  /* flag specifying the continuation parameter   */
  int    cont_type2;  /* flag specifying the second (free) parameter  */
  int    num_steps;   /* maximum number of continuation steps to take */
  double aggr;        /* step size control parameter                  */
  double step_size;   /* initial continuation step size               */
};
extern struct Loca_Struct Loca;
extern int Nnodes_box;
extern int Nnodes_per_proc;
extern int Nunk_per_node;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
int solve_continuation(double **xx,double **xx2);
double get_init_param_value(int cont_type,int);
#if !defined(_CON_CONST_H_)
double get_init_param_value(int cont_type,int);
#endif
double get_init_param_value(int cont_type,int Loca_contID);
#if !defined(_CON_CONST_H_)
typedef struct general_info_struct general_info_struct;
struct general_info_struct {
   int    method;       /* Solution strategy: see method choices above      */
   double param;        /* value of continuation parameter                  */
   double *x;           /* current solution vector                          */
   double perturb;      /* Perturbation magnitude                           */
   int    numUnks;      /* Num unknowns on this Proc (including externals)  */
   int    numOwnedUnks; /* Num unknowns on this Proc (NOT incl. externals)  */
   int    printproc;    /* Logical indicating if this Proc prints to stdout */
   int    nv_restart;   /* Restarted null vector flag                       */
   int    nv_save;      /* Null vector save flag                            */
};
typedef struct stepping_info_struct stepping_info_struct;
struct stepping_info_struct {
   double first_step;  /* Initial step size                                 */
   int    base_step;   /* Number of the first step (for output)             */
   int    max_steps;   /* Maximum # of continuation steps                   */
   int    last_step;   /* Last step flag                                    */
   double max_param;   /* parameter value to stop at                        */
   double max_delta_p; /* continuation parameter step limit                 */
   double min_delta_p; /* continuation parameter step limit                 */
   double step_ctrl;   /* step aggressiveness -- 0.0 for constant step      */
   int    max_newton_its;/* Max # Newton steps, used only for step control  */
};
typedef struct arclength_info_struct arclength_info_struct;
struct arclength_info_struct {
   double dp_ds2_goal; /* Square of target dp_ds value for rescaling
                          (desired solution contribution to arc length)     */
   double dp_ds_max;   /* High dp_ds value at which to rescale              */
   double tang_exp;    /* Power to which tang_factor is raised              */
   double tang_step_limit; /* Minimum value of tang_factor between steps    */
};
typedef struct turning_point_info_struct turning_point_info_struct;
struct turning_point_info_struct {
   double bif_param;          /* Initial guess of bifurcation parameter     */
   double *nv;                /* Restarted null vector (read in from file)  */
};
typedef struct pitchfork_info_struct pitchfork_info_struct;
struct pitchfork_info_struct {
   double bif_param;          /* Initial guess of bifurcation parameter     */
   double *psi;               /* Antisymmetry vector (also init guess for   */
                              /* the null vector, y_vec                     */
};
typedef struct hopf_info_struct hopf_info_struct;
struct hopf_info_struct {
   double bif_param;          /* Initial guess of bifurcation parameter     */
   double omega;              /* Initial guess of Hopf frequency            */
   double *y_vec;             /* Initial guess of null vector (real)        */
   double *z_vec;             /* Initial guess of null vector (imag)        */
   int    mass_flag;
};
typedef struct phase_transition_info_struct phase_transition_info_struct;
struct phase_transition_info_struct {
   double bif_param;          /* Initial guess of bifurcation parameter     */
   double *x2;                /* Initial guess of second_solution_vector    */
};
typedef struct manifold_info_struct manifold_info_struct;
#endif
#if !(defined(LOCA_MF)) && !defined(_CON_CONST_H_)
#define  MFNKMatrix int
#define  MFNVector  int
#define  MFNSpace   int
#endif
#if !defined(_CON_CONST_H_)
struct manifold_info_struct {
   int k;                     /* Manifold Dimension */
   double *param_vec;         /* Parameter Vector */
   double *param_lower_bounds;/* Parameter lower bounds */
   double *param_upper_bounds;/* Parameter upper bounds */
   double *param_steps;       /* Parameter step guess */
   MFNKMatrix phi;            /* Matrix of null vectors */
   MFNVector u0;              /* Previous solution */
   MFNVector u;               /* Current solution */
   MFNSpace space;            /* Pointer to MF's LOCA NSpace */
};
typedef struct eigen_info_struct eigen_info_struct;
struct eigen_info_struct {
   int    Num_Eigenvalues;    /* Number of Eigenvalues to Calculate          */
   int    Num_Eigenvectors;   /* Number of Eigenvectors to Write             */
   int    sort;               /* Flag to sort eigenvalues by real part       */
   double Shift_Point[3];     /* Point for shift and invert (Real, Imaginary)*/
                              /* and for cayley delta                        */
   int    Arnoldi;            /* Arnoldi Space Size                          */
   double Residual_Tol[2];    /* Convergence Tolerance for the Eigenvalue    *
                               * Residual Equation, and linear solver tol    */
   int    Max_Iter;           /* Maximum number of iterations of eigensolver */
   int    Every_n_Steps;      /* Allow for eigenvalue calc every n steps along
                                 a continuation run                          */
};
typedef struct private_info_struct private_info_struct;
struct private_info_struct {
   int mass_x;       /* flag that turns on dM/dx in komplex solves           */
   int mass_param;   /* flag that turns on dM/d(param) in komplex solves     */
   int first_iter;   /* flag for first Newton iter of each solve             */
   int step_num;     /* Current continuation step number                     */
   int nstep;        /* Current step number (for output)                     */
   double param_old; /* old value of continuation parameter                  */
   double arc_step;  /* step size of arc length variable                     */
   double dp_ds;     /* derivative of parameter w.r.t. arc length            */
   double *x_old;    /* previous solution vector                             */
   double *x_tang;   /* tangent to the solution vector w.r.t. parameter      */
   double *scale_vec;/* scaling vector for better arclength conditioning     */
};
struct con_struct {
   struct general_info_struct general_info;
   struct stepping_info_struct stepping_info;
   struct arclength_info_struct arclength_info;
   struct turning_point_info_struct turning_point_info;
   struct pitchfork_info_struct pitchfork_info;
   struct hopf_info_struct hopf_info;
   struct phase_transition_info_struct phase_transition_info;
   struct manifold_info_struct manifold_info;
   struct eigen_info_struct eigen_info;

   struct private_info_struct private_info;
};
#endif
