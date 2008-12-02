/* This file was automatically generated.  Do not edit! */
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(_CON_CONST_H_)
#define MANIFOLD_CONTINUATION          8
#endif
#if !defined(LOCA_LSA_ONLY) && !defined(_CON_CONST_H_)
#define LOCA_LSA_ONLY                  9
#endif
#if !defined(_CON_CONST_H_)
#define FIRST_ORDER_CONTINUATION       1
#define ZERO_ORDER_CONTINUATION        0
#define PHASE_TRANSITION_CONTINUATION  6
#define HOPF_CONTINUATION              5
#define PITCHFORK_CONTINUATION         4
#define TURNING_POINT_CONTINUATION     3
#define ARC_LENGTH_CONTINUATION        2
#endif
int continuation_hook(double *x,double *delta_x,void *con_void,double reltol,double abstol);
#if !defined(_CON_CONST_H_)
#define  MATRIX_ONLY     101
#define  RECOVER_MATRIX  104
#define  RHS_MATRIX_SAVE 103
#define  CONT_TANGENT      300
#endif
void vec_init(double *u);
void vec_init(double *u);
double scaled_dp(double *x,double *y);
double scaled_dp(double *x,double *y);
#if defined(LOCA_MF) && !defined(_CON_CONST_H_)
#include <mf.h>
#endif
#if defined(LOCA_MF)
void calc_rhs_multi(double *param_vec,int i,int k,double *phix,double *rhs,double *x,double perturb,int numUnks);
#endif
#if !defined(_CON_CONST_H_)
void calc_rhs_multi(double *param_vec,int i,int k,double *phix,double *rhs,double *x,double perturb,int numUnks);
#define  RHS_ONLY        100
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
typedef struct con_struct con_struct;
#endif
int manifold_alg(double *x,double *delta_x,struct con_struct *con,double reltol,double abstol);
#if !defined(_CON_CONST_H_)
#define  RHS_MATRIX      102
#endif
void matrix_residual_fill_conwrap(double *x,double *rhs,int matflag);
#if !defined(_CON_CONST_H_)
void matrix_residual_fill_conwrap(double *x,double *rhs,int matflag);
#endif
double free_energy_diff_conwrap(double *x,double *x2);
#if !defined(_CON_CONST_H_)
double free_energy_diff_conwrap(double *x,double *x2);
#endif
int phase_transition_alg(double *x,double *delta_x,struct con_struct *con,double reltol,double abstol);
#if !defined(_CON_CONST_H_)
#define  OLD_JACOBIAN_DESTROY        203
#define  HP_CONT_DMDPARAM  307
#define  HP_CONT_DMDX      306
#define  NEW_JACOBIAN                200
#endif
int komplex_linear_solver_conwrap(double *c,double *d,int jac_flag,double *omega,double *tmp);
#if !defined(_CON_CONST_H_)
int komplex_linear_solver_conwrap(double *x,double *y,int jac_flag,double *omega,double *tmp);
#define  HP_CONT_SOL3      305
#endif
void matvec_mult_conwrap(double *x,double *y);
#if !defined(_CON_CONST_H_)
void matvec_mult_conwrap(double *x,double *y);
#endif
void mass_matvec_mult_conwrap(double *x,double *y);
#if !defined(_CON_CONST_H_)
void mass_matvec_mult_conwrap(double *x,double *y);
#endif
void mass_matrix_fill_conwrap(double *x,double *rhs);
#if !defined(_CON_CONST_H_)
void mass_matrix_fill_conwrap(double *x,double *rhs);
#endif
double dp(double *x,double *y);
double dp(double *x,double *y);
#if !defined(_CON_CONST_H_)
typedef struct hopf_info_struct hopf_info_struct;
struct hopf_info_struct {
   double bif_param;          /* Initial guess of bifurcation parameter     */
   double omega;              /* Initial guess of Hopf frequency            */
   double *y_vec;             /* Initial guess of null vector (real)        */
   double *z_vec;             /* Initial guess of null vector (imag)        */
   int    mass_flag;
};
#endif
int hopf_alg(double *x,double *delta_x,struct con_struct *con,double reltol,double abstol);
#if !defined(_CON_CONST_H_)
#define  OLD_JACOBIAN                201
#endif
void vec_copy(double *dx,double *dy);
void vec_copy(double *dx,double *dy);
double ip(double *x,double *y);
double ip(double *x,double *y);
int pitchfork_alg(double *x,double *delta_x,struct con_struct *con,double reltol,double abstol);
double gsum_double_conwrap(double sum);
double gsum_double_conwrap(double sum);
#if !defined(_CON_CONST_H_)
double gsum_double_conwrap(double sum);
#endif
void assign_bif_parameter_conwrap(double tp_param);
#if !defined(_CON_CONST_H_)
void assign_bif_parameter_conwrap(double bif_param);
#define  TP_CONT_SOL4      304
#define  SAME_BUT_UNSCALED_JACOBIAN  202
#endif
double null_vector_resid(double r_val,double i_val,double *r_vec,double *i_vec,int mm_flag);
double null_vector_resid(double r_val,double i_val,double *r_vec,double *i_vec,int mm_flag);
#if !defined(_CON_CONST_H_)
#define  TP_CONT_SOL3      303
#endif
double ltransnorm(double *x,double *scale_vec);
double ltransnorm(double *x,double *y);
void calc_scale_vec_conwrap(double *x,double *scale_vec,int numUnks);
#if !defined(_CON_CONST_H_)
void calc_scale_vec_conwrap(double *x,double *scale_vec,int numUnks);
#define  TP_CONT_SOL2      302
#endif
int turning_point_alg(double *x,double *delta_x,struct con_struct *con,double reltol,double abstol);
#if !defined(_CON_CONST_H_)
int turning_point_alg(double *x,double *delta_x,struct con_struct *con,double reltol,double abstol);
#endif
#include <stdio.h>
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
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
void free_vec(double **ptr);
void free_vec(double **ptr);
void assign_parameter_conwrap(double param);
#if !defined(_CON_CONST_H_)
void assign_parameter_conwrap(double param);
#define  CHECK_JACOBIAN              204
#endif
int linear_solver_conwrap(double *x,int jac_flag,double *tmp);
#if !defined(_CON_CONST_H_)
int linear_solver_conwrap(double *x,int jac_flag,double *tmp);
#define  ARC_CONT_SOL2     301
#endif
void calc_rhs_continuation(int rhs_type,double *x,double *resid_vector,double *ab_vec,double *scale_vec,double *x_tmp,double param,double perturb,double *r_vec,int numUnks,int numOwnedUnks);
#if !defined(_CON_CONST_H_)
void calc_rhs_continuation(int rhs_type,double *x,double *a,double *x_dot,double *scale_vec,double *x_tmp,double con_par,double perturb,double *r_vec,int num_total_unknowns,int num_owned_unks);
#endif
double *alloc_vec();
double *alloc_vec();
#if !defined(_CON_CONST_H_)
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
typedef struct phase_transition_info_struct phase_transition_info_struct;
struct phase_transition_info_struct {
   double bif_param;          /* Initial guess of bifurcation parameter     */
   double *x2;                /* Initial guess of second_solution_vector    */
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
int arc_length_bordering_alg(double *x,double *delta_x,struct con_struct *con,double reltol,double abstol);
#if !defined(_CON_CONST_H_)
int arc_length_bordering_alg(double *x,double *delta_x,struct con_struct *con,double reltol,double abstol);
#endif
extern int AGS_option;
double scaled_dot_prod(double *,double *,double *,int);
double scaled_dot_prod(double *x,double *y,double *scale_vec,int n);
double scaled_dot_prod(double *x,double *y,double *sc,int n);
