/* This file was automatically generated.  Do not edit! */
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if defined(LOCA_MF) && !defined(_CON_CONST_H_)
#include <mf.h>
#endif
double gsum_double_conwrap(double sum);
double gsum_double_conwrap(double sum);
#if !defined(_CON_CONST_H_)
double gsum_double_conwrap(double sum);
#endif
void calc_scale_vec_conwrap(double *x,double *scale_vec,int numUnks);
#if !defined(_CON_CONST_H_)
void calc_scale_vec_conwrap(double *x,double *scale_vec,int numUnks);
#endif
void free_vec(double **ptr);
void free_vec(double **ptr);
#if !defined(_CON_CONST_H_)
typedef struct con_struct con_struct;
#endif
void solution_output_conwrap(int num_soln_flag,double *x,double param,double *x2,double param2,double *x3,double param3,int step_num,int num_its,struct con_struct *con);
#if !defined(_CON_CONST_H_)
void solution_output_conwrap(int type,double *x,double param,double *x2,double param2,double *x3,double param3,int stp,int nits,struct con_struct *con);
#define  NEW_JACOBIAN                200
#endif
int linear_solver_conwrap(double *x,int jac_flag,double *tmp);
#if !defined(_CON_CONST_H_)
int linear_solver_conwrap(double *x,int jac_flag,double *tmp);
#define  CONT_TANGENT      300
#endif
void calc_rhs_continuation(int rhs_type,double *x,double *resid_vector,double *ab_vec,double *scale_vec,double *x_tmp,double param,double perturb,double *r_vec,int numUnks,int numOwnedUnks);
#if !defined(_CON_CONST_H_)
void calc_rhs_continuation(int rhs_type,double *x,double *a,double *x_dot,double *scale_vec,double *x_tmp,double con_par,double perturb,double *r_vec,int num_total_unknowns,int num_owned_unks);
#endif
int nonlinear_solver_conwrap(double *x,void *con_ptr,int step_num,double lambda,double delta_s,void *aux_info);
#if !defined(_CON_CONST_H_)
int nonlinear_solver_conwrap(double *x,void *con,int step_num,double lambda,double delta_s,void *aux_info);
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
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
void perturb_solution_conwrap(double *x,double *x_old,double *scale_vec,int numOwnedUnks);
#if !defined(_CON_CONST_H_)
void perturb_solution_conwrap(double *x,double *p,double *s,int n);
#endif
void assign_bif_parameter_conwrap(double tp_param);
#if !defined(_CON_CONST_H_)
void assign_bif_parameter_conwrap(double bif_param);
#endif
void assign_parameter_conwrap(double param);
#if !defined(_CON_CONST_H_)
void assign_parameter_conwrap(double param);
#define PHASE_TRANSITION_CONTINUATION  6
#define HOPF_CONTINUATION              5
#define PITCHFORK_CONTINUATION         4
#define TURNING_POINT_CONTINUATION     3
#define ARC_LENGTH_CONTINUATION        2
#define FIRST_ORDER_CONTINUATION       1
#define ZERO_ORDER_CONTINUATION        0
#endif
void vec_copy(double *dx,double *dy);
void vec_copy(double *dx,double *dy);
double *alloc_vec();
double *alloc_vec();
#if !defined(_CON_CONST_H_)
#define MANIFOLD_CONTINUATION          8
#endif
void initialize_util_routines(int n_o,int n_t);
void initialize_util_routines(int n_o,int n_t);
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
typedef struct arclength_info_struct arclength_info_struct;
struct arclength_info_struct {
   double dp_ds2_goal; /* Square of target dp_ds value for rescaling
                          (desired solution contribution to arc length)     */
   double dp_ds_max;   /* High dp_ds value at which to rescale              */
   double tang_exp;    /* Power to which tang_factor is raised              */
   double tang_step_limit; /* Minimum value of tang_factor between steps    */
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
int con_lib(struct con_struct *con,void *aux_info);
#endif
int con_lib(struct con_struct *con,void *aux_info);
double scaled_dot_prod(double *,double *,double *,int);
double scaled_dot_prod(double *x,double *y,double *scale_vec,int n);
double scaled_dot_prod(double *x,double *y,double *sc,int n);
#if !defined(_CON_CONST_H_)
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
double solution_scale(struct con_struct *con,struct arc_scale_struct *arc);
