/* This file was automatically generated.  Do not edit! */
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
extern int *L2B_node;
void fill_test(double **x,int flag);
#define NODAL_FLAG -999
void fill_resid_and_matrix(double **x,int iter,int resid_only_flag,int unk_flag);
extern int Nnodes;
extern int Lseg_densities;
double gsum_double(double c);
double gmin_double(double c);
extern double Min_update_frac;
extern int *Pol_Sym_Seg;
extern int *Pol_Sym;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern int Max_Newton_iter;
void fix_symmetries(double **x);
int update_solution(double **x,double **delta_x,int iter);
extern double Newton_abs_tol,Newton_rel_tol;
int continuation_hook_conwrap(double **xx,double **delta_xx,void *con_ptr,double reltol,double abstol);
extern double Time_linsolver_av;
extern double Time_linsolver_first;
#if defined(NUMERICAL_JACOBIAN)
void do_numerical_jacobian(double **x);
#endif
extern double Time_manager_av;
extern double Time_manager_first;
void print_resid_norm(int iter);
#define NO_SCREEN    2 
extern double Time_fill_av;
extern double Time_fill_first;
void fill_resid_and_matrix_control(double **x,int iter,int resid_only_flag);
extern int *B2L_node;
void box2owned(double **xBox,double **xOwned);
int newton_solver(double **x,void *con_ptr);
int solve_continuation(double **xx,double **xx2);
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
#define BINODAL_FLAG  -1325  /* to let initial guess routine know we need to fill X2 */
extern int Lbinodal;
void print_profile_box(double **x,char *outfile);
#define VERBOSE      3 
extern int Iwrite;
extern int Iguess1;
void set_initial_guess(int iguess,double **xOwned);
extern int *List_coarse_nodes;
extern int Nnodes_coarse_loc;
extern int *B2G_node;
extern int Nnodes_box;
extern int *L2G_node;
extern int Nnodes_per_proc;
extern int Mesh_coarsening;
#define WTC          2
extern int Type_attr;
#define CAVWTC     6
#define BONDWTC       7
extern int Ndim;
extern int Nrho_bar_s;
#define NEQ_TYPE       9
extern int Phys2Unk_first[NEQ_TYPE];
#define HSRHOBAR       4
#define POISSON        3
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int Nlists_HW;
extern int Lhard_surf;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Type_func;
void safe_free(void **ptr);
void safe_free(void **ptr);
extern int Ncomp;
extern int Ngeqn_tot;
typedef struct Aztec_Struct Aztec_Struct;
struct Aztec_Struct {
  int    options[AZ_OPTIONS_SIZE]; /* Array used to select solver options.  */
  double params[AZ_PARAMS_SIZE];   /* User selected solver paramters.       */
#ifdef DONE_WITH_THESE
  int    proc_config[AZ_PROC_SIZE];/* Processor information.                */
  int    *data_org;                /* Array to specify data layout          */
  double status[AZ_STATUS_SIZE];   /* Information returned from AZ_solve(). */
  int    *update;                  /* vector elements updated on this node. */
  int    *external;                /* vector elements needed by this node.  */
  int    *update_index;            /* ordering of update[] and external[]   */
  int    *extern_index;            /* locally on this processor.            */
  int    *bindx;                   /* Sparse matrix to be solved is stored  */
  double *val;                     /* in these MSR arrays.                  */
  int    N_update;                 /* # of unknowns updated on this node    */
  int    nonzeros;                 /* # of nonzeros in sparse matrix        */
#endif
};
extern struct Aztec_Struct Aztec;
extern void *LinProbMgr_manager;
#define NCOMP_MAX 5
extern int Geqn_start[NCOMP_MAX];
#define CMS_G          2 
#define CMS_FIELD      1
#define DENSITY        0
#define NMER_MAX     100
extern int Unk2Phys[3 *NCOMP_MAX+NMER_MAX+NMER_MAX *NMER_MAX+13];
extern int Nunk_per_node;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
#define NONE       -1
#define NONE      -1
#define NONE -1
#define NONE        -1
extern int Type_coul;
#define CMS          0
extern int Type_poly;
extern int L_Schur;
int solve_problem(double **x,double **x2);
