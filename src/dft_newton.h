/* This file was automatically generated.  Do not edit! */
void fill_test(double **x,int flag);
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
#include "az_aztec_defs.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define SCREEN_DEBUG_RESID 2
#define SCREEN_VERBOSE     3 
#define NO_SCREEN    4 
extern int *L2G_node;
int find_length_of_file(char *filename);
#define CALC_AND_FILL_RESID_ONLY  4
extern int Nnodes;
extern int Lseg_densities;
#define HSRHOBAR       2
#define CAVWTC         4
#define BONDWTC        5
#define WJDC_FIELD     8
#define CMS_FIELD      7
#define NDIM_MAX  3
extern int Nodes_x[NDIM_MAX];
#define LAST_NODE_RESTART    4
extern int Type_bc[NDIM_MAX][2];
extern int Ndim;
void node_to_ijk(int node,int *ijk);
extern int *B2G_node;
extern double NL_update_scalingParam;
extern int *Pol_Sym_Seg;
#define WTC          2
extern int Type_poly;
#define DENSITY        0
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
extern int *Pol_Sym;
#define G_CHAIN       11 
#define NCOMP_MAX 5
#define NMER_MAX     200
extern int Unk2Phys[3 *NCOMP_MAX+2 *NMER_MAX+NMER_MAX *NMER_MAX+13];
double gsum_double(double c);
#define SCREEN_BASIC       1
double gmin_double(double c);
int update_solution_new(double **x,double **delta_x,int iter);
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern int Max_NL_iter;
void fix_symmetries(double **x);
int update_solution(double **x,double **delta_x,int iter);
extern double NL_abs_tol,NL_rel_tol;
int continuation_hook_conwrap(double **xx,double **delta_xx,void *con_ptr,double reltol,double abstol);
extern double Time_linsolver_av;
extern double Time_linsolver_first;
extern double Time_manager_av;
extern double Time_manager_first;
void print_resid_norm(int iter);
#define SCREEN_ERRORS_ONLY  0 
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
extern double Time_fill_av;
extern double Time_fill_first;
double fill_resid_and_matrix_control(double **x,int iter,int resid_only_flag);
#define FILENAME_LENGTH 300
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int *B2L_node;
extern int Nnodes_box;
void box2owned(double **xBox,double **xOwned);
extern void *ParameterList_list;
void safe_free(void **ptr);
void safe_free(void **ptr);
extern double Time_NLSolve;
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
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define NEWTON_NOX            1
#define BINODAL_FLAG  -1325  /* to let initial guess routine know we need to fill X2 */
extern int Lbinodal;
void print_profile_box(double **x,char *outfile);
#define FILES_DEBUG        2
extern int Iwrite_files;
extern void *LinProbMgr_manager;
extern double Time_InitGuess;
extern int Iguess;
void set_initial_guess(int guess_type,double **xOwned);
extern int *L2B_node;
#define PICNEWTON_NOX         5
#define PICNEWTON_BUILT_IN    4
extern int NL_Solver;
extern int Nnodes_per_proc;
extern int Nunk_per_node;
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
int solve_problem(double **x,double **x2);
void do_numerical_jacobian(double **);
void do_numerical_jacobian(double **x);
