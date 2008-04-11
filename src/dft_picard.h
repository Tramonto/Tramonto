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
#define NCOMP_MAX 5
#define NMER_MAX     100
extern int Unk2Phys[3 *NCOMP_MAX+NMER_MAX+NMER_MAX *NMER_MAX+13];
#define NDIM_MAX  3
extern int Nodes_x[NDIM_MAX];
#define LAST_NODE_RESTART    4
extern int Type_bc[NDIM_MAX][2];
extern int Ndim;
void node_to_ijk(int node,int *ijk);
extern int *L2G_node;
extern int *B2L_node;
extern double NL_update_scalingParam;
double gsum_double(double c);
double fill_resid_and_matrix_control(double **x,int iter,int resid_only_flag);
double load_WJDC_density(int iunk,int loc_inode,int inode_box,double **x,int resid_only_flag);
double load_CMS_density(int iunk,int loc_inode,int inode_box,double **x,int resid_only_flag);
#define INIT_GUESS_FLAG  2
typedef struct RB_Struct RB_Struct;
double load_euler_lagrange(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,struct RB_Struct *dphi_drb,int mesh_coarsen_flag_i,int resid_only_flag);
#define NEQ_TYPE       11 
extern int Phys2Unk_first[NEQ_TYPE];
void node_box_to_ijk_box(int node_box,int *ijk_box);
extern int Ncomp;
extern int Nseg_tot;
extern int Lseg_densities;
void FMT1stDeriv_switch(double **x,struct RB_Struct *dphi_drb);
#define NONE       -1
#define NONE      -1
#define NONE        -1
#define NONE        -1
extern int Type_func;
struct RB_Struct {
  double    S0;      /*   1/(4*pi*Ri*Ri) * Delta_fn   */
  double    S1;      /*   1/(4*pi*Ri)    * Delta_fn   */
  double    S2;      /*                    Delta_fn   */
  double    S3;      /*                    Theta_fn   */
  double    V1[NDIM_MAX];      /*  1/(4*pi*Ri) * unit_vec * Delta_Fn   */
  double    V2[NDIM_MAX];      /*                unit_vec * Delta_Fn   */
};
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern int Max_NL_iter;
#define YW_DENS        10       /* densities for Yethiraj-Woodward polymer DFTs */
void calc_init_polymer_G_wjdc(double **xInBox);
void calc_init_polymer_G_CMS(double **xInBox);
#define G_CHAIN        9 
void calc_init_CMSfield(double **xInBox);
#define CMS_FIELD      7
void calc_init_WJDC_field(double **xInBox);
#define WJDC_FIELD     8
void calc_init_BondWTC(double **xInBox);
#define BONDWTC        5
void calc_init_Xi_cavWTC(double **xInBox);
#define CAVWTC         4
#define DIFFUSION      6
#define POISSON        1
void calc_init_rho_bar(double **xInBox);
#define HSRHOBAR       2
void calc_init_mf_attract(double **xInBox);
extern int Phys2Nunk[NEQ_TYPE];
#define MF_EQ          3
#define DENSITY        0
void fix_symmetries(double **x);
int update_solution_picard(double **x,double **delta_x,int iter);
extern double NL_abs_tol,NL_rel_tol;
int continuation_hook_conwrap(double **xx,double **delta_xx,void *con_ptr,double reltol,double abstol);
void calc_density_next_iter_WJDC(double **xInBox);
void calc_density_next_iter_CMS(double **xInBox);
#define CMS          0
void calc_density_next_iter_HSperturb(double **xInBox);
#define WJDC         3
extern int Type_poly;
extern int L_HSperturbation;
void print_resid_norm_picard(double **x,int iter);
#define NO_SCREEN    2 
extern int Nnodes_box;
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
int picard_solver(double **x,void *con_ptr);
void safe_free(void **ptr);
void safe_free(void **ptr);
#define BINODAL_FLAG  -1325  /* to let initial guess routine know we need to fill X2 */
extern int Lbinodal;
void print_profile_box(double **x,char *outfile);
#define VERBOSE      3 
extern int Iwrite;
void communicate_to_fill_in_box_values(double **xInBox);
extern int *L2B_node;
extern int Iguess1;
void set_initial_guess(int iguess,double **xOwned);
extern int Nnodes_per_proc;
extern int Nunk_per_node;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
int solve_problem_picard(double **x,double **x2);
