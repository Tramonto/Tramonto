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
#include "az_aztec_defs.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
extern int *L2G_node;
double load_lambda_field(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
#define SCF_CONSTR	   9
double load_SCF_field(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
#define SCF_FIELD	  10
double load_WJDC_Geqns(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
double load_SCF_Geqns(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
double load_CMS_Geqns(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
#define G_CHAIN       11 
typedef struct RB_Struct RB_Struct;
double load_WJDC_field(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,struct RB_Struct *dphi_drb,int mesh_coarsen_flag_i,int resid_only_flag);
#define WJDC_FIELD     8
double load_CMS_field(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
#define CMS_FIELD      7
double load_bond_wtc(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
#define BONDWTC        5
double load_cavity_wtc(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
double load_nonlinear_transport_eqn(int iunk,int loc_inode,int inode_box,int *ijk_box,double **x,int resid_only_flag);
double load_linear_transport_eqn(int iunk,int loc_inode,int inode_box,int *ijk_box,double **x,int resid_only_flag);
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int Linear_transport;
#define DIFFUSION      6
double load_poisson_control(int iunk,int loc_inode,int inode_box,int *ijk_box,double **x,int resid_only_flag);
#define POISSON        1
#define THETA_CR_DATA         4
#define THETA_PAIRPOT_RCUT    2
double load_mean_field(int sten_type,int iunk,int loc_inode,int icomp,int izone,int *ijk_box,double **x,int resid_only_flag);
#define MF_EQ          3
double load_rho_bar_v(double **x,int iunk,int loc_inode,int inode_box,int izone,int *ijk_box,int resid_only_flag);
#define DELTA_FN_R            0
extern int Nrho_bar_s;
#define THETA_FN_R            1
double load_rho_bar_s(int sten_type,double **x,int iunk,int loc_inode,int inode_box,int izone,int *ijk_box,int resid_only_flag);
#define HSRHOBAR       2
#define SCFT         6	
double load_SCF_density(int iunk,int loc_inode,int inode_box,double **x,int resid_only_flag);
double load_WJDC_density(int iunk,int loc_inode,int inode_box,double **x,int resid_only_flag);
double load_CMS_density(int iunk,int loc_inode,int inode_box,double **x,int resid_only_flag);
double load_euler_lagrange(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,struct RB_Struct *dphi_drb,int mesh_coarsen_flag_i,int resid_only_flag);
#define WJDC2        4 
#define WJDC         3
extern int L_HSperturbation;
#define DENSITY        0
#define NCOMP_MAX 5
#define NMER_MAX     200
extern int Unk2Phys[3 *NCOMP_MAX+2 *NMER_MAX+NMER_MAX *NMER_MAX+13];
extern int **Zero_density_TF;
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
extern int Ncomp;
typedef struct Stencil_Struct Stencil_Struct;
extern struct Stencil_Struct ***Stencil;
#define NDIM_MAX  3
struct Stencil_Struct {
  int        Length;      /* Number of nodes that interact with current 
                             node through this stencil                    */
  int      **Offset;      /* 2D array to be set to size [Length][Ndim] that
                             gives relative position of interacting node  */
  double    *Weight;      /* 1D array of size [Length] that gives weight
                             of interacting node to total stencil         */
  double   **HW_Weight;   /* 2D array of size [Length][Nnodes_per_el] that
                             holds the weights based on what element they
                             are being contributed from. Only used for Hard
                             Walls when stencil point is a boundary node  */
};
#define SCREEN_NONE       -1 
double y_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
#define CAVWTC         4
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
#define DELTA_FN_BOND         6
double grafted_int(int sten_type,int itype_mer,int *ijk_box,int izone,int unk_G,double **x);
int node_to_node_box(int inode);
int position_to_node(double *NodePos);
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern double *Poly_graft_dist;
#define NWALL_MAX 600 
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern int WallType[NWALL_MAX];
extern int Graft_wall[NCOMP_MAX];
void node_to_position(int inode,double *NodePos);
extern int **Nodes_2_boundary_wall;
extern int Type_mer[NCOMP_MAX][NMER_MAX];
extern int ***Bonds;
extern int *Nbonds_SegAll;
extern int Grafted[NCOMP_MAX];
void safe_free(void **ptr);
void safe_free(void **ptr);
void print_residuals(int loc_inode,int iunk,double *resid_unk);
#define SCREEN_DEBUG_RESID 2
double load_standard_node(int loc_inode,int inode_box,int *ijk_box,int iunk,double **x,struct RB_Struct *dphi_drb,double *resid_unk,int mesh_coarsen_flag_i,int resid_only_flag);
double load_coarse_node_Ndim(int loc_inode,int inode_box,int iunk,double **x,int resid_only_flag);
#define FLAG_PBELEC -777
#define FLAG_BULK   -888
void load_coarse_node_1dim(int loc_inode,int inode_box,int *ijk_box,int iunk,double **x,int resid_only_flag);
#define FLAG_1DBC   -999
extern int *Mesh_coarsen_flag;
extern int L1D_bc;
extern int Nwall_type;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Mesh_coarsening;
void node_box_to_ijk_box(int node_box,int *ijk_box);
void calc_Gsum(double **x);
#define WJDC3        5 
#define CMS          0
double gsum_double(double c);
extern int Nnodes_per_el_V;
extern double Vol_el;
extern int **Nel_hit2;
extern int *L2B_node;
extern int Nnodes_per_proc;
extern double Gsum[NCOMP_MAX];
extern int ***Poly_to_Unk;
extern int Geqn_start[NCOMP_MAX];
extern int Nmer[NCOMP_MAX];
extern int Npol_comp;
extern double Size_x[NDIM_MAX];
extern int Ndim;
#define CMS_SCFT     1
extern int Type_poly;
#define NODAL_FLAG -999
void print_profile_box(double **x,char *outfile);
#define FILES_DEBUG        2
extern int Nnodes;
extern double **Array_test;
#define FILES_DEBUG_MATRIX 3 
extern int Iwrite_files;
extern int Nunk_per_node;
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
#define SCREEN_VERBOSE     3 
extern int Iwrite_screen;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
struct RB_Struct {
  double    S0;      /*   1/(4*pi*Ri*Ri) * Delta_fn   */
  double    S1;      /*   1/(4*pi*Ri)    * Delta_fn   */
  double    S2;      /*                    Delta_fn   */
  double    S3;      /*                    Theta_fn   */
  double    V1[NDIM_MAX];      /*  1/(4*pi*Ri) * unit_vec * Delta_Fn   */
  double    V2[NDIM_MAX];      /*                unit_vec * Delta_Fn   */
};
double fill_resid_and_matrix(double **x,struct RB_Struct *dphi_drb,int iter,int resid_only_flag,int unk_flag);
