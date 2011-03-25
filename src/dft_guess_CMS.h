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
#define INIT_GUESS_FLAG  2
double load_Chain_Geqns(int func_type_field,int Njacobian_types,int Njacobian_sums,void(*funcArray_Jac[3])(int,int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **),int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
extern int *Pol_Sym;
extern int ***Poly_to_Unk;
extern int *Unk_to_Bond;
extern int *Unk_to_Poly;
double CMS_Resid_Bulk_GCHAIN(int iunk,int pol_num,int jseg,int unk_B,int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x);
double CMS_Resid_GCHAIN(int iunk,int pol_num,int jseg,int unk_B,int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x);
#define NBOND_MAX 4
#define NMER_MAX     200
void calc_init_polymer_G_CMS(double **xInBox,double **xOwned);
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
double HW_boundary_weight(int icomp,int ilist,double *hw_weight,int inode_box,int *reflect_flag);
extern int **Nodes_2_boundary_wall;
extern int Lhard_surf;
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
typedef struct Stencil_Struct Stencil_Struct;
extern struct Stencil_Struct ***Stencil;
extern int Nlists_HW;
#define NCOMP_MAX 5
extern double Bond_ff[NCOMP_MAX][NCOMP_MAX];
extern double *Poly_graft_dist;
#define NDIM_MAX  3
#define NWALL_MAX 600 
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern int *L2G_node;
#define SCF_FIELD	  10
extern int ***Bonds;
extern int Type_mer[NCOMP_MAX][NMER_MAX];
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int **Nbond;
extern int Nmer[NCOMP_MAX];
extern int Geqn_start[NCOMP_MAX];
extern int Npol_comp;
#define DELTA_FN_BOND         6
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
void setup_polymer_G(double **xInBox,double **xOwned);
extern double G_CMS_b[NMER_MAX *NBOND_MAX];
#define G_CHAIN       11 
extern int Unk2Comp[NMER_MAX];
extern int *Unk_to_Seg;
extern int Nbonds;
void setup_polymer_G_newCMS(double **xOwned);
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern int *B2L_node;
extern int Nnodes_box;
void setup_polymer_rho(double **xInBox,double **xOwned,int guess_type);
#define DENSITY_MIN  1.e-20
#define CMS_SCFT     1
#define CMS          0
extern int Type_poly;
#define NSTEPS_MAX 10
extern double Rho_step[NCOMP_MAX][NSTEPS_MAX];
extern double ***Rism_cr;
extern double Xend_step[NSTEPS_MAX];
extern double Xstart_step[NSTEPS_MAX];
extern int Orientation_step[NSTEPS_MAX];
extern int Nsteps;
void node_to_position(int inode,double *NodePos);
extern int *B2G_node;
#define STEP_PROFILE     2
#define CONST_RHO        0 
void node_box_to_ijk_box(int node_box,int *ijk_box);
void setup_polymer_simple(double **xInBox,int guess_type);
#define THETA_CR_DATA         4
double int_stencil_CMSField(double **x,int inode_box,int iunk,int sten_type);
extern double **Vext;
extern int **Zero_density_TF;
extern void *LinProbMgr_manager;
void calc_init_CMSfield(double **xInBox,double **xOwned);
extern int *L2B_node;
extern double Rho_b[NCOMP_MAX];
extern double VEXT_MAX;
#define CMS_FIELD      7
#define DENSITY        0
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
extern int Nmissing_densities;
#define RESTART_FEWERCOMP  4
extern int Restart;
extern int Ncomp;
extern int Nnodes_per_proc;
void setup_polymer_field(double **xInBox,double **xOwned,int guess_type);
