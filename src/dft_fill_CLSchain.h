/* This file was automatically generated.  Do not edit! */
double HW_boundary_weight(int icomp,int ilist,double *hw_weight,int inode_box,int *reflect_flag);
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
extern int **Nodes_2_boundary_wall;
extern int Lhard_surf;
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
extern int Ncomp;
typedef struct Stencil_Struct Stencil_Struct;
extern struct Stencil_Struct ***Stencil;
#define NCOMP_MAX 5
extern int Geqn_start[NCOMP_MAX];
extern int ***Poly_to_Unk;
extern int **Nbond;
extern int Nlists_HW;
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
double load_polymer_recursion(int sten_type,int func_type_field,int Njacobian_types,int Njacobian_sums,void(*funcArray_Jac[3])(int,int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **),int iunk,int loc_inode,int inode_box,int unk_B,int itype_mer,int izone,int *ijk_box,double **x,int resid_only_flag);
double dy_dxi3_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
double dy_dxi2_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
double y_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
#define CAVWTC         6
double fill_zero_value(int iunk,int loc_inode,int inode_box,double **x,int resid_only_flag);
extern double VEXT_MAX;
extern double **Vext;
extern int **Zero_density_TF;
#define DELTA_FN_BOND         6
extern int ***Bonds;
extern int *Unk_to_Bond;
extern int *Unk_to_Seg;
extern int *Unk_to_Poly;
#define G_CHAIN        2 
double load_Chain_Geqns(int func_type_field,int Njacobian_types,int Njacobian_sums,void(*funcArray_Jac[3])(int,int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **),int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
extern int *Pol_Sym;
extern int **Poly_to_Unk_SegAll;
extern int *Nbonds_SegAll;
#define NMER_MAX     100
extern int SegChain2SegAll[NCOMP_MAX][NMER_MAX];
extern int Type_mer[NCOMP_MAX][NMER_MAX];
extern int Lseg_densities;
extern int Nmer[NCOMP_MAX];
#define NBLOCK_MAX   5
extern int Nmer_t[NCOMP_MAX][NBLOCK_MAX];
#define CMS          0
#define DENSITY        0
#define NEQ_TYPE       10 
extern int Phys2Unk_first[NEQ_TYPE];
extern int Unk2Comp[NMER_MAX];
#define WJDC         3
extern int Type_poly;
extern void *LinProbMgr_manager;
double resid_and_Jac_ChainDensity(int func_type,double **x,int iunk,int unk_B,int loc_inode,int inode_box,int resid_only_flag,double(*fp_prefactor)(int));
