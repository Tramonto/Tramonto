/* This file was automatically generated.  Do not edit! */
double HW_boundary_weight(int icomp,int ilist,double *hw_weight,int inode_box,int *reflect_flag);
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(HAS_VALUES_H)
#include <values.h>
#include <unistd.h>
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
double load_polymer_recursion(int sten_type,int iunk,int loc_inode,int inode_box,int unk_B,int itype_mer,int izone,int *ijk_box,double **x);
extern int *B2G_node;
void node_to_position(int inode,double *NodePos);
extern int ***Bonds;
#define DELTA_FN       0
#define POLYMER_GAUSS  5
#define NSTEN        8
extern int Sten_Type[NSTEN];
extern int *Unk_to_Bond;
extern int *Unk_to_Seg;
extern int *Unk_to_Poly;
double load_CMS_Geqns(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
extern int *Pol_Sym;
#define CMS_G          2 
extern int ***Poly_to_Unk;
#define NCOMP_MAX 5
extern int Geqn_start[NCOMP_MAX];
extern int **Nbond;
#define NMER_MAX     40
extern int Type_mer[NCOMP_MAX][NMER_MAX];
extern int Nmer[NCOMP_MAX];
extern double Rho_b[NCOMP_MAX];
#define NBLOCK_MAX   5
extern int Nmer_t[NCOMP_MAX][NBLOCK_MAX];
#define DENSITY        0
double load_CMS_density(int iunk,int loc_inode,int inode_box,double **x,int resid_only_flag);
extern double Charge_f[NCOMP_MAX];
#define POISSON        3
#define NONE       -1
#define NONE      -1
#define NONE -1
extern int Type_coul;
extern void *LinProbMgr_manager;
#define POLYMER_CR     4
double load_mean_field(int sten_type,int iunk,int loc_inode,int icomp,int izone,int *ijk_box,double **x,int resid_only_flag);
double fill_zero_value(int iunk,int loc_inode,int inode_box,double **x);
extern double VEXT_MAX;
extern double **Vext;
extern int **Zero_density_TF;
#define CMS_FIELD      1
#define NEQ_TYPE       8
extern int Phys2Unk_first[NEQ_TYPE];
double load_CMS_field(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
