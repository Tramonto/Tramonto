/* This file was automatically generated.  Do not edit! */
double resid_and_Jac_sten_fill(int sten_type,double **x,int iunk,int junk,int icomp,int jcomp,int loc_inode,int inode_box,int izone,int *ijk_box,int resid_only_flag,int jzone_flag,double(*fp_prefactor)(int,int,int *),double(*fp_resid)(int,int,double **),double(*fp_jacobian)(int,int,double **));
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
extern int **Zero_density_TF;
extern void *LinProbMgr_manager;
double HW_boundary_weight(int icomp,int ilist,double *hw_weight,int inode_box,int *reflect_flag);
extern int **Nodes_2_boundary_wall;
extern int Lhard_surf;
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
#define NCOMP_MAX 5
extern double Rho_b[NCOMP_MAX];
typedef struct Stencil_Struct Stencil_Struct;
extern struct Stencil_Struct ***Stencil;
#define CMS_FIELD      1
#define THETA_CR_DATA         4
#define THETA_CR_RPM_MSA      3
#define WTC          2
#define WJDC_FIELD     8
#define WJDC         3
extern int Type_poly;
#define MF_EQ          9
#define MF_VARIABLE  2
extern int Type_attr;
#define THETA_PAIRPOT_RCUT    2
extern int Nlists_HW;
#define DENSITY        0
#define NEQ_TYPE       11 
extern int Phys2Unk_first[NEQ_TYPE];
extern int *Pol_Sym_Seg;
#define NMER_MAX     100
extern int Unk2Comp[NMER_MAX];
extern int Ncomp;
extern int Nseg_tot;
extern int Lseg_densities;
#define FLAG_PBELEC -777
#define FLAG_BULK   -888
int find_jzone(int izone,int inode_box);
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Ndim;
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
double resid_and_Jac_sten_fill_sum_Ncomp(int sten_type,double **x,int iunk,int loc_inode,int inode_box,int izone,int *ijk_box,int resid_only_flag,int jzone_flag,double(*fp_prefactor)(int,int,int *),double(*fp_resid)(int,int,double **),double(*fp_jacobian)(int,int,double **));
