/* This file was automatically generated.  Do not edit! */
double prefactor_cavity_wtc(int iunk,int icomp,int *offset);
double resid_and_Jac_sten_fill_sum_Ncomp(int sten_type,double **x,int iunk,int loc_inode,int inode_box,int izone,int *ijk_box,int resid_only_flag,int jzone_flag,double(*fp_prefactor)(int,int,int *),double(*fp_resid)(int,int,double **),double(*fp_jacobian)(int,int,double **));
double load_cavity_wtc(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
double jac_rho_bar(int junk,int jnode_box,double **x);
double resid_rho_bar(int junk,int jnode_box,double **x);
double resid_and_Jac_sten_fill(int sten_type,double **x,int iunk,int junk,int icomp,int jcomp,int loc_inode,int inode_box,int izone,int *ijk_box,int resid_only_flag,int jzone_flag,double(*fp_prefactor)(int,int,int *),double(*fp_resid)(int,int,double **),double(*fp_jacobian)(int,int,double **));
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
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int *BondAll_to_ibond;
extern int *BondAll_to_isegAll;
double load_bond_wtc(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
double d2y_dxi3_dxi2(double sigma_1,double sigma_2,double xi_2,double xi_3);
double d2y_dxi3_sq(double sigma_1,double sigma_2,double xi_2,double xi_3);
double d2y_dxi2_sq(double sigma_1,double sigma_2,double xi_2,double xi_3);
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
#define PI    M_PI
extern int Nseg_tot;
#define THETA_FN_SIG   6
double load_polyTC_cavityEL(int iunk,int loc_inode,int inode_box,int icomp,int izone,int *ijk_box,double **x,int resid_only_flag);
double constant_boundary(int iunk,int jnode_box);
double HW_boundary_weight(int icomp,int ilist,double *hw_weight,int inode_box,int *reflect_flag);
extern int **Nodes_2_boundary_wall;
extern int Lhard_surf;
extern int **Zero_density_TF;
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
extern int Ncomp;
#define DELTA_FN_BOND  7
typedef struct Stencil_Struct Stencil_Struct;
extern struct Stencil_Struct ***Stencil;
extern int Nlists_HW;
extern int *Pol_Sym_Seg;
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
double load_polyTC_bondEL(int iunk,int loc_inode,int inode_box,int icomp,int izone,int *ijk_box,double **x,int resid_only_flag);
extern void *LinProbMgr_manager;
#define NCOMP_MAX 5
extern double Fac_overlap[NCOMP_MAX][NCOMP_MAX];
#define BONDWTC       7
extern int *Pol_Sym;
extern int **Poly_to_Unk_SegAll;
double dy_dxi3_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
double dy_dxi2_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
double y_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
#define CAVWTC     6
#define NMER_MAX     100
extern int Unk2Comp[NMER_MAX];
extern int **Bonds_SegAll;
extern int *Nbonds_SegAll;
#define DENSITY        0
#define NEQ_TYPE       8
extern int Phys2Unk_first[NEQ_TYPE];
double load_polyTC_diagEL(int iunk,int loc_inode,int inode_box,int icomp,int izone,int *ijk_box,double **x,int resid_only_flag);
