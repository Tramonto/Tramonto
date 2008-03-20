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
extern int **Bonds_SegAll;
extern int *BondAll_to_ibond;
extern int *BondAll_to_isegAll;
#define BONDWTC        7
#define DELTA_FN_BOND         6
double int_stencil_BondWTC(double **x,int inode_box,int iunk);
double prefactor_cavity_wtc(int iunk,int icomp,int *offset);
#define THETA_FN_SIG          5
double int_stencil_CAV(double **x,int inode_box,int iunk);
double prefactor_rho_bar_v(int iunk,int jcomp,int *offset);
double prefactor_rho_bar_s(int iunk,int jcomp,int *offset);
extern int Nrho_bar_s;
#define DELTA_FN_R            0
#define THETA_FN_R            1
#define HSRHOBAR       4
double int_stencil_HSFMT(double **x,int inode_box,int iunk);
#define NCOMP_MAX 5
extern double Rho_b[NCOMP_MAX];
double int_stencil_CMSField(double **x,int inode_box,int iunk,int sten_type);
double constant_boundary(int iunk,int jnode_box);
extern int Nnodes_box;
double HW_boundary_weight(int icomp,int ilist,double *hw_weight,int inode_box,int *reflect_flag);
extern int **Nodes_2_boundary_wall;
extern int Lhard_surf;
extern int **Zero_density_TF;
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
extern int Nlists_HW;
#define NEQ_TYPE       11 
extern int Phys2Unk_last[NEQ_TYPE];
#define DENSITY        0
extern int Phys2Unk_first[NEQ_TYPE];
#define NMER_MAX     100
extern int Unk2Comp[NMER_MAX];
extern int Lseg_densities;
void node_box_to_ijk_box(int node_box,int *ijk_box);
void node_box_to_ijk_box(int node_box,int *ijk_box);
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define NDIM_MAX  3
double int_stencil(double **x,int inode_box,int iunk,int sten_type);
extern double Esize_x[NDIM_MAX];
extern int Ndim;
extern int Ncomp;
typedef struct Stencil_Struct Stencil_Struct;
extern struct Stencil_Struct ***Stencil;
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
double int_stencil_bulk(int sten_type,int icomp,int jcomp,double(*fp_integrand)(double,int,int));
