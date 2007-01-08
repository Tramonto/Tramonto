/* This file was automatically generated.  Do not edit! */
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
extern int *B2L_node;
double HW_boundary_weight(int icomp,int ilist,double *hw_weight,int inode_box,int *reflect_flag);
extern int **Nodes_2_boundary_wall;
extern int Lhard_surf;
extern int *L2G_node;
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
typedef struct Stencil_Struct Stencil_Struct;
extern struct Stencil_Struct ***Stencil;
extern int Nlists_HW;
extern int ***Bonds;
#define NCOMP_MAX 5
#define NMER_MAX     100
extern int Type_mer[NCOMP_MAX][NMER_MAX];
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int **Nbond;
extern int Nmer[NCOMP_MAX];
extern int Geqn_start[NCOMP_MAX];
extern int Npol_comp;
#define DELTA_FN       0
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
#define NDIM_MAX  3
void setup_polymer_G(double **xOwned);
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
void setup_polymer_rho(double **xOwned,int iguess);
#define DENSITY_MIN  1.e-20
#define WTC          3
#define NONE       -1
#define NONE      -1
#define NONE -1
extern int Type_poly;
extern double **Vext;
#define NSTEPS_MAX 10
extern double Rho_step[NCOMP_MAX][NSTEPS_MAX];
extern double ***Rism_cr;
extern double Xend_step[NSTEPS_MAX];
extern double Xstart_step[NSTEPS_MAX];
extern int Orientation_step[NSTEPS_MAX];
extern int Nsteps;
void node_to_position(int inode,double *NodePos);
extern int *B2G_node;
#define STEP_PROFILE     3
#define CONST_RHO       -3 
extern int **Zero_density_TF;
void node_box_to_ijk_box(int node_box,int *ijk_box);
extern int *L2B_node;
void setup_polymer_simple(double **xOwned,int iguess);
extern double Rho_b[NCOMP_MAX];
extern double VEXT_MAX;
#define CMS_FIELD      1
#define DENSITY        0
#define NEQ_TYPE       8
extern int Phys2Unk_first[NEQ_TYPE];
extern int Ncomp;
extern int Nnodes_per_proc;
void setup_polymer_field(double **xOwned,int iguess);
