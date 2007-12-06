/* This file was automatically generated.  Do not edit! */
double constant_boundary(int iunk,int jnode_box);
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
extern int Nnodes_box;
double HW_boundary_weight(int icomp,int ilist,double *hw_weight,int inode_box,int *reflect_flag);
extern int **Nodes_2_boundary_wall;
extern int Lhard_surf;
extern int **Zero_density_TF;
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
extern int Nlists_HW;
#define NEQ_TYPE       10 
extern int Phys2Unk_last[NEQ_TYPE];
#define DENSITY        0
extern int Phys2Unk_first[NEQ_TYPE];
#define NMER_MAX     100
extern int Unk2Comp[NMER_MAX];
extern int Lseg_densities;
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
