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
double load_poisson_control(int iunk,int loc_inode,int inode_box,int *ijk_box,double **x,int resid_only_flag);
void node_box_to_ijk_box(int node_box,int *ijk_box);
extern int *L2B_node;
void calc_init_elec_pot(double **xInBox,double **xOwned);
#define NDIM_MAX  3
extern double Size_x[NDIM_MAX];
extern double Elec_pot_RTF;
extern double Elec_pot_LBB;
extern int Grad_dim;
extern double Esize_x[NDIM_MAX];
void node_to_ijk(int node,int *ijk);
extern int *L2G_node;
#define LINEAR           5
#define UNIFORM_INTERFACE  0
extern int Type_interface;
#define POISSON        1
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
extern int Nnodes_per_proc;
void setup_elec_pot(double **xOwned,int guess_type);
