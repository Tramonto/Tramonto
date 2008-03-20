/* This file was automatically generated.  Do not edit! */
double int_stencil_HSFMT(double **x,int inode_box,int iunk);
void calc_init_rho_bar(double **xInBox);
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
extern double Rhobar_b[10];
#define NDIM_MAX  3
extern double Size_x[NDIM_MAX];
extern double Rhobar_b_RTF[10];
extern double Rhobar_b_LBB[10];
extern int Grad_dim;
extern double Esize_x[NDIM_MAX];
void node_to_ijk(int node,int *ijk);
void node_to_ijk(int node,int *ijk);
extern int *B2G_node;
extern int Iliq_vap;
extern int Nwall;
extern int Lsteady_state;
#define HSRHOBAR       2
#define NEQ_TYPE       11 
extern int Phys2Unk_first[NEQ_TYPE];
extern int Nrho_bar;
extern int *L2B_node;
extern int Nnodes_per_proc;
void setup_rho_bar(double **xInBox);
