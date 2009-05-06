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
extern double X_const_mu;
extern int Grad_dim;
#define UNIFORM_INTERFACE  0
double constant_boundary(int iunk,int jnode_box);
#define NDIM_MAX  3
extern double Size_x[NDIM_MAX];
extern double Esize_x[NDIM_MAX];
void node_to_ijk(int node,int *ijk);
void setup_step_2consts(double **xInBox);
extern double **Vext;
extern int *B2L_node;
#define DIFFUSIVE_INTERFACE 1
extern int Type_interface;
#define NCOMP_MAX 5
#define NSTEPS_MAX 10
extern double Rho_step[NCOMP_MAX][NSTEPS_MAX];
extern double Xend_step[NSTEPS_MAX];
extern double Xstart_step[NSTEPS_MAX];
extern int Orientation_step[NSTEPS_MAX];
extern int Nsteps;
void node_to_position(int inode,double *NodePos);
extern int *B2G_node;
#define NMER_MAX     100
extern int Unk2Comp[NMER_MAX];
extern int **Zero_density_TF;
#define DENSITY        0
#define NEQ_TYPE       13 
extern int Phys2Unk_first[NEQ_TYPE];
extern int Nnodes_box;
void translate_xInBox_to_xOwned(double **xInBox,double **xOwned);
void setup_linear_profile(double **xInBox);
#define LINEAR           5
void setup_stepped_profile(double **xInBox);
#define STEP_PROFILE     2
void setup_exp_density(double **xInBox,double *rho,int nloop,int index);
#define EXP_RHO          1
extern int Ncomp;
extern double Rho_b[NCOMP_MAX];
extern int Nseg_tot;
extern double Rho_seg_b[NMER_MAX];
void setup_const_density(double **xInBox,double *rho,int nloop,int index);
extern int Lseg_densities;
#define CONST_RHO        0 
void setup_density(double **xInBox,double **xOwned,int iguess);
