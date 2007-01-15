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
extern double X_const_mu;
extern int Grad_dim;
double constant_boundary(int iunk,int jnode_box);
#define NDIM_MAX  3
extern double Size_x[NDIM_MAX];
extern double Esize_x[NDIM_MAX];
void node_to_ijk(int node,int *ijk);
extern int *B2G_node;
void setup_step_2consts(double **xOwned);
extern double **Vext;
#define NCOMP_MAX 5
#define NSTEPS_MAX 10
extern double Rho_step[NCOMP_MAX][NSTEPS_MAX];
extern double Xend_step[NSTEPS_MAX];
extern double Xstart_step[NSTEPS_MAX];
extern int Orientation_step[NSTEPS_MAX];
extern int Nsteps;
void node_to_position(int inode,double *NodePos);
extern int *L2G_node;
extern double Betamu[NCOMP_MAX];
#define DIFFUSION      5
extern int Lsteady_state;
#define NMER_MAX     100
extern int Unk2Comp[NMER_MAX];
extern int **Zero_density_TF;
#define DENSITY        0
#define NEQ_TYPE       8
extern int Phys2Unk_first[NEQ_TYPE];
extern int *L2B_node;
extern int Nnodes_per_proc;
void setup_linear_profile(double **xOwned);
#define LINEAR           8
#define CHOP_RHO_V       6
void chop_profile(double **xOwned,int iguess);
#define CHOP_RHO_L       5
void setup_stepped_profile(double **xOwned);
#define STEP_PROFILE     3
#define EXP_RHO_L        1
#define EXP_RHO_V        2
void setup_exp_density(double **xOwned,double *rho,int nloop,int index);
#define EXP_RHO          0
#define CONST_RHO_L     -2 
extern double Rho_coex[2];
#define CONST_RHO_V     -1 
extern int Ncomp;
extern double Rho_b[NCOMP_MAX];
extern int Nseg_tot;
extern double Rho_seg_b[NMER_MAX];
void setup_const_density(double **xOwned,double *rho,int nloop,int index);
extern int Lseg_densities;
#define CONST_RHO       -3 
void setup_density(double **xOwned,int iguess);
