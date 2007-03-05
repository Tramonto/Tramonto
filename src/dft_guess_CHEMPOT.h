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
extern double VEXT_MAX;
#define NCOMP_MAX 5
extern double Betamu_RTF[NCOMP_MAX];
extern double Betamu_LBB[NCOMP_MAX];
#define POISSON        3
extern double Charge_f[NCOMP_MAX];
#define DENSITY        0
extern int Ipot_ff_c;
extern int **Zero_density_TF;
#define DIFFUSION      5
#define NEQ_TYPE       8
extern int Phys2Unk_first[NEQ_TYPE];
extern int Ncomp;
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
void node_to_ijk(int node,int *ijk);
extern int *B2G_node;
extern int *L2B_node;
extern int Nnodes_per_proc;
extern double X_const_mu;
extern int Grad_dim;
extern double Size_x[NDIM_MAX];
void setup_chem_pot(double **xOwned);
