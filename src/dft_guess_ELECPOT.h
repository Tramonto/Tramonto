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
#define NDIM_MAX  3
extern double Size_x[NDIM_MAX];
extern double Elec_pot_RTF;
extern double Elec_pot_LBB;
extern int Grad_dim;
extern double Esize_x[NDIM_MAX];
void node_to_ijk(int node,int *ijk);
extern int *B2G_node;
extern int *L2B_node;
#define LINEAR           8
extern int Lsteady_state;
#define POISSON        3
#define NEQ_TYPE       9
extern int Phys2Unk_first[NEQ_TYPE];
extern int Nnodes_per_proc;
void setup_elec_pot(double **xOwned,int iguess);
