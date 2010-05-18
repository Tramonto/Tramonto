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
#define THETA_PAIRPOT_RCUT    2
double int_stencil(double **x,int inode_box,int iunk,int sten_type);
extern int *L2B_node;
void calc_init_mf_attract(double **xInBox,double **xOwned);
#define NCOMP_MAX 5
extern double Rho_b[NCOMP_MAX];
extern double Avdw[NCOMP_MAX][NCOMP_MAX];
#define PHASE_INTERFACE 2
#define DIFFUSIVE_INTERFACE 1
extern int Type_interface;
#define MF_EQ          3
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
extern int Nmissing_densities;
#define RESTART_FEWERCOMP  4
extern int Restart;
extern int Ncomp;
extern int Nnodes_per_proc;
void setup_mf_attract(double **xOwned);
