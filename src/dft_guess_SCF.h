/* This file was automatically generated.  Do not edit! */
void calc_init_lambda(double **xInBox);
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
extern double Rho_t;
#define NCOMP_MAX 5
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
#define NEQ_TYPE       13 
extern int Restart_field[NEQ_TYPE];
extern double **Vext;
extern int **Zero_density_TF;
#define SCF_CONSTR	   9
void calc_init_SCFfield(double **xInBox);
extern double Rho_b[NCOMP_MAX];
extern double VEXT_MAX;
#define SCF_FIELD	  10
#define DENSITY        0
extern int Phys2Unk_first[NEQ_TYPE];
extern int Ncomp;
extern int *L2B_node;
extern int Nnodes_per_proc;
void setup_polymer_SCF_field(double **xInBox,int iguess);
