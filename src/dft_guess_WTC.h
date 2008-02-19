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
#define NMER_MAX     100
extern double BondWTC_b[NMER_MAX *NMER_MAX];
#define BONDWTC        7
extern int Nbonds;
void setup_BondWTC(double **xOwned);
extern double Xi_cav_b[4];
#define NEQ_TYPE       11 
extern int Phys2Unk_first[NEQ_TYPE];
#define CAVWTC         6
extern int Phys2Nunk[NEQ_TYPE];
extern int Nnodes_per_proc;
void setup_Xi_cavWTC(double **xOwned);
