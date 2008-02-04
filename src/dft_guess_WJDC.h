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
#define NBOND_MAX 4
extern double G_WJDC_b[NMER_MAX *NBOND_MAX];
#define G_CHAIN        2 
extern int Nbonds;
void setup_polymer_G_wjdc(double **xOwned);
extern double Field_WJDC_b[NMER_MAX];
#define WJDC_FIELD     8
#define NEQ_TYPE       10 
extern int Phys2Unk_first[NEQ_TYPE];
extern int Nnodes_per_proc;
extern int Nseg_tot;
extern int Ncomp;
void setup_polymer_field_wjdc(double **xOwned);
