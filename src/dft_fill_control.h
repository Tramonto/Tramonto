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
extern int Nunk_per_node;
#define NODAL_FLAG -999
void fill_resid_and_matrix(double **x,int iter,int resid_only_flag,int unk_flag);
#define MATRIX_FILL_NODAL 1   /* set to zero for physics based ordering */
void fill_resid_and_matrix_control(double **x,int iter,int resid_only_flag);
