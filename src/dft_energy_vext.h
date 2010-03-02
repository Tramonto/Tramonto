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
extern double **Vext_coul;
double integrand_vext_elec_freen(int iunk,int inode_box,double **x);
extern int *B2L_node;
extern double **Vext;
#define DENSITY_MIN  1.e-20
#define NMER_MAX     200
extern int Unk2Comp[NMER_MAX];
extern int Lseg_densities;
#define DENSITY        0
#define NEQ_TYPE       13 
extern int Phys2Unk_first[NEQ_TYPE];
double integrand_vext_freen(int iunk,int inode_box,double **x);
