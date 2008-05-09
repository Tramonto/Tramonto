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
extern double Field_WJDC_b[NMER_MAX];
extern double Rho_seg_b[NMER_MAX];
double integrand_WJDC_freen_bulk(int iunk,int inode_box,double **x);
#define NCOMP_MAX 5
extern double Scale_fac_WJDC[NCOMP_MAX][NCOMP_MAX];
extern int **Bonds_SegAll;
extern int *Nbonds_SegAll;
extern int SegAll_to_Poly[NMER_MAX];
#define WJDC_FIELD     8
extern int Unk2Comp[NMER_MAX];
#define DENSITY        0
#define NEQ_TYPE       11 
extern int Phys2Unk_first[NEQ_TYPE];
double integrand_WJDC_freen(int iunk,int inode_box,double **x);
