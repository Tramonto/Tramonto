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
#include "az_aztec_defs.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define NMER_MAX     200
extern double BondWTC_b[NMER_MAX *NMER_MAX];
extern double Xi_cav_b[4];
extern double Rho_seg_b[NMER_MAX];
double integrand_WTC_freen_bulk(int iunk,int inode_box,double **x);
#define NCOMP_MAX 5
extern double Fac_overlap[NCOMP_MAX][NCOMP_MAX];
#define BONDWTC        5
extern int **Poly_to_Unk_SegAll;
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
double y_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
extern int **Bonds_SegAll;
extern int *Nbonds_SegAll;
#define CAVWTC         4
#define DENSITY        0
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
extern int Unk2Comp[NMER_MAX];
double integrand_WTC_freen(int iunk,int inode_box,double **x);
