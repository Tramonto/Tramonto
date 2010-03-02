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
#define WTC          2
extern int Type_poly;
#define NCOMP_MAX 5
extern double Rho_b_RTF[NCOMP_MAX];
#define NMER_MAX     200
extern double Rho_seg_RTF[NMER_MAX];
extern double Betamu_RTF[NCOMP_MAX];
double integrand_mu_freen_bulk(int iunk,int inode_box,double **x);
#define DENSITY_MIN  1.e-20
extern double Rho_b[NCOMP_MAX];
extern double Charge_f[NCOMP_MAX];
extern double Betamu[NCOMP_MAX];
extern double Rho_seg_b[NMER_MAX];
#define NONE       -1
#define NONE      -1
#define NONE        -1
#define NONE        -1
extern int Type_coul;
extern double Betamu_seg[NMER_MAX];
extern int Unk2Comp[NMER_MAX];
extern int Lseg_densities;
#define DIFFUSION      6
#define DIFFUSIVE_INTERFACE 1
extern int Type_interface;
#define DENSITY        0
#define NEQ_TYPE       13 
extern int Phys2Unk_first[NEQ_TYPE];
double integrand_mu_freen(int iunk,int inode_box,double **x);
