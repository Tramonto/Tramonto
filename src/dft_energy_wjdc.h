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
#define NCOMP_MAX 5
extern int Icomp_to_polID[NCOMP_MAX];
extern int Grafted[NCOMP_MAX];
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Grafted_Logical;
double integrand_WJDCcomp_freen_bulk(int iunk,int inode_box,double **x);
#define NMER_MAX     200
extern double Betamu_chain[NMER_MAX];
#define DIFFUSION      6
#define DIFFUSIVE_INTERFACE 1
extern int Type_interface;
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
extern int **Poly_to_Unk_SegAll;
#define G_CHAIN       11 
extern int SegChain2SegAll[NCOMP_MAX][NMER_MAX];
extern int Type_mer[NCOMP_MAX][NMER_MAX];
extern int **Nseg_type_pol;
extern int Ncomp;
extern int Nmer[NCOMP_MAX];
#define NBLOCK_MAX   20 
extern int Nmer_t[NCOMP_MAX][NBLOCK_MAX];
double integrand_WJDCcomp_freen(int iunk,int inode_box,double **x);
extern double Field_WJDC_b[NMER_MAX];
extern double Rho_seg_b[NMER_MAX];
double integrand_WJDC_freen_bulk(int iunk,int inode_box,double **x);
extern double Scale_fac_WJDC[NCOMP_MAX][NCOMP_MAX];
extern int **Bonds_SegAll;
extern int *Nbonds_SegAll;
extern int SegAll_to_Poly[NMER_MAX];
#define WJDC         3
#define WJDC_FIELD     8
#define WJDC3        5 
#define WJDC2        4 
extern int Type_poly;
extern int Unk2Comp[NMER_MAX];
#define DENSITY        0
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
double integrand_WJDC_freen(int iunk,int inode_box,double **x);
extern double *Gsum_graft;
extern double *Total_area_graft;
extern double Rho_g[NCOMP_MAX];
extern int Grafted_SegID[NCOMP_MAX];
#define GRAFT_DENSITY 1


