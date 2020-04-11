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
extern int **Poly_to_Unk_SegAll;
extern int *Nbonds_SegAll;
#define NCOMP_MAX 6
extern int Icomp_to_polID[NCOMP_MAX];
#define VERBOSE      3 
extern int Iwrite;
void chempot_chain_wjdc(double *rho,double *betamu_chain,double *field_WJDC,double *g_WJDC);
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
double y_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
extern int Geqn_start[NCOMP_MAX];
extern int ***Poly_to_Unk;
extern int **Nbond;
extern int ***Bonds;
extern int *Unk_to_Bond;
#define NMER_MAX     200
extern int SegChain2SegAll[NCOMP_MAX][NMER_MAX];
extern int *Unk_to_Seg;
extern int *Unk_to_Poly;
extern int Nbonds;
extern double Scale_fac_WJDC[NCOMP_MAX][NCOMP_MAX];
double chain_term(int kseg,int kcomp,double *rho_seg,double *xi_cav);
extern double Avdw[NCOMP_MAX][NCOMP_MAX];
extern int Ncomp;
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_attr;
extern double Fac_overlap_hs[NCOMP_MAX];
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
#define PI    3.141592653589793238462643383279502884197169399375
extern double HS_diam[NCOMP_MAX];
extern int Unk2Comp[NMER_MAX];
extern int Grafted[NCOMP_MAX];
extern int Grafted_Logical;
extern int SegAll_to_Poly[NMER_MAX];
extern int Nseg_tot;
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
#define NBOND_MAX 4
extern double G_WJDC_RTF[NMER_MAX *NBOND_MAX];
extern double Field_WJDC_RTF[NMER_MAX];
extern double Xi_cav_RTF[4];
extern double Rho_seg_RTF[NMER_MAX];
extern double Rho_b_RTF[NCOMP_MAX];
extern double Dphi_Drhobar_RTF[10];
extern double G_WJDC_LBB[NMER_MAX *NBOND_MAX];
extern double Field_WJDC_LBB[NMER_MAX];
extern double Xi_cav_LBB[4];
extern double Rho_seg_LBB[NMER_MAX];
extern double Rho_b_LBB[NCOMP_MAX];
extern double Dphi_Drhobar_LBB[10];
extern double G_WJDC_b[NMER_MAX *NBOND_MAX];
extern double Field_WJDC_b[NMER_MAX];
extern double Xi_cav_b[4];
extern double Rho_seg_b[NMER_MAX];
extern double Rho_b[NCOMP_MAX];
extern double Dphi_Drhobar_b[10];
void compute_bulk_nonlocal_wjdc_properties(char *file_echoinput,double *dphi_drhobar,double *rho,double *rho_seg,double *xi_cav,double *field_WJDC,double *g_WJDC);
#define UNIFORM_INTERFACE  0
extern int Type_interface;
void WJDC_thermo_precalc(char *file_echoinput);
