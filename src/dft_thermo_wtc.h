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
#define NCOMP_MAX 5
#define NMER_MAX     100
extern int SegChain2SegAll[NCOMP_MAX][NMER_MAX];
extern int Type_mer[NCOMP_MAX][NMER_MAX];
extern int *Unk_to_Bond;
extern int ***Bonds;
extern int *Unk_to_Poly;
extern int *Unk_to_Seg;
extern double BondWTC_RTF[NMER_MAX *NMER_MAX];
extern double BondWTC_LBB[NMER_MAX *NMER_MAX];
#define NO_BOND_PAIR -962.0
extern double BondWTC_b[NMER_MAX *NMER_MAX];
extern int Nbonds;
#define WTC          2
extern int Type_poly;
#define VERBOSE      3 
extern int Iwrite;
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
extern double Xi_cav_RTF[4];
extern double Xi_cav_LBB[4];
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
extern double Fac_overlap_hs[NCOMP_MAX];
extern double Bond_ff[NCOMP_MAX][NCOMP_MAX];
extern int Ncomp;
#define PI    M_PI
double dy_dxi3_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
double dy_dxi2_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
extern double Fac_overlap[NCOMP_MAX][NCOMP_MAX];
extern double Xi_cav_b[4];
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
double y_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
extern double Betamu_seg[NMER_MAX];
extern double Betamu_wtc[NMER_MAX];
void chempot_WTC(double *rho_seg,double *betamu);
double chain_term(int kseg,int kcomp,double *rho_seg);
extern int **Bonds_SegAll;
extern int *Nbonds_SegAll;
double pressure_WTC(double *rho_seg);
void compute_bulk_nonlocal_wtc_properties(char *output_file1);
extern double Rho_b_RTF[NCOMP_MAX];
extern double Rho_seg_RTF[NMER_MAX];
extern double Rho_b_LBB[NCOMP_MAX];
extern double Rho_seg_LBB[NMER_MAX];
extern int Lsteady_state;
#define NBLOCK_MAX   5
extern int Nmer_t_total[NBLOCK_MAX];
extern int Unk2Comp[NMER_MAX];
extern double Rho_b[NCOMP_MAX];
extern double Rho_seg_b[NMER_MAX];
extern int Nseg_tot;
void WTC_overlap();
void WTC_thermo_precalc(char *output_file1);
