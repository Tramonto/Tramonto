/* This file was automatically generated.  Do not edit! */
double dmu_drho_hs_PY(double *rho);
double dp_drho_hs_PY(double *rho);
void FMT1stDerivBulk_switch(double *n,double *inv_n3,double *dphi_drb);
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
#define VERBOSE      3 
extern int Iwrite;
void dphi_drb_bulk(double *rhobar,double *dphi_drb);
#define NCOMP_MAX 5
extern double Rho_b[NCOMP_MAX];
extern double Rho_b_RTF[NCOMP_MAX];
extern double Rho_b_LBB[NCOMP_MAX];
#define UNIFORM_INTERFACE  0
extern int Type_interface;
#define NMER_MAX     100
extern double Rho_seg_b[NMER_MAX];
void rhobar_icomp(double rho,int icomp,double *rhobar);
extern int Unk2Comp[NMER_MAX];
extern int Nseg_tot;
extern int Lseg_densities;
extern double Dphi_Drhobar_RTF[10];
extern double Dphi_Drhobar_LBB[10];
extern double Dphi_Drhobar_b[10];
extern double Rhobar_b_RTF[10];
extern double Rhobar_b_LBB[10];
extern double Rhobar_b[10];
extern int Nrho_bar;
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
void chempot_PY_hs(double *rho);
double pressure_PY_hs(double *rho);
extern double Betamu_hs_ex[NCOMP_MAX];
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
void chempot_FMT_hs(double *dphi_drhobar);
double phispt_switch(double *n);
void solutionVec_to_nOrdering(double *rhoBar_SVOrdering,double *n);
extern int Ndim;
extern int Nrho_bar_s;
#define NDIM_MAX  3
double pressure_FMT_hs(double *rhobar,double *dphi_drhobar);
extern int Type_pairPot;
double pairPot_switch(double r,double param1,double param2,double param3,double param4,int typePairPot);
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
double integrand_BH(double r,int icomp);
#define BH_DIAM             1
extern int Type_hsdiam;
void calc_HS_diams();
extern double Inv_4pirsq[NCOMP_MAX];
extern double Inv_4pir[NCOMP_MAX];
extern double HS_diam[NCOMP_MAX];
extern double Inv_rad[NCOMP_MAX];
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
#define PI    M_PI
extern double Inv_4pi;
void compute_bulk_FMT_properties(char *output_file1);
void calc_InvR_params();
extern double Fac_overlap_hs[NCOMP_MAX];
extern int Ncomp;
#define WTC          2
extern int Type_poly;
void HS_thermo_precalc(char *output_file1);
