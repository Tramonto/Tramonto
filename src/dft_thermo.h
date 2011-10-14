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
extern int Physics_scaling;
#define NCOMP_MAX 5
extern double Betamu_id[NCOMP_MAX];
#define NMER_MAX     200
#define NBOND_MAX 4
extern double G_WJDC_b[NMER_MAX *NBOND_MAX];
extern double Field_WJDC_b[NMER_MAX];
void print_to_file_comp(FILE *fp,int icomp,double val,char *var_label,int first);
void print_to_screen_comp(int icomp,double val,char *var_label);
extern double Betamu_wtc_RTF[NMER_MAX];
extern double Betamu_seg_RTF[NMER_MAX];
extern int SegChain2SegAll[NCOMP_MAX][NMER_MAX];
extern int Nmer[NCOMP_MAX];
extern double Betamu_wtc_LBB[NMER_MAX];
extern double Betamu_seg[NMER_MAX];
extern double Betamu_seg_LBB[NMER_MAX];
extern int Nseg_tot;
void chempot_WTC(double *rho_seg,double *betamu,double *xi_cav);
extern double Betamu[NCOMP_MAX];
extern double *Deltac_b;
void chempot_ELEC_MSA_GENERAL(double *rho);
void chempot_ELEC_MSA_RPM(double *rho);
#define DELTAC_GENERAL 2
#define DELTAC_RPM     1 
extern double Elec_pot_RTF;
extern double Elec_pot_LBB;
extern double Charge_f[NCOMP_MAX];
extern double Betamu_att[NCOMP_MAX];
void chempot_att(double *rho);
extern double Betamu_hs_ex[NCOMP_MAX];
void chempot_FMT_hs(double *dphi_drhobar);
extern double Betamu_RTF[NCOMP_MAX];
extern double Betamu_LBB[NCOMP_MAX];
void chempot_ideal_gas(double *rho,double *betamu);
extern double Betamu_chain[NMER_MAX];
#define PHASE_INTERFACE 2
extern double G_WJDC_RTF[NMER_MAX *NBOND_MAX];
extern double Field_WJDC_RTF[NMER_MAX];
extern double Betamu_chain_RTF[NMER_MAX];
extern double G_WJDC_LBB[NMER_MAX *NBOND_MAX];
extern double Field_WJDC_LBB[NMER_MAX];
extern double Betamu_chain_LBB[NMER_MAX];
void chempot_chain_wjdc(double *rho,double *betamu_chain,double *field_WJDC,double *g_WJDC);
extern double Xi_cav_b[4];
extern double Dphi_Drhobar_b[10];
extern double Rhobar_b[10];
extern double Rho_b[NCOMP_MAX];
extern double Rho_seg_b[NMER_MAX];
extern double Betap;
void print_to_file(FILE *fp,double val,char *var_label,int first);
void print_to_screen(double val,char *var_label);
extern double Xi_cav_RTF[4];
extern double Xi_cav_LBB[4];
double pressure_WTC(double *rho_seg,double *xi_cav);
extern int Type_coul;
double pressure_att(double *rho);
extern double Dphi_Drhobar_RTF[10];
extern double Rhobar_b_RTF[10];
extern double Dphi_Drhobar_LBB[10];
extern double Rhobar_b_LBB[10];
double pressure_FMT_hs(double *rhobar,double *dphi_drhobar);
extern double Rho_b_RTF[NCOMP_MAX];
extern double Rho_b_LBB[NCOMP_MAX];
extern double Rho_seg_RTF[NMER_MAX];
extern double Betap_RTF;
extern double Rho_seg_LBB[NMER_MAX];
double pressure_ideal_gas(double *rho);
extern double Betap_LBB;
extern int Lseg_densities;
#define UNIFORM_INTERFACE  0
extern int Type_interface;
void calc_pressure(char *output_file1,int iwrite);
void calc_chempot(char *output_file1,int iwrite);
void WJDC_thermo_precalc(char *output_file1);
extern double Scale_fac_WJDC[NCOMP_MAX][NCOMP_MAX];
extern int Npol_comp;
extern int Ncomp;
void ATT_thermo_precalc();
extern int Type_attr;
void HS_thermo_precalc(char *output_file1);
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_func;
void WTC_thermo_precalc(char *output_file1);
#define WJDC3        5 
#define WJDC2        4 
#define WJDC         3
#define WTC          2
extern int Type_poly;
extern int L_HSperturbation;
#define NO_SCREEN    4 
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
void thermodynamics(char *output_file1,int iwrite);
