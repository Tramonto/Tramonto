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
#define NMER_MAX     40
extern double Rho_seg_b[NMER_MAX];
#define NCOMP_MAX 5
extern double Betamu_id[NCOMP_MAX];
extern double Betamu[NCOMP_MAX];
void print_to_file_comp(FILE *fp,int icomp,double val,char *var_label,int first);
void print_to_screen_comp(int icomp,double val,char *var_label);
extern double Betamu_wtc_RTF[NMER_MAX];
extern double Betamu_seg_RTF[NMER_MAX];
extern double Rho_seg_RTF[NMER_MAX];
extern double Betamu_wtc_LBB[NMER_MAX];
extern double Betamu_seg[NMER_MAX];
extern double Betamu_seg_LBB[NMER_MAX];
extern int Nseg_tot;
extern double Rho_seg_LBB[NMER_MAX];
void chempot_WTC(double *rho_seg,double *betamu);
extern double *Deltac_b;
void chempot_ELEC_MSA(double *rho);
#define THETA_CHARGE   3
#define NSTEN        8
extern int Sten_Type[NSTEN];
extern double Elec_pot_RTF;
extern double Elec_pot_LBB;
extern double Charge_f[NCOMP_MAX];
extern double Betamu_att[NCOMP_MAX];
void chempot_att(double *rho);
extern double Betamu_hs_ex[NCOMP_MAX];
void chempot_FMT_hs(double *rho);
extern double Betamu_RTF[NCOMP_MAX];
extern double Betamu_LBB[NCOMP_MAX];
void chempot_ideal_gas(double *rho,double *betamu);
extern int Ncomp;
double pressure_PY_hs(double *rho);
extern double Rho_b[NCOMP_MAX];
extern double Betap;
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
void print_to_file(FILE *fp,double val,char *var_label,int first);
void print_to_screen(double val,char *var_label);
extern int Type_coul;
double pressure_att(double *rho);
double pressure_FMT_hs(double *rho);
extern double Rho_b_RTF[NCOMP_MAX];
extern double Betap_RTF;
extern double Rho_b_LBB[NCOMP_MAX];
double pressure_ideal_gas(double *rho);
extern double Betap_LBB;
extern int Lsteady_state;
void calc_chempot(char *output_file1);
void calc_pressure(char *output_file1);
void ATT_thermo_precalc();
extern int Type_attr;
void HS_thermo_precalc(char *output_file1);
extern int Type_func;
void WTC_thermo_precalc(char *output_file1);
#define WTC          3
#define NONE       -1
#define NONE      -1
#define NONE -1
extern int Type_poly;
#define NO_SCREEN    2 
extern int Iwrite;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
void thermodynamics(char *output_file1);
