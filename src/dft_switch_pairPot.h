/* This file was automatically generated.  Do not edit! */
double uSW_Integral(double r,int i,int j);
double uEXP_Integral(double r,int i,int j);
double urNandYUKAWA_Integral(double r,int i,int j);
double ur18andYUKAWA_Integral(double r,int i,int j);
double ur12andYUKAWA_Integral(double r,int i,int j);
double uLJandYUKAWA_Integral(double r,int i,int j);
double uYUKAWA_Integral(double r,int i,int j);
double uCOULOMB_Integral(double r,int i,int j);
double uLJ12_6_Integral(double r,int i,int j);
double pairPot_integral_switch(double r,int icomp,int jcomp,int typePairPot);
double uSW_ATT_noCS(double r,int i,int j);
double uEXP_ATT_noCS(double r,int i,int j);
double urNandYUKAWA_ATT_noCS(double r,int i,int j);
double ur18andYUKAWA_ATT_noCS(double r,int i,int j);
double ur12andYUKAWA_ATT_noCS(double r,int i,int j);
double uLJandYUKAWA_ATT_noCS(double r,int i,int j);
double uYUKAWA_ATT_noCS(double r,int i,int j);
double uCOULOMB_ATT_noCS(double r,int i,int j);
double uLJ12_6_ATT_noCS(double r,int i,int j);
double pairPot_ATT_noCS_switch(double r,int icomp,int jcomp,int typePairPot);
double uSW_ATT_CS(double r,int i,int j);
double uEXP_ATT_CS(double r,int i,int j);
double urNandYUKAWA_ATT_CS(double r,int i,int j);
double ur18andYUKAWA_ATT_CS(double r,int i,int j);
double ur12andYUKAWA_ATT_CS(double r,int i,int j);
double uLJandYUKAWA_ATT_CS(double r,int i,int j);
double uYUKAWA_ATT_CS(double r,int i,int j);
double uCOULOMB_ATT_CnoS(double r,int i,int j);
double uCOULOMB_ATT_CS(double r,int i,int j);
double uLJ12_6_ATT_CS(double r,int i,int j);
double pairPot_ATT_CS_switch(double r,int icomp,int jcomp,int typePairPot);
void uSW_InnerCore(int i,int j,double *rCore_left,double *rCore_right,double *epsCore);
void uEXP_InnerCore(int i,int j,double *rCore_left,double *rCore_right,double *epsCore);
void urNandYUKAWA_InnerCore(int i,int j,double *rCore_left,double *rCore_right,double *epsCore);
void ur18andYUKAWA_InnerCore(int i,int j,double *rCore_left,double *rCore_right,double *epsCore);
void ur12andYUKAWA_InnerCore(int i,int j,double *rCore_left,double *rCore_right,double *epsCore);
void uLJandYUKAWA_InnerCore(int i,int j,double *rCore_left,double *rCore_right,double *epsCore);
void uYUKAWA_InnerCore(int i,int j,double *rCore_left,double *rCore_right,double *epsCore);
void uCOULOMB_InnerCore(int i,int j,double *rCore_left,double *rCore_right,double *epsCore);
void uLJ12_6_InnerCore(int i,int j,double *rCore_left,double *rCore_right,double *epsCore);
void pairPot_InnerCore_switch(int icomp,int jcomp,int typePairPot,double *rCore_left,double *rCore_right,double *epsCore);
double uSW_DERIV1D(double r,double x,double sigma,double eps,double rcut);
double uEXP_DERIV1D(double r,double x,double sigma,double eps,double rcut,double yukawaK);
double urNandYUKAWA_DERIV1D(double r,double x,double sigma,double eps,double rcut,double yukawaK,double AYukawa,double npow);
double ur18andYUKAWA_DERIV1D(double r,double x,double sigma,double eps,double rcut,double yukawaK,double AYukawa);
double ur12andYUKAWA_DERIV1D(double r,double x,double sigma,double eps,double rcut,double yukawaK,double AYukawa);
double uLJandYUKAWA_DERIV1D(double r,double x,double sigma,double eps,double rcut,double yukawaK,double Ayukawa);
double uYUKAWA_DERIV1D(double r,double x,double sigma,double eps,double rcut,double yukawaK);
double uCOULOMB_DERIV1D(double r,double x,double z1,double z2);
double uCOULOMB_CS_DERIV1D(double r,double x,double z1,double z2,double rcut);
double uLJ12_6_DERIV1D(double r,double x,double sigma,double eps,double rcut);
double pairPot_deriv_switch(double r,double x,double param1,double param2,double param3,double param4,double param5,double param6,int typePairPot);
void uSW_setparams(int context,int i,int j,double *param1,double *param2,double *param3);
void uEXP_CS_setparams(int context,int i,int j,double *param1,double *param2,double *param3,double *param4);
void urNandYUKAWA_CS_setparams(int context,int i,int j,double *param1,double *param2,double *param3,double *param4,double *param5,double *param6);
void ur18andYUKAWA_CS_setparams(int context,int i,int j,double *param1,double *param2,double *param3,double *param4,double *param5);
void ur12andYUKAWA_CS_setparams(int context,int i,int j,double *param1,double *param2,double *param3,double *param4,double *param5);
void uLJandYUKAWA_CS_setparams(int context,int i,int j,double *param1,double *param2,double *param3,double *param4,double *param5);
void uYUKAWA_CS_setparams(int context,int i,int j,double *param1,double *param2,double *param3,double *param4);
void uCOULOMB_setparams(int context,int i,int j,double *param1,double *param2,double *param3);
void uCOULOMB_CS_setparams(int context,int i,int j,double *param1,double *param2,double *param3);
void uLJ12_6_CS_setparams(int context,int i,int j,double *param1,double *param2,double *param3);
void pairPotparams_switch(int typePairPot,int context,int i,int j,double *param1,double *param2,double *param3,double *param4,double *param5,double *param6);
double uSW(double r,double sigma,double eps,double rcut);
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
#define PAIR_SW		      5
double uEXP_CS(double r,double sigma,double eps,double rcut,double yukawaK);
#define PAIR_EXP_CS	      4
double urNandYUKAWA_CS(double r,double sigma,double eps,double rcut,double yukawaK,double AYukawa,double npow);
#define PAIR_rNandYUKAWA_CS   9
double ur18andYUKAWA_CS(double r,double sigma,double eps,double rcut,double yukawaK,double AYukawa);
#define PAIR_r18andYUKAWA_CS  8
double ur12andYUKAWA_CS(double r,double sigma,double eps,double rcut,double yukawaK,double AYukawa);
#define PAIR_r12andYUKAWA_CS  7
double uLJandYUKAWA_CS(double r,double sigma,double eps,double rcut,double yukawaK,double Ayukawa);
#define PAIR_LJandYUKAWA_CS   6
double uYUKAWA_CS(double r,double sigma,double eps,double rcut,double yukawaK);
#define PAIR_YUKAWA_CS        3
double uCOULOMB(double r,double z1,double z2);
#define PAIR_COULOMB          2
double uCOULOMB_CS(double r,double z1,double z2,double rcut);
#define PAIR_COULOMB_CS       1
double uLJ12_6_CS(double r,double sigma,double eps,double rcut);
#define PAIR_LJ12_6_CS        0
double pairPot_switch(double r,double param1,double param2,double param3,double param4,double param5,double param6,int typePairPot);
