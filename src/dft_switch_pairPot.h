/* This file was automatically generated.  Do not edit! */
double uCOULOMB_Integral(double r,int i,int j);
double uLJ12_6_Integral(double r,int i,int j);
double pairPot_integral_switch(double r,int icomp,int jcomp,int typePairPot);
double uCOULOMB_ATT_noCS(double r,int i,int j);
double uLJ12_6_ATT_noCS(double r,int i,int j);
double pairPot_ATT_noCS_switch(double r,int icomp,int jcomp,int typePairPot);
double uCOULOMB_ATT_CnoS(double r,int i,int j);
double uCOULOMB_ATT_CS(double r,int i,int j);
double uLJ12_6_ATT_CS(double r,int i,int j);
double pairPot_ATT_CS_switch(double r,int icomp,int jcomp,int typePairPot);
double uCOULOMB_DERIV1D(double r,double x,double z1,double z2);
double uCOULOMB_CS_DERIV1D(double r,double x,double z1,double z2,double rcut);
double uLJ12_6_DERIV1D(double r,double x,double sigma,double eps,double rcut);
double pairPot_deriv_switch(double r,double x,double param1,double param2,double param3,int typePairPot);
double uCOULOMB(double r,double z1,double z2);
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
#define PAIR_COULOMB       2
double uCOULOMB_CS(double r,double z1,double z2,double rcut);
#define PAIR_COULOMB_CS    1
double uLJ12_6_CS(double r,double sigma,double eps,double rcut);
#define PAIR_LJ12_6_CS     0
double pairPot_switch(double r,double param1,double param2,double param3,int typePairPot);
