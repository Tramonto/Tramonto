/* This file was automatically generated.  Do not edit! */
double pairPot_ATT_CS_switch(double r,int icomp,int jcomp,int typePairPot);
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
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
double pairPot_switch(double r,double param1,double param2,double param3,double param4,double param5,double param6,int typePairPot);
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
void print_potentials_fluid(int type_pairPot,int icomp,int jcomp);
#define EXTENDED     2
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
#define VERBOSE      3 
extern int Iwrite;
double pairPot_find_r_ZeroCut(int i,int j,double param1,double param2,double param3,double param4,double param5,double param6,int typePairPot);
extern double Rzero_ff[NCOMP_MAX][NCOMP_MAX];
double pairPot_find_rmin(int i,int j,double param1,double param2,double param3,double param4,double param5,double param6,int typePairPot);
extern double Rmin_ff[NCOMP_MAX][NCOMP_MAX];
#define FLUID_FLUID 0
extern int Type_pairPot;
void pairPotparams_switch(int typePairPot,int context,int i,int j,double *param1,double *param2,double *param3,double *param4,double *param5,double *param6);
extern int Ncomp;
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_attr;
void pot_parameters(char *output_file1);
extern int Mix_type;
void setup_pairPotentials(char *output_file1);
