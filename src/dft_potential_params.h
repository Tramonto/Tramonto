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
#define NWALL_MAX_TYPE 50 
extern double YukawaK_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Cut_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Eps_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Sigma_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Cut_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double YukawaK_w[NWALL_MAX_TYPE];
extern double YukawaK_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern int Type_uwwPot;
extern int Vext_PotentialID[NWALL_MAX_TYPE];
extern double Eps_w[NWALL_MAX_TYPE];
extern double Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Sigma_w[NWALL_MAX_TYPE];
extern double Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern int Nwall_type;
extern double YukawaK_ff[NCOMP_MAX][NCOMP_MAX];
#define PAIR_r18andYUKAWA_CS  8
#define PAIR_r12andYUKAWA_CS  7
#define PAIR_LJandYUKAWA_CS   6
#define PAIR_EXP_CS	      4
#define PAIR_YUKAWA_CS        3
extern int Type_pairPot;
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
extern double Bond_ff[NCOMP_MAX][NCOMP_MAX];
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_poly;
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
#define IDEAL_GAS    0
extern int Ipot_ff_n;
extern int Ncomp;
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
#define FILES_DEBUG        2
extern int Iwrite_files;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
void pot_parameters(char *file_echoinput);
