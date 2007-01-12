/* This file was automatically generated.  Do not edit! */
int stencil_deltaLogical(int sten);
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
#define DELTA_FN_BOND  7
#define THETA_FN_SIG   6
#define WTC          3
#define THETA_CHARGE   3
#define DELTAC     1 
extern int Type_coul;
extern int Type_attr;
#define THETA_FN       1
#define NONE       -1
#define NONE      -1
#define NONE -1
extern int Type_func;
#define POLYMER_CR     4
#define U_ATTRACT      2
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
#define DELTA_FN       0
#define CMS_SCFT     2
#define CMS          0
extern int Type_poly;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define NSTEN        8
extern int Sten_Type[NSTEN];
void setup_stencil_logicals();
