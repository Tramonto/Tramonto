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
extern double Dielec_X;
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int Ndim;
#define NWALL_MAX_TYPE 20 
extern double Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define COULOMB      1
extern int Ipot_ff_c;
void surface_atoms_inSurfaceTest(int iwall,int iwall_type,double *fluid_testpos,double **wall_pos,double dist_adjustments,int flag_X_to_center,double *delr_vext,double *delr_zone,int *logical_inwall,int *logical_nearWallDielec);
