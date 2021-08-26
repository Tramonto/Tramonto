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
#define NCOMP_MAX 6
extern double Betamu_att[NCOMP_MAX];
extern double Rho_b[NCOMP_MAX];
extern double Rho_b_RTF[NCOMP_MAX];
#define NMER_MAX     200
extern double Rho_seg_b[NMER_MAX];
extern double Rho_seg_RTF[NMER_MAX];
#define UNIFORM_INTERFACE  0
extern int Type_interface;
extern int Unk2Comp[NMER_MAX];
extern int Lseg_densities;
extern int Icomp_to_polID[NCOMP_MAX];
extern int Grafted[NCOMP_MAX];
#define WJDC3        5 
extern int Type_poly;
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Grafted_Logical;
#define DENSITY        0
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
double integrand_att_freen_bulk(int iunk,int inode_box,double **x);
#define THETA_PAIRPOT_RCUT    2
double int_stencil(double **x,int inode_box,int iunk,int sten_type);
double integrand_att_freen(int iunk,int inode_box,double **x);
