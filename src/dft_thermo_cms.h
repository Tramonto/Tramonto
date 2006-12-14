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
#define CMS_SCFT     2
double uLJatt_n(double r,int i,int j);
#define NCOMP_MAX 5
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern double Cr_rad[NCOMP_MAX][NCOMP_MAX];
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
#define LJ12_6       2
extern int Ipot_ff_n;
extern double ***Rism_cr;
extern int Last_nz_cr;
#define N_NZCR_MAX   200   /* maximum # of non-zero's in direct correlation fn */
extern double Cr_rad_hs[NCOMP_MAX][NCOMP_MAX];
extern double Deltar_cr;
extern char Cr_file4[40];
extern char Cr_file3[40];
extern char Cr_file2[40];
extern char Cr_file[40];
#define CMS          0
extern int Type_poly;
#define VERBOSE      3 
extern int Iwrite;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern double Cr_break[2];
extern double Crfac;
extern int Ncr_files;
extern double Rho_b[NCOMP_MAX];
extern int Ncomp;
void setup_polymer_cr();
