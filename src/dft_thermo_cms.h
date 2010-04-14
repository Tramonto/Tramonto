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
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
#define NCOMP_MAX 5
extern int Geqn_start[NCOMP_MAX];
extern int ***Poly_to_Unk;
extern int **Nbond;
extern int ***Bonds;
extern int *Unk_to_Bond;
#define NMER_MAX     200
extern int SegChain2SegAll[NCOMP_MAX][NMER_MAX];
extern int *Unk_to_Seg;
extern int *Unk_to_Poly;
extern int Nbonds;
extern double Avdw[NCOMP_MAX][NCOMP_MAX];
extern int Unk2Comp[NMER_MAX];
extern int SegAll_to_Poly[NMER_MAX];
extern int Nseg_tot;
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
#define NBOND_MAX 4
extern int Type_pairPot;
double pairPot_ATT_CS_switch(double r,int icomp,int jcomp,int typePairPot);
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern double Cr_rad[NCOMP_MAX][NCOMP_MAX];
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_attr;
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
extern int Ncomp;
void setup_polymer_cr();
extern double G_CMS_RTF[NMER_MAX *NBOND_MAX];
extern double Field_CMS_RTF[NMER_MAX];
extern double Rho_b_RTF[NCOMP_MAX];
extern double G_CMS_LBB[NMER_MAX *NBOND_MAX];
extern double Field_CMS_LBB[NMER_MAX];
extern double Rho_b_LBB[NCOMP_MAX];
extern double G_CMS_b[NMER_MAX *NBOND_MAX];
extern double Field_CMS_b[NMER_MAX];
extern double Rho_b[NCOMP_MAX];
void compute_bulk_nonlocal_cms_properties(char *output_file1,double *rho,double *field_CMS,double *g_CMS);
#define UNIFORM_INTERFACE  0
extern int Type_interface;
void CMS_thermo_precalc(char *output_file1);
