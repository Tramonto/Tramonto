/* This file was automatically generated.  Do not edit! */
void print_cont_variable_archived_plugin(int cont_type,FILE *fp,int Loca_contID);
void print_cont_type_archived_plugin(int cont_type,FILE *fp,int Loca_contID);
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
#define CMS_SCFT     1
void calc_InvR_params();
void calc_HS_diams();
#define BH_DIAM             1
extern int Type_hsdiam;
void scale_vext_epswf(double ratio,int icomp,int iwall);
#define NWALL_MAX 600 
extern int WallType[NWALL_MAX];
extern int Nwall;
void setup_pairPotentials(char *file_echoinput);
extern int Nwall_type;
extern int Ncomp;
#define FILES_BASIC        0
void thermodynamics(char *file_echoinput,int iwrite_screen,int iwrite_files);
void recalculate_stencils();
void setup_polymer_cr();
#define CMS          0
#define NWALL_MAX_TYPE 10 
#define NCOMP_MAX 5
void assign_param_archived_plugin(int cont_type,int Loca_contID,double param);
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
extern double Crfac;
#define CONT_CRFAC              108  /* continuous mixing of two cr files */
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
#define CONT_EPSFF_ALL		107
#define CONT_EPSWF_SOME 	106
extern double Eps_wf[NCOMP_MAX][NWALL_MAX_TYPE];
#define CONT_EPSWF_ALL	        105
extern double Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Eps_w[NWALL_MAX_TYPE];
extern int Mix_type;
#define CONT_EPSW_ALL		104
#define CONT_RHO_ALL		103
#define CONT_RHO_CONST_XSOLV    102
#define CONT_RHO_CONST_RHOTOT58 101
#define NMER_MAX     200
extern double Rho_seg_b[NMER_MAX];
extern int SegAll_to_Poly[NMER_MAX];
extern int Nseg_tot;
#define NCONT_MAX          2 /* the maximum number of solutions possible for use with Loca */
extern int Cont_ID[NCONT_MAX][2];
extern double Rho_b[NCOMP_MAX];
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_poly;
#define CONT_LOG_RHO_I          100
double get_init_param_archived_plugin(int cont_type,int Loca_contID);
