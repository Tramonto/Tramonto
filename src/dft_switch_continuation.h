/* This file was automatically generated.  Do not edit! */
void print_cont_variable_user_plugin(int cont_type,FILE *fp,int Loca_contID);
void print_cont_variable_archived_plugin(int cont_type,FILE *fp,int Loca_contID);
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
#define NWALL_MAX_TYPE 50 
extern double WallParam[NWALL_MAX_TYPE];
#define NDIM_MAX  3
extern double Size_x[NDIM_MAX];
#define NWALL_MAX 600 
extern double WallPos[NDIM_MAX][NWALL_MAX];
void print_cont_variable(int cont_type,FILE *fp,int Loca_contID);
void print_cont_type_user_plugin(int cont_type,FILE *fp,int Loca_contID);
void print_cont_type_archived_plugin(int cont_type,FILE *fp,int Loca_contID);
extern int Ndim;
#define REFLECT              2
extern int Type_bc[NDIM_MAX][2];
extern int Plane_new_nodes;
#define SWITCH_SURFACE_SEP   0
extern int Print_mesh_switch;
void print_cont_type(int cont_type,FILE *fp,int Loca_contID);
void thermodynamics(char *output_file1);
typedef struct Loca_Struct Loca_Struct;
struct Loca_Struct {
  int    method;      /* Continuation method                          */
  int    cont_type1;  /* flag specifying the continuation parameter   */
  int    cont_type2;  /* flag specifying the second (free) parameter  */
  int    num_steps;   /* maximum number of continuation steps to take */
  double aggr;        /* step size control parameter                  */
  double step_size;   /* initial continuation step size               */
};
extern struct Loca_Struct Loca;
void assign_param_user_plugin(int cont_type,int Loca_contID,double param);
void assign_param_archived_plugin(int cont_type,int Loca_contID,double param);
extern double **Vext;
extern int Nnodes_per_proc;
void scale_elec_param(double ratio);
#define BH_DIAM             1
extern int Type_hsdiam;
extern int WallType[NWALL_MAX];
void scale_vext_epswf(double ratio,int icomp,int iwall);
#define NMER_MAX     200
extern int Unk2Comp[NMER_MAX];
extern int Ntype_mer;
void scale_vext_temp(double ratio);
extern int Nwall;
void recalculate_stencils();
void setup_polymer_cr();
#define CMS          0
void calc_InvR_params();
void calc_HS_diams();
extern int Type_func;
void setup_pairPotentials(char *output_file1);
#define VEXT_HARD        1
extern int Ipot_wf_n[NWALL_MAX_TYPE];
extern int Nwall_type;
extern int Ncomp;
extern int Type_attr;
extern double Temp_elec;
#define COULOMB      1
extern int Ipot_ff_c;
#define NCOMP_MAX 5
void assign_parameter_tramonto(int cont_type,double param,int Loca_contID);
double get_init_param_user_plugin(int cont_type,int Loca_contID);
double get_init_param_archived_plugin(int cont_type,int Loca_contID);
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
#define CONT_SIGMAFF_IJ    10   /* Fluid-Fluid Interaction Diameter for IJ term */
extern double **Vext_membrane;
#define CONT_SEMIPERM_IJ   9  /* Vext_membrane */
#define CONT_ELECPARAM_ALL 8  /* Charged surface params */
extern double Elec_param_w[NWALL_MAX];
#define CONT_ELECPARAM_I   7  /* Charged surface params */
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
#define CONT_EPSFF_IJ      6   /* Fluid-Fluid Energy Params for IJ term */
extern double Eps_wf[NCOMP_MAX][NWALL_MAX_TYPE];
#define CONT_EPSWF_IJ      5    /* Wall-Fluid Energy Params for IJ term */
extern double Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Eps_w[NWALL_MAX_TYPE];
extern int Mix_type;
#define CONT_EPSW_I        4    /* Wall-Wall Energy Params for wall I */
extern double Betamu[NCOMP_MAX];
extern double Betamu_chain[NMER_MAX];
#define WJDC3        5 
#define WJDC2        4 
#define WJDC         3
#define CONT_BETAMU_I      3  /* Vary chemical potential for species I */
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
#define CONT_RHO_I         2
extern double Temp;
#define CONT_TEMP          1   /* State Parameters */
#define CONT_MESH          0   /* mesh size */
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if defined(LOCA_MF) && !defined(_CON_CONST_H_)
#include <mf.h>
#endif
double get_init_param_value(int cont_type,int);
#if !defined(_CON_CONST_H_)
double get_init_param_value(int cont_type,int);
#endif
double get_init_param_value(int cont_type,int Loca_contID);
