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
extern int *L2B_node;
extern int Nnodes_per_proc;
#define SCREEN_VERBOSE     3 
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
void safe_free(void **ptr);
void safe_free(void **ptr);
void check_zero_densities_owned(double **xOwned);
void calc_Gsum_new(double **x);
void chop_profile(double **xInBox,int guess_type);
#define RESTART_STEP       2
void setup_polymer_G_wjdc(double **xOwned);
void calc_init_polymer_G_wjdc(double **xInBox,double **xOwned);
#define WJDC2        4 
#define WJDC         3
void calc_init_polymer_G_SCF(double **xInBox,double **xOwned);
void setup_polymer_G(double **xInBox,double **xOwned);
void calc_init_polymer_G_CMS(double **xInBox,double **xOwned);
#define CALC_RHOBAR_AND_G 3
#define G_CHAIN       11 
void calc_init_lambda(double **xInBox,double **xOwned);
#define SCF_CONSTR	   9
void setup_polymer_SCF_field(double **xInBox,double **xOwned,int guess_type);
void calc_init_SCFfield(double **xInBox,double **xOwned);
#define SCF_FIELD	  10
void setup_polymer_field(double **xInBox,double **xOwned,int guess_type);
void calc_init_CMSfield(double **xInBox,double **xOwned);
#define CMS_FIELD      7
void setup_polymer_field_wjdc(double **xOwned);
void calc_init_WJDC_field(double **xInBox,double **xOwned);
#define WJDC_FIELD     8
void setup_BondWTC(double **xOwned);
void calc_init_BondWTC(double **xInBox,double **xOwned);
#define BONDWTC        5
void setup_Xi_cavWTC(double **xOwned);
void calc_init_Xi_cavWTC(double **xInBox,double **xOwned);
#define CAVWTC         4
void calc_init_chem_pot(double **xInBox,double **xOwned);
void setup_chem_pot(double **xOwned);
#define DIFFUSION      6
void calc_init_elec_pot(double **xInBox,double **xOwned);
void setup_elec_pot(double **xOwned,int guess_type);
#define CALC_ALL_FIELDS   1
#define POISSON        1
void setup_rho_bar(double **xOwned);
void calc_init_rho_bar(double **xInBox,double **xOwned);
#define HSRHOBAR       2
void setup_mf_attract(double **xOwned);
void calc_init_mf_attract(double **xInBox,double **xOwned);
#define BULK              0
extern int Iguess_fields;
#define RESTART_DENSONLY   3
#define MF_EQ          3
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
void setup_density(double **xInBox,double **xOwned,int guess_type);
void setup_polymer_rho(double **xInBox,double **xOwned,int guess_type);
#define CMS_SCFT     1
#define CMS          0
#define RESTART_FEWERCOMP  4
#define NEQ_TYPE       12 
extern int Phys2Nunk[NEQ_TYPE];
#define DENSITY        0
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
void translate_xInBox_to_xOwned(double **xInBox,double **xOwned);
void guess_restart_from_files(int start_no_info,int guess_type,double **xInBox);
extern int Imain_loop;
#define NORESTART          0
extern int Restart;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Restart_field[NEQ_TYPE];
extern int Nnodes_box;
extern int Nnodes_box_extra;
extern int Nunk_per_node;
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
extern int Grafted_Logical;
#define WJDC3        5 
extern int Type_poly;
void set_initial_guess(int guess_type,double **xOwned);
