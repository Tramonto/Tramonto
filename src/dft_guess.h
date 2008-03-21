/* This file was automatically generated.  Do not edit! */
void communicate_profile(double *x_new,double **xInBox);
void collect_x_old(double **x);
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
extern int Nnodes;
extern double *X_old;
extern int *L2B_node;
extern int Nnodes_per_proc;
void safe_free(void **ptr);
void safe_free(void **ptr);
void translate_xInBox_to_xOwned(double **xInBox,double **xOwned);
void check_zero_densities(double **xInBox);
void chop_profile(double **xInBox,int iguess);
#define YW_DENS        10       /* densities for Yethiraj-Woodward polymer DFTs */
void setup_polymer_G_wjdc(double **xInBox);
#define WJDC         3
void setup_polymer_G(double **xInBox);
void setup_polymer_G(double **xInBox);
#define G_CHAIN        9 
void setup_polymer_field(double **xInBox,int iguess);
void setup_polymer_field(double **xInBox,int iguess);
#define CMS_FIELD      7
void communicate_to_fill_in_box_values(double **xInBox);
void setup_polymer_field_wjdc(double **xInBox);
#define WJDC_FIELD     8
void setup_BondWTC(double **xInBox);
#define BONDWTC        5
void setup_Xi_cavWTC(double **xInBox);
#define CAVWTC         4
void setup_chem_pot(double **xInBox);
#define DIFFUSION      6
void setup_elec_pot(double **xInBox,int iguess);
#define POISSON        1
void setup_rho_bar(double **xInBox);
#define HSRHOBAR       2
void setup_mf_attract(double **xInBox);
#define MF_EQ          3
void setup_density(double **xInBox,int iguess);
void setup_polymer_rho(double **xInBox,int iguess);
void setup_polymer_rho(double **xInBox,int iguess);
#define CMS_SCFT     1
#define CMS          0
extern int Type_poly;
#define NEQ_TYPE       11 
extern int Phys2Nunk[NEQ_TYPE];
#define DENSITY        0
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
void guess_restart_from_files(int start_no_info,int iguess,double **xInBox);
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Restart_field[NEQ_TYPE];
extern int Imain_loop;
extern int Restart;
#define VERBOSE      3 
extern int Iwrite;
extern int Nnodes_box;
extern int Nunk_per_node;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
void set_initial_guess(int iguess,double **xOwned);
