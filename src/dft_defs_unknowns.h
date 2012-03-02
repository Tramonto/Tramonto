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
#define NEQ_TYPE       12 
extern int Constant_row_flag[NEQ_TYPE];
#define SCREEN_VERBOSE     3 
#define NCOMP_MAX 5
extern int Geqn_start[NCOMP_MAX];
#define NO_UNK        -888
extern int Nunk_per_node;
#define NMER_MAX     200
extern int Unk2Phys[3 *NCOMP_MAX+2 *NMER_MAX+NMER_MAX *NMER_MAX+13];
extern int Phys2Unk_last[NEQ_TYPE];
extern int Phys2Unk_first[NEQ_TYPE];
#define SCF_CONSTR	   9
#define SCF_FIELD	  10
extern int Ngeqn_tot;
#define G_CHAIN       11 
#define WJDC_FIELD     8
#define CMS_FIELD      7
extern int Nbonds;
extern int Nrho_bar_bond;
#define BONDWTC        5
#define WJDC3        5 
extern int Nrho_bar_cavity;
#define CAVWTC         4
extern int Npol_comp;
#define DIFFUSIVE_INTERFACE 1
extern int Type_interface;
extern int Ndiffusion;
#define DIFFUSION      6
extern int Type_coul;
extern int Npoisson;
#define POISSON        1
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_attr;
extern int ATTInA22Block;
extern int Nmf_eqns;
#define MF_EQ          3
extern int Nrho_bar_s;
extern int Ndim;
#define IDEAL_GAS    0
extern int Ipot_ff_n;
extern int Nrho_bar;
#define HSRHOBAR       2
extern int Ncomp;
extern int Nseg_tot;
extern int Phys2Nunk[NEQ_TYPE];
#define DENSITY        0
extern int L_HSperturbation;
#define SCFT         6	
#define CMS_SCFT     1
#define CMS          0
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
extern int Lseg_densities;
#define WJDC2        4 
#define WJDC         3
#define WTC          2
extern int Type_poly;
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
#define FILES_DEBUG        2
extern int Iwrite_files;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
void setup_nunk_per_node(char *file_echoinput);
void setup_matrix_constant_blocks();
void setup_matrix_constant_blocks();
