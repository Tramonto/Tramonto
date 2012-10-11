/* This file was automatically generated.  Do not edit! */
void safe_free(void **ptr);
void safe_free(void **ptr);
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
#define NCOMP_MAX 5
extern int Geqn_start[NCOMP_MAX];
extern int Ngeqn_tot;
#define SCREEN_VERBOSE     3 
extern int Nseg_type[NCOMP_MAX];
extern int Icomp_to_polID[NCOMP_MAX];
#define NMER_MAX     200
extern int Type_mer[NCOMP_MAX][NMER_MAX];
extern int Unk2Comp[NMER_MAX];
#define WJDC2        4 
#define WJDC         3
extern int SegChain2SegAll[NCOMP_MAX][NMER_MAX];
#define WTC          2
extern int SegAll_to_Poly[NMER_MAX];
extern int Nmer_comp[NCOMP_MAX];
extern int Nseg_tot;
extern int Nbonds;
extern int Ncomp;
extern int **Nseg_type_pol;
extern int *BondAll_to_ibond;
extern int *BondAll_to_isegAll;
extern int *Pol_Sym_Seg;
extern int *Pol_Sym;
extern int **Poly_to_Unk_SegAll;
extern int ***Poly_to_Unk;
extern int **Bonds_SegAll;
extern int *Nbonds_SegAll;
extern int *Unk_to_Bond;
extern int *Unk_to_Seg;
extern int *Unk_to_Poly;
void setup_chain_indexing_arrays(int nseg,int nmer_max,FILE *fpecho);
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int Grafted_SegID[NCOMP_MAX];
extern int Grafted[NCOMP_MAX];
#define WJDC3        5 
extern int Type_poly;
#define FILES_DEBUG        2
extern int Iwrite_files;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Nbond_max;
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
void setup_chain_linear_symmetric(FILE *fpecho);
#define LIN_POLY_SYM 2
void setup_chain_linear(FILE *fpecho);
#define LIN_POLY 1
void setup_chain_from_file(FILE *fpecho,char *poly_file);
#define POLY_ARCH_FILE 0
extern int ***pol_sym_tmp;
#define NBOND_MAX 4
extern int ***Bonds;
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
extern int **Nbond;
extern int Nmer[NCOMP_MAX];
extern int Npol_comp;
#define SET_IN_GUI 3
extern int Type_poly_arch;
void setup_chain_architecture(char *poly_file,FILE *fpecho);
