/* This file was automatically generated.  Do not edit! */
void print_cont_variable_user_plugin(int cont_type,FILE *fp,int Loca_contID);
void print_cont_type_user_plugin(int cont_type,FILE *fp,int Loca_contID);
void adjust_dep_params(int cont_type,int Loca_contID,double param_old,double param_new,char *file_echoinput);
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
#define FILENAME_LENGTH 300
extern char EchoInputFile_array[FILENAME_LENGTH];
void assign_param_user_plugin(int cont_type,int Loca_contID,double param,char *EchoInputFile_array);
double get_init_param_user_plugin(int cont_type,int Loca_contID);
