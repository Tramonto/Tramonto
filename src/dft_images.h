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
#define NWALL_MAX 600 
extern int WallType[NWALL_MAX];
#define NWALL_MAX_TYPE 50 
extern double WallParam[NWALL_MAX_TYPE];
#define NDIM_MAX  3
extern double WallPos[NDIM_MAX][NWALL_MAX];
void find_images2(int idim,double cut,int *image,double **image_pos,double *node_image,int iwall,int iside);
void find_images_coulomb(int idim,int *image,double **image_pos,double *node_image);
void find_images_1D(int idim,double cut,int *image,double **image_pos,double *node_image,double *node_ref);
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern double Esize_x[NDIM_MAX];
extern double Size_x[NDIM_MAX];
#define IN_WALL             -1
#define REFLECT              2
#define PERIODIC             1
extern int Type_bc[NDIM_MAX][2];
extern int Ndim;
void find_images(int idim,double cut,int *image,double **image_pos,double *node_image,double *node_ref);
