/* This file was automatically generated.  Do not edit! */
int surface_angleCutout3D_cyl(int iwall,int iwall_type,double *fluidEl_center);
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
#define PI    3.141592653589793238462643383279502884197169399375
#define NDIM_MAX  3
#define NWALL_MAX 600 
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern int Ndim;
typedef struct SurfaceGeom_Struct SurfaceGeom_Struct;
extern struct SurfaceGeom_Struct *SGeom;
extern int Lwedge_cutout[NWALL_MAX];
#define NWALL_MAX_TYPE 50 
extern int Lperiodic_overlay[NWALL_MAX_TYPE];
extern int Nperiodic_overlay[NWALL_MAX_TYPE];
#define NPERIODIC_MAX 4
struct SurfaceGeom_Struct {
  int       surfaceTypeID;    /* ID of the type of surface */
  int       orientation;  /* orientation of the surface */
  double    *halfwidth;   /* planar surface params given in halfwidth */
  double    radius;       /* radius of spherical or cylindrical surface */
  double    halflength;   /* length of finite cylinders and pores */
  double    radius2;      /* a second radius for tapered pores or cylinders */
  int       Lwedge_cutout;    /* TRUE or FALSE for wedge cutout from basic surfac */
  double    angle_wedge_start;    /* angle as measured from x0 axis */
  double    angle_wedge_end;    /* angle as measured from x0 axis */
  int       Lrough_surface;    /* TRUE or FALSE for surface roughness */
  double    roughness;          /* maximum roughness amplitude */
  double    roughness_length;    /* lengthscale for the rougness */
  int       Lperiodic_overlay;    /* TRUE or FALSE for periodic function added to surface */
  int       Nperiodic_overlay;     /* The number of periodic functions to apply */
  double    orientation_periodic[NPERIODIC_MAX];    /* maximum amplitude for a cosine wave superimposed on a cylinder */
  double    amplitude[NPERIODIC_MAX];    /* maximum amplitude for a cosine wave superimposed on a cylinder */
  double    wavelength[NPERIODIC_MAX];    /* desired wavelength of cosine wave superimposed on a cylinder */
  double    origin_PeriodicFunc[NPERIODIC_MAX];     /* The origin of periodic functions to apply */
  int    *ReflectionsAreIndependent;  /* TRUE or FALSE for treating special boundary conditions */
};
int surface_angleCutout2D(int iwall,int iwall_type,double *fluidEl_center);
