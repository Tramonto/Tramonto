/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

/*
 *  FILE: dft_vext.c
 *
 *  This file contains routines to set up an external field for a 
 *  variety of cases.
 */

#include <stdio.h>
#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

/* Prototypes for functions found in this file */
double get_wt_from_sten(double,double, double, double,
				    int, double *, double *);
double get_wt_from_sten_coul3D(double,double, double);

#define FLAG_LJ 0
#define FLAG_COULOMB 1
 

/***********************************************************************/
/*integrate_potential: In this routine we use gauss quadrature
                    (same logic as in stencil routine) to calculate
                    the potential energy at a give fluid position
                    due to a given wall element.*/
double integrate_potential(int flag, double param1, double param2, double param3, int ngp, int ngpu,
                      double *gp, double *gpu, double *gw, double *gwu,
                      double *node_pos, double *node_pos_f)
{ 
   int ig,jg,kg,jwall_type=0; 
   double weight, vext=0.0, radius, point[3],cut,z1,z2,sigma,eps;

   if (flag == FLAG_LJ){
     sigma=param1;
     eps=param2;
     cut=param3;
   }
   else if (flag == FLAG_COULOMB){
     z1=param1;
     z2=param2;
     cut=param3;
   }

   switch(Ndim){
       case 1:

         for (ig=0; ig < ngp; ig++) {
           point[0] = node_pos[0] + gp[ig] * Esize_x[0];
           radius = fabs(node_pos_f[0]-point[0]);

           if (flag==FLAG_LJ){
              radius /= cut;
              weight = get_wt_from_sten(radius, sigma, eps, cut, ngpu, gpu, gwu);
           }
           else{
              printf("Error - we can't compute coulomb external fields or wall-wall potentials\n");
              printf("        in 1 dimension...would require Ewald summation.\n");
              exit(-1);
           }
           
           vext += weight * gw[ig] * Vol_el;

         }
         break;

       case 2:
         for (ig=0; ig < ngp; ig++) {
           point[0] = node_pos[0] + gp[ig] * Esize_x[0];

           for (jg=0; jg < ngp; jg++) {
             point[1] = node_pos[1] + gp[jg] * Esize_x[1];

             radius = sqrt( (node_pos_f[0] - point[0])*
                            (node_pos_f[0] - point[0]) + 
                            (node_pos_f[1] - point[1])*
                            (node_pos_f[1] - point[1]) );
              if (flag==FLAG_LJ){
                    radius /= cut;
                    weight = get_wt_from_sten(radius, sigma,eps,cut, ngpu, gpu, gwu);
              }
              else {
                 printf("Error - we can't compute coulomb external fields or wall-wall potentials\n");
                 printf("        in 2 dimensions...would require Ewald summation.\n");
                 exit(-1);
              }

             vext += weight * gw[ig] * gw[jg] *Vol_el;
           }
         }
         break;


       case 3:
         for (ig=0; ig < ngp; ig++) {
           point[0] = node_pos[0] + gp[ig] * Esize_x[0];

           for (jg=0; jg < ngp; jg++) {
             point[1] = node_pos[1] + gp[jg] * Esize_x[1];

             for (kg=0; kg < ngp; kg++){
                point[2] = node_pos[2] + gp[kg] * Esize_x[2];
                radius = sqrt( (node_pos_f[0] - point[0])*
                            (node_pos_f[0] - point[0]) + 
                            (node_pos_f[1] - point[1])*
                            (node_pos_f[1] - point[1]) +
                            (node_pos_f[2] - point[2])*
                            (node_pos_f[2] - point[2]) );

                 if (flag==FLAG_LJ){
                            radius /= cut;
                            weight = get_wt_from_sten(radius, sigma,eps,cut,ngpu, gpu, gwu); 
                 }
                 else   weight = get_wt_from_sten_coul3D(radius, z1,z2);

                vext += weight * gw[ig] * gw[jg] * gw[kg] /* * Vol_el*/;
             }
           }
         }
         break;
   }            /* end of Ndim switch */
   return vext;
}
/***********************************************************************/
/*get_wt_from_sten: here we do the integrations out of the plane*/

double get_wt_from_sten(double r,double sigma, double eps, double rcut,
				    int ngpu, double *gpu, double *gwu)
{
  double temp, zmax, z, rho;
  int i;

  if (r >= 1.0) return(0.0);
  if (Ndim == 1) {
     temp = 0.0;
     zmax = sqrt(1.0 - r*r);
     for (i=0; i < ngpu; i++) {
        z = zmax * gpu[i];
        rho = sqrt(r*r + z*z) * rcut;
        temp += gwu[i] * z * uLJ12_6_cut(rho, sigma,eps,rcut);
     }
     return(2.0 * PI * temp * rcut *rcut * zmax);
  }
  else if (Ndim == 2) {
     temp = 0.0;
     zmax = sqrt(1 - r*r);
     for (i=0; i < ngpu; i++) {
        z = zmax * gpu[i];
        rho = sqrt(r*r + z*z) * rcut;
        temp += gwu[i] * uLJ12_6_cut(rho,sigma,eps,rcut);
     }
     return(2.0 * temp * rcut * zmax);
  }
  else {
    rho = r * rcut;
    temp = uLJ12_6_cut(rho,sigma,eps,rcut);
    return(temp);
  }

}
/***********************************************************************/
/*get_wt_from_sten_coul: here we do the integrations out of the plane*/

double get_wt_from_sten_coul3D(double r,double z1, double z2)
{
  return (uCOULOMB(r,z1,z2));
}
/******************************************************************************/
/* uLJ12_6_cut: The external field contributions in 3-dimensions                   */

double uLJ12_6_cut(double r,double sigma, double eps, double rcut)
{
  double u;
  
  if (r <= rcut) {
     u = ( POW_DOUBLE_INT(sigma/r,12) - POW_DOUBLE_INT(sigma/rcut,12) ) -
            ( POW_DOUBLE_INT(sigma/r,6) - POW_DOUBLE_INT(sigma/rcut,6) );
  }
  else u = 0.0;
  return (4.0*eps*u);
}
/*******************************************************************************/
/* uderiv_LJ12_6: The external field contributions in 3-dimensions              */

double uderiv_LJ12_6(double r,double x,double sigma, double eps, double rcut)
{
  double vext_dash;
  
  if (r <= rcut) {
     vext_dash = (1.0/sigma) * (
            -12.*x*POW_DOUBLE_INT(sigma/r,14) + 6.*x*POW_DOUBLE_INT(sigma/r,8) );
  }
  else vext_dash = 0.0;
  return (4.0*eps*vext_dash);
}
/******************************************************************************/
/* Vext_1D:  given a wall-fluid point separation, calculate the 
           9-3 Lennard Jones potential (note that prefactors are
           calculated in calling routine)                           */
double Vext_1D(double x,int icomp, int iwall_type)
{
  double vext,xmin;

  if (x <= Cut_wf[icomp][iwall_type]) {


     /* cut and shifted 9-3 LJ potential */
     vext = (2.0*PI/3.0)*(POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type],3))* (
       (2.0/15.0) * ( POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/x,9) -
       POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/Cut_wf[icomp][iwall_type],9) ) -
       ( POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/x,3) -
       POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/Cut_wf[icomp][iwall_type],3) ));
    
      /* repulsive 9 potential (no cut and shift) */
      /*vext = (2.0*PI/3.0)*( (2.0/15.0)*POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type],12)/ POW_DOUBLE_INT(x,9));*/

      /* 9-3 potential (no cut and shift) */
     /*    vext = (2.0*PI/3.0)*( (2.0/15.0)*POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type],12)/ POW_DOUBLE_INT(x,9)
                    -POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type],6)/POW_DOUBLE_INT(x,3));*/

       /* exponential attractive potential */
       /*vext = -exp(-x/Sigma_wf[icomp][iwall_type]);*/

  }
  else vext = 0.0;
  return vext;
}

/******************************************************************************/
/* Vext_1D_dash:  given a wall-fluid point separation, calculate the derivative of
               the 9-3 Lennard Jones potential (note that prefactors are
               calculated in calling routine)                           */
/* ALF: modified to have min. start at wall surface */
double Vext_1D_dash(double x,int icomp, int iwall_type)
{
  double vdash,xmin;
  xmin = pow(2./5.,1./6.);
  if (x <= Cut_wf[icomp][iwall_type]) {
/*    x += xmin;*/

     /* LJ 9-3 potential (with or without cut and shift) */
     vdash = (2.0*PI/3.0)*(POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type],2))* (
             -(9.0/5.0) * POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/x,10)
             +(9.0/2.0) * POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/x,4) );

      /* LJ 9 repulsive potential */
/*      vdash = (1.0/Sigma_wf[icomp][iwall_type])* (
    -(9.0/5.0) * POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/x,10) );*/

     /* exponential attractive potential */
     /*vext_dash =  exp(-x/Sigma_wf[icomp][iwall_type])/
                          (x*Sigma_wf[icomp][iwall_type]);*/

  }
  else vdash = 0.0;
  return vdash;
}
/******************************************************************************/
/* uLJ_wp: The external field contributions in 3-dimensions                   */

double uLJ_wp(double r,int icomp, int iwall_type)
{
  double vext;
 
  if (r <= Cut_wf[icomp][iwall_type]) {
  /* vext = -1./r + 1./Cut_wf[icomp][iwall_type];*/
  /* vext = ( POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/r,3) -
            POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/Cut_wf[icomp][iwall_type],3) ) -
            ( POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/r,6) -
            POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/Cut_wf[icomp][iwall_type],6) );  */
     vext = ( POW_DOUBLE_INT(r/Cut_wf[icomp][iwall_type],2) - 1.);
  }
  else vext = 0.0;
  return (4.0*Eps_wf[icomp][iwall_type]*vext);
}
/******************************************************************************/
/* uCOULOMB: The wall-wall LJ contributions in 3-dimensions                   */

double uCOULOMB(double r,double z1,double z2)
{
  double u;

  u = z1*z2/r;
  /*return (4.0*PI*u/Temp_elec);*/
  return (u/Temp_elec);
}
/****************************************************************************/
/* uLJatt_n:  given r12, calculate the attractive part of a cut and
           shifted 12-6 LJ potential. */

double uLJatt_n(double r,int i, int j)
{
  double uatt,r_min,sigma2,sigma6;
  double r_inv,r2_inv,r6_inv,r12_inv;
  double rc_inv,rc2_inv,rc6_inv,rc12_inv;

  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;

  rc_inv   = 1.0/Cut_ff[i][j];
  rc2_inv  = rc_inv*rc_inv;
  rc6_inv  = rc2_inv*rc2_inv*rc2_inv;
  rc12_inv = rc6_inv*rc6_inv;


  if (r <= Cut_ff[i][j]) {
     r_min = Sigma_ff[i][j] * pow(2.0,1.0/6.0);
     if (r < r_min) r = r_min;

     r_inv = 1.0/r;
     r2_inv  = r_inv*r_inv;
     r6_inv  = r2_inv*r2_inv*r2_inv;
     r12_inv = r6_inv*r6_inv;

     uatt = 4.0 * Eps_ff[i][j]* sigma6 * (
               sigma6*(r12_inv - rc12_inv)
                    - (r6_inv  - rc6_inv ) );
  }
  else uatt = 0.0;

  return uatt;
}
/****************************************************************************/
/* uLJatt_n_int:  given r_12, calculate the integral of the
                  attractive part of
                  12-6 LJ potential (not cut and shifted). */

double uLJatt_n_int(double r,int i, int j)
{
  double uatt_int, sigma2,sigma6;
  double r_inv,r3_inv,r9_inv;

  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;

  r_inv = 1.0/r;

  r3_inv  = r_inv*r_inv*r_inv;
  r9_inv  = r3_inv*r3_inv*r3_inv;

  uatt_int = 16 * PI * Eps_ff[i][j]* sigma6 * (
             - sigma6*r9_inv/9.0  + r3_inv/3.0 );

  return uatt_int;
}
/****************************************************************************/
/* uLJatt_n_noshift:  calculate the attractive part of
                  12-6 LJ potential at the minimum. */

double uLJatt_n_noshift(double r,int i, int j)
{
  double uatt,r_min,sigma2,sigma6;
  double r_inv,r2_inv,r6_inv,r12_inv;

  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;

  r_inv = 1.0/r;

  r2_inv  = r_inv*r_inv;
  r6_inv  = r2_inv*r2_inv*r2_inv;
  r12_inv = r6_inv*r6_inv;

  r_min = Sigma_ff[i][j] * pow(2.0,1.0/6.0);
  if (r < r_min) r = r_min;

  uatt = 4.0 * Eps_ff[i][j]* sigma6 * (
            sigma6*r12_inv  - r6_inv);

  return uatt;
}
/****************************************************************************/
/* uCOULOMB_att:  given r12, calculate the "attractive" part of a cut and
           shifted Coulomb potential. */

double uCOULOMB_att(double r,int i, int j)
{
  double uatt,r_min;

  if (r <= Cut_ff[i][j]) {
     r_min = Sigma_ff[i][j];
     if (r < r_min) r = r_min;

    /* uatt = 4.*PI*Charge_f[i]*Charge_f[j]/(r*Temp_elec);*/
     uatt = Charge_f[i]*Charge_f[j]/(r*Temp_elec);
  }
  else uatt = 0.0;

  return uatt;
}
/****************************************************************************/
/* uCOULOMB_att_int:  given r_12, calculate the integral of the
                  "attractive" part of COULOMB potential (no cut-shift). */

double uCOULOMB_att_int(double r,int i, int j)
{
  double uCOUL_int;

  /*uCOUL_int = 16.*PI*PI*Charge_f[i]*Charge_f[j]*r*r/Temp_elec;*/
  uCOUL_int = 4.*PI*Charge_f[i]*Charge_f[j]*r*r/Temp_elec;

  return uCOUL_int;
}


