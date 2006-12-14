/*
//@HEADER
// ******************************************************************** 
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/*
 *  FILE: dft_potentials_vext.c
 *
 *  This file contains external field implementations.  New functions added
 *  here should be given suitable keyword names and added to the define
 *  list in dft_globals_const.h as well as in the selector cases in
 *  dft_vext1D.c
 *
 */

#include "dft_potentials_vext.h"

/********************************************************************************************/
/* Vext_LJ9_3_CS: cut and shifted 9-3 LJ potential */
double Vext_LJ9_3_CS(double x,int icomp, int iwall_type)
{
   double vext,prefac;
   prefac = Eps_wf[icomp][iwall_type]*Rho_w[iwall_type];


   if (x <= Cut_wf[icomp][iwall_type]) 
        vext = prefac*(2.0*PI/3.0)*(POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type],3))* (
               (2.0/15.0) * ( POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/x,9) -
               POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/Cut_wf[icomp][iwall_type],9) ) -
               ( POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/x,3) -
               POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/Cut_wf[icomp][iwall_type],3) ));
   else vext=0.0;

   return(vext);
}
/********************************************************************************************/
/* Vextderiv_LJ9_3: derivative of LJ 9-3 potential (with or without cut and shift) */
double Vextderiv_LJ9_3(double x,int icomp, int iwall_type)
{
   double vdash,prefac;
   prefac = Eps_wf[icomp][iwall_type]*Rho_w[iwall_type];

  if (x <= Cut_wf[icomp][iwall_type])
        vdash = prefac*(2.0*PI/3.0)*(POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type],2))* (
                -(9.0/5.0) * POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/x,10)
                +(9.0/2.0) * POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/x,4) );
   else vdash=0.0;

   return(vdash);
}
/********************************************************************************************/
/* Vext_LJ9_3_v2_CS: cut and shifted 9-3 LJ potential */
double Vext_LJ9_3_v2_CS(double x,int icomp, int iwall_type)
{
   double vext,prefac;
   prefac = Eps_wf[icomp][iwall_type]*Rho_w[iwall_type];

   if (x <= Cut_wf[icomp][iwall_type]) 
        vext =  prefac* ((2.0/15.0) * ( POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/x,9) -
           POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/Cut_wf[icomp][iwall_type],9) ) -
           ( POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/x,3) -
           POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/Cut_wf[icomp][iwall_type],3) ));
   else vext=0.0;

   return(vext);
}
/********************************************************************************************/
/* Vextderiv_LJ9_3_v2: derivative of LJ 9-3 potential (with or without cut and shift) */
double Vextderiv_LJ9_3_v2(double x,int icomp, int iwall_type)
{
   double vdash,prefac;
   prefac = Eps_wf[icomp][iwall_type]*Rho_w[iwall_type];

  if (x <= Cut_wf[icomp][iwall_type])
        vdash = prefac*(POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type],-1))* (
                -(9.0/5.0) * POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/x,10)
                +(9.0/2.0) * POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/x,4) );
   else vdash=0.0;

   return(vdash);
}
/********************************************************************************************/
/* Vext_LJ9_3_noCS: 9-3 potential (no cut and shift) */
double Vext_LJ9_3_noCS(double x,int icomp, int iwall_type)
{
  double vext,prefac;
  prefac = Eps_wf[icomp][iwall_type]*Rho_w[iwall_type];

  if (x <= Cut_wf[icomp][iwall_type]) 
       vext = prefac*(2.0*PI/3.0)*( (2.0/15.0)*POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type],12)/ POW_DOUBLE_INT(x,9)
                           -POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type],6)/POW_DOUBLE_INT(x,3));
  else vext=0.0;

  return(vext);
}
/********************************************************************************************/
/* Vext_LJ9_3_shiftX_CS: cut and shifted 9-3 LJ potential */
double Vext_LJ9_3_shiftX_CS(double x,int icomp, int iwall_type)
{
   double vext, delta, vext1d_xmax, vext1d_xmin,prefac;
   prefac = Eps_wf[icomp][iwall_type]*Rho_w[iwall_type];

   delta       = 0.5*(Sigma_ff[icomp][icomp]-1);
   vext1d_xmax = Cut_wf[icomp][iwall_type]+delta;
   vext1d_xmin = vext1d_xmax-0.5;

   if       (x <= vext1d_xmin)  vext=VEXT_MAX;
   else if  (x <= vext1d_xmax) 
        vext = prefac*( (2.0/15.0) * ( POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/(x-delta),9) -
          POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/Cut_wf[icomp][iwall_type],9) ) -
          ( POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/(x-delta),3) -
          POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/Cut_wf[icomp][iwall_type],3) ));
   else vext=0.0;

   return(vext);
}
/********************************************************************************************/
/* Vext_REPULSIVE9_noCS: LJ 9 repulsive potential */
double Vext_REPULSIVE9_noCS(double x,int icomp, int iwall_type)
{
  double vext,prefac;
  prefac = Eps_wf[icomp][iwall_type]*Rho_w[iwall_type];

  if (x <= Cut_wf[icomp][iwall_type]) 
      vext = prefac*(2.0*PI/3.0)*( (2.0/15.0)*POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type],12)/ POW_DOUBLE_INT(x,9));
  else vext=0.0;

  return(vext);
}
/********************************************************************************************/
/* Vextderiv_REPULSIVE9: derivative of LJ 9 repulsive potential */
double Vextderiv_REPULSIVE9(double x,int icomp, int iwall_type)
{
  double vdash,prefac;
  prefac = Eps_wf[icomp][iwall_type]*Rho_w[iwall_type];

  if (x <= Cut_wf[icomp][iwall_type]) 
    vdash = prefac*(1.0/Sigma_wf[icomp][iwall_type])* (
           -(9.0/5.0) * POW_DOUBLE_INT(Sigma_wf[icomp][iwall_type]/x,10) );
  else vdash=0.0;

  return(vdash);
}
/********************************************************************************************/
/* Vext_EXP_ATT_noCS: exponential attractive potential */
double Vext_EXP_ATT_noCS(double x,int icomp, int iwall_type)
{
  double vext,prefac;
  prefac = Eps_wf[icomp][iwall_type]*Rho_w[iwall_type];

  if (x <= Cut_wf[icomp][iwall_type]) vext = -prefac*exp(-x/Sigma_wf[icomp][iwall_type]);
  else vext=0.0;

  return(vext);
}
/********************************************************************************************/
/* Vextderiv_EXP_ATT: derivative of exponential attractive potential */
double Vextderiv_EXP_ATT(double x,int icomp, int iwall_type)
{
  double vdash,prefac;
  prefac = Eps_wf[icomp][iwall_type]*Rho_w[iwall_type];

  if (x <= Cut_wf[icomp][iwall_type]) 
       vdash =  prefac*exp(-x/Sigma_wf[icomp][iwall_type])/
                     (x*Sigma_wf[icomp][iwall_type]);
  else vdash=0.0;

  return(vdash);
}
/********************************************************************************************/
