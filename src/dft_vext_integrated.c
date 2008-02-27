/*
//@HEADER
// ******************************************************************** 
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/*
 *  FILE: dft_vext_integrated.c
 *
 *  This file contains routines to integrate a pair potential through the elements of an
 *  arbitrary surface in order to construct an external field.
 */

#include "dft_vext_integrated.h"

/***********************************************************************/
/*integrate_potential: In this routine we use gauss quadrature
                    (same logic as in stencil routine) to calculate
                    the potential energy at a given fluid position
                    due to a given wall element.*/
double integrate_potential(double param1, double param2, double param3, double param4, int ngp, int ngpu,
                      double *gp, double *gpu, double *gw, double *gwu,
                      double *node_pos, double *node_pos_f)
{ 
   int ig,jg,kg,jwall_type=0; 
   double weight, vext=0.0, radius, point[3],cut,z1,z2,eps;

   cut=param3;
   /* note that the parameters that come to this file are as follows*/
   /**** LJ type potential ****/
   /*  param1=sigma           */
   /*  param2=eps             */
   /**** COULOMB paramters ****/
   /*  param1=z1           */
   /*  param2=z2;             */
   /***************************/

   switch(Ndim){
       case 1:

         for (ig=0; ig < ngp; ig++) {
           point[0] = node_pos[0] + gp[ig] * Esize_x[0];
           radius = fabs(node_pos_f[0]-point[0]);

           if (Type_vext3D==PAIR_COULOMB){
              printf("Error - we can't compute coulomb external fields or wall-wall potentials\n");
              printf("        in 1 dimension...would require Ewald summation.\n");
              exit(-1);
           }
           else{
              radius /= cut;
              weight = get_wt_from_sten(radius, param1, param2, param3, param4, ngpu, gpu, gwu);
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
              if(Type_vext3D==PAIR_COULOMB){
                 printf("Error - we can't compute coulomb external fields or wall-wall potentials\n");
                 printf("        in 2 dimensions...would require Ewald summation.\n");
                 exit(-1);
              }
              else{
                 radius /= cut;
                 weight = get_wt_from_sten(radius, param1,param2,param3, param4,ngpu, gpu, gwu);
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
		radius /= cut;

                weight = get_wt_from_sten(radius,param1,param2,param3,param4,ngpu, gpu, gwu); 

                vext += weight * gw[ig] * gw[jg] * gw[kg] * Vol_el;
             }
           }
         }
         break;
   }            /* end of Ndim switch */
   return vext;
}
/***********************************************************************/
/*get_wt_from_sten: here we do the integrations out of the plane*/

double get_wt_from_sten(double r,double param1, double param2, double param3,
				    double param4, int ngpu, double *gpu, double *gwu)
{
  double temp, zmax, z, rho,rcut;
  int i;
  rcut=param3;

  /* if r is larger than the cutoff, return 0
      assumes r = radius/rcut */
  if (r >= 1.0) return(0.0);

  if (Ndim == 1) {
     temp = 0.0;
     zmax = sqrt(1.0 - r*r);
     for (i=0; i < ngpu; i++) {
        z = zmax * gpu[i];
        rho = sqrt(r*r + z*z) * rcut;
        temp += gwu[i] * z * pairPot_switch(rho, param1,param2,param3,param4,Type_vext3D);
     }
     return(2.0 * PI * temp * rcut *rcut * zmax);
  }
  else if (Ndim == 2) {
     temp = 0.0;
     zmax = sqrt(1 - r*r);
     for (i=0; i < ngpu; i++) {
        z = zmax * gpu[i];
        rho = sqrt(r*r + z*z) * rcut;
        temp += gwu[i] * pairPot_switch(rho,param1,param2,param3,param4,Type_vext3D);
     }
     return(2.0 * temp * rcut * zmax);
  }
  else {
    rho = r * rcut;
    temp = pairPot_switch(rho,param1,param2,param3,param4,Type_vext3D);
    return(temp);
  }
}
/***********************************************************************/
