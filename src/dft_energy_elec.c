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

#include "dft_energy_elec.h"

#define BFD 0
#define FFD 1
#define CFD 2

double calc_deriv_epot(int,int,int*,double **);

/* The following routines were implemented to compute the Tang and Davis 
version of the electrostatic free energy */
/****************************************************************************/
double integrand_elec_PB_freen(int iunk,int inode_box, double **x)
{
     double integrand,rho_i,psi_i;
     int icomp,psiunk;

     if (Lseg_densities) icomp = Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
     else icomp = iunk-Phys2Unk_first[DENSITY];
     
     rho_i = x[iunk][inode_box];
     psiunk = Phys2Unk_first[POISSON];
     psi_i = x[psiunk][inode_box];

     if (rho_i > 0.) integrand = 0.5*rho_i*Charge_f[icomp]*psi_i;
     else integrand=0.0;
     return(integrand);
}
/****************************************************************************/
double integrand_elec_MSAcorr_freen(int iunk,int inode_box, double **x)
{
     double integrand,rho_i;

     rho_i = x[iunk][inode_box];
     integrand=0.0;


     if (Type_coul==DELTAC_RPM){
           if (rho_i > 0.0) integrand = -0.5*rho_i*int_stencil(x,inode_box,iunk,THETA_CR_RPM_MSA);
     }
     else if (Type_coul==DELTAC_GENERAL){
           if (rho_i > 0.0) integrand = -0.5*rho_i*int_stencil(x,inode_box,iunk,THETA_CR_GENERAL_MSA);
     }

     return(integrand);
}
/****************************************************************************/
double integrand_elec_MSAcorr_freen_bulk(int iunk,int inode_box, double **x)
{
     double integrand,rho_bulk;
     int icomp,iseg;

     if (Grafted_Logical==FALSE || (Type_poly==WJDC3 && 
           Grafted[Icomp_to_polID[iunk-Phys2Unk_first[DENSITY]]]==FALSE)){


     if (Lseg_densities){
         iseg = iunk-Phys2Unk_first[DENSITY];
         icomp = Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
         if(Type_interface!=UNIFORM_INTERFACE) rho_bulk = Rho_seg_RTF[iseg];
         else              rho_bulk = Rho_seg_b[iseg];
     }
     else{
           icomp = iunk-Phys2Unk_first[DENSITY];
           if(Type_interface!=UNIFORM_INTERFACE) rho_bulk = Rho_b_RTF[icomp];
           else              rho_bulk = Rho_b[icomp];
     }

     integrand = -0.5*rho_bulk*Deltac_b[icomp];
     }
     else integrand=0.0;
     return(integrand);
}
/****************************************************************************/
/* End of the Tang and Davis calculations for electrostatic free energy */


/* The Following Routines are Implemented to compute electrostatic free energy
   via the Reiner-Radke Route */
/****************************************************************************/
double integrand_maxwell_stress_freen(int iunk,int inode_box, double **x)
{
     double integrand,psi_i,sum_stress,deriv;
     int psiunk,ijk[3],int_type[3],idim;

     psiunk = Phys2Unk_first[POISSON];
     psi_i = x[psiunk][inode_box];

     node_to_ijk(B2G_node[inode_box],ijk);
     sum_stress=0.0;

     for (idim=0;idim<Ndim;idim++){

        /* identify type of finite difference to use */
        int_type[idim]=CFD; 
        if (ijk[idim] == Nodes_x[idim]-1) int_type[idim]=BFD;
        else if (ijk[idim] == 0) int_type[idim]=FFD;

        deriv=calc_deriv_epot(idim,inode_box,int_type,x);
        sum_stress += (deriv*deriv);
     }
     
     integrand=-(Temp_elec/(8.0*PI))*sum_stress;
     return(integrand);
}
/****************************************************************************
 * integrand_surface_charge:  In this routine we calculate the surface       *
 *                           integral that gives the free energy due to     *
 *                           charging up the surfaces (see the first        *
 *                           term in Eq. 6 of Reiner and Radke 1991         *
 *                           Note that we compute surface charge based on   *
 *                           the gradient of the electrostatic field rather *
 *                           than using the surface charge arrays so that   *
 *                           this code will work for constant surface charge*
 *                           or constant surface potential boundaries       */
double integrand_surface_charge(int iunk,int inode_box,int iwall,double **x)
{
  int loc_inode,idim,iel_w,surf_norm,ilist,int_type[3];
  double integrand,charge_i,prefac,deriv;

  integrand=0.0;
  iunk = Phys2Unk_first[POISSON];
  loc_inode=B2L_node[inode_box];
  int_type[0]=int_type[1]=int_type[2]=CFD;

  ilist = Nlists_HW - 1;
  if (Type_bc_elec[WallType[iwall]] != CONST_POTENTIAL){
  for (iel_w=0; iel_w<Nelems_S[ilist][loc_inode]; iel_w++){

              surf_norm = Surf_normal[ilist][loc_inode][iel_w];
              idim = abs(surf_norm) - 1;

              prefac = (double)(surf_norm/abs(surf_norm))*Area_surf_el[idim]/Nelems_S[ilist][loc_inode];

              if (surf_norm <0) int_type[idim]=BFD;
              else              int_type[idim]=FFD;
              deriv=calc_deriv_epot(idim,inode_box,int_type,x);

              charge_i = prefac*(Temp_elec/(4.0*PI))*deriv;
              integrand -= (charge_i*x[iunk][inode_box]);

  } /* end of surface element loop */

/*         charge_i = 0.0;
         for (idim=0; idim<Ndim; idim++){
             charge_i -= Charge_w_sum_els[loc_inode][idim]*Area_surf_el[idim];
         }
 */
  }
  else integrand=0.0;

  return integrand;
  
}
/****************************************************************************/
/* calc_deriv_epot : calculate a derivative of the electric potential!!   */

double calc_deriv_epot(int idim,int inode0,int *int_type, double **x)
{
   int inode1,inode2,offset1[3],offset2[3],
       ijk_box[3],reflect_flag[3],iunk,kdim;
   double deriv=0.0;

   node_box_to_ijk_box(inode0,ijk_box);

   for (kdim=0;kdim<Ndim;kdim++){
       offset1[kdim]=0; offset2[kdim]=0;
   }

   switch(int_type[idim])
   {
      case CFD:
        offset1[idim] = -1; offset2[idim] =  1; break;
      case FFD:
        offset1[idim] =  1; offset2[idim] =  2; break;
      case BFD:
        offset1[idim] = -1; offset2[idim] = -2; break;
   }

   inode1 = offset_to_node_box(ijk_box,offset1,reflect_flag);
   inode2 = offset_to_node_box(ijk_box,offset2,reflect_flag);
   iunk = Phys2Unk_first[POISSON];


   switch(int_type[idim]){
      case CFD: deriv =    x[iunk][inode2] - x[iunk][inode1];break;
      case FFD: deriv = -3*x[iunk][inode0] + 4*x[iunk][inode1] - x[iunk][inode2];break;
      case BFD: deriv =  3*x[iunk][inode0] - 4*x[iunk][inode1] + x[iunk][inode2];break;
   }
   deriv /= (2*Esize_x[idim]);

   return deriv;
}
/****************************************************************************/

