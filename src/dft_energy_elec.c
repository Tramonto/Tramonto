/*
//@HEADER
// ******************************************************************** 
// Copyright (2006) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000, there is a non-exclusive license for use of this
// work by or on behalf of the U.S. Government. Export of this program
// may require a license from the United States Government.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// ********************************************************************
//@HEADER
*/

#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

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

     if (Type_poly==WTC) icomp = Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
     else icomp = iunk-Phys2Unk_first[DENSITY];
     
     rho_i = x[iunk][inode_box];
     psiunk = Phys2Unk_first[POISSON];
     psi_i = x[psiunk][inode_box];

     
     if (rho_i > 0.) integrand = 0.5*rho_i*Charge_f[icomp]*psi_i;
     return(integrand);
}
/****************************************************************************/
double integrand_elec_MSAcorr_freen(int iunk,int inode_box, double **x)
{
     double integrand,rho_i;

     rho_i = x[iunk][inode_box];

     int_stencil(x,inode_box,iunk,THETA_CHARGE);
     if (rho_i > 0.0) integrand = -0.5*rho_i*Temporary_sum;
     return(integrand);
}
/****************************************************************************/
double integrand_elec_MSAcorr_freen_bulk(int iunk,int inode_box, double **x)
{
     double integrand,rho_bulk;
     int icomp,iseg;

     if (Type_poly==WTC){
         iseg = iunk-Phys2Unk_first[DENSITY];
         icomp = Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
         if(Lsteady_state) rho_bulk = Rho_seg_RTF[iseg];
         else              rho_bulk = Rho_seg_b[iseg];
     }
     else{
           icomp = iunk-Phys2Unk_first[DENSITY];
           if(Lsteady_state) rho_bulk = Rho_b_RTF[icomp];
           else              rho_bulk = Rho_b[icomp];
     }

     integrand = -0.5*rho_bulk*Deltac_b[icomp];
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
        if (ijk[0] == Nodes_x[0]-1) int_type[idim]=BFD;
        else if (ijk[0] == 0) int_type[idim]=FFD;

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
double integrand_surface_charge(int iunk,int inode_box,double **x)
{
  int loc_inode,idim,iel_w,surf_norm,ilist,int_type[3];
  double integrand,charge_i,prefac,deriv;

  iunk = Phys2Unk_first[POISSON];
  loc_inode=B2L_node[inode_box];
  int_type[0]=int_type[1]=int_type[2]=CFD;

  ilist = Nlists_HW - 1;
  for (iel_w=0; iel_w<Nelems_S[ilist][loc_inode]; iel_w++){

              surf_norm = Surf_normal[ilist][loc_inode][iel_w];
              idim = abs(surf_norm) - 1;

              prefac = (double)(surf_norm/abs(surf_norm))*Area_surf_el[idim]/Nelems_S[ilist][loc_inode];

              if (surf_norm <0) int_type[idim]=BFD;
              else              int_type[idim]=FFD;
              deriv=calc_deriv_epot(idim,inode_box,int_type,x);

              charge_i = prefac*(Temp_elec/(8.0*PI))*deriv;

 printf("iel_w=%d Nelems_S=%d idim=%d  deriv=%9.6f  prefac=%9.6f other=%9.6f \n",iel_w,Nelems_S[ilist][loc_inode],idim,deriv,prefac,Temp_elec/(8.*PI));
          
  } /* end of surface element loop */

/*         charge_i = 0.0;
         for (idim=0; idim<Ndim; idim++){
             charge_i -= Charge_w_sum_els[loc_inode][idim]*Area_surf_el[idim];
         }
 */
  integrand = -(charge_i*x[iunk][inode_box]);

  return integrand;
  
}
/****************************************************************************/
/* calc_deriv_epot : calculate a derivative of the electric potential!!   */

double calc_deriv_epot(int idim,int inode0,int *int_type, double **x)
{
   int inode1,inode2,offset1[3],offset2[3],
       ijk_box[3],reflect_flag[3],iunk;
   double deriv=0.0;

   node_box_to_ijk_box(inode0,ijk_box);

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

