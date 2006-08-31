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

/*
 *  FILE: dft_energy_archive.c
 *
 *  This file contains routines that are not currently used by the code, but
 *  were previously implemented in order to attempt to compute free energies
 *  of charged systems.  These routines have not been eliminated yet because
 *  this issue is not completely resolved as of 10/2005. LJDF
 *
 */

#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

#define BFD 0
#define FFD 1
#define CFD 2

double free_energy_charging_up(double **);
double energy_elec(double **,double *);
double charge_stress(double **,double *);
double calc_deriv_e(int,int,int,int *,double **,int);
double calc_deriv2(int,int,int,double **);
double energy_elec_vext_vol(double **);

/*******************************************************************************
 energy_elec_vext_vol(x); compute the wall-fluid electrostatic  contributions
                            to the free energy. This is 1/2 times the integral
                            over the fixed charge distribution times the 
                            electrostatic potential. */

double energy_elec_vext_vol(double **x)
{
  int loc_inode, iunk, inode_box, inode, iel_box,ijk_box[3],jln,ielement,reflect_flag[3];
  double charge_at_node,energy;

  reflect_flag[0]=reflect_flag[1]=reflect_flag[2]=0;

 /* The only tricky issue here is that the volumetric charge distribution
    is assigned on an element by element basis while the electrostatic
    potentials are defined on the nodes of the calculation.  So, we
    compute the average charge at each node for this
    integral. */

    energy = 0.0;

    for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
       inode = L2G_node[loc_inode];
       inode_box = L2B_node[loc_inode];
   
      iunk = Phys2Unk_first[POISSON];
    
       charge_at_node = 0.0; 
       for (jln=0; jln< Nnodes_per_el_V; jln++) { 
          ielement = node_to_elem(inode,jln,reflect_flag);
          iel_box = el_to_el_box(ielement); 
          if (iel_box > 0){
             charge_at_node += Charge_vol_els[iel_box]/Nnodes_per_el_V;
          }
       }
      energy += 0.5*charge_at_node*x[iunk][inode_box]*Vol_el;

    }
    printf("PROC=%d: RETURNING energy=%9.6f \n",Proc,energy);
    return (energy*Vol_el);
}

/*******************************************************************************
 * energy_elec: For charged surfaces, find the electrostatic contribution   *
 *             to the force on each wall.                                  */

double energy_elec(double **x, double *sum3)
{
   int iunk,loc_inode,iwall,idim=0,ilist,
     iel_w,inode,surf_norm,ijk[3],inode_box;
   int wall_count[NWALL_MAX],blocked,nblock=0;
   double prefac,nodepos[3],geom_factor[3],deriv_x[3];
   double sum=0.0,wall_avg_psi[NWALL_MAX],
          sblock=0.0,sopen=0.0;
   double deriv_x_A[3],deriv_x_B[3];
   *sum3 = 0.;

   for (iwall=0; iwall<Nwall; iwall++){
        wall_avg_psi[iwall] = 0.;
        wall_count[iwall] = 0;
   }

   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){

       inode_box = L2B_node[loc_inode];
       inode = L2G_node[loc_inode];
       node_to_ijk(inode,ijk);

       for (ilist=0; ilist<Nlists_HW;ilist++) {
          iwall = Nodes_2_boundary_wall[ilist][inode_box];
          if (iwall !=-1 ) break;
       }
       /* iwall = Nodes_2_boundary_wall[Nlists_HW-1][inode_box];
          iwall = Nodes_2_boundary_wall[0][inode_box];*/

       if (iwall != -1){       /* sum up the boundary element contributions */

          node_to_position(inode,nodepos); 
/*          calc_geom_factor(iwall,nodepos,geom_factor);*/

          geom_factor[0]=geom_factor[1]=geom_factor[2]=1.0;

	  iunk = Phys2Unk_first[POISSON];

          for (idim=0; idim<Ndim; idim++)  deriv_x[idim] = 0.0;

            for (iel_w=0; iel_w<Nelems_S[ilist][loc_inode]; iel_w++){
              surf_norm = Surf_normal[ilist][loc_inode][iel_w];
              idim = abs(surf_norm) - 1;

              prefac = (double)(surf_norm/abs(surf_norm))*Area_surf_el[idim]/
                                                           (Nnodes_per_el_S);
 
              if (surf_norm < 0) {
                 deriv_x_A[idim] = calc_deriv_e(idim,inode_box,BFD,&blocked,x,ilist);
                 deriv_x_B[idim] = calc_deriv_e(idim,inode_box,FFD,&blocked,x,ilist);
              }
              else{               
                  deriv_x_A[idim] = calc_deriv_e(idim,inode_box,FFD,&blocked,x,ilist);
                  deriv_x_B[idim] = calc_deriv_e(idim,inode_box,BFD,&blocked,x,ilist);
              }

/*              printf("iel_w: %d  surf_norm: %d  idim: %d\n",iel_w,surf_norm,idim);*/

              sum += x[iunk][inode_box]*(deriv_x_A[idim]-deriv_x_B[idim])*prefac*geom_factor[idim];

              if (blocked){
                 sblock += deriv_x[idim]*prefac;
                 nblock++;
              }
              else
                 sopen += deriv_x[idim]*prefac;

              if (Type_bc_elec[WallType[iwall]]==2){
                       wall_avg_psi[iwall] += x[iunk][inode_box];
                       wall_count[iwall]++;
              }

          } /* end of surface element loop */ 
       }    /* end of test for boundary node */

       /* now calculate the volume term !! */

  /*   for (loc_node_el=0; loc_node_el< Nnodes_per_el_V ; loc_node_el++){

          iel = node_to_elem(inode, loc_node_el, reflect_flag);
          iel_box = el_to_el_box(iel);

*          if (iel >= 0 && Wall_elems[Nlists_HW-1][iel_box] != -1){*
            if (iel >= 0 && Wall_elemsF[0][iel_box] != -1){
*            if (iel >= 0 ){*

              for (idim=0; idim<Ndim; idim++){
                 deriv2_x[idim] = 0.0;
                 switch(Ndim)
                 {
                  case 3:
                    if ( 
                         (idim == 0 && (loc_node_el==0 || loc_node_el == 2
                                   || loc_node_el==4 || loc_node_el==6 )  
                                    && ijk[0] < Nodes_x[0]-1 )              ||
                         (idim == 1 && (loc_node_el==0 || loc_node_el == 1
                                   || loc_node_el==4 || loc_node_el==5 ) 
                                    && ijk[1] < Nodes_x[1]-1 )              ||
                         (idim == 2 && (loc_node_el==0 || loc_node_el == 1
                                   || loc_node_el==2 || loc_node_el==3 ) 
                                    && ijk[2] < Nodes_x[2]-1 )               
                         ){
                         deriv2_x[idim] = calc_deriv2(idim,inode_box,FFD,x);
                         sum2 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                    }
                    break;
                  case 2:
                    if ( 
                         (idim == 0 && (loc_node_el==0 || loc_node_el == 2 ) 
                                    && ijk[0] < Nodes_x[0]-1               ) ||
                         (idim == 1 && (loc_node_el==0 || loc_node_el == 1 ) 
                                    && ijk[1] < Nodes_x[1]-1               )
                         ){
                         deriv2_x[idim] = calc_deriv2(idim,inode_box,FFD,x);
                         sum2 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                    }
                    break;
                  case 1:
                    if ( loc_node_el==0 && ijk[0] < Nodes_x[0]-1 ){
                         deriv2_x[idim] = calc_deriv2(idim,inode_box,FFD,x);
                         sum2 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                    }
                    break;
                 }
              }   * end of idim loop * 
          }       * end of test for a wall element * 
       }          * end of loop over elements that include inode */


    }       /* end of loop over nodes on this processor */
    sum *= -Temp_elec/(8.0*PI);

/*     calculate charging free energy without derivative */
    for (iwall=0; iwall<Nwall; iwall++){
       if (Type_bc_elec[WallType[iwall]]==2) *sum3 = 0.5*wall_avg_psi[iwall]*
                    Area_surf_el[idim]*Elec_param_w[iwall]/ Nnodes_per_el_S;
           /*printf("wall_count %d    Nnodes_per_el_S %d   S_area_tot %6.3f\n",
             wall_count[iwall],Nnodes_per_el_S,S_area_tot[iwall]);*/
 /*  if (Proc==0 && nblock != 0)
      printf("Avg deriv: blocked %9.6f     open %9.6f\n",sblock/nblock,
               sopen/(wall_count[iwall]-nblock));*/
    }
/*if (Proc==0) printf("Surface integral: %9.6f from grad_psi.n    %9.6f from q\n",
          sum,*sum3);*/

/*  if (sum != 0.0 || sum2 != 0.0) 
    printf ("Proc: %d  surface term: %9.6f  volume term: %9.6f\n",
             Proc,sum,sum2*Temp_elec/(8.0*PI));*/

/*  if (Type_bc_elec[WallType[iwall]]==2) sum = *sum3;*/
    return (sum);
}
/*******************************************************************************
 * charge_stress: For charged surfaces, find the contribution of the stress *
 *             tensor to the local pressure.                                */

double charge_stress(double **x,double *sum2)
{
   int loc_inode,iwall,idim,
       inode,ijk[3],inode_box,iel_box;
   int loc_node_el,iel,reflect_flag[3]={FALSE,FALSE,FALSE};
   double deriv2_x[3];
   double sum3=0.0;
   *sum2=0.0;

   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){

       inode_box = L2B_node[loc_inode];
       inode = L2G_node[loc_inode];
       node_to_ijk(inode,ijk);

       iwall = Nodes_2_boundary_wall[Nlists_HW-1][inode_box];

       /* now calculate the volume term !! */

       for (loc_node_el=0; loc_node_el< Nnodes_per_el_V ; loc_node_el++){

          iel = node_to_elem(inode, loc_node_el, reflect_flag);
          iel_box = el_to_el_box(iel);

          if (iel >= 0){

              for (idim=0; idim<Ndim; idim++){
                 deriv2_x[idim] = 0.0;
                 switch(Ndim)
                 {
                  case 3:
                    if ( 
                         (idim == 0 && (loc_node_el==0 || loc_node_el == 2
                                   || loc_node_el==4 || loc_node_el==6 )  
                                    && ijk[0] < Nodes_x[0]-1 )              ||
                         (idim == 1 && (loc_node_el==0 || loc_node_el == 1
                                   || loc_node_el==4 || loc_node_el==5 ) 
                                    && ijk[1] < Nodes_x[1]-1 )              ||
                         (idim == 2 && (loc_node_el==0 || loc_node_el == 1
                                   || loc_node_el==2 || loc_node_el==3 ) 
                                    && ijk[2] < Nodes_x[2]-1 )               
                         ){
                         deriv2_x[idim] = calc_deriv2(idim,inode_box,FFD,x);
                         sum3 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                         if (Wall_elems[Nlists_HW-1][iel_box] != -1)
                           *sum2 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                    }
                    break;
                  case 2:
                    if ( 
                         (idim == 0 && (loc_node_el==0 || loc_node_el == 2 ) 
                                    && ijk[0] < Nodes_x[0]-1               ) ||
                         (idim == 1 && (loc_node_el==0 || loc_node_el == 1 ) 
                                    && ijk[1] < Nodes_x[1]-1               )
                         ){
                         deriv2_x[idim] = calc_deriv2(idim,inode_box,FFD,x);
                         sum3 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                         if (Wall_elems[Nlists_HW-1][iel_box] != -1)
                           *sum2 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                    }
                    break;
                  case 1:
                    if ( loc_node_el==0 && ijk[0] < Nodes_x[0]-1 ){
                         deriv2_x[idim] = calc_deriv2(idim,inode_box,FFD,x);
                         sum3 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                         if (Wall_elems[Nlists_HW-1][iel_box] != -1)
                           *sum2 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                    }
                    break;
                 }
              }   /* end of idim loop */

          }       /* end of test for a wall element */
       }          /* end of loop over elements that include inode */


    }       /* end of loop over nodes on this processor */

    *sum2 *= Temp_elec/(8.0*PI);
    sum3 *= Temp_elec/(8.0*PI);
/*  if (Proc==0) printf ("integral of (grad_phi)^2 over wall: %9.6f   total volume: %9.6f\n",*sum2,sum3);*/

    return (sum3);
}
/****************************************************************************
 * free_energy_charging_up:  In this routine we calculate the surface       *
 *                           integral that gives the free energy due to     *
 *                           charging up the surfaces                       */
double free_energy_charging_up(double **x)
{
  int loc_inode,idim,inode,iunk, inode_box;
  double psi_i,charging_int,charge_i=0.0;
  int i,iel_w,surf_norm,ilist;
  double psi_deriv[2][2],prefac,dphi_term;

  charging_int = 0.0;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
      inode = L2G_node[loc_inode];
      inode_box = L2B_node[loc_inode];

      if  (Nodes_2_boundary_wall[Nlists_HW-1][node_to_node_box(inode)] != -1){

         iunk = Phys2Unk_first[POISSON];

         ilist = Nlists_HW - 1;
         for (iel_w=0; iel_w<Nelems_S[ilist][loc_inode]; iel_w++){
              
              surf_norm = Surf_normal[ilist][loc_inode][iel_w];
              idim = abs(surf_norm) - 1;

              prefac = (double)(surf_norm/abs(surf_norm))*Area_surf_el[idim];

              if (prefac > 0.0) {  /* RHS*/
                   i = 1;
                   psi_deriv[i][0] = (
                               -    x[iunk][inode_box+2] 
                               + 4.*x[iunk][inode_box+1] 
                               - 3.*x[iunk][inode_box]         )/(2.*Esize_x[0]);
                   psi_deriv[i][1] = (x[iunk][inode_box] -
                                      x[iunk][inode_box-1])/Esize_x[0];
              }
              else {              /* LHS */
                    i = 0;
                    psi_deriv[i][0] = -(
                                  3.*x[iunk][inode_box] 
                                - 4.*x[iunk][inode_box-1]
                                +    x[iunk][inode_box-2])
                                    /(2.0*Esize_x[0]);
                    psi_deriv[i][1] = (x[iunk][inode_box+1] 
                                                      - x[iunk][inode_box])/Esize_x[0];
              }
              charge_i = (Temp_elec/(8.0*PI))*psi_deriv[i][0]; 


          } /* end of surface element loop */ 

/*         charge_i = 0.0;
         for (idim=0; idim<Ndim; idim++){
             charge_i -= Charge_w_sum_els[loc_inode][idim]*Area_surf_el[idim];
         }
 */
 

         psi_i = x[iunk][inode_box];
         charging_int -= (psi_i * charge_i);
/*         printf("loc_inode: %d  charge_i: %9.6f   psi_i: %9.6f  charging_int: %9.6f\n",
                 loc_inode,charge_i,psi_i,charging_int);   */

      }   /* end of test for boundary node */
  }       /* end of loc_inode loop */

/*  printf("PI=%9.6f\n",PI);   */
  dphi_term = (Temp_elec/(8.0*PI))*WallParam[0]
                  *psi_deriv[1][1]*psi_deriv[1][1];
  printf("dphi_term: %9.6f  charging_int: %9.6f \n",dphi_term,charging_int);
  charging_int += dphi_term;

  return charging_int;

  /* collate all the charging_int values from the different processors */
}
/****************************************************************************
 * calc_deriv_e : calculate a derivative of the electric potential!!   */

double calc_deriv_e(int idim,int inode0,int flag,int *blocked, double **x, 
                  int ilist)
{
   int inode1,inode2,offset1[3],offset2[3],
       iwall1,iwall2,
       jdim,ijk_box[3],reflect_flag[3],iunk;
   double deriv=0.0;

   node_box_to_ijk_box(inode0,ijk_box);

   for (jdim=0; jdim<Ndim; jdim++) {
       offset1[jdim] = 0;
       offset2[jdim] = 0;
   }

   switch(flag)
   {
      case CFD:
        offset1[idim] = -1;
        offset2[idim] =  1;
        break;

      case FFD:
        offset1[idim] =  1;
        offset2[idim] =  2;
        break;

      case BFD:
        offset1[idim] = -1;
        offset2[idim] = -2;
        break;
   }

   inode1 = offset_to_node_box(ijk_box,offset1,reflect_flag);
   inode2 = offset_to_node_box(ijk_box,offset2,reflect_flag);
   iunk = Phys2Unk_first[POISSON];

   iwall1 = Nodes_2_boundary_wall[ilist][inode1];
   iwall2 = Nodes_2_boundary_wall[ilist][inode2];

   if (iwall1 == -1 && iwall2 == -1) {
      *blocked = FALSE;
      switch(flag)
      {
         case CFD:
           deriv =    x[iunk][inode2] - x[iunk][inode1];break;
         case FFD:
           deriv = -3*x[iunk][inode0] + 4*x[iunk][inode1] - x[iunk][inode2];break;
         case BFD:
           deriv =  3*x[iunk][inode0] - 4*x[iunk][inode1] + x[iunk][inode2];break;
      }
      deriv /= (2.0*Esize_x[idim]);
   }
   else{
      *blocked = TRUE;
      switch(flag)
      {
         case FFD:
           deriv = x[iunk][inode1] - x[iunk][inode0]; break;
         case BFD:
           deriv = x[iunk][inode0] - x[iunk][inode1]; break;
      }
      deriv /= (Esize_x[idim]);
   }

   return deriv;
}
/***************************************************************************
 * calc_deriv2 : calculate a derivative of the electric potential (first order)!!   */

double calc_deriv2(int idim,int inode0,int flag,double **x)
{
   int inode1,offset[3], iunk,
       jdim,ijk_box[3],reflect_flag[3];
   double deriv=0.0;


   node_box_to_ijk_box(inode0,ijk_box);

   for (jdim=0; jdim<Ndim; jdim++) offset[jdim] = 0;

   switch(flag)
   {
      case FFD:
        offset[idim] =  1;
        break;

      case BFD:
        offset[idim] = -1;
        break;
   }

   inode1 = offset_to_node_box(ijk_box,offset,reflect_flag);
   iunk = Phys2Unk_first[POISSON];

   switch(flag)
   {
      case FFD:
        deriv = x[iunk][inode1] - x[iunk][inode0]; break;
      case BFD:
        deriv = x[iunk][inode0] - x[iunk][inode1];  break;
   }
   deriv /= (Esize_x[idim]);

   return deriv;
}
/*****************************************************************************/

