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

/* ---------------------------------------------------------
dft_out_utils.c:

Here are some routines that are used generically in the postprocessing
phase of Tramonto - function pointers specify the particulars of integrands.
------------------------------------------------------------*/

#include "dft_out_utils.h"

/*****************************************************************************************/
double integrateInSpace(double(*fp_integrand)(int,int,double **),int iunk,
                                                int **nelhit,double **x,double *profile)
{

   double sum,sum_i,integrand;
   int loc_inode,inode_box,iwall;
   double area;

   sum_i=0.0,sum=0.0;
   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
      inode_box = L2B_node[loc_inode];
      integrand = (*fp_integrand)(iunk,inode_box,x);
      sum_i += integrand*nelhit[iunk-Phys2Unk_first[DENSITY]][inode_box]*Vol_el/((double)Nnodes_per_el_V);
      if (profile != NULL) profile[loc_inode]+=integrand;
  }       /* end of loc_inode loop */

  sum=gsum_double(sum_i);

  if (Lper_area && Area>0.0) {
      sum /= Area;
      sum *= (Fac_vol/Fac_area);
  }
  else sum *= Fac_vol;

  return(sum);
}
/*****************************************************************************************************/
double integrateInSpace_SumInComp(double(*fp_integrand)(int,int,double**),
                                                int **nelhit,double **x,double *profile)
{

   double sum,sum_i,integrand;
   int loc_inode,inode_box,iloop,nloop,iunk;

   sum_i=0.0,sum=0.0;

   nloop=Ncomp;
   if (Lseg_densities) nloop=Nseg_tot;

   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
    inode_box = L2B_node[loc_inode];
    for (iloop=0; iloop<nloop; iloop++){

         iunk = Phys2Unk_first[DENSITY]+iloop;
         integrand = (*fp_integrand)(iunk,inode_box,x);
         sum_i += integrand*
                   nelhit[iloop][inode_box]*Vol_el/((double)Nnodes_per_el_V);

         if (profile != NULL) profile[loc_inode]+=integrand;
    }   
  }       /* end of loc_inode loop */

  sum=gsum_double(sum_i);

  if (Lper_area && Area>0.0) {
      sum /= Area;
      sum *= (Fac_vol/Fac_area);
  }
  else{
      sum*= Fac_vol;
  }
  return(sum);
}
/*****************************************************************************************************/
double integrateOverSurface(double(*fp_integrand)(int,int,int,double **), int iunk, double **x, double *profile)
{

  double sum,sum_i,integrand,fac;
  int loc_inode,inode_box,iwall,ijk[3],idim;

  sum_i=0.0,sum=0.0;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){

      inode_box = L2B_node[loc_inode];     
      iwall=Nodes_2_boundary_wall[Nlists_HW-1][inode_box];
      if  (iwall != -1){ /* identify surface nodes */

         integrand = (*fp_integrand)(iunk,inode_box,iwall,x);

         fac=1.0;
         node_to_ijk(L2G_node[loc_inode],ijk);
         for (idim=0;idim<Ndim;idim++){ 
            if ((ijk[idim]==0 || ijk[idim]==Nodes_x[idim]-1)  && Type_bc[idim][1] != PERIODIC ) fac*=0.5;
         } 
         sum_i += integrand*fac;
         if (profile != NULL) profile[loc_inode]+=integrand;
      }       /* end of loc_inode loop */
  }

  sum=gsum_double(sum_i);

  if (Lper_area && Area>0.0) {
      sum /= Area;
      sum *= (Fac_vol/Fac_area);
  }
  else sum *= Fac_vol;

  return(sum);
}
/*****************************************************************************************/
/*setup_domain_multipliers: Here compute area and Fac_vol and Fac_area that
  are needed for all integrals */
void setup_domain_multipliers()
{
  int idim,iwall;

  /* compute Fac_vol and Fac_area */
   Fac_area = 1.0;
   Fac_vol = 1.0;
   for (idim = 0; idim<Ndim; idim++) {
       if (!(Type_bc[idim][0] == REFLECT && Type_bc[idim][1] == REFLECT) && Lcount_reflect){

       if (Type_bc[idim][0] == REFLECT){

          Fac_vol *= 2.0;
          if (WallPos[idim][0] == -0.5*Size_x[idim]) Fac_area *= 2.0;
       }
       else if (Type_bc[idim][1] == REFLECT) {

          Fac_vol *= 2.0;
          if (WallPos[idim][0] == 0.5*Size_x[idim]) Fac_area *= 2.0;
       }
       }
   }


  /* compute surface area */
  Area = 0.0;
  if (Nwall == 0) Area = 1.0;
  else{
     if (Nlink == Nwall) Area = S_area_tot[Nlists_HW-1][0];
     else
        for (iwall=0; iwall<Nwall; iwall++){
           if (Link[iwall]==0)
           Area += S_area_tot[Nlists_HW-1][iwall];
        }
  }

  return;
}
/***************************************************************************/
/* setup_integrals:  here we store arrays of elements hit per node and
                     the list propertiese for post processing integrated
                     parameters (adsorption, free energy etc.). */
void setup_integrals()
{
  int loc_inode, icomp, iel, iunk, inode_box, nel_hit,inode,iel_box,
      nel_hit2,ilist,idim,ielement,semiperm,iwall,jcomp;
  int reflect_flag[3],ijk[3],i,nloop,iloop;
  
  
  for (idim=0; idim<Ndim; idim++) reflect_flag[idim]=FALSE;
  
  if (Lseg_densities) nloop=Nseg_tot;
  else                nloop=Ncomp;

  Nel_hit = (int **) array_alloc (2,nloop, Nnodes_box, sizeof(int));
  Nel_hit2 = (int **) array_alloc (2,nloop, Nnodes_box, sizeof(int));
  
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
      inode = L2G_node[loc_inode];
      inode_box = L2B_node[loc_inode];
      node_to_ijk(inode,ijk);
      
      for (iloop=0; iloop<nloop; iloop++){
      
         if (Lseg_densities) icomp=Unk2Comp[iloop];
         else                icomp=iloop;
         
         if (Lhard_surf){
            if (Nlists_HW==1 || Nlists_HW==2) List[0] = 0;
            else                              List[0] = icomp;
            List[1] = Nlists_HW-1;            
         }  
         else{
             List[0]=0;
             List[1]=0;
         }

         if (Lhard_surf && Nwall>0) Imax = 2;
         else Imax = 1;

         /* define integration rules to density discontinuities */
         nel_hit = Nnodes_per_el_V;
         for (iel=0; iel<Nnodes_per_el_V; iel++){
             ielement = node_to_elem(inode,iel,reflect_flag);
             iel_box = el_to_el_box(ielement);
             if (ielement == -1) nel_hit--;
             else if (ielement != -2){
                iwall =  Wall_elems[List[1]][iel_box];
                semiperm=FALSE;
                for (jcomp=0; jcomp<Ncomp; jcomp++)
                  if (iwall>=0 && Lsemiperm[WallType[iwall]][jcomp]) semiperm=TRUE;
                if (iwall !=-1 && !semiperm) nel_hit--;
             }
         }

         /* define integration rules to core surfaces */
         nel_hit2 = Nnodes_per_el_V;
         for (iel=0; iel<Nnodes_per_el_V; iel++){
             ielement = node_to_elem(inode,iel,reflect_flag);
             iel_box = el_to_el_box(ielement);
             if (ielement == -1) nel_hit2--;
             else if (ielement != -2){
                iwall =  Wall_elems[List[0]][iel_box];
                semiperm=FALSE;
                for (jcomp=0; jcomp<Ncomp; jcomp++)
                  if (iwall>=0 && Lsemiperm[WallType[iwall]][jcomp]) semiperm=TRUE;
                if (iwall!=-1 && !semiperm) nel_hit2--;
             }
         }

         for (idim=0; idim<Ndim; idim++){
            if ( (Type_bc[idim][0] == REFLECT || Type_bc[idim][0] == IN_BULK || Type_bc[idim][0]==LAST_NODE)
                                                          &&  ijk[idim] == 0)  {
                  nel_hit /= 2;
                  nel_hit2 /= 2;
            }
            if ( (Type_bc[idim][1] == REFLECT || Type_bc[idim][1] == IN_BULK || Type_bc[idim][0]==LAST_NODE)
                                            &&  ijk[idim] == Nodes_x[idim]-1)  {
                  nel_hit /= 2;
                  nel_hit2 /= 2;
            }
         }

         Nel_hit[iloop][inode_box]=nel_hit;
         Nel_hit2[iloop][inode_box]=nel_hit2;

      } /* end iloop loop */
    }
}
/**************************************************************************************/
