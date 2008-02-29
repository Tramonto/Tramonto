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
// but WITHOUT ANY WARRANTY; without even the implied warranty of// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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
 *  FILE: dft_linsolver_setup.c
 *
 *  This file contains logic to control which of the physics specific
 *  Schur solver approaches will be used to solve the problem at hand.
 */
#include "dft_linsolver_setup.h"

/*******************************************************************************/
void linsolver_setup_control()
{
   if (L_Schur && Type_poly == CMS)    linsolver_setup_CMSTYPE();
   else if (L_Schur && Type_poly == WJDC)   linsolver_setup_WJDCTYPE();
   else if (L_Schur && Type_func != NONE) linsolver_setup_HSTYPE();
   else {
    //   LinProbMgr_manager = dft_basic_lin_prob_mgr_create(Nunk_per_node, Aztec.options, Aztec.params, MPI_COMM_WORLD);
         LinProbMgr_manager = dft_basic_lin_prob_mgr_create(Nunk_per_node, ParameterList_list, MPI_COMM_WORLD);
   }
   return;
}
/*******************************************************************************/
void linsolver_setup_CMSTYPE()
{
  int iunk,i;
  int *geq, *ginveq, *cmseq, *densityeq;
  int count_density,count_cms_field,count_geqn,count_ginv_eqn;
  int index_save;
  int count_poisson;
  int *poissoneq;

  /* Construct dft_Linprobmgr with information on number of unknowns*/
   densityeq = (int *) array_alloc(1, Nunk_per_node, sizeof(int));
   cmseq = (int *) array_alloc(1, Nunk_per_node, sizeof(int));
   geq = (int *) array_alloc(1, Nunk_per_node, sizeof(int));
   ginveq = (int *) array_alloc(1, Nunk_per_node,  sizeof(int));
   poissoneq = (int *) array_alloc(1, Nunk_per_node, sizeof(int));
   count_poisson = 0;
   
   count_density=count_cms_field=count_geqn=count_ginv_eqn=0;
   for (iunk=0;iunk<Nunk_per_node;iunk++){
     switch(Unk2Phys[iunk]){
     case DENSITY:
       densityeq[count_density++]=iunk; break; 
     case CMS_FIELD:                  
       cmseq[count_cms_field++]=iunk; break; 
     case G_CHAIN:                  
       if ((iunk-Geqn_start[0])%2 == 0)
         geq[count_geqn++]=iunk; 
       else
         ginveq[count_ginv_eqn++]=iunk; 
       break;
     case POISSON:
       poissoneq[count_poisson++]=iunk; break;
     default:
        printf("ERROR: every unknown should be linked to a physics type and added to id lists for solver iunk=%d\n",iunk);
        exit(-1);
        break;
     } 
   }
   /* now invert the order of the ginverse equations ! */
   for (i=0;i<count_ginv_eqn/2;i++){
     index_save = ginveq[i];
     ginveq[i] = ginveq[count_ginv_eqn-1-i];
     ginveq[count_ginv_eqn-1-i]=index_save;
   }
   // LinProbMgr_manager = dft_poly_lin_prob_mgr_create(Nunk_per_node, Aztec.options, Aztec.params, MPI_COMM_WORLD);
   LinProbMgr_manager = dft_poly_lin_prob_mgr_create(Nunk_per_node, ParameterList_list, MPI_COMM_WORLD);
   dft_poly_lin_prob_mgr_setgequationids(LinProbMgr_manager, Ngeqn_tot/2, geq);
   dft_poly_lin_prob_mgr_setginvequationids(LinProbMgr_manager, Ngeqn_tot/2, ginveq);
   dft_poly_lin_prob_mgr_setcmsequationids(LinProbMgr_manager, Ncomp, cmseq);
   dft_poly_lin_prob_mgr_setdensityequationids(LinProbMgr_manager, Ncomp, densityeq);
   dft_poly_lin_prob_mgr_setpoissonequationids(LinProbMgr_manager, count_poisson, poissoneq);
   /*dft_poly_lin_prob_mgr_setfieldondensityislinear(LinProbMgr_manager,TRUE);*/
   safe_free((void *) &densityeq);
   safe_free((void *) &cmseq);
   safe_free((void *) &geq);
   safe_free((void *) &ginveq);
   safe_free((void *) &poissoneq);
}
/*******************************************************************************/
void linsolver_setup_HSTYPE()
{
  int iunk,i;
  int *densityeq, *indnonlocaleq, *depnonlocaleq;
  int count_density, count_indnonlocal,count_depnonlocal;
  int one_particle_size;

  /* Construct dft_Linprobmgr with information on number of unknowns*/

   densityeq = (int *) array_alloc(1, Nunk_per_node, sizeof(int));
   indnonlocaleq = (int *) array_alloc(1, Nunk_per_node, sizeof(int));
   depnonlocaleq = (int *) array_alloc(1, Nunk_per_node, sizeof(int));

   count_density=count_indnonlocal=count_depnonlocal=0;
   one_particle_size=FALSE;
   if ((Lhard_surf && Nlists_HW == 2) || (!Lhard_surf && Nlists_HW == 1)) one_particle_size=TRUE;
   for (iunk=0;iunk<Nunk_per_node;iunk++){
     switch(Unk2Phys[iunk]){
     case POISSON:
     case DENSITY:
       densityeq[count_density++]=iunk; break; 
     case HSRHOBAR:                  
       if (one_particle_size && (
             (iunk > Phys2Unk_first[HSRHOBAR]+1 && 
              iunk < Phys2Unk_first[HSRHOBAR]+Nrho_bar_s)  ||
              iunk >= Phys2Unk_first[HSRHOBAR]+Nrho_bar_s+Ndim) )
         depnonlocaleq[count_depnonlocal++]=iunk; 
       else
         indnonlocaleq[count_indnonlocal++]=iunk; 
       break;
     case BONDWTC:
/*       if (Pol_Sym[iunk-Phys2Unk_first[BONDWTC]] == -1)
          indnonlocaleq[count_indnonlocal++]=iunk; 
       else 
          depnonlocaleq[count_depnonlocal++]=iunk; 
       break;*/
     case CAVWTC:
       indnonlocaleq[count_indnonlocal++]=iunk; break;   
      default:
        printf("ERROR: every unknown should be linked to a physics type and added to id lists for solver iunk=%d\n",iunk);
        exit(-1);
        break;
     } 
   }
   // LinProbMgr_manager = dft_hardsphere_lin_prob_mgr_create(Nunk_per_node, Aztec.options, Aztec.params, MPI_COMM_WORLD);
   LinProbMgr_manager = dft_hardsphere_lin_prob_mgr_create(Nunk_per_node, ParameterList_list, MPI_COMM_WORLD);
   dft_hardsphere_lin_prob_mgr_setindnonlocalequationids(LinProbMgr_manager, count_indnonlocal, indnonlocaleq);
   dft_hardsphere_lin_prob_mgr_setdepnonlocalequationids(LinProbMgr_manager, count_depnonlocal, depnonlocaleq);
   dft_hardsphere_lin_prob_mgr_setdensityequationids(LinProbMgr_manager, count_density, densityeq);
   if (Type_attr != NONE || Type_poly==WTC || Mesh_coarsening || Type_coul != NONE)
                dft_hardsphere_lin_prob_mgr_seta22blockisdiagonal(LinProbMgr_manager, FALSE);
   else         dft_hardsphere_lin_prob_mgr_seta22blockisdiagonal(LinProbMgr_manager, TRUE);
   safe_free((void *) &densityeq);
   safe_free((void *) &indnonlocaleq);
   safe_free((void *) &depnonlocaleq);
   return;
}
/*******************************************************************************/
void linsolver_setup_WJDCTYPE()
{
  int iunk,i;
  double **xOwned, **x2Owned;
  int *geq, *ginveq, *wjdceq, *densityeq, *indnonlocaleq, *depnonlocaleq;
  int count_density,count_cms_field,count_geqn,count_ginv_eqn;
  int count_indnonlocal,count_depnonlocal,index_save;
  int one_particle_size;
  int count_poisson;
  int *poissoneq;

  /* Construct dft_Linprobmgr with information on number of unknowns*/
  densityeq = (int *) array_alloc(1, Nunk_per_node, sizeof(int));
  indnonlocaleq = (int *) array_alloc(1, Nunk_per_node, sizeof(int));
  depnonlocaleq = (int *) array_alloc(1, Nunk_per_node, sizeof(int));
  wjdceq = (int *) array_alloc(1, Nunk_per_node, sizeof(int));
  geq = (int *) array_alloc(1, Nunk_per_node, sizeof(int));
  ginveq = (int *) array_alloc(1, Nunk_per_node,  sizeof(int));
  poissoneq = (int *) array_alloc(1, Nunk_per_node, sizeof(int));
   
  count_density=count_cms_field=count_geqn=count_ginv_eqn=0;
  count_indnonlocal=count_depnonlocal=0;
  count_poisson = 0;

  one_particle_size=FALSE;
  if ((Lhard_surf && Nlists_HW == 2) || (!Lhard_surf && Nlists_HW == 1)) one_particle_size=TRUE;

  for (iunk=0;iunk<Nunk_per_node;iunk++){
     switch(Unk2Phys[iunk]){
     case DENSITY:
       densityeq[count_density++]=iunk; break; 
     case HSRHOBAR:
       if (one_particle_size && (
             (iunk > Phys2Unk_first[HSRHOBAR]+1 &&
              iunk < Phys2Unk_first[HSRHOBAR]+Nrho_bar_s)  ||
              iunk >= Phys2Unk_first[HSRHOBAR]+Nrho_bar_s+Ndim) )
         depnonlocaleq[count_depnonlocal++]=iunk;
       else
         indnonlocaleq[count_indnonlocal++]=iunk;
       break;
     case CAVWTC:
       indnonlocaleq[count_indnonlocal++]=iunk; break;
     case WJDC_FIELD:                  
       wjdceq[count_cms_field++]=iunk; break; 
     case G_CHAIN:                  
       if ((iunk-Geqn_start[0])%2 == 0)
         geq[count_geqn++]=iunk; 
       else
         ginveq[count_ginv_eqn++]=iunk; 
       break;
     case POISSON:
       poissoneq[count_poisson++]=iunk; break; 
     default:
        printf("ERROR: every unknown should be linked to a physics type and added to id lists for solver iunk=%d\n",iunk);
        exit(-1);
        break;
     } 
  }
  /* now invert the order of the ginverse equations ! */
  for (i=0;i<count_ginv_eqn/2;i++){
     index_save = ginveq[i];
     ginveq[i] = ginveq[count_ginv_eqn-1-i];
     ginveq[count_ginv_eqn-1-i]=index_save;
  }


   // LinProbMgr_manager = dft_wjdc_lin_prob_mgr_create(Nunk_per_node, Aztec.options, Aztec.params, MPI_COMM_WORLD);

/*  Turn on these calls when the routines exist */
/*
   LinProbMgr_manager = dft_wjdctype_lin_prob_mgr_create(Nunk_per_node, ParameterList_list, MPI_COMM_WORLD);
  dft_wjdctype_lin_prob_mgr_setindnonlocalequationids(LinProbMgr_manager, count_indnonlocal, indnonlocaleq);
  dft_wjdctype_lin_prob_mgr_setdepnonlocalequationids(LinProbMgr_manager, count_depnonlocal, depnonlocaleq);
  dft_wjdctype_lin_prob_mgr_setgequationids(LinProbMgr_manager, Ngeqn_tot/2, geq);
  dft_wjdctype_lin_prob_mgr_setginvequationids(LinProbMgr_manager, Ngeqn_tot/2, ginveq);
  dft_wjdctype_lin_prob_mgr_setwjdcequationids(LinProbMgr_manager, Ncomp, wjdceq);
  dft_wjdctype_lin_prob_mgr_setdensityequationids(LinProbMgr_manager, Ncomp, densityeq);
  dft_wjdctype_lin_prob_mgr_setpoissonequationids(LinProbMgr_manager, count_poisson, poissoneq);
*/

   /*dft_poly_lin_prob_mgr_setfieldondensityislinear(LinProbMgr_manager,TRUE);*/
  safe_free((void *) &densityeq);
  safe_free((void *) &indnonlocaleq);
  safe_free((void *) &depnonlocaleq);
  safe_free((void *) &wjdceq);
  safe_free((void *) &geq);
  safe_free((void *) &ginveq);
  safe_free((void *) &poissoneq);
  return;
}
/*******************************************************************************/



