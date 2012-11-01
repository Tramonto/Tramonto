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
 *  FILE: dft_vext.c
 *
 *  This file contains routines to set up an external field for a 
 *  variety of cases.
 */

#include "dft_vext.h"

/******************************************************************************/
void setup_external_field_n( int **nelems_w_per_w, int ***elems_w_per_w) 

/* Set up an external field for the Classical Fluids Density Functional Theory
 *
 *    Authors:  Laura Frink, 9225
 */

{
   /* Local variable declarations */
   
   char *yo = "external_field_neutral";
   int idim,icomp,iwall,iunk,loc_inode,ilist;
   int **elems_w_per_w_global, *nelems_w_per_w_global;
   int inode_box;
   /*int **nodes_vext_max,*nnodes_vext_max;*/
   /*int i,ijk[3],inode;*/
   double *vext_tmpOwned, *vext_tmpBox;
   double t1=0.0,t_zeroTF=0.0,t_vext=0.0;
   double t_tot_min,t_zero_min,t_vext_min;
   double t_tot_max,t_zero_max,t_vext_max;
   int count,count_max=0;
  
  /********************** BEGIN EXECUTION ************************************/

  if (Proc==0 &&Iwrite==VERBOSE) printf("\n %s: Setting up External Field ... \n",yo);
  t1 -= MPI_Wtime();
  /*
   * Allocate and zero the arrays we will calculate here
   */

   Vext = (double **) array_alloc (2,Nnodes_per_proc,Ncomp, sizeof(double));
   Vext_set = (double **) array_alloc (2,Nnodes_per_proc,Ncomp, sizeof(double));

   /* set up the Vext_dash array if we will need it */
   Lvext_dash=FALSE;
   for (iwall=0;iwall<Nwall;iwall++){
      if (Ipot_wf_n[WallType[iwall]] != VEXT_NONE && Ipot_wf_n[WallType[iwall]] != VEXT_HARD &&
         Ipot_wf_n[WallType[iwall]] != VEXT_3D_INTEGRATED) Lvext_dash=TRUE;
   }
   if (Nwall>0 && Lvext_dash) 
         Vext_dash =  (double ***) array_alloc (3, Nnodes_per_proc,Ncomp*Nwall, 
                                          Ndim, sizeof(double));

  for (icomp=0; icomp<Ncomp; icomp++){
     for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){

        Vext[loc_inode][icomp] = UNINIT_VEC;

        if (Lvext_dash){
        for (iwall=0; iwall<Nwall; iwall++){
          for (idim=0; idim<Ndim; idim++){
             iunk = iwall*Ncomp + icomp;
             Vext_dash[loc_inode][iunk][idim] = UNINIT_VEC;
          }
        }
        }

     }
  }

  /* 
   * For each desired case, write a new subroutine, 
   * add the case as a choice to the input file, 
   * and add a case statement here.    
   */

  setup_semiperm(nelems_w_per_w,elems_w_per_w);
  setup_vext_max();
  setup_zero();


  t_vext -=MPI_Wtime();
  for (iwall=0;iwall<Nwall_Images;iwall++){
     switch(Ipot_wf_n[WallType_Images[iwall]])
     {
         
       case VEXT_NONE: 
       case VEXT_HARD : break;

       case VEXT_DIST_TO_SURF:              /* 1-dimensional LJ potential */
       case VEXT_DIST_TO_CENTER:            /* 12-6 LJ potentials for atoms : a 3D potential for 3D systems */
         setup_vext_XRSurf(iwall);
         break;

       case VEXT_3D_INTEGRATED:             /* a more proper treatment of 2D or 3D systems where the 12-6 LJ 
                                        potential is explicitly integratedover unusual geometries */
         elems_w_per_w_global = (int **) array_alloc(1,Nlists_HW, sizeof(int *));
         nelems_w_per_w_global = (int *) array_alloc(1,Nlists_HW, sizeof(int));
         comm_wall_els(iwall,nelems_w_per_w,elems_w_per_w, nelems_w_per_w_global, elems_w_per_w_global);
/*         if (RealWall_Images[iwall]<Nwall){*/
               setup_integrated_LJ_walls(iwall,nelems_w_per_w_global, elems_w_per_w_global);

               for (ilist=0; ilist<Nlists_HW; ilist++){
                   safe_free((void *) &elems_w_per_w_global[ilist]);
                }
                safe_free((void *) &elems_w_per_w_global);
                safe_free((void *) &nelems_w_per_w_global);
         /*}*/
       
         break;

       default:
           printf ("ERROR:no function set up for the chosen external field\n");
           printf ("iwall=%d  Type=%d Ipot_wf_n=%d\n",RealWall_Images[iwall],WallType[iwall],Ipot_wf_n[WallType[iwall]]);
           exit (-1);
     }
  }
  t_vext +=MPI_Wtime();

/*  Now set up Zero_density_TF everywhere in the box based on
    the external field array */

  vext_tmpOwned = (double *) array_alloc (1,Nnodes_per_proc,sizeof(double));
  vext_tmpBox = (double *) array_alloc (1,Nnodes_box,sizeof(double));
  if (Nwall > 0){
     for (icomp=0; icomp<Ncomp; icomp++){
         for (loc_inode=0;loc_inode<Nnodes_per_proc;loc_inode++){
              vext_tmpOwned[loc_inode]=Vext[loc_inode][icomp];
              vext_tmpBox[L2B_node[loc_inode]]=vext_tmpOwned[loc_inode];
         }
         (void) dft_linprobmgr_importnodalr2c(LinProbMgr_manager,vext_tmpOwned,vext_tmpBox);
         for (inode_box=0;inode_box<Nnodes_box;inode_box++){
              if (vext_tmpBox[inode_box] >=VEXT_MAX) Zero_density_TF[inode_box][icomp]=TRUE;
              else Zero_density_TF[inode_box][icomp]=FALSE;
         }
     }
  }
  safe_free((void *) &vext_tmpOwned);
  safe_free((void *) &vext_tmpBox);

  /*if (Nwall > 0) {

     t_zeroTF -= MPI_Wtime();
     nodes_vext_max = (int **) array_alloc (1,Ncomp,sizeof(int *));
     nnodes_vext_max = (int *) array_alloc (1,Ncomp,sizeof(int));

     comm_vext_max(nnodes_vext_max,nodes_vext_max);

     for (icomp=0; icomp<Ncomp; icomp++){
        for (i=0; i<nnodes_vext_max[icomp]; i++){
            inode = nodes_vext_max[icomp][i];
            node_to_ijk(inode,ijk);
            inode_box = node_to_node_box_no_bound(inode);
            if (inode_box >=0){
               Zero_density_TF[inode_box][icomp] = TRUE;
            }
        }
     }

     for (icomp=0; icomp<Ncomp; icomp++)
           safe_free((void *) &nodes_vext_max[icomp]);

     safe_free((void *) &nodes_vext_max);
     safe_free((void *) &nnodes_vext_max);

     t_zeroTF += MPI_Wtime();
  }*/

/* Finally, we may still have discrepancies in the Zero_TF array
   on different processors when there are semi-permeable
   membranes present. */

   /*if (Num_Proc>1) correct_zeroTF_array();*/

  safe_free((void *) &Vext_set);

  t1 += MPI_Wtime();

  t_zero_max = gmax_double(t_zeroTF);
  t_vext_max = gmax_double(t_vext);
  t_tot_max = gmax_double(t1);

  t_zero_min = gmin_double(t_zeroTF);
  t_vext_min = gmin_double(t_vext);
  t_tot_min = gmin_double(t1);

  count=0;
  for (inode_box=0; inode_box<Nnodes_box; inode_box++)
     if (B2L_node[inode_box] >=0){
        if (!Zero_density_TF[inode_box][0]) count++;
     }
  count_max = gsum_int(count);


  if (Proc==0&&Iwrite==VERBOSE) {
    printf("---------------------------------------------------------------\n");
    printf("---TIMINGS IN EXTERNAL FIELD CALC-----\n");
    printf("\t\t\tMAX \t\tMIN .... all in seconds\n");
    printf("total  \t\t\t %g \t\t%g\n",t_tot_max,t_tot_min);
    printf("vext calc \t\t %g \t\t%g\n",t_vext_max,t_vext_min);
    printf("Zero_density_TF \t\t %g \t\t%g\n",t_zero_max,t_zero_min);
    printf("Zero_density_TF : %d nonzeros of %d Nodes \n",count_max,Nnodes);
    printf("---------------------------------------------------------------\n");
  }
  return;
}
/*****************************************************************************/
/* setup_zero: Set the external field to zero everywhere in the fluid */
void setup_zero()
{
   int loc_inode, icomp, iwall, iunk, idim,inode_box;
   for (icomp=0; icomp<Ncomp; icomp++) {
      
      for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {
         inode_box=L2B_node[loc_inode];

          if (!Zero_density_TF[inode_box][icomp]){
            Vext[loc_inode][icomp] = 0.0;
            if (Lvext_dash){
            for (iwall=0; iwall<Nwall; iwall++) {
                iunk = iwall*Ncomp + icomp;
                for (idim=0; idim<Ndim; idim++) 
                    Vext_dash[loc_inode][iunk][idim] = 0.0;
             }  /* end of loop over walls*/
             }
          }
          else if (Vext[loc_inode][icomp] != VEXT_MAX){
             Zero_density_TF[inode_box][icomp]=FALSE;
          }
      }      /* end of loop over local nodes */
   }         /* end of loop over components */
}
/******************************************************************************/
/* setup_vext_max:  Set the external field to the maximum value
in the first 1/2 sigma from the wall for each component.                      */
void setup_vext_max()
{
  int ilist, loc_inode, icomp, iwall, iunk, idim ,inode_box;
  int inode;
  double xpos[3];

  for (icomp=0; icomp<Ncomp; icomp++){ 
     if (Nlists_HW == 1 || Nlists_HW == 2) ilist = 0;
     else                                  ilist = icomp;

     for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {
         inode_box=L2B_node[loc_inode];
         inode=L2G_node[loc_inode];
         node_to_position(inode,xpos);
         if (Zero_density_TF[inode_box][icomp]) {
            Vext[loc_inode][icomp] = Vext_set[loc_inode][icomp];
            if (Lvext_dash){
            for (iwall=0; iwall<Nwall; iwall++) {
                iunk = iwall*Ncomp + icomp;
                for (idim=0; idim<Ndim; idim++) Vext_dash[loc_inode][iunk][idim] = 0.0;
            }  /* end of loop over iwall */
            }
         }
         if (Grafted[Icomp_to_polID[icomp]]!=FALSE && icomp==Grafted_TypeID[Icomp_to_polID[icomp]] && (
             Nodes_2_boundary_wall[ilist][inode_box]<0 || 
             WallType[Nodes_2_boundary_wall[ilist][inode_box]] != Graft_wall[Icomp_to_polID[icomp]])){
             Vext[loc_inode][icomp] = Vext_set[loc_inode][icomp];
             Zero_density_TF[inode_box][icomp] = TRUE;
             if (Lvext_dash){
             for (iwall=0; iwall<Nwall; iwall++) {
                iunk = iwall*Ncomp + icomp;
                for (idim=0; idim<Ndim; idim++) Vext_dash[loc_inode][iunk][idim] = 0.0;
             }   /*end of loop over iwall */
             }
         }
     }        /* end of loop over nodes */
  }           /* end of loop over components */
}
/******************************************************************************/
/*setup_semiperm: In this routine we set up the set point energies
               for local nodes that are inside surfaces ( may or may
               not be semipermeable.  If not Vext_set=VEXT_MAX).*/
void setup_semiperm(int **nelems_w_per_w, int ***elems_w_per_w)
{
  int iwall,ilist,iel_box,wall_el,inode_box,ijk_box[3],ijk_box_0[3];
  int loc_inode,icomp,i,idim;
  int flag,ijk[3];
  
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) 
     for (icomp=0; icomp<Ncomp; icomp++) Vext_set[loc_inode][icomp]=VEXT_MAX;

  for (iwall=0; iwall<Nwall; iwall++)
    for (icomp=0; icomp<Ncomp; icomp++) {
      if (Lsemiperm[WallType[iwall]][icomp]==TRUE){
        if (Nlists_HW==1 || Nlists_HW==2) ilist = 0;
        else                             ilist = icomp;
	
        for(wall_el=0; wall_el<nelems_w_per_w[ilist][iwall]; wall_el++){
	  iel_box = elems_w_per_w[ilist][wall_el][iwall];
	  inode_box = element_box_to_node_box(iel_box);
	  node_box_to_ijk_box(inode_box,ijk_box_0);
	  
	  /*loop through all nodes on this element starting with LLB corner*/
	  for (i=0; i<Nnodes_per_el_V; i++){
	    for (idim=0; idim<Ndim; idim++) ijk_box[idim]=ijk_box_0[idim];
	    if (i>3) ijk_box[2]=ijk_box_0[2]+1;
	    if (i==2 || i==3 || i==6 || i==7) ijk_box[1]=ijk_box_0[1]+1;
	    if (i==1 || i==3 || i==5 || i==7) ijk_box[0]=ijk_box_0[0]+1;
            flag=TRUE;
	    for (idim=0; idim<Ndim; idim++) {
                  ijk[idim]=ijk_box[idim]+Min_IJK_box[idim];
                  if (ijk[idim] >= Nodes_x[idim]) flag=FALSE;
            }
            if (flag){
	    inode_box = ijk_box_to_node_box(ijk_box);
	    loc_inode=B2L_node[inode_box];
			if(loc_inode >=0)
					Vext_set[loc_inode][icomp] = Vext_membrane[WallType[iwall]][icomp];
            }
	  }
        }
      }
    }
   return;
}
/******************************************************************************/
/* setup_vext_XRSurf:  this option works for cylindrical or spherical surfaces
                        where a 1D potential, Vext(r) is set up such that r=0 is
                        taken at the surface of the 'wall' rather than at its center. */
void setup_vext_XRSurf(int iwall)
{
  int iwall_type,icomp,loc_inode,ijk_box[3],ijk_box_tmp[3],iunk, inode_box_tmp, idim;
  int flag_on_cutoff,logical_setup_vext;
  double xORr, fluid_pos[3], sign, vtmpUP, vtmpDOWN, xORr_tmp;
  double param1,param2,param3,param4,param5,param6;

  iwall_type=WallType_Images[iwall];
    for (icomp=0; icomp<Ncomp;icomp++){
       logical_setup_vext=TRUE;
       if (Grafted_Logical==TRUE && Type_poly==WJDC3 &&
           Grafted[Icomp_to_polID[icomp]]!=FALSE && icomp==Grafted_TypeID[Icomp_to_polID[icomp]] && 
           WallType[RealWall_Images[iwall]] == Graft_wall[Icomp_to_polID[icomp]]){ logical_setup_vext=FALSE;
       }
     
       if (logical_setup_vext==TRUE){ 
       if (Type_vext[iwall_type]==VEXT_PAIR_POTENTIAL) pairPotparams_switch(Vext_PotentialID[iwall_type],WALL_FLUID, icomp,iwall,
                                                       &param1,&param2,&param3,&param4,&param5,&param6);
       for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {
           flag_on_cutoff=FALSE;
           if (Zero_density_TF[L2B_node[loc_inode]][icomp]==FALSE) {

              node_box_to_ijk_box(L2B_node[loc_inode],ijk_box);
              node_to_position(L2G_node[loc_inode],fluid_pos);

/*            Code for specific potential - need to think about generalizing for semipermeable membranes. */   
/*              if (Lsemiperm[iwall_type][icomp] && 
                  r<1.1224620*Sigma_wf[icomp][iwall_type]){
                  Vext_set[loc_inode][icomp] = Vext_membrane[WallType[iwall]][icomp];
              }*/

              xORr=X_wall[L2B_node[loc_inode]][iwall];
              if (fabs(Vext[loc_inode][icomp]-Vext_set[loc_inode][icomp])>1.e-6){
                 if (xORr>0.001){

                     if (Type_vext[iwall_type]==VEXT_PAIR_POTENTIAL) {
                          Vext[loc_inode][icomp] += pairPot_switch(xORr,param1,param2,param3,param4,param5,param6,Vext_PotentialID[iwall_type]);
                     }
                     else Vext[loc_inode][icomp]+= Vext_1D(xORr,icomp,iwall_type);
                     if (fabs(xORr-Cut_wf[icomp][iwall_type])<1.e-6) flag_on_cutoff=TRUE;
                 }
                 else Vext[loc_inode][icomp] = Vext_set[loc_inode][icomp]; 
              } 
              if (Vext[loc_inode][icomp] >= Vext_set[loc_inode][icomp]) {
                    Vext[loc_inode][icomp] = Vext_set[loc_inode][icomp]; 
              }
        
              for (idim=0;idim<Ndim; idim++) ijk_box_tmp[idim]=ijk_box[idim];
    
              for (idim=0;idim<Ndim; idim++){
                 xORr_tmp=Xwall_delUP[L2B_node[loc_inode]][iwall][idim];
                 if (xORr>0.001){
                     if (Type_vext[iwall_type]==VEXT_PAIR_POTENTIAL) {
                          vtmpUP = pairPot_switch(xORr_tmp,param1,param2,param3,param4,param5,param6,Vext_PotentialID[iwall_type]);
                      }
                      else vtmpUP = Vext_1D(xORr_tmp,icomp,iwall_type);
                 }
                 else            vtmpUP=Vext_set[loc_inode][icomp];

                 xORr_tmp=Xwall_delDOWN[L2B_node[loc_inode]][iwall][idim];
                 if (xORr>0.001){ 
                     if (Type_vext[iwall_type]==VEXT_PAIR_POTENTIAL) {
                          vtmpDOWN = pairPot_switch(xORr_tmp,param1,param2,param3,param4,param5,param6,Vext_PotentialID[iwall_type]);
                      }
                      else {vtmpDOWN = Vext_1D(xORr_tmp,icomp,iwall_type);}
                 }
                 else vtmpDOWN=Vext_set[loc_inode][icomp];
 
                 iunk=RealWall_Images[iwall]*Ncomp+icomp;
                 if (WallPos_Images[iwall][idim] > fluid_pos[idim]) sign=-1.0;
                 else sign=1.0;

                 if (Vext[loc_inode][icomp]==Vext_set[loc_inode][icomp]){
                       Vext_dash[loc_inode][iunk][idim]=0.0;
                 }
                 else{
                     if ((fabs(vtmpUP-Vext_set[loc_inode][icomp])<1.e-6 && fabs(vtmpDOWN-Vext_set[loc_inode][icomp])>1.e-6) ||
                         (fabs(vtmpUP)<1.e-6 && flag_on_cutoff==TRUE)){
                        Vext_dash[loc_inode][iunk][idim]+=sign*(Vext[loc_inode][icomp]-vtmpDOWN)/VDASH_DELTA;
                     }
                     else if ((fabs(vtmpDOWN-Vext_set[loc_inode][icomp])<1.e-6 && fabs(vtmpUP-Vext_set[loc_inode][icomp])>1.e-6) ||
                         (fabs(vtmpDOWN<1.e-6) && flag_on_cutoff==TRUE)){
                        Vext_dash[loc_inode][iunk][idim]+=sign*(vtmpUP-Vext[loc_inode][icomp])/VDASH_DELTA;
                     }
                     else{
                        Vext_dash[loc_inode][iunk][idim]+=sign*(vtmpUP-vtmpDOWN)/(2.0*VDASH_DELTA);
                     }
                 }
              }

       }  /*zero check */
    } /*loc_inode*/
    }
   } /*icomp */
   return;
}

/***************************************************************************/
/* setup_integrated_LJ_walls:  In this routine we assume that materials properties
                            are constant in the surfaces of interest, and that
                            surface-fluid interactions are described by 12-6
                            Lennard-Jones potentials.*/
void setup_integrated_LJ_walls(int iwall, int *nelems_w_per_w,int **elems_w_per_w)
{
   int icomp,iwall_type,ilist,idim, i,inode_box,ijk[3];
   int loc_inode,inode_llb,inode,iel,wall_el;
   int max_els,image,image_x,image_xy,iel_y,iel_z;
   int ngp,ngpu,ngp1,ngpu1,ngp2,ngpu2,ngp3,ngpu3;
   double *gp,*gw,*gpu,*gwu;
   double gp1[12],gw1[12],gpu1[40],gwu1[40];
   double gp2[12],gw2[12],gpu2[12],gwu2[12];
   double gp3[12],gw3[12],gpu3[12],gwu3[12];
   double max_cut,**image_pos,node_pos[3],node_pos_w[3],
          node_pos_f[3],node_pos_w2[3],vext,r_center_sq;
   double param1,param2,param3,param4,param5,param6;
   int logical_setup_vext;

   iwall_type = WallType_Images[iwall];
   max_cut = 0.0;
   for (icomp=0; icomp<Ncomp; icomp++)
         if (max_cut < Cut_wf[icomp][iwall_type]) 
                 max_cut = Cut_wf[icomp][iwall_type];
  
   max_els = (int)(POW_DOUBLE_INT(max_cut,Ndim)/Vol_el) + 1;
   image_pos = (double **) array_alloc (2, max_els, Ndim, sizeof(double));

   ilist = Nlists_HW-1; /* only add the contributions of the solid*/

  /* put 20 or fewer gauss quadrature points in each element */
   gp = (double *) array_alloc (1, 20, sizeof(double));
   gw = (double *) array_alloc (1, 20, sizeof(double));
   gpu = (double *) array_alloc (1, 20, sizeof(double));
   gwu = (double *) array_alloc (1, 20, sizeof(double));

  ngp1  = 6;  ngp2 = 3; ngp3  = 3;
  ngpu1 = 20; ngpu2 = 12; ngpu3 = 6;

  set_gauss_quad(ngp1, gp1, gw1); 
  set_gauss_quad(ngp2, gp2, gw2); 
  set_gauss_quad(ngp3, gp3, gw3); 

  set_gauss_quad(ngpu1, gpu1, gwu1);
  set_gauss_quad(ngpu2, gpu2, gwu2);
  set_gauss_quad(ngpu3, gpu3, gwu3);

    for (wall_el=0; wall_el< nelems_w_per_w[ilist]; wall_el++){
      iel = elems_w_per_w[ilist][wall_el];
      inode_llb = element_to_node(iel);

      image= 0;
      node_to_position (inode_llb,node_pos);
      for (idim=0; idim<Ndim; idim++) {
        node_pos_w[idim] = node_pos[idim] + 0.5*Esize_x[idim];
        image_pos[image][idim] = node_pos_w[idim];
      }

      for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {

        inode_box = L2B_node[loc_inode];
        node_to_ijk(B2G_node[inode_box],ijk);

        for (icomp=0; icomp<Ncomp; icomp++) {
          logical_setup_vext=TRUE;
          if (Grafted_Logical==TRUE && Type_poly==WJDC3 &&
              Grafted[Icomp_to_polID[icomp]]!=FALSE && 
              icomp==Grafted_TypeID[Icomp_to_polID[icomp]] && 
              WallType[RealWall_Images[iwall]] == Graft_wall[Icomp_to_polID[icomp]]){ 

              logical_setup_vext=FALSE;
           }
     
       if (!Zero_density_TF[inode_box][icomp] && logical_setup_vext==TRUE) {

           image = 1;
           inode = L2G_node[loc_inode];
           node_to_position(inode,node_pos_f);

         if (Vext[loc_inode][icomp] < Vext_set[loc_inode][icomp] ||
             Vext[loc_inode][icomp]==0.0) {

        find_images(0,Cut_wf[icomp][iwall_type],&image,image_pos,
                       node_pos_w,node_pos_f);


           if (Ndim > 1) {
             image_x = image;
             for (iel_y=0; iel_y<image_x; iel_y++){
                for (idim=0; idim<Ndim; idim++) 
                    node_pos_w2[idim] = image_pos[iel_y][idim];
                find_images(1,Cut_wf[icomp][iwall_type],&image,image_pos,
                           node_pos_w2,node_pos_f);
             }
          }

          if (Ndim == 3) {
             image_xy = image;
             for (iel_z=0; iel_z<image_xy; iel_z++){
                for (idim=0; idim<Ndim; idim++) 
                    node_pos_w2[idim] = image_pos[iel_z][idim];
                find_images(2,Cut_wf[icomp][iwall_type],&image,image_pos,
                           node_pos_w2,node_pos_f);
             }
          }

          /* 
           *  now that we have all the images associated with
           *  this particular surface element for icomp...calculate
           *     the external field !!
           */
          
          for (i=0; i<image; i++){
             r_center_sq = 0.0;
             for (idim=0; idim<Ndim; idim++) {
                 /* adjust to llb node of element because gauss points are on interval [0,1] */
                 node_pos_w2[idim] = image_pos[i][idim]-0.5*Esize_x[idim];
                 r_center_sq += (node_pos_w2[idim]-node_pos_f[idim])*
                                (node_pos_w2[idim]-node_pos_f[idim]);
             }
             if (Lsemiperm[iwall_type][icomp] && r_center_sq<Sigma_wf[icomp][iwall_type]){
                 Vext_set[loc_inode][icomp] = Vext_membrane[WallType[iwall]][icomp];
             }
             if (r_center_sq < 4.0) {
                ngp = ngp1; ngpu=ngpu1; 
                gp = &gp1[0]; gpu = &gpu1[0]; 
                gw = &gw1[0]; gwu = &gwu1[0];
             }
             else if (r_center_sq < 16.0){
                ngp = ngp2; ngpu=ngpu2; 
                gp = &gp2[0]; gpu = &gpu2[0]; 
                gw = &gw2[0]; gwu = &gwu2[0];
             }
             else{
                ngp = ngp3; ngpu=ngpu3; 
                gp = &gp3[0]; gpu = &gpu3[0]; 
                gw = &gw3[0]; gwu = &gwu3[0];
             }


             pairPotparams_switch(Vext_PotentialID[iwall_type],WALL_FLUID,icomp,RealWall_Images[iwall],
                                  &param1,&param2,&param3,&param4,&param5,&param6);
             vext = integrate_potential(Vext_PotentialID[iwall_type],param1,param2,param3,param4,param5,param6,
                       ngp, ngpu, gp, gpu, gw, gwu, node_pos_w2, node_pos_f);

	     
             Vext[loc_inode][icomp] += vext;

          }  /* end of loop over the surface elements and its images */
             if (Vext[loc_inode][icomp] >= Vext_set[loc_inode][icomp]) 
                   Vext[loc_inode][icomp] = Vext_set[loc_inode][icomp];
         }


         }
       }   /* end of icomp loop */
      
      }       /* end of fluid node loop */
    }         /* end of wall element loop */
  safe_free((void *) &image_pos);
  return;
}
/*******************************************************************************/
/*comm_wall_els: In this routine we collect all the wall elements
               on to every processor for calculation of the 
               external fields. */
void comm_wall_els(int iwall,int **nelems_w_per_w, int ***elems_w_per_w,
                   int *nelems_w_per_w_global, int **elems_w_per_w_global)
{
  int ilist,iel_box,iel,wall_el,inode_box;
  int n_loc,n_tot;
  int *wall_elem_tmp;

  /* set up for broadcasting the wall elements */
  /* NOTE:  elems_w_per_w_global has different order of arguments
     then elems_w_per_w  */


     for (ilist=0; ilist<Nlists_HW; ilist++){

        wall_elem_tmp = (int *) array_alloc (1,nelems_w_per_w[ilist][iwall],sizeof(int));

        n_loc = 0;
        for(wall_el=0; wall_el<nelems_w_per_w[ilist][iwall]; wall_el++){
            iel_box = elems_w_per_w[ilist][wall_el][iwall];
            inode_box = element_box_to_node_box(iel_box);
            if (B2L_node[inode_box] != -1){
               iel = el_box_to_el(iel_box);
               wall_elem_tmp[n_loc++] = iel;
            }
        }
       
        /*
         * get the total number of wall elements 
         * on all the processors. 
         */
        n_tot = gsum_int(n_loc);

        /*
         * allocate arrays that contain global element numbers
         * for all wall elements in the domain.
         */
        elems_w_per_w_global[ilist] = (int *) array_alloc (1,n_tot,sizeof(int));
        nelems_w_per_w_global[ilist] = n_tot;


        comm_loc_to_glob_vec(&n_loc, wall_elem_tmp, 
                              elems_w_per_w_global[ilist]);
        safe_free((void *) &wall_elem_tmp);

     }
return;
}
/*******************************************************************************/
/* correct_zeroTF_array:  This is needed because the Vext arrays
   are based on local coordinates and the Zero_TF are based on box
   coordinates.  So, we need to gather the local Zero_TF
   values into a global array, and then redistribute the correct
   values into the original box based array. */
void correct_zeroTF_array()
{
  int *index,*unk_loc,*unk,*unk_global;
  int inode,inode_box,loc_inode,icomp;

  for (icomp=0; icomp<Ncomp; icomp++) {
 
     
    unk = (int *) array_alloc (1, Nnodes, sizeof(int));
    if (Proc==0){
         index = (int *) array_alloc (1, Nnodes, sizeof(int));
         unk_global = (int *) array_alloc (1, Nnodes, sizeof(int));
    }

    /* collect the global indices from all processors */
     MPI_Gatherv(L2G_node,Nnodes_per_proc,MPI_INT,
              index,Comm_node_proc,Comm_offset_node,
              MPI_INT,0,MPI_COMM_WORLD);

     unk_loc = (int *) array_alloc (1, Nnodes_per_proc, sizeof(int));
     for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ ){
         inode_box = L2B_node[loc_inode];
         unk_loc[loc_inode] = Zero_density_TF[inode_box][icomp];
     }

     /* collect the unknowns from all the processors */

     MPI_Gatherv(unk_loc,Nnodes_per_proc,MPI_INT,
              unk_global,Comm_node_proc,Comm_offset_node,
              MPI_INT,0,MPI_COMM_WORLD);
     safe_free((void *) &unk_loc);

     if (Proc == 0){
        for (inode=0; inode<Nnodes; inode++){
            unk[index[inode]] = unk_global[inode];
        }
        safe_free((void *) &unk_global);
        safe_free((void *) &index);
     }
     MPI_Bcast(unk,Nnodes,MPI_INT,0,MPI_COMM_WORLD);

     for (inode_box=0; inode_box<Nnodes_box; inode_box++)
        Zero_density_TF[inode_box][icomp]=unk[B2G_node[inode_box]];

     safe_free((void *) &unk);
  }
}
/*******************************************************************************/
/*comm_vext_max: In this routine we collect all the nodes where 
               Vext = VEXT_MAX and Zero_density_TF = FALSE */
void comm_vext_max(int *nnodes_vext_max, int **nodes_vext_max)
{
  int inode_box,loc_inode,icomp,inode,ijk[3];
  int n_loc,n_tot;
  int *vext_max_nodes_tmp;

  /* set up for broadcasting the wall elements */
  /* NOTE:  elems_w_per_w_global has different order of arguments
     then elems_w_per_w  */


  for (icomp=0; icomp<Ncomp; icomp++){

     vext_max_nodes_tmp = (int *) array_alloc (1,Nnodes_per_proc,sizeof(int));

     n_loc = 0;
     for(loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
         inode_box = L2B_node[loc_inode];
         if ( Vext[loc_inode][icomp]==VEXT_MAX  &&
             !Zero_density_TF[inode_box][icomp]) {
             inode = L2G_node[loc_inode];
             vext_max_nodes_tmp[n_loc++] = inode;
             node_to_ijk(inode,ijk);
         }
     }
       
     /*
      * get the total number of wall elements 
      * on all the processors. 
      */
     n_tot = gsum_int(n_loc);
 
     /*
      * allocate arrays that contain global element numbers
      * for all wall elements in the domain.
      */

     if (n_tot > 0) {
     nodes_vext_max[icomp] = (int *) array_alloc (1,n_tot,sizeof(int));

     nnodes_vext_max[icomp] = n_tot;
     comm_loc_to_glob_vec(&n_loc, vext_max_nodes_tmp, 
                               nodes_vext_max[icomp]);
     }
     else{
        nodes_vext_max[icomp] = (int *) array_alloc (1,1,sizeof(int));
        nnodes_vext_max[icomp] = n_tot;
     }
     safe_free((void *) &vext_max_nodes_tmp);

  }
  return;
}
/***********************************************************************/
/*comm_loc_to_glob_vec: core communication routine for gathering up a
          global vector that has different lengths on each
          processor. */
void comm_loc_to_glob_vec(int *n_loc, int *in_loc_vec, 
                                    int *out_glob_vec)
{
  int ierror,iproc;
  int *range, *displs;
    

  /* allocate array that tells how long each incoming unit of 
     information will be and fill this array with n_loc from
     each processor */

  range = (int *) array_alloc (1,Num_Proc,sizeof(int));

  ierror = MPI_Allgather(n_loc,1,MPI_INT,range,1,
                                 MPI_INT,MPI_COMM_WORLD);
  /* gather up a global array */

  displs = (int *) array_alloc (1,Num_Proc,sizeof(int));
  displs[0] = 0;
  for (iproc=1; iproc<Num_Proc; iproc++){
     displs[iproc] = displs[iproc-1] + range[iproc-1];
  }

  ierror = MPI_Allgatherv(in_loc_vec,*n_loc,MPI_INT,
                          out_glob_vec, range,
                          displs,MPI_INT,MPI_COMM_WORLD);
  safe_free((void *) &range);
  safe_free((void *) &displs);
  return;
}
/***************************************************************/
/* read_external_field_n ..... read in the external field from
  the file dft_vext.dat */
void read_external_field_n()
{
   int idim, i,ijk[3],icomp,inode,loc_inode,c,inode_box;
   double **vext_tmp,**vext_static_tmp,pos_old0[3],pos_old,vtmp;
   FILE *fp,*fp2;

   Vext = (double **) array_alloc (2,Nnodes_per_proc,Ncomp, sizeof(double));
   vext_tmp = (double **) array_alloc (2,Nnodes,Ncomp, sizeof(double));

   if (Restart_Vext == READ_VEXT_STATIC){
       Vext_static = (double **) array_alloc (2,Nnodes_per_proc,Ncomp, sizeof(double));
       vext_static_tmp = (double **) array_alloc (2,Nnodes,Ncomp, sizeof(double));
   }

   /* first read in the external field components found in the default file */
 
   if (Proc==0){
      if (Iwrite==VERBOSE) printf("setting up the external field from the file %s\n",vext_file_array);
      if( (fp=fopen(vext_file_array,"r"))==NULL){
          printf("the file %s does not exist\n",vext_file_array);
      }

      for (i=0; i<Nnodes; i++) {
         for (idim=0;idim<Ndim;idim++)         {
             fscanf(fp,"%lf",&pos_old);
             if (i==0) pos_old0[idim]=pos_old;
             pos_old -= pos_old0[idim];
             ijk[idim] = round_to_int(pos_old/Esize_x[idim]);
         }
         inode=ijk_to_node(ijk);

             /* now read in the external field params */
         for (icomp=0; icomp<Ncomp; icomp++) fscanf(fp,"%lf",&vext_tmp[inode][icomp]);
         while ((c=getc(fp)) != EOF && c !='\n') ;
      }
      fclose (fp);
   }

   /* now if necessary read in the external field components found in file 2 */
   if (Restart_Vext == READ_VEXT_SUMTWO || Restart_Vext == READ_VEXT_STATIC){
      if (Proc==0){
         if (Iwrite==VERBOSE) printf("setting up the external field from the file %s\n",vext_file2_array);
         if( (fp2=fopen(vext_file2_array,"r"))==NULL){
             printf("the file %s does not exist\n",vext_file2_array);
         }

         for (i=0; i<Nnodes; i++) {
            for (idim=0;idim<Ndim;idim++)         {
                fscanf(fp2,"%lf",&pos_old);
                if (i==0) pos_old0[idim]=pos_old;
                pos_old -= pos_old0[idim];
                ijk[idim] = round_to_int(pos_old/Esize_x[idim]);
            }
            inode=ijk_to_node(ijk);
   
                /* now read in the external field params */
            for (icomp=0; icomp<Ncomp; icomp++){ 
               fscanf(fp2,"%lf",&vtmp);
               if (vext_tmp[inode][icomp] >= VEXT_MAX || vtmp >= VEXT_MAX){
                  vext_tmp[inode][icomp] = VEXT_MAX;
               }
               else{ vext_tmp[inode][icomp] += vtmp; }
               if (Restart_Vext == READ_VEXT_STATIC) vext_static_tmp[inode][icomp] = vtmp;
            }
            while ((c=getc(fp2)) != EOF && c !='\n') ;
         }
         fclose (fp2);
      }
    }

    MPI_Bcast(*vext_tmp,Nnodes*Ncomp,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {
        inode=L2G_node[loc_inode];
        for (icomp=0;icomp<Ncomp;icomp++){
             Vext[loc_inode][icomp]=vext_tmp[inode][icomp];
                 /* now if we need to keep track of a static part, fill up the Vext_static array */
             if (Restart_Vext == READ_VEXT_STATIC) Vext_static[loc_inode][icomp]=vext_static_tmp[inode][icomp];
        }
    }

    /* finally take care of Zero_density_TF */
    for (inode_box=0; inode_box<Nnodes_box; inode_box++) {
       inode=B2G_node[inode_box];
       for (icomp=0;icomp<Ncomp;icomp++){
          if (vext_tmp[inode][icomp] >= VEXT_MAX) Zero_density_TF[inode_box][icomp]=TRUE;
          else                                    Zero_density_TF[inode_box][icomp]=FALSE;
       }
    }
    safe_free((void *) &vext_tmp);
    if (Restart_Vext == READ_VEXT_STATIC) safe_free((void *) &vext_static_tmp);
    return;
}
/***************************************************************/
