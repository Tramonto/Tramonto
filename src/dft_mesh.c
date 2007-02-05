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
 *  FILE: dft_mesh.c
 *
 *  This routine does all the mesh setup chores that are independent of 
 *  the surface geometry.  The primary output are the number of nodes 
 *  and elements (total, wall, fluid, and boundary).  
 *
 *  Finally, some aztec parameters are set up here.
 *  
 */

#include "dft_mesh.h"

/********************** BEGIN EXECUTION ************************************/
void set_up_mesh (char *output_file1,char *output_file2)
{
 /* Local variable declarations */
  char *yo = "set_up_mesh";
  double t1;
  FILE *fp1=NULL;
  int i,inode,ijk[3],flag,idim,coarse_fac,count,print_flag;
  int N_update=0; /*local variables to replace AztecStruct global ones of same name*/
  int *update=NULL;

  if (Proc==0 && Iwrite != NO_SCREEN){
       printf("\n-------------------------------------------------------------------\n");
       printf("%s: Setting up the mesh ... \n",yo);
  }
  t1 = MPI_Wtime();
  if (Proc==0) {
    if( (fp1 = fopen(output_file1,"a+")) == NULL) {
      printf("Can't open file %s\n", output_file1);
      exit(1);
    }
  }
  
  /*
   * set up all the basic domain parameters for this run:
   * Nnodes,Nelements,Nodes_x[],Elements_x[],Nodes_plane,
   * Elements_plane, Vol_el, Area_surf_el,Nunk_per_node,
   * Nunknowns. 
   */
  setup_basic_domain(fp1);

  /* 
   * Initialize Aztec settings. We may read in 
   * some of these eventually 
   */
  initialize_Aztec(&N_update, &update);
  Nnodes_per_proc = N_update / Nunk_per_node;

  /* 
   * Do load balancing assuming equal weights
   * on all the nodes in the domain. The load balance
   * routine resets the value of Nodes_per_proc.
   */
  load_balance(0, NULL, &N_update, &update);

/*set up arrays on processor 0 that are needed for communications
  involved with printing arrays */
  Comm_node_proc = (int *) array_alloc (1, Num_Proc, sizeof(int));
  Comm_unk_proc = (int *) array_alloc (1, Num_Proc, sizeof(int));
  Comm_offset_node = (int *) array_alloc (1, Num_Proc, sizeof(int));
  Comm_offset_unk = (int *) array_alloc (1, Num_Proc, sizeof(int));
  for (i=0; i<Num_Proc; i++)
     Comm_node_proc[i]=Comm_unk_proc[i]=Comm_offset_node[i]=Comm_offset_unk[i] = 0;

  MPI_Gather(&Nnodes_per_proc,1,MPI_INT,
             Comm_node_proc,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc == 0){
     for (i=0; i<Num_Proc; i++)
       Comm_unk_proc[i] = Comm_node_proc[i]*Nunk_per_node;

     Comm_offset_node[0] = 0; Comm_offset_unk[0] = 0;
     for (i=1; i<Num_Proc; i++){
        Comm_offset_node[i] = Comm_offset_node[i-1] + Comm_node_proc[i-1];
        Comm_offset_unk[i]  = Comm_offset_unk[i-1]  + Comm_unk_proc[i-1];
     }
  }

  /*
   * Set up all mesh variables based on the initial
   * load balance.
   */

  print_flag=TRUE;
  if (L1D_bc || Load_Bal_Flag == LB_WEIGHTS || (Load_Bal_Flag >= LB_TIMINGS
				      && Mesh_coarsening != FALSE) ) print_flag=FALSE;

  control_mesh(fp1,output_file2,print_flag, update);

  if (L1D_bc || Load_Bal_Flag == LB_WEIGHTS || (Load_Bal_Flag >= LB_TIMINGS
				      && Mesh_coarsening != FALSE) ) {
    print_flag=TRUE; 
     /*
      * Now redo load balancing based on mesh info
      * and then free all the old global mesh arrays.
      */
    load_balance(1,NULL, &N_update, &update);

    if (Iwrite==VERBOSE) printf("Proc: %d Nodes_per_proc: %d\n",Proc,Nnodes_per_proc);

    free_mesh_arrays();

      MPI_Gather(&Nnodes_per_proc,1,MPI_INT,
             Comm_node_proc,1,MPI_INT,0,MPI_COMM_WORLD);

      if (Proc == 0){
         for (i=0; i<Num_Proc; i++)
           Comm_unk_proc[i] = Comm_node_proc[i]*Nunk_per_node;

         Comm_offset_node[0] = 0; Comm_offset_unk[0] = 0;
         for (i=1; i<Num_Proc; i++){
            Comm_offset_node[i] = Comm_offset_node[i-1] + Comm_node_proc[i-1];
            Comm_offset_unk[i]  = Comm_offset_unk[i-1]  + Comm_unk_proc[i-1];
         }
      }

     /*
      * And redo the mesh variables based on the new
      * load balancing.
      */
     if (Proc==0 && Iwrite == VERBOSE) printf("\n%s: Setting up the mesh again after load balance... \n",yo);
     control_mesh(fp1,output_file2,print_flag, update);
  }
    /* now reset the Mesh_coarsen_flag on the one row of nodes where
       we will do the computation beyond the 1D boundary */
    if (L1D_bc){
       count=0;
       for (i=0; i < Nnodes_box; i++) {
         inode = B2G_node[i];
         node_to_ijk(inode, ijk);
         flag=TRUE;
         for (idim=0; idim<Ndim; idim++) 
            if(Mesh_coarsen_flag[i] != FLAG_1DBC || (ijk[idim] != 0 && idim != Grad_dim)) flag=FALSE;
         if (flag){
              count++;
              Mesh_coarsen_flag[i] = Nodes_to_zone[i];
              if (Mesh_coarsening != FALSE && Nodes_to_zone[i] > 0){
                      /* reset to negative flag if residual is not even to be set */
                 coarse_fac = POW_INT(2,Nodes_to_zone[i]);
                 if      (ijk[0]%coarse_fac) Mesh_coarsen_flag[i] = -1;
                 else if (ijk[1]%coarse_fac) Mesh_coarsen_flag[i] = -2;
                 else if (ijk[2]%coarse_fac) Mesh_coarsen_flag[i] = -3;
              }
          }
       }
    }

  safe_free((void *) &update);
  
  if (Proc==0) fclose(fp1);
  if (Proc==0 && Iwrite !=NO_SCREEN){
       printf("mesh set up took %g secs\n",MPI_Wtime()-t1);
       printf("-------------------------------------------------------------------\n");
  }

  return;
}
/********************************************************/


/*************************************************
free_mesh_arrays:  here just free all the
      arrays that are set up in the mesh routine. */
void free_mesh_arrays(void)
{
   int flag,i;
   /*mesh arrays*/
   safe_free((void *) &Nodes_2_boundary_wall);
   safe_free((void *) &Wall_elems);
   safe_free((void *) &Nodes_to_zone);
   safe_free((void *) &Nwall_owners);
   safe_free((void *) &Wall_owners);

   /*if (Mesh_coarsening != FALSE || L1D_bc)*/ safe_free((void *) &Mesh_coarsen_flag);

   if (Ipot_ff_c == COULOMB) safe_free((void *) &Dielec);

   /* external field arrays */
   safe_free((void *) &Vext);
   if (Lvext_dash && Restart != 4) safe_free((void *) &Vext_dash);
   safe_free((void *) &Zero_density_TF);
   flag=FALSE;
   for (i=0;i<Nwall_type;i++)if (Ipot_wf_n[i]==VEXT_1D_XMIN) flag=TRUE;
   if (flag){
       safe_free((void *) &X_wall);
       safe_free((void *) &X_wall2);
   }
   if (Restart != 4){
       safe_free((void *) &Uww);
       safe_free((void *) &Uww_link);
   }

   /* solution arrays*/
   safe_free((void *) &B2G_unk);
   safe_free((void *) &B2G_node);
   safe_free((void *) &L2B_node);
   safe_free((void *) &B2L_node);
   safe_free((void *) &L2G_node);

   return;
}
/************************************************/



/********************************************************
control_mesh: this routine calls other functions
   involved in setting up the mesh.                     */

void control_mesh(FILE *fp1,char *output_file2,int print_flag, int *update)
{
  int  **elems_f,  
    ***elems_w_per_w;

  int  *nelems_f, 
       **nelems_w_per_w;
  int  ***el_type;

  int  *elem_zones;

  int icomp,ilist,iwall,inode;
  double sigma_test;
  int flag,iel_box,i,nwall_max,idim;
  int inode_box,iel,reflect_flag[3],j,k; 

  int count_zero,count_zero_all,count_coarse_resid,count_coarse_r_all,
      count_coarse_jac,count_coarse_jac_all;
   

  char *output_TF="dft_zeroTF.dat";

  int *comm_icount_proc, *comm_offset_icount;

  int ntot_per_list,ntot_per_list_all_procs,*el_tmp,*elems_w_per_w_proc_0_tmp,ncount;
  
  reflect_flag[0]=reflect_flag[1]=reflect_flag[2] = FALSE;

  /*
   * set up all the basic parameters for local box units:
   * Min_IJK_box[], Max_IJK_box[],
   * Nnodes_box,Nelements_box,Nodes_x_box[],
   * Elements_x_box[],Nodes_plane_box, Elements_plane_box,
   * Nunknowns_box. 
   */
  setup_basic_box(fp1, update);

  /*
   * Set up all arrays for surface geometry, external fields, and
   * surface charges
   */

  if (Nwall != 0) {

     /* figure out how many lists we need to work with - this 
        depends on the presence of discontinuities in the
        density profile. */

     Nlists_HW = 2;
     if (Lhard_surf && Ipot_ff_n != IDEAL_GAS){
        if (Ncomp > 1){
           sigma_test = Sigma_ff[0][0];
           flag = FALSE;
           for (icomp=1; icomp<Ncomp; icomp++){
               if (Sigma_ff[icomp][icomp] != sigma_test) flag = TRUE;
           }
           if (flag == TRUE) Nlists_HW = Ncomp + 1;
        }
     }


     /*
      * SET UP SOME ARRAYS 
      */
 
     elems_f = (int **) array_alloc (2, Nlists_HW, Nelements_box, sizeof(int));

     nelems_f = (int *) array_alloc (1, Nlists_HW, sizeof(int));

     elems_w_per_w =(int ***)array_alloc(3, Nlists_HW,Nelements_box,Nwall,sizeof(int));
     nelems_w_per_w  =(int **) array_alloc(2, Nlists_HW, Nwall, sizeof(int));

     elem_zones = (int *) array_alloc (1, Nelements_box, sizeof(int));
     for (iel_box=0; iel_box<Nelements_box; iel_box++) elem_zones[iel_box] = 0;

     Wall_elems = (int **) array_alloc(2, Nlists_HW, Nelements_box, sizeof(int));
     el_type = (int ***) array_alloc(3, Nwall,Nlists_HW, Nelements_box, sizeof(int));
     if (Ipot_ff_c == COULOMB){
       Dielec = (double *) array_alloc (1, Nelements_box, sizeof(double));
       for (iel=0; iel<Nelements_box; iel++) Dielec[iel]=Dielec_bulk;
     }

     Nwall_owners = (int **) array_alloc (2, Nlists_HW, Nelements_box, sizeof(int));

     nwall_max = Nlink;
     for (idim=0; idim<Ndim; idim++){
       if (Type_bc[idim][0] == PERIODIC || Type_bc[idim][0] == REFLECT ||
           Type_bc[idim][1] == PERIODIC || Type_bc[idim][1] == REFLECT ) nwall_max *= 3;
     }
     Wall_owners = (int ***) array_alloc (3, Nlists_HW, Nelements_box, nwall_max,sizeof(int));
     for (i=0; i<Nlists_HW; i++){
        for (j=0; j<Nelements_box; j++) {
            Nwall_owners[i][j] = 0;
            for (k=0; k<nwall_max; k++) Wall_owners[i][j][k] = 0;
        }
     }



     /*										 
      * Separate wall elements from fluid elements, and assign fluid elements to 
      *  the different quadrature zones.  (surface geometry dependent routines) 
      *  For now, all the processers do this step for all the nodes in the domain !!
      */
     if (Nwall != 0) {
          /* initialize some arrays for use in setting up Nodes_2_boundary_wall array
             later - these arrays are filled in in setup_surface */

         Index_wall_nodes = (int *) array_alloc (1, Nnodes_box, sizeof(int));
         Nodes_wall_box = (int *) array_alloc (1, Nnodes_box, sizeof(int));
         Nwall_touch_node = (int *) array_alloc (1, Nnodes_box, sizeof(int));
         Wall_touch_node = (int **) array_alloc (2, Nnodes_box, Nwall*Nlists_HW, sizeof(int));
         List_wall_node = (int **) array_alloc (2, Nnodes_box, Nwall*Nlists_HW, sizeof(int));
         Nnodes_wall_box=0;
         for (inode_box=0; inode_box<Nnodes_box; inode_box++) Nwall_touch_node[inode_box]=0;
         for (inode_box=0; inode_box<Nnodes_box; inode_box++) Index_wall_nodes[inode_box]=-1;

         setup_surface(fp1,nelems_f, nelems_w_per_w, elems_f,
                                   elems_w_per_w,elem_zones,
                                   el_type);

         /* gather the elems_w_per_w array into a global list of wall elements on processor zero */
         /* first gather and sum the nelems_w_per_w on processor zero -- this gives a maximum possible number of wall elements */
         /* to do --- every processor - loop through list - convert box unit elements to global elements */
         /* pass the global wall element array to processor zero.  For each entry in the passed array - proc 0 must look 
            for matches with other processor and discard duplicates - alternatively, we could allocate a global array, copy
            all versions into this global array and then count the wall element entries. */


/*          for (ilist=0;ilist<Nlists_HW;ilist++){
             Vol_in_surfs[ilist]=0.0;

             ntot_per_list=0;
             for (iwall=0;iwall<Nwall;iwall++){
                 ntot_per_list+=nelems_w_per_w[ilist][iwall];
             }

             ntot_per_list_all_procs=gsum_int(ntot_per_list);
             if (Proc==0){
                  elems_w_per_w_proc_0_tmp  =(int *) array_alloc(1, ntot_per_list_all_procs, sizeof(int));
             }
             el_tmp  =(int *) array_alloc(1, ntot_per_list, sizeof(int));

             j=0;
             for (iwall=0;iwall<Nwall;iwall++) 
                for (i=0;i<nelems_w_per_w[ilist][iwall];i++){
                     el_tmp[j]=el_box_to_el(elems_w_per_w[ilist][i][iwall]);
                     j++;
                }

             comm_icount_proc = (int *) array_alloc (1, Num_Proc, sizeof(int));
             comm_offset_icount = (int *) array_alloc (1, Num_Proc, sizeof(int));
 
             MPI_Gather(&j,1,MPI_INT,comm_icount_proc,1,MPI_INT,0,MPI_COMM_WORLD);

             if (Proc == 0){
                comm_offset_icount[0] = 0;
                for (i=1; i<Num_Proc; i++)
                   comm_offset_icount[i] = comm_offset_icount[i-1] + comm_icount_proc[i-1];
             }   

             MPI_Gatherv(el_tmp,ntot_per_list,MPI_INT,
                         elems_w_per_w_proc_0_tmp,comm_icount_proc,comm_offset_icount,
                         MPI_INT,0,MPI_COMM_WORLD);

             safe_free((void *) &comm_icount_proc);     
             safe_free((void *) &comm_offset_icount);

             if (Proc==0){
                 ncount=ntot_per_list_all_procs; *(start out assuming all entries are unique)*
                 for (i=1;i<ntot_per_list_all_procs;i++){
                    j=0;
                                                          * increment j if no match is found *
                    while (j<i && elems_w_per_w_proc_0_tmp[i]!=elems_w_per_w_proc_0_tmp[j]){j++;} 
                    if (j != i) ncount--;  * deincrement the counter if there was a match *
                 }
                 safe_free((void *) &elems_w_per_w_proc_0_tmp);
             }
             safe_free((void *) &el_tmp);
             if (Proc==0){
                  Vol_in_surfs[ilist]+=(double)ncount*Vol_el;
                  if (Iwrite != NO_SCREEN) printf("total volume in surfaces for ilist=%d is %9.6f (%d elements)\n",
                         ilist,Vol_in_surfs[ilist],ncount);
             }
          }
*/

     }

/*     if (Num_Proc>1) MPI_Barrier(MPI_COMM_WORLD);*/

     if (Imain_loop == 0 && Proc ==0){
        fprintf (fp1,"\n---------------------------------------------------------------\n");
        fprintf (fp1, " Have set up elements for the  selected surface geometry\n");
        fprintf (fp1, " For this problem, Nlists_HW= %d \n",Nlists_HW);
        for (ilist=0; ilist<Nlists_HW; ilist++){
           fprintf (fp1,"ilist: %d\n",ilist);
           fprintf (fp1,"\t nelems_f[ilist]: %d Proc %d\n",nelems_f[ilist],Proc);
           for (iwall=0; iwall<Nwall; iwall++)
                 fprintf (fp1,"\t iwall: %d \t nelems_w_per_w[ilist][iwall]: %d\n",
                                         iwall, nelems_w_per_w[ilist][iwall]);
        }
        if (Proc==0) 
        fprintf (fp1,"---------------------------------------------------------------\n");
        /*fprintf (fp1,"%d Quadrature zones have been assigned\n",Nzone);*/
     }

    /* 
     * From the elem_zones array, assemble an array of Nodes_to_zones. 
     */

     zones_el_to_nodes(elem_zones);

     /* set mesh coarsening flag for residual zones */

     /*if (Mesh_coarsening !=FALSE || L1D_bc)*/ set_mesh_coarsen_flag();

     safe_free((void *) &elem_zones);

    /* 
     * From the Elem arrays, assemble lists of wall, fluid, and Boundary nodes.
     * Each processor sets up the arrays that pertain to the nodes it owns only.
     * The global arrays set up are Nodes_2_boundary_wall and Wall_elems.
     */
     setup_zeroTF_and_Node2bound_new (fp1, el_type);

     /* some temporary list arrays can be trashed now */
     safe_free((void *) &elems_f);
     safe_free((void *) &nelems_f);
     safe_free((void *) &el_type);
     safe_free((void *) &Nodes_wall_box);
     safe_free((void *) &Nwall_touch_node);
     safe_free((void *) &Index_wall_nodes);
     safe_free((void *) &Wall_touch_node);
     safe_free((void *) &List_wall_node);

    /*
     * Set up the neutral part of the external field.
     */
 
     if (Restart == 4) { 
                       read_external_field_n(output_file2);
                       read_zero_density_TF(output_TF);
     }
     else              setup_external_field_n(nelems_w_per_w,elems_w_per_w);
     if (Nwall > 1 && Lprint_pmf) setup_wall_wall_potentials(nelems_w_per_w,elems_w_per_w);


     /* the temporary list arrays can be trashed now */
     safe_free((void *) &elems_w_per_w);
     safe_free((void *) &nelems_w_per_w);

     /* finally print mesh info as we have set up the mesh */
     for (icomp=0; icomp<Ncomp; icomp++){
       count_zero=0;
       count_coarse_resid=0;
       count_coarse_jac=0;

        for (i=0; i<Nnodes_per_proc; i++){
     
            if (Zero_density_TF[L2B_node[i]][icomp]) count_zero++;
            if (Mesh_coarsening !=FALSE || L1D_bc) if (Mesh_coarsen_flag[L2B_node[i]] < 0) count_coarse_resid++;
            if (Nodes_to_zone[L2B_node[i]] > 0 || Coarser_jac>0) count_coarse_jac++;
        }
        count_zero_all=gsum_int(count_zero);
        if (icomp==0){
          count_coarse_r_all=gsum_int(count_coarse_resid);
          count_coarse_jac_all=gsum_int(count_coarse_jac);
          if (Proc==0 && Iwrite != NO_SCREEN && print_flag) {
              printf("**************************************************************\n");
              printf("..............MESH SUMMARY..........\n");     
              printf("Total number of nodes in calculation = \t %d\n",Nnodes);
              printf("Number of coarsened residual nodes = \t %d \t or %g percent\n",
                        count_coarse_r_all,100.*count_coarse_r_all/Nnodes);
              printf("Number of coarsened jacobian nodes = \t %d \t  or %g percent \n",
                        count_coarse_jac_all,100.*count_coarse_jac_all/Nnodes);
           
              printf("--------------------------------------------------------------\n");
          }
        }
        if (Proc==0 && Iwrite != NO_SCREEN && print_flag){
             printf("Number of zero density nodes for icomp: %d = \t %d \t or %g percent \n",
                  icomp,count_zero_all,100.*count_zero_all/Nnodes);
             if (icomp==Ncomp-1) 
              printf("**************************************************************\n");
        }
     }

  }
  else {    /* if zero walls */
     Nlists_HW = 1;
     Nodes_to_zone = (int *) array_alloc (1, Nnodes_box, sizeof(int));
     Wall_elems = (int **) array_alloc(2, Nlists_HW, Nelements_box, sizeof(int));
     Nodes_2_boundary_wall =(int **)array_alloc(2, Nlists_HW, Nnodes_box, sizeof(int));
     Zero_density_TF = (int **) array_alloc (2, Nnodes_box,Ncomp+1,sizeof(int));
     Dielec = (double *) array_alloc (1, Nelements_box, sizeof(double));

     for (inode = 0; inode < Nnodes_box; inode++){
        Nodes_to_zone[inode] = 0;
        Nodes_2_boundary_wall[0][inode] = -1;
        for (icomp=0; icomp<Ncomp+1; icomp++)
            Zero_density_TF[inode][icomp] = FALSE;
     }
     for (iel_box=0; iel_box<Nelements_box; iel_box++) {
           Wall_elems[0][iel_box] = -1;
           if (Ipot_ff_c == COULOMB) Dielec[iel_box] = Dielec_bulk;
     }

      
     setup_external_field_n(NULL,NULL);

  }

  if (Lsteady_state && Ndim==1){
      Area_IC = (double *) array_alloc (1, Nnodes_box, sizeof(double));
      setup_area_IC();
  }

  return;
}
/************************************************************************/
/*setup_basic_domain: Here we set up the basic parameters of the calculation.
                 Number of nodes and elements etc.*/
void setup_basic_domain(FILE *fp1)
{
  int idim,jdim;

  Nnodes = 1;
  Nelements = 1;
  if (Proc==0 && Imain_loop == 0) {
     fprintf(fp1,"\n--------------------------------------------------------------\n");
     fprintf(fp1,"\n idim \t Nodes_x[idim] \t Nnodes \t Elements_x[idim] \t Nelements \n");
  }
  for (idim = 0; idim < Ndim; idim++) {
      if (round_to_int(fmod(Size_x[idim],Esize_x[idim])) == 0){


         Nodes_x[idim] = round_to_int(Size_x[idim]/Esize_x[idim] + 1.);
         Elements_x[idim] = Nodes_x[idim] - 1;
         if (Type_bc[idim][0] == PERIODIC &&
             Type_bc[idim][1] == PERIODIC) Nodes_x[idim]--;

         Nnodes = Nnodes*Nodes_x[idim];
         Nelements = Nelements*Elements_x[idim];
         if (Proc==0 && Imain_loop == 0) fprintf(fp1,"   %d \t    %d \t \t %d \t \t \t %d \t \t %d\n",
                       idim,Nodes_x[idim],Nnodes,Elements_x[idim],Nelements);
      }
      else{
         if (Proc==0 && Imain_loop == 0){
            fprintf(fp1,"\n \nERROR: Esize_x and Size_x are not commensurate\n");
            fprintf(fp1,"Esize_x: %lf\tSize_x: %lf\n",Esize_x[idim],Size_x[idim]);
            fprintf(fp1,"idim: %d \t fmod(Size_x[idim],Esize_x[idim]): %g\n",
                    idim,fmod(Size_x[idim],Esize_x[idim]) );
         }
         exit(-1);
      }
  } 
  for (idim = Ndim; idim < 3; idim++) Nodes_x[idim] = 1;

  if (Ndim == 3) {
      Nodes_plane = Nodes_x[0]*Nodes_x[1];
      Elements_plane = Elements_x[0]*Elements_x[1];
      if (Proc==0 && Imain_loop == 0) fprintf (fp1,"\nNodes_plane: %d \t Elements_plane: %d \n",
                                                Nodes_plane,Elements_plane);
  }

  /* a couple of constants we'll need later - the volume of an
     element - for volume integrals, and the area of a surface
     element - for surface integrals.  In the latter case, the
     surface areas are stored by dimension as our grid can be
     sized differently in the different directions!! */

  Vol_el = 1.0;
  for (idim=0; idim<Ndim; idim++) Vol_el *= Esize_x[idim];

  for (idim=0; idim<Ndim; idim++){
     Area_surf_el[idim]=1.0;
     for (jdim=0; jdim<Ndim; jdim++){
         if (jdim != idim) Area_surf_el[idim] *= Esize_x[jdim];
     }
  }

  /*calculated total unknowns in the problem of interest ! */ 
  Nunknowns = Nnodes * Nunk_per_node;

  return;
}
/******************************************************************
setup_basic_box: here we set up the box coordinates for
                 parallel computations. */
void setup_basic_box(FILE *fp1, int *update)
{
  int idim,inode_box,i_box,loc_inode,inode,icomp;
  double max_cut=0.0;
  int ijk_1D=0;
  FILE *fp2;

  if (Iwrite==VERBOSE) {
    if( (fp2 = fopen("proc_mesh.dat","w+")) == NULL) {
      printf("Can't open file proc_mesh.dat\n");
      exit(1);
    }
  }
  if (Ipot_ff_n > HARD_SPHERE) 
    {
     for (icomp=0; icomp<Ncomp; icomp++)
       {
         if (Cut_ff[icomp][icomp] > max_cut) 
	   {
	     max_cut = Cut_ff[icomp][icomp];
	   }
       }
    }
  else 
    {
      max_cut = 0.5;
    }

  if (L1D_bc) ijk_1D = (int) X_1D_bc/Esize_x[Grad_dim]+1;


  for (idim=0; idim<Ndim; idim++){
      Min_IJK_box[idim] = Min_IJK[idim] - Max_sten_length[idim];
      Max_IJK_box[idim] = Max_IJK[idim] + Max_sten_length[idim];

      if (L1D_bc && idim!= Grad_dim){
            if (Min_IJK[Grad_dim] <= ijk_1D || 
               Max_IJK[Grad_dim] >= Nodes_x[Grad_dim]-ijk_1D) Min_IJK_box[idim] = 0;
      }

     /* set Pflag to catch cases where (1) the cut-offs are so large
        that the box coordinates will require the entire domain; or
        (2) the processor owns the entire range of x[idim].
        In these cases, we can apply the old type of periodic 
        boundary conditions and not extend the box nodes beyond the
        boundaries. */


      Pflag[idim] = FALSE;
      if (Type_bc[idim][0]==PERIODIC && (
          max_cut >0.5*Size_x[idim] ||
          (Min_IJK[idim]==0 && Max_IJK[idim]==Nodes_x[idim]-1))) 
                                               Pflag[idim]=TRUE;

      if ((Type_bc[idim][0] != PERIODIC && Min_IJK_box[idim]<0) ||Pflag[idim]) 
          Min_IJK_box[idim] = 0;
      if ((Type_bc[idim][1] != PERIODIC && Max_IJK_box[idim] >= Nodes_x[idim]) ||Pflag[idim])
          Max_IJK_box[idim] = Nodes_x[idim] - 1;
 
  }

  Nnodes_box = 1;
  Nelements_box = 1;
  if (Num_Proc > 1) MPI_Barrier(MPI_COMM_WORLD);
  if (Imain_loop == 0 && Proc==0) {
    fprintf(fp1,"\n-------------------------------------------------------\n");
    fprintf(fp1,"\n \t idim \t Nodes_x[idim] \t Nnodes \t Elements_x[idim] \t Nelements...box units\n");
  }
  for (idim = 0; idim < Ndim; idim++) {
    Nodes_x_box[idim] = Max_IJK_box[idim] - Min_IJK_box[idim] + 1;
    if (Pflag[idim] && Type_bc[idim][0]==PERIODIC && 
	Max_IJK_box[idim]==Nodes_x[idim]-1) Elements_x_box[idim] = Nodes_x_box[idim];
    else Elements_x_box[idim] = Nodes_x_box[idim] - 1;
    
    Nnodes_box = Nnodes_box*Nodes_x_box[idim];
    Nelements_box = Nelements_box*Elements_x_box[idim];
    if (Imain_loop == 0 && Proc==0) 
      fprintf(fp1,"  Proc: %d   %d \t    %d \t \t %d \t \t \t %d \t \t %d\n",
	      Proc,idim,Nodes_x_box[idim],Nnodes_box,
	      Elements_x_box[idim],Nelements_box);
  } 
  for (idim = Ndim; idim < 3; idim++) Nodes_x_box[idim] = 1;


  if (Ndim == 3) {
    Nodes_plane_box = Nodes_x_box[0]*Nodes_x_box[1];
    Elements_plane_box = Elements_x_box[0]*Elements_x_box[1];
    if (Imain_loop == 0 && Proc==0) 
      fprintf (fp1,"\n Proc: %d  Nodes_plane_box: %d \t Elements_plane_box: %d \n",
	       Proc, Nodes_plane_box,Elements_plane_box);
  }

  Nunknowns_box = Nnodes_box * Nunk_per_node;
  if (Proc == 0) fprintf (fp1," Proc: %d  Nunknowns_box: %d \n",
                                           Proc, Nunknowns_box);

  if (Iwrite==VERBOSE){
     for (idim=0; idim<Ndim; idim++) fprintf(fp2, "Proc=%d Pflag[%d]=%d\n",Proc,idim,Pflag[idim]);
     fprintf(fp2,"Proc: %d  Nnodes_per_proc: %d  Nnodes_box: %d",Proc,Nnodes_per_proc,Nnodes_box); 
     for (idim=0; idim<Ndim; idim++) fprintf(fp2,"Min_IJK[%d]: %d  Max_IJK[%d]:%d ",idim,Min_IJK[idim],idim,Max_IJK[idim]); 
     for (idim=0; idim<Ndim; idim++) 
          fprintf(fp2,"Min_IJK_box[%d]: %d  Max_IJK_box[%d]:%d \n",idim,Min_IJK_box[idim],idim,Max_IJK_box[idim]); 
  }

  /*
   * Set up flags indicating if multiple box nodes map to the same
   * global node, which can happen in periodic BC. Last spot in
   * array is true if any of the others is true.
   */

  Non_unique_G2B[0] = Non_unique_G2B[1] = 
    Non_unique_G2B[2] = Non_unique_G2B[3] = 0;
  for (idim=0; idim<Ndim; idim++) {
   if (Nodes_x_box[idim] > Nodes_x[idim]) {
      Non_unique_G2B[idim] = Nodes_x_box[idim] / Nodes_x[idim];
      Non_unique_G2B[3] = TRUE;
   }
  }


/*  if (Num_Proc>1) MPI_Barrier(MPI_COMM_WORLD);
  if (Proc==0) printf("Proc #: | Min_IJK_box  Max_IJK_box | x Ndim\n");
  if (Proc%(int)(sqrt(Num_Proc)) == 0 || Non_unique_G2B[Ndim]) {
    printf("Proc %d: ",Proc);
    for (idim=0; idim<Ndim; idim++)
      printf("| %4d - %4d |",Min_IJK_box[idim], Max_IJK_box[idim]);
    if (Non_unique_G2B[3]) {
      printf("  Non_unique_G2B:");
      for (idim=0; idim<Ndim; idim++)
        printf("  %d",Non_unique_G2B[idim]);
    }
    printf("\n");
  }
*/

  /* set up B2G_node and B2G_unk mapping from box to global */
  
  B2G_node = (int *) array_alloc(1, Nnodes_box, sizeof(int));
  B2G_unk  = (int *) array_alloc(1, Nunknowns_box, sizeof(int));
  for (inode_box=0; inode_box <Nnodes_box; inode_box++){
     B2G_node[inode_box] = node_box_to_node(inode_box);
  }
  for (i_box=0; i_box <Nunknowns_box; i_box++){
     B2G_unk[i_box] = unk_box_to_unk(i_box);
  }

  /* set up L2B_node from local to box */

  L2B_node = (int *) array_alloc(1, Nnodes_per_proc, sizeof(int));
  L2G_node = (int *) array_alloc(1, Nnodes_per_proc, sizeof(int));
  B2L_node = (int *) array_alloc(1, Nnodes_box, sizeof(int));
  for (inode_box=0; inode_box<Nnodes_box; inode_box++) {
    B2L_node[inode_box] = -1;
  }
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {
    if (MATRIX_FILL_NODAL) {
      inode = update[Nunk_per_node * loc_inode] / Nunk_per_node;
    }
    else inode = update[loc_inode];
    inode_box = node_to_node_box(inode);
    L2B_node[loc_inode] = inode_box;
    B2L_node[inode_box] = loc_inode;
    L2G_node[loc_inode] = inode;
  }
  if (Iwrite==VERBOSE) fclose(fp2);
  return;
}
/*****************************************************************************
boundary_setup: This routine sets up the boundaries and surface charges.
                These calls were removed from mesh_setup to be called
                after load balancing, since some arrays are local to a Proc.*/
void boundary_setup(char *output_file1)
{
  FILE *fp1=NULL;

  if (Proc==0) {
    if( (fp1 = fopen(output_file1,"a+")) == NULL) {
      printf("Can't open file %s\n", output_file1);
      exit(1);
    }
  }
  if (Nwall != 0) {

    /*
     * Determine the local boundary properties.  These include the
     * unit normals, the total areas of each surface, and the surface
     * area in each dimension on each processor. 
     */
     boundary_properties(fp1);
  }

  /*
   * Now that we have number the surface elements, it's time to assign      
   * a local surface charge to each surface element if we are doing charged 
   * surfaces with constant surface charge boundary conditions !!!          
   */
  if(Ipot_wf_c == COULOMB) {
     Vol_charge_flag = FALSE;
     Surf_charge_flag = FALSE;
     setup_surface_charge(fp1);
  }

  if (Proc==0) fclose(fp1);
}
/*****************************************************************************
boundary_free: This routine frees the arrays formed in boundary_setup       */
void boundary_free(void)
{

  if (Nwall != 0) {

    /*
     * free arrays allocated in boundary_properties
     */

      safe_free((void *) &Nelems_S);
      safe_free((void *) &Surf_normal);
      safe_free((void *) &Surf_elem_type);
      safe_free((void *) &Surf_elem_to_wall);
      safe_free((void *) &S_area);
      safe_free((void *) &S_area_tot);


     /*
      * Free arrays in setup_surface_charge
      */

     if(Ipot_wf_c){
          if (Vol_charge_flag) safe_free((void *) &Charge_vol_els);
          if (Surf_charge_flag) safe_free((void *) &Charge_w_sum_els);
     }

  }
}

/*****************************************************************************
setup_zeroTF_and_Node2bound: This routine take the information about which elements
                    are wall and fluid elements, and translates that into
                    lists of wall nodes and fluid nodes for each wall list */
void setup_zeroTF_and_Node2bound (FILE *fp1,int ***el_type)
{
 int reflect_flag[3], **n_fluid_els, iwall,
      iel, idim, ilist, loc_node_el,
      inode, inode_box, iel_box, icomp, iside,
      periodic_flag, wall_flag;
 int **n_el_in_box, **countw;
 int test,all_equal,ijk[3], reflect,count,*els_one_owner,i;
 int dim_reflect[3],ndim_reflect,wall_chk; 
 int link_chk,ilink;
 int iw,jw,nonat_count,*nonatomic_walls;

 Nodes_2_boundary_wall = (int **) array_alloc(2, Nlists_HW, Nnodes_box, sizeof(int));
 Zero_density_TF = (int **) array_alloc (2, Nnodes_box,Ncomp+1,sizeof(int));

 /* zero the arrays that will be set up here */

 for (idim=0; idim<Ndim; idim++) reflect_flag[idim] = FALSE;

 for (inode_box=0; inode_box<Nnodes_box; inode_box++){
    for (ilist=0; ilist<Nlists_HW; ilist++) {
            Nodes_2_boundary_wall[ilist][inode_box] = -1;
    }
    for (icomp=0; icomp<Ncomp; icomp++) {
            Zero_density_TF[inode_box][icomp] = FALSE;
    }
    if (Lhard_surf) 
       Zero_density_TF[inode_box][Nlists_HW-1] = FALSE;
 }

   /* special case of zeolites for M.Mitchell 
     for (inode_box=0; inode_box<Nnodes_box; inode_box++) {
        node_box_to_ijk_box(inode_box,ijk_box);
        flag=TRUE;
        for (idim=0; idim<Ndim; idim++){
           ijk[idim] = ijk_box[idim] + Min_IJK_box[idim];
           if (ijk[idim]<0 || ijk[idim] >= Nodes_x[idim]) flag=FALSE;
        }
        if (flag){
           rsq = 0.0;
           for (idim=0; idim<Ndim; idim++) { 
              rsq += (ijk[idim]*Esize_x[idim]-0.5*Size_x[idim])*
                     (ijk[idim]*Esize_x[idim]-0.5*Size_x[idim]);
           }
           if (rsq > 0.5*Size_x[0]*0.5*Size_x[0]) {
              for (icomp=0; icomp<Ncomp; icomp++)
                    Zero_density_TF[inode_box][icomp] = TRUE;
           }
        }
     }*/

   /* build an array of surfaces that are _NOT_ atoms.  We need
      this because every 3d system seems to have O(1000) atoms,
      and the remainder of this piece of code is very very 
      slow.  If it could be sped up, we could remove this logic
      (and the atomic cases in dft_vext.c).   If we are doing 
     walls based on atomic centers (LJ or HS) potentials -
     make all adjustments to Zero_TF in the vext routine */

    nonatomic_walls = (int *) array_alloc (1, Nwall, sizeof(int));
    nonat_count=0;
    for (iwall=0; iwall<Nwall;iwall++){
      if (Surface_type[WallType[iwall]] != atomic_centers || Nwall<1000){
          nonatomic_walls[nonat_count]=iwall;
          nonat_count++;
      }
    }

 n_fluid_els = (int **) array_alloc (2, nonat_count,Nlists_HW, sizeof(int));
 n_el_in_box = (int **) array_alloc (2, nonat_count,Nlists_HW, sizeof(int));
 els_one_owner = (int *) array_alloc (1, Nnodes_per_el_V, sizeof(int));
 countw = (int **) array_alloc (2, nonat_count,Nlists_HW,sizeof(int));

  if (nonat_count > 0 || Nwall<0){
printf("setting up problem with old ZeroTF routine\n");
 /* 
  * loop over all nodes in the box coordinates of this processor
  * to set up Nodes_2_boundary_wall array 
  */
 for (inode_box=0; inode_box<Nnodes_box; inode_box++) {

    inode = B2G_node[inode_box];
    node_to_ijk(inode,ijk);

    for (iwall=0; iwall< nonat_count; iwall++)
       for (ilist=0; ilist<Nlists_HW; ilist++){
          n_fluid_els[iwall][ilist] = 0;
          n_el_in_box[iwall][ilist] = 0;
       }

    for (iwall=0; iwall<nonat_count; iwall++){
       for (ilist=0; ilist<Nlists_HW; ilist++) countw[iwall][ilist]=0;
    }

    for (loc_node_el=0; loc_node_el<Nnodes_per_el_V ; loc_node_el++){
   
       iel = node_to_elem_return_dim(inode, loc_node_el,reflect_flag,
                                     &idim,&iside,&periodic_flag);
       if (iel >= 0) 
	 {
	   iel_box = el_to_el_box(iel); 
	 }
       else 
	 {
	   iel_box = 0;
	 }

       for (iw=0; iw<nonat_count; iw++){
       iwall=nonatomic_walls[iw];
       for (ilist=0; ilist<Nlists_HW; ilist++){
          if (iel >= 0){

            if (iel_box >= 0) {
               n_el_in_box[iw][ilist]++;
               if (el_type[iwall][ilist][iel_box] == FLUID_EL)
                                   n_fluid_els[iw][ilist]++;
               wall_flag = FALSE;
               for (jw=0; jw<nonat_count; jw++)
                   if (el_type[nonatomic_walls[jw]][ilist][iel_box] == WALL_EL) wall_flag = TRUE;
               if (wall_flag) countw[iw][ilist]++;
            }

          }
          else if (iel == -2){
             n_fluid_els[iw][ilist]++;
             n_el_in_box[iw][ilist]++;  /*i.e. we know what is beyond this boundary */
          }
          else if (iel == -1){
             /* find if iwall touches boundary */
             if (Touch_domain_boundary[iwall][ilist][idim][iside]){
                n_el_in_box[iw][ilist]++;  
                countw[iw][ilist]++;
             }
          }

       }  /* end of loop over lists */
       }  /* end of loop over walls */
    }     /* end of loop over elements that touch a given node */


    /* if any one wall owns all the elements around a given 
       node, the number of fluid elements is zero for all
       walls. */
    for (iw=0; iw<nonat_count; iw++)
       for (ilist=0; ilist<Nlists_HW; ilist++) 
          if (countw[iw][ilist]==Nnodes_per_el_V) {
              for (jw=0; jw<nonat_count; jw++)
                        n_fluid_els[jw][ilist]=0;
          }


/*    for (iw=0; iw<Nwall; iw++){
    for (ilist=0; ilist<Nlists_HW; ilist++){
        printf("inode_box: %d  iwall: %d  ilist: %d  n_fluid_els: %d\n",
                  inode_box,nonatomic_walls[iw],ilist,n_fluid_els[iw][ilist]);
        if (n_fluid_els[iw][ilist] != 0)
           printf("iwall %d   ilist %d   n_fluid: %d   n_boc: %d \n",
                    nonatomic_walls[iw],ilist,n_fluid_els[iw][ilist],n_el_in_box[iw][ilist]);
    }}*/

    /* Now set up the Nodes_2_boundary_wall array */
    for (iw=0; iw<nonat_count; iw++){
    for (ilist=0; ilist<Nlists_HW; ilist++){

        if (n_fluid_els[iw][ilist] == 0){
            if ( Nlists_HW == Ncomp+1)
                Zero_density_TF[inode_box][ilist] = TRUE;
            else if (ilist == 0){
              for (icomp=0; icomp<Ncomp; icomp++) 
                Zero_density_TF[inode_box][icomp] = TRUE;
            }
        }

        else if (n_fluid_els[iw][ilist] != n_el_in_box[iw][ilist]) {
            Nodes_2_boundary_wall[ilist][inode_box] = nonatomic_walls[iw]; 
        }
    }
    }

    /*
     * for hard sphere cases, we need to double check for 1 sigma
     *  surface separations where the density is nonzero at a point
     *  where all surrounding elements are wall elements. 
     */

    if (Lhard_surf){
       for (ilist=0; ilist<Nlists_HW; ilist++){
          if (Zero_density_TF[inode_box][ilist]){

               /*
                * all wall elements around this point 
                * check to see if there are at least two
                * elements that have only one owner, and
                * that those owners are different. If so,
                * set Zero_den_TF = FALSE  and set Nodes_2_boundary
                * wall to -2 to indicate that there are more
                * than 1 wall involved. 
                */


              count = 0;
              for (loc_node_el=0; loc_node_el<Nnodes_per_el_V; loc_node_el++){
  
                 iel = node_to_elem(inode, loc_node_el,reflect_flag);
                 if (iel >= 0) {
                   iel_box = el_to_el_box(iel); 
                   if (Nwall_owners[ilist][iel_box] == 1){
                      els_one_owner[count++] = iel_box;
                   } 
                 }
               } /* end loop over elements touching inode_box */

               all_equal = FALSE;
               if (count >= 2) {

                  /* test to see if all wall owners are the same*/
                  all_equal = TRUE;
                  test =  Wall_owners[ilist][els_one_owner[0]][0];
                  for (i=1; i<count; i++){
                      if (Wall_owners[ilist][els_one_owner[i]][0] != test)
                                                    all_equal = FALSE;
                  }
                  if (all_equal == FALSE){
                      Nodes_2_boundary_wall[ilist][inode_box] = -2;
                      if (ilist < Nlists_HW-1)
                      Zero_density_TF[inode_box][ilist] = FALSE;
                  }
               }

               /* 
                * If on a reflective boundary, the wall owners need
                * not be different for this to be a wall-wall boundary.
                * only condition is that there is only one owner on
                * at least one of the elements.  However, need to be
                * careful about reflections that go through the center
                * of a wall .... these do not result in new wall-wall
                * boundary nodes !!!
                */  
               node_to_ijk(inode,ijk);
               reflect = FALSE;
               ndim_reflect=0;
               for (idim=0; idim<Ndim; idim++){
                   if((ijk[idim] == 0 && Type_bc[idim][0] == REFLECT) ||
                      (ijk[idim] == Nodes_x[idim]-1 && Type_bc[idim][1] == REFLECT)) {
                       reflect = TRUE;
                       dim_reflect[ndim_reflect++] = idim;
                   }
               }

               /* assume that this node is not on centerline of a surface to start */
               test=FALSE;
               if (all_equal == TRUE && reflect){
                   link_chk = Wall_owners[ilist][els_one_owner[0]][0];
                   for (i=0; i<ndim_reflect; i++){
                     idim = dim_reflect[i];
                     if(Xtest_reflect_TF[link_chk][idim]){
                      for (ilink=0; ilink<Nwall_this_link[link_chk]; ilink++){
                         wall_chk=Link_list[link_chk][ilink];
                         if (ijk[idim]==0 && Type_bc[idim][0] == REFLECT &&
                             fabs(WallPos[idim][wall_chk]+0.5*Size_x[idim]) > 0.00000001) 
                             test = TRUE;
                         else if (ijk[idim] == Nodes_x[idim]-1 && Type_bc[idim][1] == REFLECT &&
                             fabs(0.5*Size_x[idim]-WallPos[idim][wall_chk]) > 0.00000001) 
                             test = TRUE;
                      }
                     }  
                   }

                  if (test == TRUE) { /* node proves to _not_ be located on a center line */
                     Nodes_2_boundary_wall[ilist][inode_box] = -2;
                     if (ilist < Nlists_HW-1) Zero_density_TF[inode_box][ilist] = FALSE;
                  }
               } 

          }  /*end of Zero_density_TF check */
       }     /* end of loop over lists */
    }        /* end of x-tra checks for the hard wall cases */

 }        /* End of loop over nodes in local box coordinates */
 }


 if (Imain_loop==0 && Proc==0 && Iwrite==VERBOSE) {
    fprintf (fp1,"\n---------------------------------------------------------------\n");
    fprintf (fp1,"Proc: %d \n",Proc);
    fprintf(fp1,"Have set up arrays Wall_elems, and Nodes_2_boundary_wall\n");

    for (ilist=0; ilist<Nlists_HW; ilist++){
         for (inode_box=0; inode_box<Nnodes_box; inode_box++){
             node_to_ijk(B2G_node[inode_box],ijk);
             if (Nodes_2_boundary_wall[ilist][inode_box] != -1)
                fprintf(fp1,"ilist: %d   inode: %d  ijk: %d iwall: %d\n",
                           ilist,B2G_node[inode_box],ijk[0],
                           Nodes_2_boundary_wall[ilist][inode_box]);
         }
    }
    fprintf (fp1,"\n---------------------------------------------------------------\n");
 } 

 safe_free((void *) &n_fluid_els);
 safe_free((void *) &n_el_in_box);
 safe_free((void *) &Touch_domain_boundary);
 safe_free((void *) &els_one_owner);
 safe_free((void *) &countw);
 safe_free((void *) &nonatomic_walls);

 return;
}
/*****************************************************************************
setup_zeroTF_and_Node2bound_new: This routine take the information about which elements
                    are wall and fluid elements, and translates that into
                    lists of wall nodes and fluid nodes for each wall list */
void setup_zeroTF_and_Node2bound_new (FILE *fp1,int ***el_type)
{
 int reflect_flag[3], *n_fluid_els, iwall,
      iel, idim, ilist, loc_node_el,
      inode,inode_box,iel_box,icomp, iside,
      periodic_flag,wall_flag,jwall;
 int n_el_in_box, *countw_per_w,*countw_per_list,*countw_per_link;
 int *counted_list,*counted_link;
 int test,all_equal,ijk[3], reflect,count,*els_one_owner,i;
 int dim_reflect[3],ndim_reflect,wall_chk; 
 int link_chk,ilink;
 int jw,index,index_w,turn_off;
 double node_pos[3];

 Nodes_2_boundary_wall = (int **) array_alloc(2, Nlists_HW, Nnodes_box, sizeof(int));
 Zero_density_TF = (int **) array_alloc (2, Nnodes_box,Ncomp+1,sizeof(int));

 /* zero the arrays that will be set up here */

 for (idim=0; idim<Ndim; idim++) reflect_flag[idim] = FALSE;

 for (inode_box=0; inode_box<Nnodes_box; inode_box++){
    for (ilist=0; ilist<Nlists_HW; ilist++) {
            Nodes_2_boundary_wall[ilist][inode_box] = -1;
    }
    for (icomp=0; icomp<Ncomp+1; icomp++) {
            Zero_density_TF[inode_box][icomp] = FALSE;
    }
 }

   /* special case of zeolites for M.Mitchell 
     for (inode_box=0; inode_box<Nnodes_box; inode_box++) {
        node_box_to_ijk_box(inode_box,ijk_box);
        flag=TRUE;
        for (idim=0; idim<Ndim; idim++){
           ijk[idim] = ijk_box[idim] + Min_IJK_box[idim];
           if (ijk[idim]<0 || ijk[idim] >= Nodes_x[idim]) flag=FALSE;
        }
        if (flag){
           rsq = 0.0;
           for (idim=0; idim<Ndim; idim++) { 
              rsq += (ijk[idim]*Esize_x[idim]-0.5*Size_x[idim])*
                     (ijk[idim]*Esize_x[idim]-0.5*Size_x[idim]);
           }
           if (rsq > 0.5*Size_x[0]*0.5*Size_x[0]) {
              for (icomp=0; icomp<Ncomp; icomp++)
                    Zero_density_TF[inode_box][icomp] = TRUE;
           }
        }
     }*/

 n_fluid_els = (int *) array_alloc (1, Nwall*Nlists_HW, sizeof(int));
 els_one_owner = (int *) array_alloc (1, Nnodes_per_el_V, sizeof(int));
/* keep track for each wall and each list */
 countw_per_w = (int *) array_alloc (1, Nwall*Nlists_HW,sizeof(int));
/* keep track for each list but be blind to the particular wall */
 countw_per_list = (int *) array_alloc (1, Nlists_HW,sizeof(int));
 counted_list = (int *) array_alloc (1, Nlists_HW,sizeof(int));
/* keep track for each linked wall and each list */
 countw_per_link = (int *) array_alloc (1, Nlink*Nlists_HW,sizeof(int));
 counted_link = (int *) array_alloc (1, Nlink*Nlists_HW,sizeof(int));

 /* 
  * loop over all nodes in the box coordinates of this processor
  * to set up Nodes_2_boundary_wall array -----
  * only loop over nodes in the box that actually do touch at least one wall !! 
  */
 for (index=0; index<Nnodes_wall_box; index++) {

    inode_box = Nodes_wall_box[index];
    inode = B2G_node[inode_box];
    node_to_ijk(inode,ijk);

    for (index_w=0; index_w<Nwall_touch_node[index]; index_w++){
       n_fluid_els[index_w] = 0;
       countw_per_w[index_w]=0;
    }
    for (ilist=0;ilist<Nlists_HW;ilist++){
         countw_per_list[ilist]=0;
         for (ilink=0;ilink<Nlink;ilink++) countw_per_link[ilist+Nlists_HW*ilink]=0;
    }
    n_el_in_box=0; 

    for (loc_node_el=0; loc_node_el<Nnodes_per_el_V ; loc_node_el++){
       iel = node_to_elem_return_dim(inode, loc_node_el,reflect_flag,
                                     &idim,&iside,&periodic_flag);
       if (iel >= 0) iel_box = el_to_el_box(iel); 
       else iel_box = 0;

       if (iel_box>=0 || iel==-2 || iel==-1) n_el_in_box++;

       for (ilist=0;ilist<Nlists_HW; ilist++){
              counted_list[ilist]=FALSE;
             for (ilink=0;ilink<Nlink;ilink++) counted_link[ilist+Nlists_HW*ilink]=FALSE;
       }

       for (index_w=0; index_w<Nwall_touch_node[index]; index_w++){
         iwall=Wall_touch_node[index][index_w];
         ilist=List_wall_node[index][index_w];

          if (iel >= 0){

            if (iel_box >= 0) {
               if (el_type[iwall][ilist][iel_box] == FLUID_EL) n_fluid_els[index_w]++;
               else{
                    countw_per_w[index_w]++;
                    if (!counted_list[ilist]){
                          countw_per_list[ilist]++;
                          counted_list[ilist]=TRUE;
                    }
                    if (!counted_link[ilist+Nlists_HW*Link[iwall]]){
                         countw_per_link[ilist+Nlists_HW*Link[iwall]]++;
                         counted_link[ilist+Nlists_HW*Link[iwall]]=TRUE;
                    }
               }

/*               wall_flag = FALSE;
               for (jw=0; jw<Nwall_touch_node[index]; jw++){
                 if (ilist==List_wall_node[index][jw]){
                   jwall=Wall_touch_node[index][jw];
                   if (el_type[jwall][ilist][iel_box] == WALL_EL) wall_flag = TRUE;
                 }
               }
               if (wall_flag) countw[index_w]++;*/
            }

          }
          else if (iel == -2){        /* this indicates a bulk fluid boundary */
             n_fluid_els[index_w]++;
          }
          else if (iel == -1){         /* indicates an semi-infinite surface boundary */
             if (Touch_domain_boundary[iwall][ilist][idim][iside]){
                countw_per_w[index_w]++;
                if (!counted_list[ilist]){
                      countw_per_list[ilist]++;
                      counted_list[ilist]=TRUE;
                }
                if (!counted_link[ilist+Nlists_HW*Link[iwall]]){
                     countw_per_link[ilist+Nlists_HW*Link[iwall]]++;
                     counted_link[ilist+Nlists_HW*Link[iwall]]=TRUE;
                }
             }
          }

       }  /* end of loop over walls */
    }     /* end of loop over elements that touch a given node */


    /* in the case where n_el_in_box==Nnodes_per_el_V...
         if (i) any one wall or (ii) any linked wall or (iii) any list has all walls
         in the absence of hard surfaces then ...
       the number of fluid elements is zero for all walls that touch this node
       and have the same list number*/
    if (n_el_in_box==Nnodes_per_el_V){
       for (index_w=0; index_w<Nwall_touch_node[index]; index_w++){
           ilist=List_wall_node[index][index_w];
           iwall=Wall_touch_node[index][index_w];
           if (countw_per_w[index_w]==Nnodes_per_el_V || 
               countw_per_list[ilist]==Nnodes_per_el_V ||
               countw_per_link[ilist+Nlists_HW*Link[iwall]]==Nnodes_per_el_V){ 
               for (jw=0; jw<Nwall_touch_node[index]; jw++) 
                    if (ilist==List_wall_node[index][jw]) n_fluid_els[jw]=0;
           }
       }
    }

    /* Now set up the Nodes_2_boundary_wall array */

    for (index_w=0; index_w<Nwall_touch_node[index]; index_w++){
        ilist=List_wall_node[index][index_w];

        if (n_fluid_els[index_w] == 0){
            if ( Nlists_HW == Ncomp+1 || ilist==Ncomp)
                Zero_density_TF[inode_box][ilist] = TRUE;
            else if (ilist == 0){
                 for (icomp=0; icomp<Ncomp; icomp++) {
                   Zero_density_TF[inode_box][icomp] = TRUE;
                 }
            }
        }

        else if ((n_fluid_els[index_w] != n_el_in_box) /*&& n_fluid_els[index_w]>0*/) {
            Nodes_2_boundary_wall[ilist][inode_box] = Wall_touch_node[index][index_w]; 
        }
    }


    /*
     * for hard sphere cases, we need to double check for 1 sigma
     *  surface separations where the density is nonzero at a point
     *  where all surrounding elements are wall elements. 
     */
/* turn this off */ turn_off=TRUE;
    if (Lhard_surf && !turn_off){
       for (ilist=0; ilist<Nlists_HW; ilist++){
          if (Zero_density_TF[inode_box][ilist]){

               /*
                * all wall elements around this point 
                * check to see if there are at least two
                * elements that have only one owner, and
                * that those owners are different. If so,
                * set Zero_den_TF = FALSE  and set Nodes_2_boundary
                * wall to -2 to indicate that there are more
                * than 1 wall involved. 
                */


              count = 0;
              for (loc_node_el=0; loc_node_el<Nnodes_per_el_V; loc_node_el++){
  
                 iel = node_to_elem(inode, loc_node_el,reflect_flag);
                 if (iel >= 0) {
                   iel_box = el_to_el_box(iel); 
                   if (Nwall_owners[ilist][iel_box] == 1){
                      els_one_owner[count++] = iel_box;
                   } 
                 }
               } /* end loop over elements touching inode_box */

               all_equal = FALSE;
               if (count >= 2) {

                  /* test to see if all wall owners are the same*/
                  all_equal = TRUE;
                  test =  Wall_owners[ilist][els_one_owner[0]][0];
                  for (i=1; i<count; i++){
                      if (Wall_owners[ilist][els_one_owner[i]][0] != test)
                                                    all_equal = FALSE;
                  }
                  if (all_equal == FALSE){
                      Nodes_2_boundary_wall[ilist][inode_box] = -2;
                      if (ilist < Nlists_HW-1)
                      Zero_density_TF[inode_box][ilist] = FALSE;
                  }
               }

               /* 
                * If on a reflective boundary, the wall owners need
                * not be different for this to be a wall-wall boundary.
                * only condition is that there is only one owner on
                * at least one of the elements.  However, need to be
                * careful about reflections that go through the center
                * of a wall .... these do not result in new wall-wall
                * boundary nodes !!!
                */  
               node_to_ijk(inode,ijk);
               reflect = FALSE;
               ndim_reflect=0;
               for (idim=0; idim<Ndim; idim++){
                   if((ijk[idim] == 0 && Type_bc[idim][0] == REFLECT) ||
                      (ijk[idim] == Nodes_x[idim]-1 && Type_bc[idim][1] == REFLECT)) {
                       reflect = TRUE;
                       dim_reflect[ndim_reflect++] = idim;
                   }
               }

               /* assume that this node is not on centerline of a surface to start */
               test=FALSE;
               if (all_equal == TRUE && reflect){
                   link_chk = Wall_owners[ilist][els_one_owner[0]][0];
                   for (i=0; i<ndim_reflect; i++){
                     idim = dim_reflect[i];
                     if(Xtest_reflect_TF[link_chk][idim]){
                      for (ilink=0; ilink<Nwall_this_link[link_chk]; ilink++){
                         wall_chk=Link_list[link_chk][ilink];
                         if (ijk[idim]==0 && Type_bc[idim][0] == REFLECT &&
                             fabs(WallPos[idim][wall_chk]+0.5*Size_x[idim]) > 0.00000001) 
                             test = TRUE;
                         else if (ijk[idim] == Nodes_x[idim]-1 && Type_bc[idim][1] == REFLECT &&
                             fabs(0.5*Size_x[idim]-WallPos[idim][wall_chk]) > 0.00000001) 
                             test = TRUE;
                      }
                     }  
                   }

                  if (test == TRUE) { /* node proves to _not_ be located on a center line */
                     Nodes_2_boundary_wall[ilist][inode_box] = -2;
                     if (ilist < Nlists_HW-1) Zero_density_TF[inode_box][ilist] = FALSE;
                  }
               } 

          }  /*end of Zero_density_TF check */
       }     /* end of loop over lists */
    }        /* end of x-tra checks for the hard wall cases */

 }        /* End of loop over nodes in local box coordinates */

 if (Imain_loop==0 && Proc==0 && Iwrite==VERBOSE) {
    fprintf (fp1,"\n---------------------------------------------------------------\n");
    fprintf (fp1,"Proc: %d \n",Proc);
    fprintf(fp1,"Have set up arrays Wall_elems, and Nodes_2_boundary_wall\n");

    for (ilist=0; ilist<Nlists_HW; ilist++){
         for (inode_box=0; inode_box<Nnodes_box; inode_box++){
             node_to_ijk(B2G_node[inode_box],ijk);
             if (Nodes_2_boundary_wall[ilist][inode_box] != -1)
                fprintf(fp1,"ilist: %d   inode: %d   iwall: %d\n",
                           ilist,B2G_node[inode_box], Nodes_2_boundary_wall[ilist][inode_box]);
         }
    }
    fprintf (fp1,"\n---------------------------------------------------------------\n");
 } 

 safe_free((void *) &n_fluid_els);
 safe_free((void *) &Touch_domain_boundary);
 safe_free((void *) &els_one_owner);
 safe_free((void *) &countw_per_w);
 safe_free((void *) &countw_per_list);
 safe_free((void *) &countw_per_link);

 return;
}
/***************************************************************************
*boundary_properties:  In this routine we identify the surface elements     *
*                      that each boundary node is a part of.  We store      *
*                      the normal to the surface element for future use.    */

void boundary_properties(FILE *fp1)
{
  int norm,ilist,loc_inode,iwall,i,j,*iel,idim,inode,iel_s,
      loc_node_el,reflect_flag[3], ijk[3],inode_box,surf_norm,
      sum[3],sum_all[3],*iel_box,flag,el,type,normal;
  double  s_area_tot_proc,esize1=1.0,esize2=1.0,sarea_sum;
  int test,all_equal,reflect,count,*els_one_owner;
  int dim_reflect[3],ndim_reflect,ielement,ielement_box; 
  int flag_fluid;
  int link_chk;

  els_one_owner = (int *) array_alloc (1, Nnodes_per_el_V, sizeof(int));
  iel =     (int *) array_alloc(1, Nnodes_per_el_V, sizeof(int));
  iel_box = (int *) array_alloc(1, Nnodes_per_el_V, sizeof(int));
  Nelems_S = (int **) array_alloc(2, Nlists_HW, Nnodes_per_proc, sizeof(int));
  Surf_normal = (int ***) array_alloc (3, Nlists_HW, Nnodes_per_proc,
                                      2*Nnodes_per_el_V,sizeof(int));
  Surf_elem_type = (int **) array_alloc (2, Nnodes_per_proc,
                                         2*Nnodes_per_el_V,sizeof(int));
  S_area = (double ***) array_alloc(3, Nlists_HW, Nwall,Ndim, sizeof(double));
  S_area_tot = (double **) array_alloc(2, Nlists_HW, Nwall,sizeof(double));
  Surf_elem_to_wall = (int ***) array_alloc (3, Nlists_HW,Nnodes_per_proc,
                                         2*Nnodes_per_el_V,sizeof(int));
  norm = Nnodes_per_el_S;

  for (ilist=0; ilist<Nlists_HW; ilist++){

     for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
        Nelems_S[ilist][loc_inode] = 0;
        for (iel_s=0; iel_s<Nnodes_per_el_S; iel_s++){
            Surf_normal[ilist][loc_inode][iel_s] = 0;
            if (ilist == Nlists_HW-1) Surf_elem_type[loc_inode][iel_s]= -1;
        }
     }
 
     for (iwall=0; iwall<Nwall; iwall++)
        for(idim=0; idim<Ndim; idim++) S_area[ilist][iwall][idim] = 0.0;
  }

  for (idim=0; idim<Ndim; idim++) reflect_flag[idim] = FALSE;

  for (ilist=0; ilist<Nlists_HW; ilist++){
      for (loc_inode=0; loc_inode<Nnodes_per_proc;loc_inode++){

          inode = L2G_node[loc_inode];
          inode_box = node_to_node_box(inode);
          flag_fluid = FALSE;
        if (Nodes_2_boundary_wall[ilist][inode_box] >= 0){   /* wall-fluid
                                                              boundary nodes */
            flag = -1;
            find_local_els(inode,iel,iel_box,flag);

            if (Ndim == 1){
                type = 0; esize1=1.0; esize2=1.0; idim=0;

                if (Wall_elems[ilist][iel_box[0]] != -1) {
                   el = 0; normal = -1; 
                }
                else{
                   el = 1; normal = 1;
                }
                surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                               normal,idim,esize1,esize2);
            }
            else{
                for (i=0; i<Ndim-1 ; i++){
                   j=4*i;
                 
                   /* check for surface element in left/right (on top elements)*/

                   if ( (iel[j] >= 0 && iel[j+1]>= 0) && (
                        (Wall_elems[ilist][iel_box[j]] != -1  && 
                         Wall_elems[ilist][iel_box[j+1]] == -1) ||
                        (Wall_elems[ilist][iel_box[j]] == -1 && 
                         Wall_elems[ilist][iel_box[j+1]] != -1)  )  ) {

                       if (j == 0) type = UP_BACK;
                       else        type = UP_FRONT;
                                      esize1 = Esize_x[1];
                       if (Ndim == 3) esize2 = Esize_x[2];
                       idim = 0;

                       if (Wall_elems[ilist][iel_box[j]] != -1){
                              el = j;   normal = -1;
                       }
                       else {  
                              el = j+1; normal =  1; 
                       }
                       surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                        normal,idim,esize1,esize2);
                   }

                   /* check for surface element in top/bottom (on right elements)*/

                   if ( (iel[j] >= 0 && iel[j+2]>= 0) && (
                        (Wall_elems[ilist][iel_box[j]] != -1  && 
                         Wall_elems[ilist][iel_box[j+2]] == -1) ||
                        (Wall_elems[ilist][iel_box[j]] == -1 && 
                         Wall_elems[ilist][iel_box[j+2]] != -1)  )  ) {

                       if (j == 0) type = RIGHT_BACK;
                       else        type = RIGHT_FRONT;
                                      esize1 = Esize_x[0];
                       if (Ndim == 3) esize2 = Esize_x[2];
                       idim = 1;

                       if (Wall_elems[ilist][iel_box[j]] != -1){
                              el = j;   normal = -2;
                       }
                       else {  
                              el = j+2; normal =  2; 
                       }
                       surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                        normal,idim,esize1,esize2);
                   }

                   /* check for surface element in top/bottom (on left elements)*/
                   if ( (iel[j+1] >= 0 && iel[j+3]>= 0) && (
                        (Wall_elems[ilist][iel_box[j+1]] != -1  && 
                         Wall_elems[ilist][iel_box[j+3]] == -1) ||
                        (Wall_elems[ilist][iel_box[j+1]] == -1 && 
                         Wall_elems[ilist][iel_box[j+3]] != -1)  )  ) {

                       if (j == 0) type = LEFT_BACK;
                       else        type = LEFT_FRONT;
                                      esize1 = Esize_x[0];
                       if (Ndim == 3) esize2 = Esize_x[2];
                       idim = 1;

                       if (Wall_elems[ilist][iel_box[j+1]] != -1){
                              el = j+1; normal = -2;
                       }
                       else {  
                              el = j+3; normal =  2; 
                       }
                       surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                        normal,idim,esize1,esize2);
                   }

                   /* check for surface element in right/left (on bottom elements)*/
                   if ( (iel[j+2] >= 0 && iel[j+3]>= 0) && (
                        (Wall_elems[ilist][iel_box[j+2]] != -1  && 
                         Wall_elems[ilist][iel_box[j+3]] == -1) ||
                        (Wall_elems[ilist][iel_box[j+2]] == -1 && 
                         Wall_elems[ilist][iel_box[j+3]] != -1)  )  ) {

                       if (j == 0) type = DOWN_BACK;
                       else        type = DOWN_FRONT;
                                      esize1 = Esize_x[1];
                       if (Ndim == 3) esize2 = Esize_x[2];
                       idim = 0;

                       if (Wall_elems[ilist][iel_box[j+2]] != -1){
                              el = j+2; normal = -1;
                       }
                       else {  
                              el = j+3; normal =  1; 
                       }
                       surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                        normal,idim,esize1,esize2);
                   }
                }
                if (Ndim == 3){
                   for (i=0; i<4 ; i++){
                      j=i+4;

                      /* check for surface element in back/front */
                      if ( (iel[i] >= 0 && iel[j]>= 0) && (
                           (Wall_elems[ilist][iel_box[i]] != -1  && 
                            Wall_elems[ilist][iel_box[j]] == -1) ||
                           (Wall_elems[ilist][iel_box[i]] == -1 && 
                            Wall_elems[ilist][iel_box[j]] != -1)  )  ) {

                          type = i+8;
                                         esize1 = Esize_x[0];
                          if (Ndim == 3) esize2 = Esize_x[1];
                          idim = 2;
   
                          if (Wall_elems[ilist][iel_box[i]] != -1){
                                 el = i; normal = -3;
                          }
                          else {  
                                 el = j; normal = 3; 
                          }
                          surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                            normal,idim,esize1,esize2);

                      }
                   }
   
                }  /* end of 3D case */
            }     /* end of 2D and 3D cases */


            /* check to see if this Wall-Fluid boundary node
               is _also_ a wall-wall boundary node !!!*/

    if (Lhard_surf){

              count = 0;
              for (loc_node_el=0; loc_node_el<Nnodes_per_el_V; loc_node_el++){
  
                 ielement = node_to_elem(inode, loc_node_el,reflect_flag);
                 if (ielement >= 0) {
                   ielement_box = el_to_el_box(ielement); 
                   if (Nwall_owners[ilist][ielement_box] == 1){
                      els_one_owner[count++] = ielement_box;
                   } 
                 }
               } /* end loop over elements touching inode_box */

               all_equal = FALSE;
               if (count >= 2) {
                  all_equal = TRUE;
                  test =  Wall_owners[ilist][els_one_owner[0]][0];
                  for (i=1; i<count; i++){
                      if (Wall_owners[ilist][els_one_owner[i]][0] != test)
                                                    all_equal = FALSE;
                  }
                  if (all_equal == FALSE) flag_fluid = TRUE;
               }

               /* 
                * If on a reflective boundary, the wall owners need
                * not be different for this to be a wall-wall boundary.
                * only condition is that there is only one owner on
                * at least one of the elements.  However, need to be
                * careful about reflections that go through the center
                * of a wall .... these do not result in new wall-wall
                * boundary nodes !!!
                */  
               node_to_ijk(inode,ijk);
               reflect = FALSE;
               ndim_reflect=0;
               for (idim=0; idim<Ndim; idim++){
                   if((ijk[idim] == 0 && Type_bc[idim][0] == REFLECT) ||
                      (ijk[idim] == Nodes_x[idim]-1 && Type_bc[idim][1] == REFLECT)) {
                       reflect = TRUE;
                       dim_reflect[ndim_reflect++] = idim;
                   }
               }

               /* assume that this node is on centerline of a surface to start */
               test = FALSE;
               if (all_equal == TRUE && reflect){
                   link_chk = Wall_owners[ilist][els_one_owner[0]][0];
                   for (i=0; i<ndim_reflect; i++){
                     idim = dim_reflect[i];
                     if (Xtest_reflect_TF[link_chk][idim]){
                         if (ijk[idim]==0 && Type_bc[idim][0] == REFLECT){
                               test = TRUE;
                         }
                         else if (ijk[idim] == Nodes_x[idim]-1 && Type_bc[idim][1] == REFLECT){
                               test = TRUE;
                         }
                     }
                   }
               } 
               if (test == TRUE) flag_fluid = TRUE;

    }        /* end of x-tra checks for the hard wall cases */

         
        } /* end of test for if this is a wall-fluid boundary node */

      /*now check for wall-wall boundary nodes */
        if (Nodes_2_boundary_wall[ilist][inode_box] == -2 || flag_fluid){ 
            flag = -2;
            find_local_els(inode,iel,iel_box,flag);

            if (Ndim == 1){
                type = 0; esize1=1.0; esize2=1.0; idim=0;
                if (iel[0] >=0){
                   el = 0; normal = -1; 
                   surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                     normal,idim,esize1,esize2);
                }
                if (iel[1] >= 0){
                   el = 1; normal = 1; 
                   surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                     normal,idim,esize1,esize2);
                }
             }
             else{
                for (i=0; i<Ndim-1 ; i++){
                   j=4*i;
                 
                   /* check for surface element in left/right (on top elements)*/

                   if ( (iel[j] >= 0 && iel[j+1]>= 0 && 
                         Nwall_owners[ilist][iel_box[j]] == 1   && 
                         Nwall_owners[ilist][iel_box[j+1]] == 1 &&
                         Wall_elems[ilist][iel_box[j]] != Wall_elems[ilist][iel_box[j+1]] ) ||

                        ( iel[j] == flag && iel[j+1]>= 0 && 
                         Nwall_owners[ilist][iel_box[j+1]] == 1 &&
                         fabs(0.5*Size_x[0] - WallPos[0][Wall_elems[ilist][iel_box[j+1]]]) 
                                                                          > 0.00000001 ) ||

                        ( iel[j] >= 0 && iel[j+1] == flag && 
                         Nwall_owners[ilist][iel_box[j]] == 1 &&
                         fabs(WallPos[0][Wall_elems[ilist][iel_box[j]]] + 0.5*Size_x[0]) 
                                                                          > 0.00000001 ) ){

                       if (j == 0) type = UP_BACK;
                       else        type = UP_FRONT;
                                      esize1 = Esize_x[1];
                       if (Ndim == 3) esize2 = Esize_x[2];
                       idim = 0;

                       if (iel[j] >= 0) {
                          el = j;   normal = -1;
                          surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                            normal,idim,esize1,esize2);
                       }
                       if (iel[j+1] >= 0){ 
                          el = j+1; normal =  1; 
                          surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                            normal,idim,esize1,esize2);
                       }
                   }

                   /* check for surface element in top/bottom (on right elements)*/
                   if ( ( iel[j] >= 0 && iel[j+2]>= 0 && 
                         Nwall_owners[ilist][iel_box[j]] == 1   &&
                         Nwall_owners[ilist][iel_box[j+2]] == 1 &&
                         Wall_elems[ilist][iel_box[j]] != Wall_elems[ilist][iel_box[j+2]] ) ||

                        ( iel[j] == flag && iel[j+2]>= 0 && 
                         Nwall_owners[ilist][iel_box[j+2]] == 1 &&
                         fabs(0.5*Size_x[1] - WallPos[1][Wall_elems[ilist][iel_box[j+2]]]) 
                                                                          > 0.00000001 ) ||

                        ( iel[j] >= 0 && iel[j+2] == flag && 
                         Nwall_owners[ilist][iel_box[j]] == 1 &&
                         fabs(WallPos[1][Wall_elems[ilist][iel_box[j]]] + 0.5*Size_x[1]) 
                                                                          > 0.00000001 ) ){

                       if (j == 0) type = RIGHT_BACK;
                       else        type = RIGHT_FRONT;
                                      esize1 = Esize_x[0];
                       if (Ndim == 3) esize2 = Esize_x[2];
                       idim = 1;

                       if (iel[j] >=0) {
                          el = j;   normal = -2;
                          surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                            normal,idim,esize1,esize2);
                       }
                       if (iel[j+2] >=0) {
                          el = j+2; normal =  2; 
                          surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                            normal,idim,esize1,esize2);
                       }
                   }

                   /* check for surface element in top/bottom (on left elements)*/
                   if ( ( iel[j+1] >= 0 && iel[j+3]>= 0 && 
                        Nwall_owners[ilist][iel_box[j+1]] == 1   &&
                        Nwall_owners[ilist][iel_box[j+3]] == 1   &&
                        Wall_elems[ilist][iel_box[j+1]] != Wall_elems[ilist][iel_box[j+3]] ) ||

                        ( iel[j+1] == flag && iel[j+3]>= 0      && 
                         Nwall_owners[ilist][iel_box[j+3]] == 1 &&
                         fabs(0.5*Size_x[1] - WallPos[1][Wall_elems[ilist][iel_box[j+3]]]) 
                                                                          > 0.00000001 ) ||

                        ( iel[j+1] >= 0 && iel[j+3] == flag     && 
                         Nwall_owners[ilist][iel_box[j+1]] == 1 &&
                         fabs(WallPos[1][Wall_elems[ilist][iel_box[j+1]]] + 0.5*Size_x[1]) 
                                                                          > 0.00000001 ) ){

                       if (j == 0) type = LEFT_BACK;
                       else        type = LEFT_FRONT;
                                      esize1 = Esize_x[0];
                       if (Ndim == 3) esize2 = Esize_x[2];
                       idim = 1;

                       if (iel[j+1] >=0) {
                          el = j+1; normal = -2;
                          surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                            normal,idim,esize1,esize2);
                       }
                       if (iel[j+3] >=0) {
                          el = j+3; normal =  2; 
                          surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                            normal,idim,esize1,esize2);
                       }
                   }

                   /* check for surface element in right/left (on bottom elements)*/
                   if ( ( iel[j+2] >= 0 && iel[j+3]>= 0 && 
                        Nwall_owners[ilist][iel_box[j+2]] == 1   &&
                        Nwall_owners[ilist][iel_box[j+3]] == 1 &&
                        Wall_elems[ilist][iel_box[j+2]] != Wall_elems[ilist][iel_box[j+3]] )||

                        ( iel[j+2] == flag && iel[j+3]>= 0 && 
                         Nwall_owners[ilist][iel_box[j+3]] == 1 &&
                         fabs(0.5*Size_x[0] - WallPos[0][Wall_elems[ilist][iel_box[j+3]]]) 
                                                                          > 0.00000001 ) ||

                        ( iel[j+2] >= 0 && iel[j+3] == flag && 
                         Nwall_owners[ilist][iel_box[j+2]] == 1 &&
                         fabs(WallPos[0][Wall_elems[ilist][iel_box[j+2]]] + 0.5*Size_x[0]) 
                                                                          > 0.00000001 ) ){

                       if (j == 0) type = DOWN_BACK;
                       else        type = DOWN_FRONT;
                                      esize1 = Esize_x[1];
                       if (Ndim == 3) esize2 = Esize_x[2];
                       idim = 0;

                       if (iel[j+2] >=0){
                          el = j+2; normal = -1;
                          surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                            normal,idim,esize1,esize2);
                       }
                       if (iel[j+3] >=0){
                          el = j+3; normal =  1; 
                          surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                            normal,idim,esize1,esize2);
                       }
                   }
                }
                if (Ndim == 3){
                   for (i=0; i<4 ; i++){
                      j=i+4;

                      /* check for surface element in back/front */
                      if ( ( iel[i] >= 0 && iel[j]>= 0 && 
                          Nwall_owners[ilist][iel_box[i]] == 1   &&
                          Nwall_owners[ilist][iel_box[j]] == 1 &&
                          Wall_elems[ilist][iel_box[i]] != Wall_elems[ilist][iel_box[j]] ) ||

                        ( iel[i] == flag && iel[j]>= 0        && 
                         Nwall_owners[ilist][iel_box[j]] == 1 &&
                         fabs(0.5*Size_x[2] - WallPos[2][Wall_elems[ilist][iel_box[j]]]) 
                                                                          > 0.00000001 ) ||

                        ( iel[i] >= 0 && iel[j] == flag && 
                         Nwall_owners[ilist][iel_box[i]] == 1 &&
                         fabs(WallPos[2][Wall_elems[ilist][iel_box[i]]] + 0.5*Size_x[2]) 
                                                                          > 0.00000001 ) ){

                         type = i+8;
                                        esize1 = Esize_x[0];
                         if (Ndim == 3) esize2 = Esize_x[1];
                         idim = 2;
  
                         if(iel[i] >= 0) {
                            el = i; normal = -3;
                            surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                            normal,idim,esize1,esize2);
                         }
                         if(iel[j] >= 0) {
                            el = j; normal = 3; 
                            surf_el_to_list(loc_inode,ilist,iel_box,el,type,
                                            normal,idim,esize1,esize2);
                         }
                      }
                   }
   
                }  /* end of 3D case */
             }     /* end of 2D and 3D cases */

        } /* end of test for if this is a wall-wall boundary node */

    }     /* end of loop over processor nodes */

    for (iwall=0; iwall<Nwall; iwall++){
       s_area_tot_proc = 0.0;
       for (idim=0; idim<Ndim; idim++) {
           s_area_tot_proc += S_area[ilist][iwall][idim];
           S_area_tot[ilist][iwall] = gsum_double(s_area_tot_proc);
           S_area[ilist][iwall][idim] = gsum_double(S_area[ilist][iwall][idim]);
       }
    }
  }             /* end of loop over lists */
  safe_free((void *) &iel);

  if (Proc == 0 && Imain_loop>=0) {
     fprintf (fp1,"\n---------------------------------------------------------------\n");
     fprintf(fp1,"Done calculating surface normals of boundary elements : Iloop: %d\n",Imain_loop);
     fprintf(fp1,"\nProc: %d\n",Proc);
  }

  for (ilist=0; ilist<Nlists_HW; ilist++){
     for (idim=0; idim<Ndim; idim++) sum[idim]=0;
     if (Proc==0 && Imain_loop>=0) fprintf (fp1,"\nilist: %d\n",ilist); 
     for (loc_inode=0; loc_inode<Nnodes_per_proc;loc_inode++){

         inode = L2G_node[loc_inode];
         inode_box = node_to_node_box(inode);
         iwall = Nodes_2_boundary_wall[ilist][inode_box];

         node_to_ijk(inode,ijk);

         if (iwall != -1){      /* proceed if this is a boundary node */
            node_to_ijk(inode,ijk);
            if (Imain_loop>=0) {
                  for (iel_s=0; iel_s<Nelems_S[ilist][loc_inode]; iel_s++){
                     surf_norm = Surf_normal[ilist][loc_inode][iel_s]; 
                     idim = abs(surf_norm) - 1;
                     sum[idim] += surf_norm/abs(surf_norm);
                  }
            }

         }
      }
      for (idim=0; idim<Ndim; idim++){
      sum_all[idim] = gsum_double(sum[idim]);
        if (Proc == 0)
          fprintf(fp1,"ilist: %d Summing surface normals idim: %d yields: %d\n", 
                                                      ilist,idim,sum_all[idim]);
      }

  }
  if (Proc==0 && Imain_loop>=0) {
     fprintf (fp1,"\n---------------------------------------------------------------\n");
     fprintf(fp1,"Done calculating surface areas.Iloop: %d\n",Imain_loop);
     fprintf(fp1,"\nProc: %d\n",Proc);
  }

  if (Proc == 0 && Imain_loop>=0) {
     for (ilist=0; ilist<Nlists_HW; ilist++){
        fprintf (fp1,"\nilist: %d\n",ilist); 
        sarea_sum=0.0;
        for (iwall=0; iwall<Nwall; iwall++){
            fprintf (fp1,"\t iwall: %d \t S_area_tot[ilist][iwall]: %9.6f\n",
                                            iwall, S_area_tot[ilist][iwall]); 
            sarea_sum += S_area_tot[ilist][iwall];
            for (idim=0; idim<Ndim; idim++)
                 fprintf (fp1,"\t\t idim: %d \t S_area[ilist][iwall][idim]: %9.6f\n",
                                                  idim, S_area[ilist][iwall][idim]); 
       }
            fprintf (fp1,"\t total of walls: %d \t sarea_sum=%9.6f\n", ilist, sarea_sum);
     }
  }
  if (Proc ==0 && Imain_loop >=0) 
       fprintf (fp1,"\n---------------------------------------------------------------\n");

  safe_free ((void *) &iel_box);
  safe_free ((void *) &els_one_owner);
  return;
}
/****************************************************************/
/*find_local_els: in this routine, find all the elements in both
    global and box coordinates that surround a given node. */
void find_local_els(int inode,int *iel, int *iel_box,int flag)
{
  int loc_node_el,idim,ijk[3],i;
  int reflect_flag[3]; 
  reflect_flag[0]=reflect_flag[1]=reflect_flag[2] = FALSE;

  for (loc_node_el=0; loc_node_el<Nnodes_per_el_V; loc_node_el++)
     iel[loc_node_el] = node_to_elem(inode, loc_node_el,reflect_flag);

  node_to_ijk(inode,ijk);

  for (idim=0; idim<Ndim; idim++){
     if ((Type_bc[idim][0] == REFLECT || Type_bc[idim][0]==LAST_NODE) && ijk[idim]==0){
        if (idim == 0){
           if      (Ndim == 1) iel[1] = flag;
           else if (Ndim == 2) iel[1] = iel[3] = flag;
           else if (Ndim == 3) iel[1] = iel[3] = iel[5] = iel[7] = flag;
        }
        else if (idim == 1){
           if      (Ndim == 2) iel[2] = iel[3] = flag;
           else if (Ndim == 3) iel[2] = iel[3] = iel[6] = iel[7] = flag;
        }
        else if (idim == 2)    iel[4] = iel[5] = iel[6] = iel[7] = flag;
     }
     else if ((Type_bc[idim][1] == REFLECT|| Type_bc[idim][0]==LAST_NODE)  && ijk[idim]==Nodes_x[idim]-1){
        if (idim == 0){
           if      (Ndim == 1) iel[0] = flag;
           else if (Ndim == 2) iel[0] = iel[2] = flag;
           else if (Ndim == 3) iel[0] = iel[2] = iel[4] = iel[6] = flag;
        }
        else if (idim == 1){
           if      (Ndim == 2) iel[0] = iel[1] = flag;
           else if (Ndim == 3) iel[0] = iel[1] = iel[4] = iel[5] = flag;
        }
        else if (idim == 2)    iel[0] = iel[1] = iel[2] = iel[3] = flag;
     }
  }

  for (i=0; i<Nnodes_per_el_V; i++){
      if (iel[i] >=0) iel_box[i] = el_to_el_box(iel[i]);
      else            iel_box[i] = flag;
  }
  return;
}
/*****************************************************************/
/*surf_el_to_list: count this as a surface element:
    add entries for normals, surface area etc. !*/

void surf_el_to_list(int loc_inode, int ilist, int *iel_box,
		     int el, int type, int normal, int idim, 
		     double esize1,double esize2)
{
    int jwall, i;
    int inode,ijk[3];
    for (i=0; i< Nwall_owners[ilist][iel_box[el]];i++) {

       if (Nwall_owners[ilist][iel_box[el]] == 1) 
             jwall = Wall_elems[ilist][iel_box[el]];
       else
             jwall = Wall_owners[ilist][iel_box[el]][i];

       if (jwall < Nwall){

          inode = L2G_node[loc_inode];
          node_to_ijk(inode,ijk);

          Surf_normal[ilist][loc_inode][Nelems_S[ilist][loc_inode]] = normal;
          if (Ndim >1) Surf_elem_type[loc_inode][Nelems_S[ilist][loc_inode]] = type;
          Surf_elem_to_wall[ilist][loc_inode][Nelems_S[ilist][loc_inode]] = jwall;
          Nelems_S[ilist][loc_inode]++;

          if (Ndim == 1) S_area[ilist][jwall][idim] += 1.0;
          else if (Ndim == 2) S_area[ilist][jwall][idim] += esize1/Nnodes_per_el_S;
          else if (Ndim == 3) S_area[ilist][jwall][idim] += esize1*esize2/Nnodes_per_el_S;
 
       }
    }
    return;
}
/*********************************************************************
setup_surface_charge: assign a surface charge to each of the surface
                      boundary nodes based on the number of surface
                      elements this node contributes to.  This routine
                      is only set up for constant charge per unit area 
                      at the moment.                               */ 

void setup_surface_charge(FILE *fp1)
{
  int iwall,iwall_test,iwall_type, loc_inode, idim, inode_box;
  int ncharge_s,ncharge_v,iel;

  /* we need different arrays to store the charge that is attributed to a 
     smeared surface charge and a local volume charge.  They are added
     to the residuals differently in dft_fill_pde.c when filling Poisson's eqn */

  ncharge_s=ncharge_v=0;
  for (iwall_type=0; iwall_type<Nwall_type; iwall_type++){
       if (Type_bc_elec[iwall_type] == CONST_CHARGE) ncharge_s++;
       else if (Type_bc_elec[iwall_type] ==ATOMIC_CHARGE) ncharge_v++;
  }

  if (ncharge_s >0){
        Charge_w_sum_els = (double **) array_alloc(2, Nnodes_per_proc, 
                                                   Ndim, sizeof(double));
        for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++)
           for (idim=0; idim<Ndim; idim++) Charge_w_sum_els[loc_inode][idim] = 0.0;
  }
  if (ncharge_v >0 || Nlocal_charge != 0){
     Charge_vol_els = (double *) array_alloc(1, Nelements_box, sizeof(double));
     for (iel=0; iel<Nelements_box; iel++)
           Charge_vol_els[iel] = 0.0;

  }


  for (iwall=0; iwall<Nwall; iwall++){
     if (Type_bc_elec[WallType[iwall]] == CONST_CHARGE){

        for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
           inode_box = L2B_node[loc_inode];
           iwall_test = Nodes_2_boundary_wall[Nlists_HW-1][inode_box];
           if (iwall_test == iwall){
                bc_setup_const_charge(iwall,loc_inode);
                Surf_charge_flag=TRUE;
           }
        }
     }
     else if (Type_bc_elec[WallType[iwall]] == ATOMIC_CHARGE){
        setup_volume_charge1(iwall);
        Vol_charge_flag = TRUE;
     }
  }

  if (Nlocal_charge != 0) {
      if (Nlocal_charge == -1) setup_linear_grad_of_charge();
      setup_volume_charge2();
      Vol_charge_flag = TRUE;
  }

  if (Iwrite == VERBOSE){
      if (Proc==0) printf ("PRINTING CHARGE DISTRIBUTIONS: ncharge_v=%d ncharge_s=%d\n",ncharge_v,ncharge_s);
      if (ncharge_v>0 || Nlocal_charge !=0) print_charge_vol(Charge_vol_els,"dft_charge_vol.dat");
      if (ncharge_s>0) print_charge_surf(Charge_w_sum_els,"dft_charge_surf.dat");
  }


        
           /********** PRINTING **********/
/*  if (Proc == 0 && Imain_loop==0) {
     fprintf (fp1,"\n---------------------------------------------------------------\n");
     fprintf(fp1,"Done setting up the surface charges.\n");
     fprintf(fp1,"\nProc: %d\n",Proc);
  }

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
      charge_TF = FALSE;
      for (idim=0; idim<Ndim; idim++) 
         if (Charge_w_sum_els[loc_inode][idim] > 0.0)charge_TF = TRUE;
    
      if (charge_TF == TRUE) {
         inode = L2G_node[loc_inode];
         if (Proc == 0 && Imain_loop==0) {
            fprintf (fp1,"inode: %d \n", inode); 
            for (idim=0; idim<Ndim; idim++)
               fprintf (fp1,"\t idim: %d \t Charge_w_sum_els[inode][idim]: %9.6f\n",
                                           idim, Charge_w_sum_els[loc_inode][idim]); 
         }
      }
  }
  if (Proc == 0 && Imain_loop==0) 
     fprintf (fp1,"\n---------------------------------------------------------------\n"); 
*/
  return;
}

/****************************************************************************
 * bc_setup_const_charge: This routine sets up the array Charge_w_sum_els   *
 *                        which is used to set the boundary condition if    *
 *                        the constant charge case is considered.           */

void bc_setup_const_charge(int iwall, int loc_inode)
{
  int iel,idim,inode,ijk[3],jdim;
  double weight;

  inode=L2G_node[loc_inode];
  node_to_ijk(inode,ijk);

  for (iel = 0; iel<Nelems_S[Nlists_HW-1][loc_inode];iel++){
      if (Surf_normal[Nlists_HW-1][loc_inode][iel] !=0){
            idim = abs(Surf_normal[Nlists_HW-1][loc_inode][iel]) - 1;
            weight=1.0;
            for (jdim=0;jdim<Ndim;jdim++) {
                if (jdim != idim){
                   if ( (ijk[jdim]==0 && Type_bc[jdim][0]>=REFLECT) || 
                        (ijk[jdim]==Nodes_x[jdim]-1 && Type_bc[jdim][1]>=REFLECT)){
                        weight*=2.0;
                   }
                }
            }
            Charge_w_sum_els[loc_inode][idim] += 
                         weight*Elec_param_w[iwall]/(double)Nnodes_per_el_S;
      }
  }

  return;
}
/******************************************************************/
/*setup_linear_grad_of_charge: in this routine we set up a linear
      charge gradient along in an indicated direction between two 
      points*/
void setup_linear_grad_of_charge(void)
{

    int i,j,idim=0;
    double charge_1,charge_2,charge_x_1[3],charge_x_2[3];
    /* first determine direction of interest */
    for (i=0; i<Ndim; i++) 
       if (Charge_x[i][0] != Charge_x[i][1]) {
          idim = i; 
          Nlocal_charge = 1+ (int)(fabs(Charge_x[i][1]-Charge_x[i][0])/Esize_x[i]);
          break;
       }
    charge_1 = Charge[0];
    charge_2 = Charge[1];
    for (i=0; i<Ndim; i++){
          charge_x_1[i] = Charge_x[i][0];
          charge_x_2[i] = Charge_x[i][1];
    }
   safe_free ((void *) &Charge);
   safe_free ((void *) &Charge_x);
   safe_free ((void *) &Charge_Diam);

   Charge_Diam = (double *) array_alloc (1, Nlocal_charge, sizeof(double));
   Charge      = (double *) array_alloc (1, Nlocal_charge, sizeof(double));
   Charge_x    = (double **) array_alloc (2, Ndim, Nlocal_charge,sizeof(double));

   
   for (i=0; i<Nlocal_charge; i++){
        for (j=0; j<Ndim; j++) Charge_x[j][i] = charge_x_1[j];

        Charge[i] = (charge_1*(Nlocal_charge-1 - i) + charge_2*i)/(Nlocal_charge-1);
        Charge_x[idim][i] = (charge_x_1[idim]*(Nlocal_charge-1 - i) 
                                              + charge_x_2[idim]*i)/(Nlocal_charge-1);
        Charge_Diam[i] = 0.0;
   }

    return;
}
/****************************************************************************/
/*setup_volume_charge1:  This routine sets up charges on a per element
                    basis for atomic surfaces where certain atom types are
                    charged. */
void setup_volume_charge1(int iwall)
{

  int nmax,idim,iel,nelems,nelems_unique, *elems;
  double r,x[3],charge_per_el, esize_x_min = 10.;

  nmax = Nwall;
  for (idim=0; idim<Ndim; idim++) esize_x_min = AZ_MIN(Esize_x[idim],esize_x_min);

  elems = (int *) array_alloc(1, Nelements_box, sizeof(int));

  
  if (Charge_type_atoms == POINT_CHARGE || 
        Sigma_ww[WallType[iwall]][WallType[iwall]]<= esize_x_min) r = 0.6 * esize_x_min;
  else if (Surface_type[WallType[iwall]] == atomic_centers )  
            r = 0.5*Sigma_ww[WallType[iwall]][WallType[iwall]];
  else      r = WallParam[WallType[iwall]];

  for (idim=0; idim<Ndim; idim++) x[idim] = WallPos[idim][iwall];

  nelems = 0;
  nelems_unique = 0;
  els_charge_spheres(r,x,&nelems,&nelems_unique,elems,Charge_type_atoms);
  nelems_unique = gsum_int(nelems_unique);

  charge_per_el = Elec_param_w[iwall]/(double) nelems_unique;

  for (iel=0; iel<nelems; iel++){
      Charge_vol_els[elems[iel]]+=charge_per_el;
  }
  if (Proc==0 && Iwrite==VERBOSE){ 
        printf("iwall: %d size=%g has a charge of: %g being smeared over %d elements\n",
                         iwall, r, Elec_param_w[iwall],nelems_unique);
  }
  safe_free((void *) &elems);
  return;
}
/****************************************************************************/
/*setup_volume_charge2:  This routine sets up charges on a per element
                    basis for cases where charge sources are placed in 
                    the volume.  The number of walls in this case may be
                    zero.  */
void setup_volume_charge2(void)
{

  int i,nmax,idim,iel,nelems,nelems_unique, *elems;
  double r,x[3],charge_per_el, esize_x_min = 10.;

  nmax = Nlocal_charge;
  for (idim=0; idim<Ndim; idim++) esize_x_min = AZ_MIN(Esize_x[idim],esize_x_min);

  elems = (int *) array_alloc(1, Nelements_box, sizeof(int));

  for (i=0; i<nmax; i++){
      if (Charge_type_local==BACKGROUND){
          charge_per_el = Charge[i]/(double) Nelements;
          for (iel=0; iel<Nelements_box; iel++) Charge_vol_els[iel] += charge_per_el;
        if (Proc==0&&Iwrite!=NO_SCREEN) printf("a uniform background charge of %g total is being smeared over every element\n",
                                                                           Charge[i]);
      }
      else {
         if (Charge_type_local==POINT_CHARGE || Charge_Diam[i] == 0.0) {
             r = 0.5*esize_x_min + 0.0001*0.5*esize_x_min;
         }
         else if (Charge_type_local== SMEAR_CHARGE){           
             r = 0.5*Charge_Diam[i];
         }
         for (idim=0; idim<Ndim; idim++) x[idim] = Charge_x[idim][i];

         nelems = 0;
         nelems_unique = 0;
         els_charge_spheres(r,x,&nelems,&nelems_unique,elems,Charge_type_local);
         nelems_unique = gsum_int(nelems_unique);

         charge_per_el = Charge[i]/(double) nelems_unique;

         for (iel=0; iel<nelems; iel++){
            Charge_vol_els[elems[iel]] += charge_per_el;
         }
         if (Proc==0 && Iwrite==VERBOSE) 
              printf(" a charge of: %g being smeared over %d elements\n", 
                                                         Charge[i],nelems_unique);
      }
  }
  safe_free((void *) &elems);
  return;
}
/****************************************************************************/
/* els_charge_spheres: Assign each element to either the fluid
                                      or the fixed surface spheres */

void els_charge_spheres(double radius,double *x,int *nelems,
                        int *nelems_unique, int *elems,int charge_type) 
{
   int iel_box,iel,inode,inode_box,idim,iel_save=0,inode_save=0;
   double node_pos[3],xtest[3];
   double r12_sq_sum,r12,dx,r12min=1000.,esize_x_min=1000.;
   int image_flag, ijk[3], ijk_box[3],image_flag_save=0;

   for (idim=0; idim<Ndim; idim++) esize_x_min = AZ_MIN(Esize_x[idim],esize_x_min);
   for (iel_box = 0; iel_box < Nelements_box; iel_box++){

       inode_box = element_box_to_node_box(iel_box);

       /* check if this box--node is outside of real domain */

       image_flag = FALSE;
       node_box_to_ijk_box(inode_box, ijk_box);
         for (idim=0; idim<Ndim; idim++) {
           ijk[idim] = ijk_box[idim] + Min_IJK_box[idim];

           if (Type_bc[idim][0] == PERIODIC && ijk[idim] < 0)
              image_flag = TRUE;
           if (Type_bc[idim][1] == PERIODIC && ijk[idim] >= Nodes_x[idim])
              image_flag = TRUE;
         }

       
       iel=el_box_to_el(iel_box);
       inode = element_to_node(iel);
       node_to_position(inode,node_pos);

       switch(Ndim)   /* xtest = position at center of element */
       {
          case 3:  xtest[2] = node_pos[2] + 0.5*Esize_x[2];
          case 2:  xtest[1] = node_pos[1] + 0.5*Esize_x[1];
          default: xtest[0] = node_pos[0] + 0.5*Esize_x[0];
       }

       r12_sq_sum = 0.0;
       for (idim = 0; idim < Ndim; idim++) {
          dx =  xtest[idim] -  x[idim];
          r12_sq_sum = r12_sq_sum + dx*dx;
       }
       r12 = sqrt(r12_sq_sum);
       if (r12<r12min){
             r12min=r12;
             iel_save=iel_box;
             inode_save=inode_box;
             image_flag_save=image_flag;
       }

       if (r12 <= radius && charge_type != POINT_CHARGE) {
           elems[(*nelems)++] = iel_box;
           if ((B2L_node[inode_box] != -1) && !image_flag) (*nelems_unique)++;
       }
    }     /* end of loop over elements       */

    /* if no elements were flagged or we are doing point charges, pile all the
       charge into one element */

    if ( (*nelems == 0 && r12min < esize_x_min) || charge_type==POINT_CHARGE){
        (*nelems)=1;
        elems[0]=iel_save;
        if ((B2L_node[inode_save] != -1) && !image_flag_save) (*nelems_unique)=1;
    }
}
/****************************************************************************/
/* zones_el_to_nodes:  Assigns a zone to each of the nodes in the problem. */

void zones_el_to_nodes(int *elem_zones)
{
  int iel_box, inode_box, idim, izone, ijk_box[3], ijk_plus[3], ijk_tmp[3];
  int nzone_max;
  Nodes_to_zone = (int *) array_alloc (1, Nnodes_box, sizeof(int));

  if (Coarser_jac != 5) nzone_max=Nzone-1;
  else                  nzone_max=Nzone-2;
  
  for (izone = nzone_max; izone>=0; izone--){
     for (iel_box = 0; iel_box < Nelements_box; iel_box++){
        if (elem_zones[iel_box] == izone) {

           inode_box = element_box_to_node_box(iel_box);
           Nodes_to_zone[inode_box] = izone;

           node_box_to_ijk_box(inode_box,ijk_box);
           for (idim=0; idim<Ndim; idim++) {
             ijk_tmp[idim] = ijk_box[idim];
             if (ijk_box[idim]+1 < Nodes_x_box[idim] || (
                 ijk_box[idim]+1 == Nodes_x_box[idim] &&
                 Max_IJK_box[idim]+1 !=Nodes_x[idim])) 
                                       ijk_plus[idim] = ijk_box[idim]+1;
             else                      ijk_plus[idim] = 0;
           }

           ijk_tmp[0] = ijk_plus[0];
           Nodes_to_zone[ijk_box_to_node_box(ijk_tmp)] = izone;

           if (Ndim > 1) {
               ijk_tmp[1] = ijk_plus[1];

               ijk_tmp[0] = ijk_box[0];
               Nodes_to_zone[ijk_box_to_node_box(ijk_tmp)] = izone;

               ijk_tmp[0] = ijk_plus[0];
               Nodes_to_zone[ijk_box_to_node_box(ijk_tmp)] = izone;
           }

           if (Ndim==3) {

               ijk_tmp[2] = ijk_plus[2];

               ijk_tmp[0] = ijk_box[0];
               ijk_tmp[1] = ijk_box[1];
               Nodes_to_zone[ijk_box_to_node_box(ijk_tmp)] = izone;

               ijk_tmp[0] = ijk_plus[0];
               Nodes_to_zone[ijk_box_to_node_box(ijk_tmp)] = izone;

               ijk_tmp[0] = ijk_box[0];
               ijk_tmp[1] = ijk_plus[1];
               Nodes_to_zone[ijk_box_to_node_box(ijk_tmp)] = izone;

               ijk_tmp[0] = ijk_plus[0];
               Nodes_to_zone[ijk_box_to_node_box(ijk_tmp)] = izone;

           }
          
       }
     }
  }
  if (Iwrite == VERBOSE) print_Nodes_to_zone(Nodes_to_zone,"dft_zones.dat");
}
/****************************************************************************/
void set_mesh_coarsen_flag(void)

/* Sets up flags for mesh coarsening. If Mesh_coarsen_flag[i] >= 0,
 * then it is the "izone" for the residual calculation. If < 0, then
 * it is the dimension which to average the neighbors values
 */

{
  /* coarse_fac is the power-of-2 coarsening of the mesh */
  int i, inode, coarse_fac, ijk[3],count,count_coarse,icomp,nodes_coarse;

  Mesh_coarsen_flag = (int *) array_alloc(1, Nnodes_box, sizeof(int));
  List_coarse_nodes = (int *) array_alloc(1, Nnodes_per_proc, sizeof(int));

  /* set mesh coarsening flag */
  count=0;count_coarse=0;

  ijk[0] = ijk[1] = ijk[2] = 0;

  for (i=0; i < Nnodes_box; i++) {

    /* Set to zone number for coarse Residual calc */

      Mesh_coarsen_flag[i] = Nodes_to_zone[i];

      inode = B2G_node[i];
      node_to_ijk(inode, ijk);

      if (Mesh_coarsening != FALSE && Nodes_to_zone[i] > 0){

         if (Mesh_coarsening == BULK_ZONE){
                Mesh_coarsen_flag[i] = FLAG_BULK;
                if (B2L_node[i] >=0) count_coarse++;
         }
         else if (Mesh_coarsening == PB_ZONE){
                Mesh_coarsen_flag[i] = FLAG_PBELEC;
                if (B2L_node[i] >=0) count_coarse++;
         }
         else{
                      /* reset to negative flag if residual is not even to be set */
         coarse_fac = POW_INT(2,Nodes_to_zone[i]);
         if      (ijk[0]%coarse_fac) Mesh_coarsen_flag[i] = -1;
         else if (ijk[1]%coarse_fac) Mesh_coarsen_flag[i] = -2;
         else if (ijk[2]%coarse_fac) Mesh_coarsen_flag[i] = -3;
         }
         if (B2L_node[i] >=0 && Mesh_coarsen_flag[i] < 0) List_coarse_nodes[count_coarse++]=B2G_node[i];
               
      }
      if (L1D_bc){
         if (ijk[Grad_dim]*Esize_x[Grad_dim] <= X_1D_bc+0.00000001 ||
             ijk[Grad_dim]*Esize_x[Grad_dim] >= Size_x[Grad_dim] - (X_1D_bc+0.00000001) ) {
             Mesh_coarsen_flag[i] = FLAG_1DBC;
             if (B2L_node[i] >=0) List_coarse_nodes[count++]=B2L_node[i];
          }
      }
  }
  if (L1D_bc){
     Nnodes_coarse_loc=count;
     nodes_coarse=gsum_int(count);
     if (Proc==0&&Iwrite!=NO_SCREEN) printf(" %d nodes of %d total will be set to the 1D boundary region\n",nodes_coarse,Nnodes);
  }
  else{
     if (Mesh_coarsening==BULK_ZONE || Mesh_coarsening==PB_ZONE) Nnodes_coarse_loc=0;
     else                                                        Nnodes_coarse_loc=count_coarse;
     nodes_coarse=gsum_int(count_coarse);
     if (Proc==0&&Iwrite!=NO_SCREEN) printf(" %d nodes of %d total will be coarsened\n",nodes_coarse,Nnodes);
  }

}
/***********************************************************
setup_area_IC:  Set up array of area as a function of of
                distance down the length of the pore.  Allows
	        some flexibility in 1D ion channel models */
void setup_area_IC(void)
{
  int inode_box, inode,idim,i;
  double nodepos[3],rad,xleft,xright;
  char *filename = "Area_IC.dat";
  FILE *fp=NULL;

  if (Iwrite == VERBOSE && Proc==0) fp=fopen(filename,"w");

  idim = Grad_dim;
  for (inode_box=0; inode_box <Nnodes_box; inode_box++){
    inode = B2G_node[inode_box];
    node_to_position(inode,nodepos); 

    switch(Geom_flag)
    {
       case OPTION_CYL: 
          Area_IC[inode_box] = PI*Pore_rad_L_IC[0]*Pore_rad_L_IC[0];
          break;

       case OPTION_VARY:
          if (nodepos[idim]<= -0.5*Size_x[idim]+ X_const_mu+ 0.00000001)  
               Area_IC[inode_box] = PI*Pore_rad_L_IC[0]*Pore_rad_L_IC[0];
          else if(nodepos[idim] >=  0.5*Size_x[idim]-X_const_mu-0.00000001)
               Area_IC[inode_box] = PI*Pore_rad_R_IC[Nseg_IC-1]
                                      *Pore_rad_R_IC[Nseg_IC-1];
          else {
             xleft = -0.5*Size_x[idim] + X_const_mu;
             for (i=0; i<Nseg_IC; i++){
                xright = xleft + Lseg_IC[i];
                if(  nodepos[idim]> xleft - 0.00000001 &&  
                     nodepos[idim]< xright+0.00000001){
                   rad = Pore_rad_L_IC[i] + (nodepos[idim]-xleft)*
                         (Pore_rad_R_IC[i]-Pore_rad_L_IC[i])/(xright-xleft);
                   Area_IC[inode_box] = PI*rad*rad;
                }
                xleft = xright;
             }
          }
          break;

       default:
          Area_IC[inode_box] = 1.0;
          break;
    }
    if (Iwrite==VERBOSE && Proc==0) 
        fprintf(fp," %d  %9.6f \n", inode_box,Area_IC[inode_box]);
  }
  if (Iwrite == VERBOSE && Proc==0) fclose(fp);
  return;
}
/****************************************************************************/
/* bc_setup_local_charge: This routine sets up the arrays Charge_w and      *
 * El_Area for charged walls where the surface charge is smeared only       *
 * locally over elements near a given charge.                               */

/*void bc_setup_local_charge(void)
*{
*  int iel,match[3],idim,i,ielmax,iwall,icharge,n_centers;
* int inode,ijk[3],**ijk_node,jnode,ijk_node_min[3];
* int count,jwall;
* double node_pos[3],rsq,***el_center,charge_area,rsqmin,center_min[3];

* n_centers = Nnodes_b[Nlists_HW-1];
* el_center = (double ***) array_alloc (3, Nwall, n_centers, Ndim, sizeof(double)); 
* ijk_node = (int **) array_alloc(2,4,3,sizeof(int));

* for (iwall=0; iwall<Nwall; iwall++) {
*    ielmax = Surf_elem_pos[Nlists_HW-1][iwall]+Surf_elem_neg[Nlists_HW-1][iwall];
*    for (iel=0; iel<ielmax; iel++) {
*         El_Area[iwall][iel] = 1.0;
*         Charge_w[iwall][iel] = 0.0;
*    }
* }
*
* for (icharge=0; icharge < Nlocal_charge; icharge++) { 
*   charge_area = 0.0;

*   rsqmin = 20.0;   * a big number .... we're looking for something small * 
*   for (iwall=0; iwall<Nwall; iwall++) { 

*      *
*      * Loop over the surface elements ... assign element areas
*      * and search for the element nearest to the charge.
*      * 
*     ielmax = Surf_elem_pos[Nlists_HW-1][iwall]+Surf_elem_neg[Nlists_HW-1][iwall];
*     for (iel=0; iel<ielmax; iel++){ 

*       if (iel < Surf_elem_pos[Nlists_HW-1][iwall]){
*           count = 0;
*           for (jwall=0; jwall<iwall; jwall++) count += Surf_elem_pos[Nlists_HW-1][jwall];
*           find_match_dim(iwall,count+iel,ijk_node,match,Surf_nodes_pos);
*       }
*       else  {
*           i = iel - Surf_elem_pos[Nlists_HW-1][iwall];
*           count = 0;
*           for (jwall=0; jwall<iwall; jwall++) count += Surf_elem_neg[Nlists_HW-1][jwall];
*           find_match_dim(iwall,count+i,ijk_node,match,Surf_nodes_neg);
*       }

*       find_lbb_node(ijk_node,node_pos);

*        * 
*        * Find the element area and the center of the surface element.
*        * And locate the element that is nearest to the point of charge
*        * 
*       rsq = 0.0;
*       for (idim=0; idim<Ndim; idim++){
*          if (match[idim]==FALSE) {
*             El_Area[iwall][iel] *= Esize_x[idim];
*             el_center[iwall][iel][idim] = node_pos[idim] + 0.5*Esize_x[idim];
*          }
*          else el_center[iwall][iel][idim] = node_pos[idim];
*          rsq += (el_center[iwall][iel][idim] - Charge_x[icharge][idim])*
*                 (el_center[iwall][iel][idim] - Charge_x[icharge][idim]);
*       }
*       if (rsq < rsqmin) {
*          for (idim=0; idim<Ndim; idim++) 
*               center_min[idim] = el_center[iwall][iel][idim];
*          rsqmin = rsq;
*       }
*     }  * end of loop over surface elements * 
*   }   * end of loop over iwall * 

*   printf("\n \n center_min: %f %f ",center_min[0],center_min[1]);

*    *
*    * Now find the total area over which the charge is to be spread.
*    * Increment counter n_charges[iwall][iel];
*    * 
*   for (iwall=0; iwall<Nwall; iwall++) { 
*     ielmax = Surf_elem_pos[Nlists_HW-1][iwall]+Surf_elem_neg[Nlists_HW-1][iwall];
*     for (iel=0; iel<ielmax; iel++){ 
*        rsq = 0.0;
*        for (idim = 0; idim<Ndim; idim++)
*           rsq +=  (el_center[iwall][iel][idim] - center_min[idim]) 
*                  *(el_center[iwall][iel][idim] - center_min[idim]);

*        printf("\n iel: %d  rsq: %f",iel,rsq);
*        if (rsq <= Charge_Diam[icharge]*Charge_Diam[icharge])
*           charge_area += El_Area[iwall][iel];

*     }  * end of loop over surface elements * 
*   }   * end of loop over iwall * 

*   printf("\n charge_area: %f",charge_area);
*   printf("\n charge_diam^sq: %f",Charge_Diam[icharge]*Charge_Diam[icharge]);

*    *
*    *Finally assign the charge to each element !!
*    * 
*   for (iwall=0; iwall<Nwall; iwall++) { 
*     ielmax = Surf_elem_pos[Nlists_HW-1][iwall]+Surf_elem_neg[Nlists_HW-1][iwall];
*     for (iel=0; iel<ielmax; iel++){ 
*        rsq = 0.0;
*        for (idim = 0; idim<Ndim; idim++)
*           rsq +=  (el_center[iwall][iel][idim] - center_min[idim]) 
*                  *(el_center[iwall][iel][idim] - center_min[idim]);

*        if (rsq <= Charge_Diam[icharge]*Charge_Diam[icharge])
*           Charge_w[iwall][iel] += (Charge_loc[icharge]/charge_area)
*                                                *El_Area[iwall][iel];
*           printf("\n iel: %d charge_w[iwall][iel]: %f",iel,
*                  Charge_w[iwall][iel]);

*     }  * end of loop over surface elements * 
*   }   * end of loop over iwall * 

* }       * end of loop over icharge * 

* safe_free((void *) &el_center);
* safe_free((void *) &ijk_node);
*  
*/
/****************************************************************************/
/* find_match_dim(void);given a wall and surface element number, this routine
                    determines in which dimensions the positions do not 
                    match.  Esize_x in these dimensions must be multiplied
                    to determine the surface area of the element.           */
/*void find_match_dim(int iwall,int i, int **ijk_node,int *match,
                                             int ***surf_nodes)
{
 int idim,inode,jnode,ijk[3];


 for (inode=0; inode<Nnodes_per_el_S; inode++){
    node_to_ijk(surf_nodes[Nlists_HW-1][i][inode],ijk);
    for (idim=0; idim<Ndim; idim++) ijk_node[inode][idim] = ijk[idim];
 }

 for (idim=0; idim<Ndim; idim++) match[idim]=TRUE;

 for (inode=0; inode<Nnodes_per_el_S-1; inode++)
    for (jnode=inode+1; jnode<Nnodes_per_el_S; jnode++)
       for (idim=0; idim<Ndim; idim++)
           if (ijk_node[inode][idim] != ijk_node[jnode][idim]) match[idim]=FALSE;
}
*/
/****************************************************************************/
/*find_lbb_node: Given a set of surface nodes ijk_node[][], return the
                 position of the one node that is in the 
                 left bottom back corner. */
/*void find_lbb_node(int **ijk_node,double *node_pos)
{
 int inode,idim,ijk_node_lbb[3],logical,jnode;

 logical = FALSE;
 for (inode=0; inode<(Nnodes_per_el_S-1); inode++)
    for (jnode=inode; jnode<(Nnodes_per_el_S); jnode++)
       for (idim=0; idim<Ndim; idim++) 
           if (abs(ijk_node[inode][idim] - ijk_node[jnode][idim]) > 1)
             logical = TRUE;

 if (logical) for (idim=0; idim<Ndim; idim++) ijk_node_lbb[idim] = -1;
 else         for (idim=0; idim<Ndim; idim++) ijk_node_lbb[idim] = 1e6;

 for (inode=0; inode<Nnodes_per_el_S; inode++)
    for (idim=0; idim<Ndim; idim++) 
        if ((ijk_node[inode][idim] < ijk_node_lbb[idim] && logical == FALSE) 
          ||( ijk_node[inode][idim]> ijk_node_lbb[idim] && logical == TRUE) )  
          ijk_node_lbb[idim] = ijk_node[inode][idim];
    
 
 inode = ijk_to_node(ijk_node_lbb);
 node_to_position(inode,node_pos);
}
*/
/****************************************************************************/
void initialize_Aztec(int* N_update, int *update[])
/*
 * This sets up various options for Aztec, the linear solver package.
 * Most are defaults. All info is stored in the structure Aztec.
 */
{
  /* Define partitioning:  matrix rows (ascending order) owned by this node */
  /* Partitioning Nnodes groups of Nunk_per_node                            */
  int flag;
 
/*  if (Load_Bal_Flag == LB_LINEAR || Ndim != 2) flag = AZ_linear;
  else                            flag = AZ_box;*/
  flag = AZ_linear;

  MY_read_update(N_update, update, Nnodes, Nodes_x, Nunk_per_node, flag);

 /*
  * Set up linear solver options. Some of these can be moved to the
  * input file for more user control.
  */

  AZ_defaults(Aztec.options, Aztec.params);

  switch (Az_solver) {
    case 1:  Aztec.options[AZ_solver]   = AZ_cg; break;
    case 2:  Aztec.options[AZ_solver]   = AZ_tfqmr; break;
    case 3:  Aztec.options[AZ_solver]   = AZ_cgs; break;
    case 4:  Aztec.options[AZ_solver]   = AZ_bicgstab; break;
    default: Aztec.options[AZ_solver]   = AZ_gmres;
  }
  switch (Az_scaling) {
    case 1:  Aztec.options[AZ_scaling]   = AZ_Jacobi; break;
    case 2:  Aztec.options[AZ_scaling]   = AZ_sym_row_sum; break;
    case -1: Aztec.options[AZ_scaling]   = AZ_none; break;
    default: Aztec.options[AZ_scaling]   = AZ_row_sum;
  }
  switch (Az_preconditioner) {
    case 1:  Aztec.options[AZ_precond]   = AZ_Jacobi; break;
    case 2:  Aztec.options[AZ_precond]   = AZ_sym_GS; break;
    case 3:  Aztec.options[AZ_precond]   = AZ_ls; break;

    case 4:  Aztec.options[AZ_precond]   = AZ_dom_decomp; 
             Aztec.options[AZ_subdomain_solve]=AZ_ilut; 
             Aztec.params[AZ_ilut_fill]  = Az_ilut_fill_param;    break;

    case 5:  Aztec.options[AZ_precond]   = AZ_dom_decomp; 
             Aztec.options[AZ_subdomain_solve]=AZ_ilut; 
             Aztec.params[AZ_ilut_fill]  = Az_ilut_fill_param;  
	     /* Parameters to improve condition numer of preconditioner */
             Aztec.params[AZ_athresh]  = 1.0e-5;  
             Aztec.params[AZ_rthresh]  = 1.01;    break;

    case -1: Aztec.options[AZ_precond]   = AZ_none; break;
    default: Aztec.options[AZ_precond]   = AZ_dom_decomp;
             Aztec.options[AZ_subdomain_solve]   = AZ_ilu;
  }
  Aztec.params[AZ_tol]  = Az_tolerance; /* This is linear solver convergence */
                                        /* criterion (often called eta_k)    */

  Aztec.options[AZ_conv]     = AZ_r0;
  if (Iwrite==NO_SCREEN)  Aztec.options[AZ_output]   = 0;  /* no output */
  else  Aztec.options[AZ_output]   = 10;  /* lots of output */
  if (Iwrite == NO_SCREEN)   Aztec.options[AZ_output]   = AZ_none;
  Aztec.options[AZ_pre_calc] = AZ_calc;
  Aztec.options[AZ_max_iter] = Max_gmres_iter;
  Aztec.options[AZ_poly_ord] = 3;
  Aztec.options[AZ_overlap]  = AZ_none;
  if (Az_kspace < 1) Aztec.options[AZ_kspace]   = Max_gmres_iter;
  else               Aztec.options[AZ_kspace]   = Az_kspace;
  Aztec.options[AZ_orthog]   = AZ_classic;
  Aztec.options[AZ_aux_vec]  = AZ_resid;

  Aztec.options[AZ_keep_info] = TRUE;

  Aztec.params[AZ_drop] = 0.0;
}
/****************************************************************************/
/****************************************************************************/

void MY_read_update(int *N_update, int *update[],
                    int N, int *nodes_x, int chunk, int input_option)

/*******************************************************************************

  This routine initializes update[] to the global indices updated by this
  processor and initializes N_update to the total number of elements to be
  updated.

  If input_option == AZ_linear Do a linear partitioning of the chunks.
     Specifically, proc 0 is assigned the first floor( (N+P-1)/P ) chunks,
     processor 1 is assigned the next floor( (N+P-2)/P ) chunks.
  If input_option == AZ_file Processor 0 reads the file '.update'.  This file
     should contain nprocs lists.  Each list consists of a number telling how
     many global indices are in the list followed by a list of global indices.
     The first list is then sent to processor 'nprocs-1', the second list is
     sent to processor 'nprocs-2', etc.
  If input_option == AZ_box
     we do a box partitioning of the unknowns (see comments below).

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N_update:        On Output, number of unknowns updated by this processor.

  update:          On Output, list of unknowns updated by this processor in
                   ascending order.

  N:               Total number of chunks to be distributed.

  chunk:           Size of each chunk to be treated as a single unit.
                   The unknowns contained in the kth chunk are given
                   by {k*chunk, k*chunk + 1, ..... , (k+1)*chunk - 1}
                   and 'N*chunk' is the total number of unknowns to be
                   distributed.

  input_option:    AZ_linear   ==> perform linear partitioning
                   AZ_file     ==> read partioning from file '.update'
                   AZ_box      ==> perform a box partitioning.

*******************************************************************************/

{

  /* local variables */

  int   t1, t2, i;
  int   ii, j;
  int   proc_x, proc_y;
  int   npts;
  int   pts_x, pts_y;
  int   px, py;
  int   start_x, start_y;
  int   end_x, end_y;
  int   pt_number;
  int   count;
  int   proc, nprocs;


  /**************************** execution begins ******************************/

  proc   = Proc;
  nprocs = Num_Proc;

  /*
   * Figure out which chunks should be assigned to this processor using a box
   * decomposition.  That is, it is assumed that all the chunks are ordered
   * naturally corresponding to an m x m x m box where m = N^(1/3).  Boxes of
   * chunks are assigned to processors.
   *
   * NOTE: it is assumed that nprocs = 2^power and that the number of chunks in
   * each direction is divisible by the number of processors in each direction.
   */
  if (input_option == AZ_box) {


    /* find approx number of nodes on a side of each proc's box */

    npts = pow( (double) N / (double) nprocs, 1.0 / (double) Ndim);

    if (Ndim == 2) {

      proc_x = nodes_x[0] / npts;
      if (proc_x == 0) proc_x = 1;
      if ((nprocs % proc_x !=0) && (nprocs % (proc_x + 1) ==0)) proc_x++;

      proc_y = nprocs / proc_x;
      proc_x = nprocs / proc_y; /* maybe extra procs can be used...*/

      pts_x = nodes_x[0] / proc_x;

      if (proc < proc_x*proc_y) {
        pts_y  = nodes_x[1] / proc_y;

        px = proc % proc_x;
        py = (proc / proc_x) % proc_y;

        start_x = px * pts_x;
        end_x   = start_x + pts_x;
        start_y = py * pts_y;
        end_y   = start_y + pts_y;

        /* divy up remainder of nodes  */
        if (nodes_x[0] % proc_x >= proc_x - px) {
           start_x += (nodes_x[0] % proc_x) - (proc_x - px);
           end_x   += (nodes_x[0] % proc_x) - (proc_x - px) + 1;
        }
        if (nodes_x[1] % proc_y >= proc_y - py) {
           start_y += (nodes_x[1] % proc_y) - (proc_y - py);
           end_y   += (nodes_x[1] % proc_y) - (proc_y - py) + 1;
        }


        *N_update = chunk * (end_x - start_x) * (end_y - start_y);

/*
printf("Proc %d nx %d ny %d px %d py %d st_x %d e_x %d st_y %d e_y %d N_update %d\n",
         proc, nodes_x[0], nodes_x[1], px, py, start_x, end_x-1, start_y, end_y-1, *N_update);
*/

        /* if (!AZ_using_fortran) */
              *update     = (int *) calloc(*N_update, sizeof(int));

        /* set update[] */

        count = 0;
        for (j = start_y; j < end_y; j++ ) {
          for (i = start_x; i < end_x; i++ ) {
            for (ii = 0; ii < chunk; ii++ ) {
              pt_number = (i + j * nodes_x[0]) * chunk + ii;
              (*update)[count++] = pt_number;
            }
          }
        }
      }
      else {
        *N_update = 0;
        *update = NULL;
      }
    }
    else {
     printf("AZ_BOX NOT YET SET FOR NDIM=3 :  %d\n", Ndim);
     exit(-1);
    }
  }
  else if (input_option == AZ_linear) {

    /*
     * Figure out which chunks should be assigned to this processor for linear
     * partitioning.  This means that processor 0 is assigned the chunks
     * approximately corresponding to 0, ... , N/nprocs and processor 1 is
     * approximately assigned the chunks 1+N/nprocs to 2*N/nprocs.
     */

    t1 = N/nprocs;
    t2 = N - t1 * nprocs;

    if ( proc >= t2) t2 += (proc * t1);
    else {
      t1++;
      t2    = proc*t1;
    }
    *N_update = t1*chunk;
    t2   *= chunk;

    /*if (!AZ_using_fortran) */
      *update = (int *) calloc(*N_update,sizeof(int));
    if ( (*update == NULL) && (*N_update != 0)) {
      (void) fprintf (stderr, "Not enough space to allocate 'update'\n");
      exit(-1);
    }

    if (MATRIX_FILL_NODAL){
        for (i = 0; i < *N_update; i++) (*update)[i] = i + t2;
    }
    else{
                 /* t1=numNodesThisProc, (i+t2/chunk)=global node number */
     for (ii = 0; ii < chunk; ii++)
       for (i = 0; i < t1; i++)
         (*update)[i + ii*t1] = (i + t2/chunk) + ii*N;
    }

  }
  else
    (void) fprintf(stderr,"Unknown input option (%d) in MY_read_update()\n",
                   input_option);

} /* MY_read_update */
