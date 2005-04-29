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
 *  FILE: dft_output.c
 *
 *  This file contains routines that post-process the density profiles.
 *
 */

#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

/*******************************************************************************
collect_x_old: This gathers all of the densities into X_old on proc 0.        */ 

void collect_x_old(double **x)
{
  int i,iunk,loc_inode, loc_i,idim,nunk_per_proc;
  int *index=NULL;
  double *unk_global, *unk_loc;

  Nodes_old = Nnodes;
  for (idim=0; idim<Ndim; idim++) Nodes_x_old[idim] = Nodes_x[idim];

  /* allocate temporary arrays */
  /* X_old gets allocated here and will be freed in dft_guess.c */

  nunk_per_proc = Nnodes_per_proc*Nunk_per_node;
  unk_loc = (double *) array_alloc (1, nunk_per_proc, sizeof(double));

  for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ )
     for (iunk=0; iunk<Nunk_per_node; iunk++){
     unk_loc[iunk+Nunk_per_node*loc_inode] = x[iunk][L2B_node[loc_inode]];  /* always use nodal ordering here */
  }

  if (Proc == 0) {
    unk_global = (double *) array_alloc (1, Nunknowns, sizeof(double));
    index = (int *) array_alloc (1, Nnodes, sizeof(int));
  }
  else {
    unk_global=NULL;
    index=NULL;
  }

  /* collect the node numbers from all the processors */

  MPI_Gatherv(L2G_node,Nnodes_per_proc,MPI_INT,
              index,Comm_node_proc,Comm_offset_node,
              MPI_INT,0,MPI_COMM_WORLD);

  /* collect the unknowns from all the processors */

  MPI_Gatherv(unk_loc,nunk_per_proc,MPI_DOUBLE,
              unk_global,Comm_unk_proc,Comm_offset_unk,
              MPI_DOUBLE,0,MPI_COMM_WORLD);
  safe_free((void *) &unk_loc);

  if (Proc == 0){
     for (i=0; i<Nnodes; i++){
        for (iunk=0; iunk<Nunk_per_node; iunk++)
           X_old[index[i]*Nunk_per_node+iunk] = unk_global[i*Nunk_per_node+iunk];
     }
     safe_free((void *) &unk_global);
     safe_free((void *) &index);
  }
  safe_free((void *) &unk_loc);

  return;
}
/*******************************************************************************
collect_vext_old: This gathers all of the external field into Vext_old on proc 0   */ 

void collect_vext_old()
{
  int i,iunk,loc_inode,nunk_per_proc;
  int *index=NULL;
  int *comm_comp_proc_tmp,*comm_offset_comp_tmp;
  double *unk_global, *unk_loc;

  /* allocate temporary arrays */
  /* X_old gets allocated here and will be freed in dft_guess.c */

  nunk_per_proc = Nnodes_per_proc*Ncomp;
  unk_loc = (double *) array_alloc (1, nunk_per_proc, sizeof(double));

  for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ )
     for (iunk=0; iunk<Ncomp; iunk++){
     unk_loc[loc_inode*Ncomp+iunk] = Vext[loc_inode][iunk];
  }


  if (Proc == 0) {
    unk_global = (double *) array_alloc (1, Nnodes*Ncomp, sizeof(double));
    index = (int *) array_alloc (1, Nnodes, sizeof(int));
  }
  else {
    unk_global=NULL;
    index=NULL;
  }

  /* collect the node numbers from all the processors */
 

  MPI_Gatherv(L2G_node,Nnodes_per_proc,MPI_INT,
              index,Comm_node_proc,Comm_offset_node,
              MPI_INT,0,MPI_COMM_WORLD);

  /* collect the unknowns from all the processors */

  
  comm_comp_proc_tmp = (int *) array_alloc (1, Num_Proc, sizeof(int));
  comm_offset_comp_tmp = (int *) array_alloc (1, Num_Proc, sizeof(int));
  for (i=0; i<Num_Proc; i++){ comm_comp_proc_tmp[i] = 0; comm_offset_comp_tmp[i] = 0;}
  if (Proc==0){
     for (i=0; i<Num_Proc; i++){
        comm_comp_proc_tmp[i] = Comm_unk_proc[i]*Ncomp/Nunk_per_node;
        comm_offset_comp_tmp[i] = Comm_offset_unk[i]*Ncomp/Nunk_per_node;
     }
  }


  MPI_Gatherv(unk_loc,nunk_per_proc,MPI_DOUBLE,
              unk_global,comm_comp_proc_tmp,comm_offset_comp_tmp,
              MPI_DOUBLE,0,MPI_COMM_WORLD);
  safe_free((void *) &unk_loc);
  safe_free((void *) &comm_comp_proc_tmp);
  safe_free((void *) &comm_offset_comp_tmp);

  if (Proc == 0){
     for (i=0; i<Nnodes; i++){
        for (iunk=0; iunk<Ncomp; iunk++)
           Vext_old[index[i]*Ncomp+iunk] = unk_global[i*Ncomp+iunk];
     }
     safe_free((void *) &unk_global);
     safe_free((void *) &index);
  }
  safe_free((void *) &unk_loc);

  return;
}
/*******************************************************************************
print_profile_box: This routine prints out the density profile. It gathers
the solution vector on Proc 0 and then prints it out using print_profile      */ 

void print_profile_box(double **x, char *outfile)
{

  if (Proc == 0){
    X_old = (double *) array_alloc (1, Nnodes*Nunk_per_node, sizeof(double));
    Vext_old = (double *) array_alloc (1, Nnodes*Ncomp, sizeof(double));
  }

  collect_x_old(x);
  collect_vext_old();

  if (Proc==0) {
     print_profile(outfile);
     safe_free((void *) &X_old);
     safe_free((void *) &Vext_old);
  }
}
/*******************************************************************************
print_profile: This routine prints out the density profile.        
this routine is only ever called by Proc 0                                    */ 

void print_profile(char *output_file4)
{
  int icomp,iunk,i,inode,ijk[3],idim,ipol,iseg,itype_mer,ibond,unk_GQ,unk_B;
  int unk_field, node_start;
  double kappa_sq,kappa,r,rsq,bondproduct,site_dens=0.,sumsegdens[NCOMP_MAX];
  char *unk_char;
  
  char gfile[20],gfile2[20];
  FILE *ifp=NULL,*fp6=NULL,*fp7=NULL;
  /* 
   *  print out the densities (and electrostatic potential)
   *  to the file dft_dens.dat or dft_dens.?.?   
   */

           /* open primary output file .... densities, electrostatic potential, and CMS fields */
     ifp = fopen(output_file4,"w");

           /* open file for CMS_G variables ... */
     if (Type_poly != -1){
       sprintf(gfile,"%sg",output_file4);
       fp6 = fopen(gfile,"w");
     } 

           /* open file for segment densities */
     if (Type_poly != NONE || Type_poly_TC){
       sprintf(gfile2,"%s_site",output_file4);
       fp7 = fopen(gfile2,"w");
     }

           /* print order of unknowns at the top of the file */
     for (i=0; i<NEQ_TYPE; i++){
         switch(i){
            case DENSITY: 
               unk_char = "DENSITY"; 
               if (Phys2Nunk[i] > 0){
                 fputs (unk_char,ifp); 
                 fprintf(ifp,"\n"); break;
               }
            case POISSON: 
               unk_char = "POISSON"; 
               if (Phys2Nunk[i] > 0){
                 fputs (unk_char,ifp); 
                 fprintf(ifp,"\n"); break;
               }
            case DIFFUSION: 
               unk_char = "CHEMPOT";
               if (Phys2Nunk[i] > 0){
                 fputs (unk_char,ifp); 
                 fprintf(ifp,"\n"); break;
               }
            case CMS_FIELD: 
               unk_char = "CMSFIELD";
               if (Phys2Nunk[i] > 0){
                 fputs (unk_char,ifp); 
                 fprintf(ifp,"\n"); break;
               }
            case RHOBAR_ROSEN: 
               unk_char="RHOBAR_ROSEN";
               if (Phys2Nunk[i] > 0){
                 fputs (unk_char,ifp); 
                 fprintf(ifp,"\n"); break;
               }
         }
     }

           /* compute DeBroglie wavelength for charged systems */
     if (Npoisson >0 ){
        kappa_sq = 0.0;
        if (Ipot_ff_c == COULOMB){
        for(icomp = 0; icomp<Ncomp; icomp++)
           kappa_sq += (4.0*PI/Temp_elec)*Rho_b[icomp]*
                        Charge_f[icomp]*Charge_f[icomp];
        }
        kappa = sqrt(kappa_sq);
     }

     for (inode=0; inode<Nnodes; inode++){
        node_to_ijk(inode,ijk);
        node_start = Nunk_per_node*inode;

                         /* print ijk coordinates of this node in the files */ 
        for (idim=0; idim<Ndim; idim++) {
                                    fprintf(ifp,"%9.6f\t ", ijk[idim]*Esize_x[idim]);
            if (Type_poly != NONE)  fprintf(fp6,"%9.6f\t ",ijk[idim]*Esize_x[idim]);
            if (Type_poly != NONE || Type_poly_TC) fprintf(fp7,"%9.6f\t ", ijk[idim]*Esize_x[idim]);
        }

        for (iunk=0; iunk<Nunk_per_node; iunk++){
            
            switch(Unk2Phys[iunk]){
               case DENSITY:   icomp = iunk-Phys2Unk_first[DENSITY]; break;
               case DIFFUSION: icomp = iunk-Phys2Unk_first[DIFFUSION]; break;
               case CMS_FIELD: icomp = iunk-Phys2Unk_first[CMS_FIELD]; break;
            }
            switch(Unk2Phys[iunk]){
                case DENSITY:
/*                fprintf(ifp,"%22.17f\t", X_old[iunk+node_start]/Rho_b[icomp]);*/
                case DIFFUSION:
/*                if (Ipot_ff_n != IDEAL_GAS)
                       fprintf(ifp,"%22.17f\t", X_old[iunk+node_start]
                            + 3.0*log(Sigma_ff[icomp][icomp]) + 1.5*log(Mass[icomp]*Temp)  );*/
                case POISSON:
                case RHOBAR_ROSEN:
                  fprintf(ifp,"%22.17f\t", X_old[iunk+node_start]);
                  break;

                case CMS_FIELD:
                   if (X_old[iunk+node_start] > 0.0 /*DENSITY_MIN*/ /*Rho_b[icomp]*exp(-VEXT_MAX)*/){
                      fprintf(ifp,"%22.17f\t", -log(X_old[iunk+node_start]));
                   }
                   else fprintf(ifp,"%22.17f\t", VEXT_MAX);
                   break;

                case CMS_G:
                   if (Type_poly == 2){
                      ipol=Unk_to_Poly[iunk-Phys2Unk_first[CMS_G]];
                      iseg=Unk_to_Seg[iunk-Phys2Unk_first[CMS_G]];
                      itype_mer=Type_mer[ipol][iseg];
                      unk_field = Phys2Unk_first[CMS_FIELD]+itype_mer;
                         fprintf(fp6,"%22.17f\t", X_old[iunk+node_start]*X_old[unk_field+node_start]);
                   }
                   else  fprintf(fp6,"%22.17f\t", X_old[iunk+node_start]);
                   break;

                case DENSITY_SEG:
                  fprintf(fp7,"%22.17f\t", X_old[iunk+node_start]);
                  break;
            }

        }    /* end loop over unknowns in the run */

                /* print the Poisson-Boltzmann solution based on the computed electrostatic field */
        if (Ipot_ff_c == 1 && !Sten_Type[POLYMER_CR]){
        for (icomp=0; icomp<Ncomp; icomp++)
          fprintf(ifp,"%20.15f\t",
                  Rho_b[icomp]*exp(-Charge_f[icomp]*X_old[Phys2Unk_first[POISSON]+node_start]
                                                              -Vext_old[inode*Ncomp+icomp]));
        }
 
               /* print segment densities for a CMS polymer run */
        if (Type_poly != NONE){
           for (icomp=0; icomp<Npol_comp; icomp++){
               sumsegdens[icomp]=0.0;
               for(iseg=0;iseg<Nmer[icomp];iseg++){
                    itype_mer=Type_mer[icomp][iseg];
                    bondproduct=1.0;
                    for(ibond=0;ibond<Nbond[icomp][iseg];ibond++){
                         unk_GQ  = Phys2Unk_first[CMS_G] + Poly_to_Unk[icomp][iseg][ibond];
                         bondproduct *= X_old[unk_GQ+node_start];
                    }  
                   unk_B=Phys2Unk_first[CMS_FIELD]+itype_mer;
                   if (Type_poly==2)
                      site_dens=bondproduct*X_old[unk_B+node_start]*Rho_b[itype_mer];
                   else
                      site_dens=bondproduct*POW_DOUBLE_INT(X_old[unk_B+node_start],-(Nbond[icomp][iseg]-1))
                                           *Rho_b[itype_mer]/Nmer_t[icomp][itype_mer];

                   sumsegdens[itype_mer]+=site_dens;
                   fprintf(fp7,"%22.17f\t", site_dens);
              }
           }
           for (itype_mer=0; itype_mer<Ntype_mer; itype_mer++) fprintf(fp7,"%22.17f\t", sumsegdens[itype_mer]);
        }
 
                /* add a carriage return to the file to start a new line */
        fprintf(ifp,"\n");
        if (Type_poly != NONE) fprintf(fp6,"\n");
        if (Type_poly != NONE || Type_poly_TC) fprintf(fp7,"\n");

                /* add some blank lines for improved graphics in 2D and 3D gnuplot */
        if (ijk[0] == Nodes_x[0]-1) fprintf(ifp,"\n");

     }    /* loop over all nodes  */

          /* close files */
     fclose(ifp);
     if (Type_poly != NONE) fclose(fp6);
     if (Type_poly != NONE || Type_poly_TC) fclose(fp7);

  return;
}
/*******************************************************************************
print_gofr: This routine prints out the density profile.        
this routine is only ever called by Proc 0                                    */ 

void print_gofr(char *output_file6)
{
  int icomp,i,inode,ijk[3],idim,nunk_print,npol=0,itype_mer,iwall,iunk;
  double kappa_sq,kappa,r,rsq;
  FILE *ifp=NULL;
  /* 
   *  print out the densities (and electrostatic potential)
   *  to the file dft_dens.dat or dft_dens.?.?   
   */

     ifp = fopen(output_file6,"w");

     if (Type_poly==NONE) nunk_print = Nunk_per_node;
     else nunk_print = 2*Ncomp;

     for (iwall=0; iwall<Nwall; iwall++){  /*compute g(r) for different atoms
                                             in one linked wall --- e.g. 
                                             could represent H2O as 3 atoms
                                             then need to compute gHH,gOO,gHO */
     for (inode=0; inode<Nnodes; inode++){
        node_to_ijk(inode,ijk);
 
        rsq=0.0; r=0.0;
        for (idim=0; idim<Ndim; idim++) {
            rsq = rsq+(ijk[idim]*Esize_x[idim]-(WallPos[idim][iwall]+0.5*Size_x[idim]))*
                      (ijk[idim]*Esize_x[idim]-(WallPos[idim][iwall]+0.5*Size_x[idim]));
        }
        if (rsq > 0.0) r=sqrt(rsq); 
        fprintf(ifp,"%9.6f\t ",r);

        for (iunk=0; iunk<Nunk_per_node; iunk++){
            if (Unk2Phys[iunk]==DENSITY){
                icomp = iunk-Phys2Unk_first[DENSITY];
                fprintf(ifp,"%22.17f\t", -log(X_old[iunk]/Rho_b[icomp]));
            }
        }

        fprintf(ifp,"\n");
        if (ijk[0] == Nodes_x[0]-1) fprintf(ifp,"\n");

     } /* loop over nodes */
     }    /* loop over all walls  */
     fclose(ifp);
  return;
}
/************************************************************************
print_zeroTF: This routine collects the zero_TF array and prints it out  */
void print_zeroTF(int **zero_TF, char *output_file)
{
  int icomp,loc_inode,inode,ijk[3],*index,idim,inode_box;
  int *unk_loc,*unk_global,**unk;
  FILE *ifp=NULL;

  if (Proc == 0) {
       ifp = fopen(output_file,"w");
       unk = (int **) array_alloc (2, Nnodes, Ncomp+1, sizeof(int));
  }

  unk_loc = (int *) array_alloc (1, Nnodes_per_proc, sizeof(int));

  for (icomp=0; icomp<Ncomp+1; icomp++){

     /*  define the local array for each component separately */

     for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ ){
         inode_box = L2B_node[loc_inode];
         unk_loc[loc_inode] = zero_TF[inode_box][icomp];
     }

     if (Proc ==0){
       index = (int *) array_alloc (1, Nnodes, sizeof(int));
       unk_global = (int *) array_alloc (1, Nnodes, sizeof(int));
     }

     /* collect the global indices from all processors */
     MPI_Gatherv(L2G_node,Nnodes_per_proc,MPI_INT,
              index,Comm_node_proc,Comm_offset_node,
              MPI_INT,0,MPI_COMM_WORLD);

     /* collect the unknowns from all the processors */

     MPI_Gatherv(unk_loc,Nnodes_per_proc,MPI_INT,
              unk_global,Comm_node_proc,Comm_offset_node,
              MPI_INT,0,MPI_COMM_WORLD);

     if (Proc == 0){
        for (inode=0; inode<Nnodes; inode++){
            unk[index[inode]][icomp] = unk_global[inode];
        }
        safe_free((void *) &unk_global);
        safe_free((void *) &index);
     }
  }
  safe_free((void *) &unk_loc);
     
  /* 
   *  now print out the Zero_density_TF array to the file dft_zeroTF.dat
   */
  if (Proc ==0){

     for (inode=0; inode<Nnodes; inode++){
        node_to_ijk(inode,ijk);
        for (idim=0; idim<Ndim; idim++)
            fprintf(ifp,"%9.6f\t ",
            ijk[idim]*Esize_x[idim]);

        for (icomp=0; icomp<Ncomp+1; icomp++){
            fprintf(ifp,"%d\t", unk[inode][icomp]);
        }

        fprintf(ifp,"\n");
        if (ijk[0] == Nodes_x[0]-1) fprintf(ifp,"\n");
     }    /* loop over all nodes  */
     fclose(ifp);
     safe_free((void *) &unk);
  }       /* end of Proc ==0 test */


  return;
}
/************************************************************************
print_Nodes_to_zone: This routine collects and prints nodes_to_zone  */
void print_Nodes_to_zone(int *node_to_zone, char *output_file)
{
  int loc_inode,inode,ijk[3],*index,idim,inode_box;
  int *unk,*unk_loc, *unk_global;
  FILE *ifp=NULL;

  if (Proc == 0){
     ifp = fopen(output_file,"w");
     unk = (int *) array_alloc (1, Nnodes, sizeof(int));
  }

  unk_loc = (int *) array_alloc (1, Nnodes_per_proc, sizeof(int));

  for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ ){
      inode_box = L2B_node[loc_inode];
      unk_loc[loc_inode] = node_to_zone[inode_box];
  }

  if (Proc ==0){
    index = (int *) array_alloc (1, Nnodes, sizeof(int));
    unk_global = (int *) array_alloc (1, Nnodes, sizeof(int));
  }

  /* collect the global indices from all processors */
  MPI_Gatherv(L2G_node,Nnodes_per_proc,MPI_INT,
           index,Comm_node_proc,Comm_offset_node,
           MPI_INT,0,MPI_COMM_WORLD);

  /* collect the unknowns from all the processors */

  MPI_Gatherv(unk_loc,Nnodes_per_proc,MPI_INT,
           unk_global,Comm_node_proc,Comm_offset_node,
           MPI_INT,0,MPI_COMM_WORLD);
  safe_free((void *) &unk_loc);

  if (Proc == 0){
     for (inode=0; inode < Nnodes; inode++ ){
        unk[index[inode]] = unk_global[inode];
     }

     if (Proc==0) safe_free((void *) &index);
     if (Proc==0) safe_free((void *) &unk_global);
  }

  /* 
   * now print the array.
   */
  if (Proc ==0){

     for (inode=0; inode<Nnodes; inode++){
        node_to_ijk(inode,ijk);
        for (idim=0; idim<Ndim; idim++)
            fprintf(ifp,"%9.6f\t ",
            (double)ijk[idim]*Esize_x[idim]);

            fprintf(ifp,"%d \t", unk[inode]);

        fprintf(ifp,"\n");
        if (ijk[0] == Nodes_x[0]-1) fprintf(ifp,"\n");
     }    /* loop over all nodes  */
     fclose(ifp);
     safe_free((void *) &unk);
  }       /* end of Proc ==0 test */

  return;
}
/************************************************************************
print_charge_surf: This routine collects and prints Charge_w_sum_els  */
void print_charge_surf(double **charge_w_sum, char *output_file)
{
  int i,icount;
  int loc_inode,inode,ijk[3],*index_loc,*index,idim;
  int reflect_flag[3];
  int *comm_icount_proc, *comm_offset_icount;
  double **unk,*unk_loc, *unk_global;
  FILE *ifp=NULL;


  if (Proc == 0){
     unk = (double **) array_alloc (2, Nnodes, Ndim, sizeof(double));
     ifp = fopen(output_file,"w");
  }
  reflect_flag[0] = reflect_flag[1] = reflect_flag[2] = FALSE;

  for (idim=0; idim<Ndim; idim++) {

  index_loc = (int *) array_alloc (1, Nnodes_per_proc, sizeof(int));
  unk_loc = (double *) array_alloc (1, Nnodes_per_proc, sizeof(double));

  icount=0;
  for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ ){
      index_loc[icount] = loc_inode;
      unk_loc[icount++] = charge_w_sum[loc_inode][idim];
  }

  if (Proc ==0){
    unk_global = (double *) array_alloc (1, Nnodes, sizeof(double));
    index = (int *) array_alloc (1, Nnodes, sizeof(int));
  }

  comm_icount_proc = (int *) array_alloc (1, Num_Proc, sizeof(int));
  comm_offset_icount = (int *) array_alloc (1, Num_Proc, sizeof(int));

  MPI_Gather(&icount,1,MPI_INT,
             comm_icount_proc,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc == 0){
     comm_offset_icount[0] = 0; 
     for (i=1; i<Num_Proc; i++){
        comm_offset_icount[i] = comm_offset_icount[i-1] + comm_icount_proc[i-1];
     }
  }

  /* collect the global indices from all processors */
  MPI_Gatherv(index_loc,Nnodes_per_proc,MPI_INT,
           index,comm_icount_proc,comm_offset_icount,
           MPI_INT,0,MPI_COMM_WORLD);
  safe_free((void *) &index_loc);

  /* collect the unknowns from all the processors */

  MPI_Gatherv(unk_loc,Nnodes_per_proc,MPI_DOUBLE,
           unk_global,comm_icount_proc,comm_offset_icount,
           MPI_DOUBLE,0,MPI_COMM_WORLD);
  safe_free((void *) &unk_loc);

  
  safe_free((void *) &comm_icount_proc);
  safe_free((void *) &comm_offset_icount);

  if (Proc == 0) {
     for (inode=0; inode<Nnodes; inode++){
         unk[index[inode]][idim] = unk_global[inode];
     }
     safe_free((void *) &unk_global);
     safe_free((void *) &index);
  }

  }


  /* 
   *  now print out the volumetric charges
   */
  if (Proc ==0){

     for (inode=0; inode<Nnodes; inode++){
        node_to_ijk(inode,ijk);
        for (idim=0; idim<Ndim; idim++)
            fprintf(ifp,"%9.6f\t ", ((double)ijk[idim])*Esize_x[idim]);

        for (idim=0; idim<Ndim; idim++)
            fprintf(ifp,"%9.6f\t", unk[inode][idim]);

        fprintf(ifp,"\n");
        if (ijk[0] == Nodes_x[0]) fprintf(ifp,"\n");
     }    /* loop over all nodes  */

     safe_free((void *) &unk);
     fclose(ifp);
  }       /* end of Proc ==0 test */
  return;
}
/************************************************************************
print_free_energy_profile: This routine collects and prints freen_profile_1D  */
void print_freen_profile_1D(double *freen_profile_1D, char *output_file)
{
  int i,iel,iel_box,icount,logical;
  int loc_inode,inode,ijk[3],*index_loc,*index,idim,inode_box;
  int reflect_flag[3];
  int *comm_icount_proc, *comm_offset_icount;
  double *unk,*unk_loc, *unk_global,charge_total;
  FILE *ifp=NULL;

  reflect_flag[0] = reflect_flag[1] = reflect_flag[2] = FALSE;

  if (Proc == 0){
     unk = (double *) array_alloc (1, Nnodes, sizeof(double));
     ifp = fopen(output_file,"w");
  }

  index_loc = (int *) array_alloc (1, Nnodes_per_proc, sizeof(int));
  unk_loc = (double *) array_alloc (1, Nnodes_per_proc, sizeof(double));

  icount=0;
  for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ ){
      inode_box = L2B_node[loc_inode];
      inode = L2G_node[loc_inode];
      index_loc[icount] = inode;
      unk_loc[icount++] = freen_profile_1D[loc_inode];
  }

  if (Proc ==0){
    unk_global = (double *) array_alloc (1, Nnodes, sizeof(double));
    index = (int *) array_alloc (1, Nnodes, sizeof(int));
  }

  comm_icount_proc = (int *) array_alloc (1, Num_Proc, sizeof(int));
  comm_offset_icount = (int *) array_alloc (1, Num_Proc, sizeof(int));

  MPI_Gather(&icount,1,MPI_INT, comm_icount_proc,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc == 0){
     comm_offset_icount[0] = 0; 
     for (i=1; i<Num_Proc; i++){
        comm_offset_icount[i] = comm_offset_icount[i-1] + comm_icount_proc[i-1];
     }
  }

  /* collect the global indices from all processors */
  MPI_Gatherv(index_loc,icount,MPI_INT,
           index,comm_icount_proc,comm_offset_icount,
           MPI_INT,0,MPI_COMM_WORLD);
  safe_free((void *) &index_loc);

  /* collect the unknowns from all the processors */

  MPI_Gatherv(unk_loc,icount,MPI_DOUBLE,
           unk_global,comm_icount_proc,comm_offset_icount,
           MPI_DOUBLE,0,MPI_COMM_WORLD);
  safe_free((void *) &unk_loc);

  
  safe_free((void *) &comm_icount_proc);
  safe_free((void *) &comm_offset_icount);

  if (Proc == 0) {
     for (inode=0; inode<Nnodes; inode++){
         unk[index[inode]] = unk_global[inode];
     }
     safe_free((void *) &unk_global);
     safe_free((void *) &index);
  }

  /* 
   *  now print out the free energy profile
   */
  if (Proc ==0){
     for (inode=0; inode<Nnodes; inode++){
        node_to_ijk(inode,ijk);
/*        for (idim=0; idim<Ndim; idim++)*/  /* only do free energy profile in 1D so far */
            idim=0;
            fprintf(ifp,"%9.6f\t ",
            ((double)ijk[idim])*Esize_x[idim]);

            fprintf(ifp,"%9.6f\t", unk[inode]);

        fprintf(ifp,"\n");
        if (ijk[0] == Nodes_x[0]-1) fprintf(ifp,"\n");
     }    /* loop over all nodes  */

     safe_free((void *) &unk);
     fclose(ifp);
  }       /* end of Proc ==0 test */
  return;
}
/************************************************************************
print_charge_vol: This routine collects and prints Charge_vol_els  */
void print_charge_vol(double *charge_els, char *output_file)
{
  int i,iel,iel_box,icount,logical;
  int loc_inode,inode,ijk[3],*index_loc,*index,idim,inode_box;
  int reflect_flag[3];
  int *comm_icount_proc, *comm_offset_icount;
  double *unk,*unk_loc, *unk_global,charge_total;
  FILE *ifp=NULL;

  reflect_flag[0] = reflect_flag[1] = reflect_flag[2] = FALSE;

  if (Proc == 0){
     unk = (double *) array_alloc (1, Nelements, sizeof(double));
     ifp = fopen(output_file,"w");
  }

  index_loc = (int *) array_alloc (1, Nnodes_per_proc, sizeof(int));
  unk_loc = (double *) array_alloc (1, Nnodes_per_proc, sizeof(double));

  icount=0;
  for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ ){
      inode_box = L2B_node[loc_inode];
      inode = L2G_node[loc_inode];
      node_to_ijk(inode,ijk);

      logical = FALSE;
      for (idim=0; idim<Ndim; idim++) 
          if (ijk[idim] == Nodes_x[idim] -1 
              && Type_bc[idim][1] != PERIODIC) logical = TRUE;

      if (!logical) {
         iel   = node_to_elem(inode,0,reflect_flag);
         iel_box = node_box_to_elem_box_reflect(inode_box,0,reflect_flag);
         index_loc[icount] = iel;
         unk_loc[icount++] = charge_els[iel_box];
      }
  }

  if (Proc ==0){
    unk_global = (double *) array_alloc (1, Nnodes, sizeof(double));
    index = (int *) array_alloc (1, Nnodes, sizeof(int));
  }

  comm_icount_proc = (int *) array_alloc (1, Num_Proc, sizeof(int));
  comm_offset_icount = (int *) array_alloc (1, Num_Proc, sizeof(int));

  MPI_Gather(&icount,1,MPI_INT, comm_icount_proc,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc == 0){
     comm_offset_icount[0] = 0; 
     for (i=1; i<Num_Proc; i++){
        comm_offset_icount[i] = comm_offset_icount[i-1] + comm_icount_proc[i-1];
     }
  }

  /* collect the global indices from all processors */
  MPI_Gatherv(index_loc,icount,MPI_INT,
           index,comm_icount_proc,comm_offset_icount,
           MPI_INT,0,MPI_COMM_WORLD);
  safe_free((void *) &index_loc);

  /* collect the unknowns from all the processors */

  MPI_Gatherv(unk_loc,icount,MPI_DOUBLE,
           unk_global,comm_icount_proc,comm_offset_icount,
           MPI_DOUBLE,0,MPI_COMM_WORLD);
  safe_free((void *) &unk_loc);

  
  safe_free((void *) &comm_icount_proc);
  safe_free((void *) &comm_offset_icount);

  if (Proc == 0) {
     charge_total=0.0;
     for (iel=0; iel<Nelements; iel++){
         unk[index[iel]] = unk_global[iel];
         charge_total+=unk_global[iel];
     }
     safe_free((void *) &unk_global);
     safe_free((void *) &index);
     printf("THE NET VOLUME CHARGE IN THE DOMAIN IS : %9.6f\n",charge_total);
  }

  /* 
   *  now print out the volumetric charges
   */
  if (Proc ==0){

     for (iel=0; iel<Nelements; iel++){
        inode = element_to_node(iel);
        node_to_ijk(inode,ijk);
        for (idim=0; idim<Ndim; idim++)
            fprintf(ifp,"%9.6f\t ",
            ((double)ijk[idim]+0.5)*Esize_x[idim]);

            fprintf(ifp,"%9.6f\t", unk[iel]);

        fprintf(ifp,"\n");
        if (ijk[0] == Nodes_x[0]-1) fprintf(ifp,"\n");
     }    /* loop over all nodes  */

     safe_free((void *) &unk);
     fclose(ifp);
  }       /* end of Proc ==0 test */
  return;
}
/************************************************************************
print_vext: This routine collects the vext array and prints it out  */
void print_vext(double **vext, char *output_file)
{
  int icomp,loc_inode,inode,ijk[3],*index,idim;
  double *unk_loc,*unk_global,**unk,rsq,r;
  FILE *ifp=NULL;

  if (Proc == 0) {
       ifp = fopen(output_file,"w");
       unk = (double **) array_alloc (2, Nnodes, Ncomp, sizeof(double));
  }

  unk_loc = (double *) array_alloc (1, Nnodes_per_proc, sizeof(double));

  for (icomp=0; icomp<Ncomp; icomp++){

     /* 
      *  first collect all the icomp unknowns on processor 0 
      */

     for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ ){
         unk_loc[loc_inode] = vext[loc_inode][icomp];
     }

     if (Proc ==0){
       index = (int *) array_alloc (1, Nnodes, sizeof(int));
       unk_global = (double *) array_alloc (1, Nnodes, sizeof(double));
     }

     /* collect the global indices from all processors */
     MPI_Gatherv(L2G_node,Nnodes_per_proc,MPI_INT,
              index,Comm_node_proc,Comm_offset_node,
              MPI_INT,0,MPI_COMM_WORLD);

     /* collect the unknowns from all the processors */

     MPI_Gatherv(unk_loc,Nnodes_per_proc,MPI_DOUBLE,
              unk_global,Comm_node_proc,Comm_offset_node,
              MPI_DOUBLE,0,MPI_COMM_WORLD);

     if (Proc == 0){
        for (inode=0; inode<Nnodes; inode++){
            unk[index[inode]][icomp] = unk_global[inode];
        }
        safe_free((void *) &unk_global);
        safe_free((void *) &index);
     }
  }
  safe_free((void *) &unk_loc);

     
  /* 
   *  now print out the densities (and electrostatic potential)
   *  to the file dft_dens.dat. 
   */
  if (Proc ==0){

     for (inode=0; inode<Nnodes; inode++){
        node_to_ijk(inode,ijk);
        r=0.0; rsq=0.0;
        for (idim=0; idim<Ndim; idim++){
            fprintf(ifp,"%9.6f\t ", ijk[idim]*Esize_x[idim]);
            rsq +=  ijk[idim]*Esize_x[idim]* ijk[idim]*Esize_x[idim];
        }
        r=sqrt(rsq);
        if (Lprint_gofr) fprintf(ifp,"%9.6f\t ",r);

        for (icomp=0; icomp<Ncomp; icomp++)
            fprintf(ifp,"%22.17f\t", unk[inode][icomp]);

        fprintf(ifp,"\n");
        if (ijk[0] == Nodes_x[0]-1) fprintf(ifp,"\n");
     }    /* loop over all nodes  */
     fclose(ifp);
     safe_free((void *) &unk);
  }       /* end of Proc ==0 test */
  return;
}
/************************************************************************
print_rho_bar: This routine collects arrays of the rho_bar type
               and prints them out  */
void print_rho_bar(struct RB_Struct *rho_bar, char *output_file)
{
  int loc_inode,inode,ijk[3],idim,inode_box,isten,imax,*index;
  double *unk_loc,*unk_global;
  struct  RB_Struct *unk = NULL;
  struct  RB_Struct *point;
  FILE *ifp=NULL;

  if (Proc == 0) {
       ifp = fopen(output_file,"w");
       unk = (struct RB_Struct *) array_alloc
             (1, Nnodes, sizeof(struct RB_Struct));
  }
  
  if (Ndim == 1)      imax = 6;
  else if (Ndim == 2) imax = 8;
  else                imax = 10;

  unk_loc = (double *) array_alloc (1, Nnodes_per_proc, sizeof(double));

  for (isten=0; isten<imax; isten++){

     /* 
      *  first collect all the icomp unknowns on processor 0 
      */

     for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ ){
         inode_box = L2B_node[loc_inode];

         if (B2L_1stencil[inode_box] >=0) {
            point = &(rho_bar[B2L_1stencil[inode_box]]);

            if      (isten == 0) unk_loc[loc_inode] = point->S0;
            else if (isten == 1) unk_loc[loc_inode] = point->S1;
            else if (isten == 2) unk_loc[loc_inode] = point->S2;
            else if (isten == 3) unk_loc[loc_inode] = point->S3;
            else if (isten == 4)  unk_loc[loc_inode] = point->V1[0];
            else if (isten == 5)  unk_loc[loc_inode] = point->V2[0];
            else if (isten == 6) unk_loc[loc_inode] = point->V1[1];
            else if (isten == 7) unk_loc[loc_inode] = point->V2[1];
            else if (isten == 8) unk_loc[loc_inode] = point->V1[2];
            else if (isten == 9) unk_loc[loc_inode] = point->V2[2];
         }
         else unk_loc[loc_inode] = 0.0; 
     }

     if (Proc ==0){
       index = (int *) array_alloc (1, Nnodes, sizeof(int));
       unk_global = (double *) array_alloc (1, Nnodes, sizeof(double));
     }

     /* collect the global indices from all processors */

     MPI_Gatherv(L2G_node,Nnodes_per_proc,MPI_INT,
              index,Comm_node_proc,Comm_offset_node,
              MPI_INT,0,MPI_COMM_WORLD);


     /* collect the unknowns from all the processors */

     MPI_Gatherv(unk_loc,Nnodes_per_proc,MPI_DOUBLE,
              unk_global,Comm_node_proc,Comm_offset_node,
              MPI_DOUBLE,0,MPI_COMM_WORLD);

     if (Proc == 0){
        for (inode=0; inode<Nnodes; inode++){
            if     (isten == 0) unk[index[inode]].S0 = unk_global[inode]; 
            else if(isten == 1) unk[index[inode]].S1 = unk_global[inode]; 
            else if(isten == 2) unk[index[inode]].S2 = unk_global[inode]; 
            else if(isten == 3) unk[index[inode]].S3 = unk_global[inode]; 
            else if(isten == 4) unk[index[inode]].V1[0] = unk_global[inode]; 
            else if(isten == 5) unk[index[inode]].V2[0] = unk_global[inode]; 
            else if(isten == 6) unk[index[inode]].V1[1] = unk_global[inode]; 
            else if(isten == 7) unk[index[inode]].V2[1] = unk_global[inode]; 
            else if(isten == 8) unk[index[inode]].V1[2] = unk_global[inode]; 
            else if(isten == 9) unk[index[inode]].V2[2] = unk_global[inode]; 
        }
        safe_free((void *) &unk_global);
        safe_free((void *) &index);
     }
  }
  safe_free((void *) &unk_loc);
     
  /* 
   *  now print out rho_bars or dphi_drho_bars)
   *  to the input file. 
   */
  if (Proc ==0){

     for (inode=0; inode<Nnodes; inode++){
        fprintf(ifp,"%d\t ",inode);

        node_to_ijk(inode,ijk);

        for (idim=0; idim<Ndim; idim++)
            fprintf(ifp,"%9.6f\t ",ijk[idim]*Esize_x[idim]);

        fprintf(ifp,"%10.6f\t %10.6f\t %10.6f\t %10.6f\t", 
                              unk[inode].S0,unk[inode].S1,
                              unk[inode].S2,unk[inode].S3);

        for (idim=0; idim<Ndim; idim++)
            fprintf(ifp,"%10.6f\t %10.6f\t", 
                    unk[inode].V1[idim],unk[inode].V2[idim]);

        fprintf(ifp,"\n");
        if (ijk[0] == Nodes_x[0]-1) fprintf(ifp,"\n");

     }    /* loop over all nodes  */
     fclose(ifp);
     safe_free((void *) &unk);
  }       /* end of Proc ==0 test */
  return;
}
/****************************************************************************/
