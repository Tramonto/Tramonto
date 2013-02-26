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
 *  FILE: dft_chain_architecture.c
 *
 *  This file contains routines to set up fluids composed of chains
 *  and initializes several arrays that are used for indexing in bonded systems.
 */

#include "dft_chain_architecture.h"

/*************************************************************************************/
void setup_chain_architecture(char *poly_file,FILE *fpecho)

{
   /* Local variable declarations */
   
   int pol_number, nseg,nmer_max;
   
   if (Type_poly_arch != SET_IN_GUI){    
       nseg=nmer_max=0;
       for (pol_number=0; pol_number<Npol_comp; pol_number++){
            nseg += Nmer[pol_number];
            if (Nmer[pol_number] > nmer_max) nmer_max=Nmer[pol_number];
       }   
       Nbond = (int **) array_alloc (2, Npol_comp,nmer_max,sizeof(int));
       Bonds = (int ***) array_alloc (3, Npol_comp,nmer_max,NBOND_MAX,sizeof(int));
       pol_sym_tmp = (int ***) array_alloc (3, Npol_comp,nmer_max,NBOND_MAX,sizeof(int));

      if (Type_poly_arch==POLY_ARCH_FILE) setup_chain_from_file(fpecho,poly_file);
      else if (Type_poly_arch==LIN_POLY) setup_chain_linear(fpecho);
      else if (Type_poly_arch==LIN_POLY_SYM) setup_chain_linear_symmetric(fpecho);
   }

   return;
}
/*************************************************************************************/
/* read in a file that contains chain architecture */
void setup_chain_from_file(FILE *fpecho, char *poly_file)
{
   int pol_number,iseg,seg_id,ibond,graft_point_indicated;
   FILE *fppoly;

   if( (fppoly  = fopen(poly_file,"r")) == NULL) {
           if (Iwrite_screen != SCREEN_NONE) printf("Can't open poly file %s and no other file was specified in GUI.\n", poly_file);
           exit(-1);
   }

   Nbond_max=0; 
   for (pol_number=0; pol_number<Npol_comp; ++pol_number){
      graft_point_indicated=FALSE;
      for (iseg=0; iseg<Nmer[pol_number]; iseg++){
           fscanf(fppoly,"%d %d",&seg_id, &Nbond[pol_number][iseg]);
           if (Nbond[pol_number][iseg]>Nbond_max) Nbond_max=Nbond[pol_number][iseg];
           if (seg_id != iseg){
               if (Iwrite_screen != SCREEN_NONE) printf("problem with seg_ids in file - there should be no skips, and every chain starts at 0\n");
               exit(-1);
           }

	for (ibond=0; ibond<Nbond[pol_number][iseg]; ibond++){
	    fscanf(fppoly,"%d  %d", &Bonds[pol_number][iseg][ibond],&pol_sym_tmp[pol_number][iseg][ibond]);
	    if (Iwrite_files == FILES_DEBUG) fprintf(fpecho,"%d  ",Bonds[pol_number][iseg][ibond]);
            if(Type_poly==WJDC3 && Grafted[pol_number] && iseg==Grafted_SegID[pol_number]){
               if (Bonds[pol_number][iseg][ibond]==-1) graft_point_indicated=TRUE;
            }
	} /* end of loop over ibond */
        if (Type_poly==WJDC3 && Grafted[pol_number] && iseg==Grafted_SegID[pol_number] && graft_point_indicated==FALSE){
             Bonds[pol_number][iseg][Nbond[pol_number][iseg]]=-1;
             pol_sym_tmp[pol_number][iseg][Nbond[pol_number][iseg]]=-1;
             Nbond[pol_number][iseg]+=1;
        }
      } /* end of loop over iseg */
   }
   if (Proc==0 && Iwrite_files==FILES_DEBUG) fclose (fppoly);
   return;
}
/*************************************************************************************/
/* automatically set up a linear chain */
void setup_chain_linear(FILE *fpecho){

   int pol_number,iseg,ibond;
   Nbond_max=2;

   for (pol_number=0; pol_number<Npol_comp; ++pol_number){
      for (iseg=0; iseg<Nmer[pol_number]; iseg++){
	Nbond[pol_number][iseg]=2;

	for (ibond=0; ibond<Nbond[pol_number][iseg]; ibond++){
          if ((iseg==0 && ibond==0) || (iseg==Nmer[pol_number]-1 && ibond==Nbond[pol_number][iseg]-1)){
               Bonds[pol_number][iseg][ibond]=-1;
               pol_sym_tmp[pol_number][iseg][ibond]=-1;
	       if (Proc==0 && Iwrite_files==FILES_DEBUG) fprintf(fpecho,"%d  ",Bonds[pol_number][iseg][ibond]);
          }
          else if (ibond==0){
               Bonds[pol_number][iseg][ibond]=iseg-1;
               pol_sym_tmp[pol_number][iseg][ibond]=-1;
	       if (Proc==0 && Iwrite_files==FILES_DEBUG) fprintf(fpecho,"%d  ",Bonds[pol_number][iseg][ibond]);
          }
          else if (ibond==1){
               Bonds[pol_number][iseg][ibond]=iseg+1;
               pol_sym_tmp[pol_number][iseg][ibond]=-1;
	       if (Proc==0 && Iwrite_files==FILES_DEBUG) fprintf(fpecho,"%d  ",Bonds[pol_number][iseg][ibond]);
          }
          else{
             printf("problem in linear chain code - can only have ibond=0 or ibond=1...ibond=%d\n",ibond);
             exit(-1);
          }

	} /* end of loop over ibond */
      } /* end of loop over iseg */
   }
   return;
}
/*************************************************************************************/
/* automatically set up a linear chain with proper symmetries in place */
void setup_chain_linear_symmetric(FILE *fpecho){

   int pol_number,iseg,iseg_sym,ibond;
   Nbond_max=2;

   for (pol_number=0; pol_number<Npol_comp; ++pol_number){
      for (iseg=0; iseg<Nmer[pol_number]; iseg++){

        Nbond[pol_number][iseg]=2;

        for (ibond=0; ibond<Nbond[pol_number][iseg]; ibond++){
          if ((iseg==0 && ibond==0) || (iseg==Nmer[pol_number]-1 && ibond==Nbond[pol_number][iseg]-1)){ 
               Bonds[pol_number][iseg][ibond]=-1;
               if (iseg==0)  pol_sym_tmp[pol_number][iseg][ibond]=-1;
               else          pol_sym_tmp[pol_number][iseg][ibond]= 0;
	       if (Proc==0 && Iwrite_files==FILES_DEBUG) fprintf(fpecho,"%d  ",Bonds[pol_number][iseg][ibond]);
          }
          else if (ibond==0){ 
               Bonds[pol_number][iseg][ibond]=iseg-1;
               if (iseg<Nmer[pol_number]/2) pol_sym_tmp[pol_number][iseg][ibond]=-1;
               else {
                  iseg_sym=(Nmer[pol_number]-1)-iseg;
                  pol_sym_tmp[pol_number][iseg][ibond]=2*iseg_sym+1;
               }
	       if (Proc==0 && Iwrite_files==FILES_DEBUG) fprintf(fpecho,"%d  ",Bonds[pol_number][iseg][ibond]);
          }
          else if (ibond==1){ 
               Bonds[pol_number][iseg][ibond]=iseg+1;
               if (iseg<Nmer[pol_number]/2) pol_sym_tmp[pol_number][iseg][ibond]=-1;
               else {
                   iseg_sym=(Nmer[pol_number]-1)-iseg;
                   pol_sym_tmp[pol_number][iseg][ibond]=2*iseg_sym;
               }
	       if (Proc==0 && Iwrite_files==FILES_DEBUG) fprintf(fpecho,"%d  ",Bonds[pol_number][iseg][ibond]);
          }
          else{
             printf("problem in linear chain code - can only have ibond=0 or ibond=1...ibond=%d\n",ibond);
             exit(-1);
          }
        } /* end of loop over ibond */
      } /* end of loop over iseg */
   }
   return;
}
/*************************************************************************************/
void setup_chain_indexing_arrays(int nseg, int nmer_max,FILE *fpecho)
{
    int *nbond_tot;
    int nbond_all,iseg,seg_tot,end_count_all,icomp,pol_number,nunk,end_count,ibond,pol_num2;
 
    /* first define some arrays we need for our calculations */

    Unk_to_Poly = (int *) array_alloc (1, NBOND_MAX*nseg, sizeof(int));
    Unk_to_Seg = (int *) array_alloc (1, NBOND_MAX*nseg, sizeof(int));
    Unk_to_Bond = (int *) array_alloc (1, NBOND_MAX*nseg, sizeof(int));
    Nbonds_SegAll = (int *) array_alloc (1, NMER_MAX,sizeof(int));
    Bonds_SegAll = (int **) array_alloc (2, NMER_MAX,NBOND_MAX,sizeof(int));
    Poly_to_Unk = (int ***) array_alloc (3, Npol_comp,nmer_max,NBOND_MAX,sizeof(int));
    Poly_to_Unk_SegAll = (int **) array_alloc (2, NMER_MAX,NBOND_MAX,sizeof(int));
    Pol_Sym = (int *) array_alloc (1, nseg*NBOND_MAX,sizeof(int));
    Pol_Sym_Seg = (int *) array_alloc (1, nseg,sizeof(int));
    BondAll_to_isegAll = (int *) array_alloc (1, nseg*NBOND_MAX,sizeof(int));
    BondAll_to_ibond = (int *) array_alloc (1, nseg*NBOND_MAX,sizeof(int));
    nbond_tot = (int *) array_alloc (1, Npol_comp, sizeof(int));
    Nseg_type_pol = (int **) array_alloc (2, Npol_comp,Ncomp,sizeof(int));

    /* now zero some counters and arrays */
    nbond_all = 0; 
    Nbonds=0;
    Nseg_tot=0;
    seg_tot=0;
    end_count_all=0;
    for (icomp=0;icomp<Ncomp;icomp++){
         Nmer_comp[icomp]=0;
       for (pol_number=0; pol_number<Npol_comp; ++pol_number){
             Nseg_type_pol[pol_number][icomp]=0;
       }
    } 

    for (pol_number=0; pol_number<Npol_comp; ++pol_number){
      Nseg_tot += Nmer[pol_number];
      if(Nseg_tot > NMER_MAX) {
      	printf("Error: too many polymer segments, must increase NMER_MAX\n");
      	exit(-1);
      }
      nbond_tot[pol_number]=0;
      nunk=0; 
      for (iseg=0; iseg<Nmer[pol_number]; iseg++){
        Pol_Sym_Seg[seg_tot]=-1;
        end_count=0;
        Nbonds_SegAll[seg_tot]=0;
        SegAll_to_Poly[seg_tot]=pol_number;

	for (ibond=0; ibond<Nbond[pol_number][iseg]; ibond++){
          if (Type_poly!=WTC || (Type_poly==WTC && Bonds[pol_number][iseg][ibond] != -1)){
                                  /* note we don't want to include ends for the WTC polymers!*/
	    Unk_to_Poly[nbond_all] = pol_number;
  	    Unk_to_Seg[nbond_all]  = iseg;
	    Unk_to_Bond[nbond_all] = ibond;
	    Poly_to_Unk[pol_number][iseg][ibond] = nunk;
	    if(Bonds[pol_number][iseg][ibond] != -1)
	      Bonds_SegAll[seg_tot][Nbonds_SegAll[seg_tot]]=Bonds[pol_number][iseg][ibond]+SegChain2SegAll[pol_number][0];
	    else
	      Bonds_SegAll[seg_tot][Nbonds_SegAll[seg_tot]]=Bonds[pol_number][iseg][ibond];
	    Poly_to_Unk_SegAll[seg_tot][Nbonds_SegAll[seg_tot]] = nbond_all;
	    if (pol_sym_tmp[pol_number][iseg][ibond] != -1) Pol_Sym[nbond_all]=pol_sym_tmp[pol_number][iseg][ibond]-end_count_all;
            else                                            Pol_Sym[nbond_all]=-1;
            BondAll_to_isegAll[nbond_all]=seg_tot;
	    BondAll_to_ibond[nbond_all]=Nbonds_SegAll[seg_tot];

            if (Pol_Sym[nbond_all]!= -1 && (Type_poly==WTC || Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3)){
                if(Pol_Sym_Seg[seg_tot]==-1){ 
                     Pol_Sym_Seg[seg_tot] = BondAll_to_isegAll[Pol_Sym[nbond_all]];  /* need to debug for chain mixtures */
                }
            }
	    nbond_all++;
	    nunk++;
            Nbonds++; 
            Nbonds_SegAll[seg_tot]++;
          }
          else if (Type_poly==WTC && Bonds[pol_number][iseg][ibond] == -1){
              end_count++;
              end_count_all++;
          }
	} /* end of loop over ibond */
	nbond_tot[pol_number] += (Nbond[pol_number][iseg]-end_count);
        Unk2Comp[seg_tot]=Type_mer[pol_number][iseg];
        Nmer_comp[Unk2Comp[seg_tot]]++;
        seg_tot++;
        Nseg_type_pol[pol_number][Type_mer[pol_number][iseg]]++;
        Icomp_to_polID[Type_mer[pol_number][iseg]]=pol_number;
      } /* end of loop over iseg */
    }



    for (icomp=0;icomp<Ncomp;icomp++) Nseg_type[icomp]=0;
    for (iseg=0;iseg<Nseg_tot;iseg++) Nseg_type[Unk2Comp[iseg]]++;
    if (Proc==0 && Iwrite_files==FILES_DEBUG){
       fprintf(fpecho,"\n********************\n BOND DETAILS \n **********************\n");
       fprintf(fpecho,"\t total number of bonds is %d\n",Nbonds);
       for (ibond=0;ibond<Nbonds; ibond++){
           fprintf(fpecho,"Unk_to_Poly[ibond=%d]=%d Unk_to_Seg[]=%d Unk_to_Bond[]=%d\n",
            ibond,Unk_to_Poly[ibond],Unk_to_Seg[ibond],Unk_to_Bond[ibond]);
       }
       for (pol_number=0; pol_number<Npol_comp; ++pol_number){
          for (iseg=0;iseg<Nmer[pol_number];iseg++){
              for (ibond=0;ibond<Nbond[pol_number][iseg];ibond++){
	        if(Bonds[pol_number][iseg][ibond] != -1)
                  fprintf(fpecho,"Poly_to_Unk[%d][%d][%d]=%d\n",
                     pol_number, iseg,ibond,Poly_to_Unk[pol_number][iseg][ibond]);
              }
          }
       }
       fprintf(fpecho,"****************\n END BOND DETAILS \n **********************\n");
    }
    if (Proc==0 && Iwrite_screen==SCREEN_VERBOSE){
       for (pol_number=0; pol_number<Npol_comp; ++pol_number){
          for (iseg=0;iseg<Nmer[pol_number];iseg++){
             printf("SegChain2SegAll[%d][%d]=%d\n", pol_number, iseg,SegChain2SegAll[pol_number][iseg]);
          }
       }
       printf("Total Number of segments in the problem=%d\n",Nseg_tot);
       for (iseg=0;iseg<Nseg_tot;iseg++){
           printf("Nbonds_SegAll[%d]=%d",iseg,Nbonds_SegAll[iseg]);
           printf("\t seg %d is bonded to...",iseg);
           for(ibond=0;ibond<Nbonds_SegAll[iseg];ibond++)
                printf("%d  ",Bonds_SegAll[iseg][ibond]); 
           printf("\n");
       }
    }
    
    if (Type_poly != WTC){  /*POLYMER INPUT FOR EITHER CMS OR WJDC FUNCTIONAL */
    /* set start value of Geqns for each of the polymers in the system.  
       It is necessary to account for Ncomp Boltz and Ncomp Rho eqns */
   
       Ngeqn_tot=0; 
       for (pol_number=0; pol_number<Npol_comp; ++pol_number){
          Geqn_start[pol_number] = 0;
          Ngeqn_tot += (nbond_tot[pol_number]);
          for (pol_num2=0; pol_num2<pol_number; pol_num2++)
              Geqn_start[pol_number] += (nbond_tot[pol_num2]);
       }
       safe_free((void *)  &nbond_tot); 
       if (Iwrite_screen==SCREEN_VERBOSE && Proc==0) printf("The total number of g equations will be %d\n",Ngeqn_tot);
       for (pol_number=0; pol_number<Npol_comp; ++pol_number)
       if (Iwrite_screen==SCREEN_VERBOSE && Proc==0) printf("The start unknown for polymer %d is %d \n",
                                                       pol_number,Geqn_start[pol_number]);
   }
   return;
}
/*************************************************************************************/
