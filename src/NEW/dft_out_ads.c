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

/***********************************************************************
 * calc_adsorption:  Calculate the excess adsorption:                  *
 *                  Gamma^ex = integral [rho-rho_b]dr                  *
 *                  for this processor.                                */

double calc_adsorption(FILE *fp,double *x,double fac_area,double fac_vol)
{
  double area,ads_icomp,ads_ex_icomp,ads_return=0.0;
  int loc_inode, icomp,loc_i,iwall;
  int i;
  static int first=TRUE;

  for (i=0; i<2; i++) {
       for (icomp=0; icomp<Ncomp; icomp++){
            Ads[icomp][i] = 0.0;
            Ads_ex[icomp][i] = 0.0;
       }
  }

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
    for (icomp=0; icomp<Ncomp; icomp++){
      for (i=0; i<Imax; i++){

         if (Sten_Type[POLYMER_CR])
           loc_i = Aztec.update_index[Ncomp+icomp+Nunk_per_node * loc_inode];
         else
           loc_i = Aztec.update_index[icomp + Nunk_per_node * loc_inode];

	 /* Note: if change def. of ads, must change polymer free energy also */
          Ads_ex[icomp][i] += (x[loc_i]*Nel_hit2[i][loc_i]-Rho_b[icomp]*Nel_hit[i][loc_i])
	  *Vol_el/((double)Nnodes_per_el_V); 
	  Ads[icomp][i] += (x[loc_i]*Nel_hit2[i][loc_i])*Vol_el/((double)Nnodes_per_el_V);

      }
    }   /* end of icomp loop */
  }       /* end of loc_inode loop */

  /* get sum of Ads[icomp] and if you're proc #0 print it !! */
  area = 0.0;
  if (Nwall == 0) area = 1.0;
  else{
     if (Nlink == Nwall) area = S_area_tot[Nlists_HW-1][0];
     else
        for (iwall=0; iwall<Nwall; iwall++){
           if (Link[iwall]==0)
           area += S_area_tot[Nlists_HW-1][iwall];
        }
  }


  for (i=0; i<Imax; i++){
     for (icomp=0; icomp<Ncomp; icomp++){
         ads_icomp = AZ_gsum_double(Ads[icomp][i],Aztec.proc_config);
         ads_ex_icomp = AZ_gsum_double(Ads_ex[icomp][i],Aztec.proc_config);
         Ads[icomp][i] = ads_icomp;
         Ads_ex[icomp][i] = ads_ex_icomp;

         if (Lper_area && area>0.0) {
             Ads[icomp][i]/= area;
             Ads[icomp][i]*= (fac_vol/fac_area);
             Ads_ex[icomp][i]/= area;
             Ads_ex[icomp][i]*= (fac_vol/fac_area);
         }
         else{
             Ads[icomp][i]*= fac_vol;
             Ads_ex[icomp][i]*= fac_vol;
         }

         if (i==Imax-1 && icomp==0) ads_return = Ads_ex[icomp][i];
     }
  }

  if (Proc == 0 && fp!=NULL) {
     if (Iwrite != NO_SCREEN){
        printf("\n Summary of Adsorption Results\n");
        printf("\t===== Integrating over domain volume =====\n");
     }
     if (Lhard_surf) i=1;
     else i=0;
     if (Iwrite != NO_SCREEN){
     for (icomp=0;icomp<Ncomp; icomp++)
        printf("\t\t Ads[icomp=%d]=%9.6f\n",icomp,Ads[icomp][i]);
     for (icomp=0;icomp<Ncomp; icomp++)
        printf("\t\t Ads_ex[icomp=%d]=%9.6f\n",icomp,Ads_ex[icomp][i]);
        printf("\t=========================================\n");
     if (Ipot_ff_n != IDEAL_GAS && Lhard_surf){
        printf("\t===== Integrating over fluid volume only =====\n");
        for (icomp=0;icomp<Ncomp; icomp++)
           printf("\t\t Ads[icomp=%d]=%9.6f\n",icomp,Ads[icomp][0]);
        for (icomp=0;icomp<Ncomp; icomp++)
           printf("\t\t Ads_ex[icomp=%d]=%9.6f\n",icomp,Ads_ex[icomp][0]);
        printf("\t=========================================\n");
     }
     }
     if (first){
       for (icomp=0;icomp<Ncomp; icomp++) fprintf(fp,"ads[%d]=%9.6f ",icomp,Ads[icomp][i]);
       for (icomp=0;icomp<Ncomp; icomp++) fprintf(fp,"ads_ex[%d]=%9.6f ",icomp,Ads_ex[icomp][i]);
       fprintf(fp," area=%9.4f fac_area=%9.4f fac_vol=%9.4f",area,fac_area,fac_vol);
       first=FALSE;
     }
     else{
       for (icomp=0;icomp<Ncomp; icomp++) fprintf(fp,"%9.6f ",Ads[icomp][i]);
       for (icomp=0;icomp<Ncomp; icomp++) fprintf(fp,"%9.6f ",Ads_ex[icomp][i]);
       fprintf(fp," %9.4f %9.4f %9.4f",area,fac_area,fac_vol);
     }
   }

  return(ads_return);
}
/***************************************************************************
 * calc_surface_charge:  Calculate the total surface charge based on the   *
 *                       principle of charge neutrality. i.e. the charge   *
 *                       on the surface = the charge in the fluid !!       */

void calc_surface_charge(FILE *fp, double *x,double fac_area,double fac_vol)
{
  double charge_sum,area;
  int loc_inode, icomp, iel, loc_i,nel_hit,idim,ilist,
      inode,ielement,iwall,ijk[3], iel_box=0;
  int reflect_flag[3],semiperm,jcomp;
  static int first=TRUE;
  

  charge_sum = 0.0;
  for (idim=0; idim<Ndim; idim++) reflect_flag[idim]=FALSE;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
      /* convert local node to global */
      inode = Aztec.update[Nunk_per_node * loc_inode] / Nunk_per_node;
      node_to_ijk(inode,ijk);

      for (icomp=0; icomp<Ncomp; icomp++){

         if (Nlists_HW==1 || Nlists_HW==2) ilist = 0;
         else                              ilist = icomp;


         nel_hit = Nnodes_per_el_V;

         for (iel=0; iel<Nnodes_per_el_V; iel++){
             ielement = node_to_elem(inode,iel,reflect_flag);
             if (ielement >= 0) iel_box = el_to_el_box(ielement);
             if (ielement != -2){
                iwall = Wall_elems[ilist][iel_box];
                semiperm=FALSE;
                for (jcomp=0; jcomp<Ncomp; jcomp++) 
                    if (iwall>=0 &&Lsemiperm[WallType[iwall]][jcomp]) semiperm=TRUE;
                if (ielement == -1 || 
                    (iwall != -1 && !semiperm )) nel_hit--;
             }
         }
         for (idim=0; idim<Ndim; idim++){
            if ( (Type_bc[idim][0] == REFLECT || Type_bc[idim][0] == IN_BULK || Type_bc[idim][0]==LAST_NODE) 
                                               &&  ijk[idim] == 0)  nel_hit /= 2;
            if ( (Type_bc[idim][1] == REFLECT || Type_bc[idim][1] == IN_BULK || Type_bc[idim][0]==LAST_NODE) 
                                 &&  ijk[idim] == Nodes_x[idim]-1)  nel_hit /= 2;
         }

         if (Type_poly==-1) loc_i = Aztec.update_index[icomp + Nunk_per_node * loc_inode];
         else               loc_i = Aztec.update_index[icomp + Ncomp + Nunk_per_node * loc_inode];

         charge_sum += (Charge_f[icomp]*x[loc_i])
                       *Vol_el*((double)nel_hit)/((double)Nnodes_per_el_V);

      }   /* end of icomp loop */
  }       /* end of loc_inode loop */

  /* get sum of charge_sum and if you're proc #0 print it !! */
  charge_fluid = AZ_gsum_double(charge_sum,Aztec.proc_config);

  area = 0.0;
  if (Nwall == 0) area = 1.0;
  else{
     if (Nlink == Nwall)
       area = S_area_tot[Nlists_HW-1][0];
     else
       for (iwall=0; iwall<Nwall; iwall++){
          if (Link[iwall]==0)
             area += S_area_tot[Nlists_HW-1][iwall];
       }
  }

  if (Lper_area && area > 0.0) {
      charge_fluid /= area;
      charge_fluid *= (fac_vol/fac_area);
  }
  else charge_fluid *= fac_vol;

  if (Proc == 0) {
     printf("\n----------------------------------------------------------\n");
     printf(" charge in fluid (after vol/area factors): %9.6f\n",charge_fluid); 
     if(first){
        fprintf(fp,"fluid charge=%9.6f   ",charge_fluid);
        first=FALSE;
     }
     else{
        fprintf(fp," %9.6f ",charge_fluid);
     }
     printf("----------------------------------------------------------\n");
  }
}
/****************************************************************************/
