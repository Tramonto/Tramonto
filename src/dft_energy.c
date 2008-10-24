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

#include "dft_energy.h"

/****************************************************************************/
double calc_free_energy(FILE *fp, double **x)
{
double omega_sum, omega_s_sum, omega_id, omega_id_b,omega_id_surf_ex,
       omega_hs,omega_hs_b,omega_hs_surf_ex,
       omega_att,omega_att_b,omega_att_surf_ex,
       omega_wtc,omega_wtc_b,omega_wtc_surf_ex,
       omega_psirho,omega_psirho_surf_ex,
       omega_MSA,omega_MSA_b,omega_MSA_surf_ex,
       omega_vext,omega_vext_surf_ex,
       omega_vext_elec,omega_vext_elec_surf_ex,
       omega_WTC,omega_WTC_b,omega_WTC_surf_ex,
       omega_WJDC,omega_WJDC_b,omega_WJDC_surf_ex,
       omega_CMS,omega_CMS_b,omega_CMS_surf_ex,
       omega_maxwell_stress,omega_surface_charge,
       omega_osmotic, omega_osmotic_b,omega_osmotic_surf_ex,
       omega_mu,omega_mu_b,omega_mu_surf_ex;
static int first=TRUE,loc_inode;
double energy,volume;
int iunk,idim;
int lfirst;

       double L,L1,L2,sum,energy_RR,energy_RR_plus,energy_RR_minus,derivative_neg,lambda;
       int icomp;


  if (Ndim==1 && Iwrite==VERBOSE){
    Integration_profile = (double *) array_alloc(1, Nnodes_per_proc, sizeof(double));
    for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
      Integration_profile[loc_inode]=0.0;
    }
  }
  else Integration_profile=NULL;

  if (!first && Proc==0 && Iwrite != NO_SCREEN) printf("---------------------- FREE ENERGY ----------------------------\n");

    omega_sum=0.0;
    omega_s_sum=0.0;

    /* lfirst: extra flag to keep track when we call this routine from somewhere 
       other than dft_out_main.c, in which case fp == NULL */
    lfirst = FALSE;
    if (Proc==0 && fp!=NULL) lfirst = TRUE;
    MPI_Bcast(&lfirst,1,MPI_INT,0,MPI_COMM_WORLD);

    if (L_HSperturbation){
                                    /* IDEAL GAS CONTRIBUTIONS */

      /* if we are calling this for the first time from the post-processing, don't
	 calculate all the integrals; otherwise do calculate */
      if(!first || !lfirst) {

       if (Type_poly != WJDC && Type_poly !=WJDC2 && Type_poly!=WJDC3){
          omega_id=omega_id_b=0.0;
          for (iunk=Phys2Unk_first[DENSITY];iunk<Phys2Unk_last[DENSITY];iunk++) {
   	     if (Lseg_densities) icomp=Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
   	     else                icomp=iunk-Phys2Unk_first[DENSITY];
             if (fabs(Charge_f[icomp])<1.e-12 || Type_coul==NONE){
                 omega_id+=integrateInSpace(&integrand_ideal_gas_freen,iunk,Nel_hit2,x,Integration_profile);

                 omega_id_b+=integrateInSpace(&integrand_ideal_gas_freen_bulk,iunk,Nel_hit,x,Integration_profile);
             }
           }
          omega_id_surf_ex = omega_id-omega_id_b;
          if (Proc==0 && Iwrite != NO_SCREEN){
                print_to_screen(omega_id,"IDEAL GAS");
                print_to_screen(omega_id_b,"BULK TERM: IDEAL GAS");
                print_to_screen(omega_id_surf_ex,"SURF.EX.: IDEAL GAS");
           }
          omega_sum += omega_id;
          omega_s_sum += omega_id_surf_ex;

                                    /* CHEMICAL POTENTIAL CONTRIBUTIONS */

          omega_mu=integrateInSpace_SumInComp(&integrand_mu_freen,Nel_hit2,x,Integration_profile);
          if (Proc==0 && Iwrite != NO_SCREEN){
              print_to_screen(omega_mu,"CHEM.POTENTIAL");
          }

          omega_mu_b=integrateInSpace_SumInComp(&integrand_mu_freen_bulk,Nel_hit,x,Integration_profile);
          omega_mu_surf_ex = omega_mu-omega_mu_b;
          if (Proc==0 && Iwrite != NO_SCREEN){
               print_to_screen(omega_mu_b,"BULK TERM: CHEM.POTENTIAL");
               print_to_screen(omega_mu_surf_ex,"SURF.EX.: CHEM.POTENTIAL");
          }
          omega_sum += omega_mu;
          omega_s_sum += omega_mu_surf_ex;
       }

                                 /* NEUTRAL EXTERNAL FIELD CONTRIBUTIONS */
       if (Nwall!=0){
          omega_vext=integrateInSpace_SumInComp(&integrand_vext_freen,Nel_hit2,x,Integration_profile);
          if (Proc==0 && Iwrite != NO_SCREEN){
                print_to_screen(omega_vext,"NEUTRAL EXT.FIELD");
          }

          omega_vext_surf_ex = omega_vext; /* note Vext=0 in the bulk */
          if (Proc==0 && Iwrite != NO_SCREEN){
                print_to_screen(omega_vext_surf_ex,"SURF.EX.: NEUTRAL EXT.FIELD");
          }
          omega_sum += omega_vext;
          omega_s_sum += omega_vext_surf_ex;
       }

                                    /* HARD SPHERE CONTRIBUTIONS */
       if (Type_func != NONE){
          omega_hs=integrateInSpace(&integrand_hs_freen,0,Nel_hit,x,Integration_profile);
          if (Proc==0 && Iwrite != NO_SCREEN){
               print_to_screen(omega_hs,"HARD SPHERE TERM");
          }

          omega_hs_b=integrateInSpace(&integrand_hs_freen_bulk,0,Nel_hit,x,Integration_profile);
          omega_hs_surf_ex = omega_hs-omega_hs_b;
          if (Proc==0 && Iwrite != NO_SCREEN){
              print_to_screen(omega_hs_surf_ex,"SURF.EX.: HARD SPHERE TERM");
          }
          omega_sum += omega_hs;
          omega_s_sum += omega_hs_surf_ex;
      }

                                    /* ATTRACTION CONTRIBUTIONS */
      if (Type_attr != NONE){
         omega_att=integrateInSpace_SumInComp(&integrand_att_freen,Nel_hit2,x,Integration_profile);
         if (Proc==0 && Iwrite != NO_SCREEN){
               print_to_screen(omega_att,"MEAN_FIELD_POTENTIAL");
         }

         omega_att_b=integrateInSpace_SumInComp(&integrand_att_freen_bulk,Nel_hit,x,Integration_profile);
         if (Proc==0 && Iwrite != NO_SCREEN){
               print_to_screen(omega_att_b,"BULK: MEAN_FIELD_POTENTIAL");
         }
         omega_att_surf_ex = omega_att-omega_att_b;
         if (Proc==0 && Iwrite != NO_SCREEN){
               print_to_screen(omega_att_surf_ex,"SURF.EX.: MEAN_FIELD_POTENTIAL");
         }

         omega_sum += omega_att;
         omega_s_sum += omega_att_surf_ex;
      }

                                    /* COULOMB SYSTEMS */
      if (Type_coul != NONE){

         /* Here we do some very simple analytical things to see if we can reproduce the sum rule of the most
            simple electrostatic systems */

            sum=0;
            for (icomp=0; icomp<Ncomp;icomp++) sum+=Charge_f[icomp]*Charge_f[icomp]*Rho_b[icomp];
            lambda = sqrt(Temp_elec/(4*PI*sum));
	    if (Proc==0 && Iwrite != NO_SCREEN) {
	      print_to_screen(lambda, "Debye wave length");
	    }
            /* printf("Debye wave length: lambda=%9.6f\n",lambda);*/

         /* Reiner and Radke method for computing the free energy of a PB electrolyte near a charged surface */

                /* Maxwell Stress Term */
         omega_maxwell_stress=integrateInSpace(&integrand_maxwell_stress_freen,0,Nel_hit,x,Integration_profile); 
         if (Proc==0 && Iwrite != NO_SCREEN){
             print_to_screen(omega_maxwell_stress,"MAXWELL STRESS TERM");
         }
         omega_sum += omega_maxwell_stress;
         omega_s_sum += omega_maxwell_stress;             /* note that the maxwell stress in the bulk is 0 */

                /* osmotic pressure contribution ...only for charged species */

         omega_osmotic=0.;
         omega_osmotic_b=0.;
         for (iunk=Phys2Unk_first[DENSITY];iunk<Phys2Unk_last[DENSITY];iunk++) {
           if (fabs(Charge_f[iunk-Phys2Unk_first[DENSITY]])>1.e-15){
              omega_osmotic-=integrateInSpace(&integrand_adsorption,iunk,Nel_hit2,x,Integration_profile);

              omega_osmotic_b-=integrateInSpace(&integrand_adsorption_bulk,iunk,Nel_hit,x,Integration_profile);
           }
         }
         omega_osmotic_surf_ex=omega_osmotic-omega_osmotic_b;

         if (Proc==0 && Iwrite != NO_SCREEN){
             print_to_screen(omega_osmotic,"OSMOTIC PRESSURE TERM");
             print_to_screen(omega_osmotic_surf_ex,"SURF.EX.: OSMOTIC PRESSURE");
         }
         omega_sum += omega_osmotic;
         omega_s_sum += omega_osmotic_surf_ex; 

                /* surface charge term */

         omega_surface_charge=integrateOverSurface(&integrand_surface_charge,iunk,x,Integration_profile);

         if (Proc==0 && Iwrite != NO_SCREEN){
               print_to_screen(omega_surface_charge,"SURFACE CHARGE TERM");
         }
         omega_sum += omega_surface_charge;
         omega_s_sum += omega_surface_charge;

                                    /* PSI-RHO TERM FROM TANG-DAVIS PAPER */
/*       omega_psirho=0.5*integrateInSpace_SumInComp(&integrand_elec_PB_freen,Nel_hit2,x,Integration_profile);
         omega_psirho_surf_ex = omega_psirho; 
         if (Proc==0 && Iwrite != NO_SCREEN){
             print_to_screen(omega_psirho,"PSI*RHO ELEC TERM");
         }

         omega_sum += omega_psirho;
         omega_s_sum += omega_psirho_surf_ex;*/

                                 /* CHARGED EXTERNAL FIELD CONTRIBUTIONS */
         /* term 2 based on Tang and Davis papers for electrostatics 
         if (Vext_coul != NULL){
         omega_vext_elec=integrateInSpace_SumInComp(&integrand_vext_elec_freen,Nel_hit2,x,Integration_profile);
         if (Proc==0 && Iwrite != NO_SCREEN){
              print_to_screen(omega_vext_elec,"CHARGED EXT.FIELD");
         }

         omega_vext_elec_surf_ex = omega_vext_elec; * note Vext=0 in the bulk *
         if (Proc==0 && Iwrite != NO_SCREEN){
              print_to_screen(omega_vext_elec_surf_ex,"SURF.EX.: CHARGED EXT.FIELD");
         }
         omega_sum += omega_vext_elec;
         omega_s_sum += omega_vext_elec_surf_ex;
         }
         */

         /* term 3 based on Tang and Davis papers for electrostatics 
         if (Type_coul == DELTAC){     * MSA CORRECTIONS FOR ELECTROLYTES *
            omega_MSA=integrateInSpace_SumInComp(&integrand_elec_MSAcorr_freen,Nel_hit2,x,Integration_profile);
            if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(omega_MSA,"MSA CORRECTIONS");
  
            omega_MSA_b=integrateInSpace_SumInComp(&integrand_elec_MSAcorr_freen_bulk,Nel_hit,x,Integration_profile);
            omega_MSA_surf_ex = omega_MSA-omega_MSA_b;
            if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(omega_MSA_surf_ex,"SURF.EX.: MSA CORRECTIONS");
  
            omega_sum += omega_MSA;
            omega_s_sum += omega_MSA_surf_ex;
         }
         */
      }

                                    /* WTC CONTRIBUTIONS */
       if (Type_poly == WTC){
          omega_WTC=integrateInSpace_SumInComp(&integrand_WTC_freen,Nel_hit2,x,Integration_profile);
          if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(omega_WTC,"WTC BONDS");

          omega_WTC_b=integrateInSpace_SumInComp(&integrand_WTC_freen_bulk,Nel_hit,x,Integration_profile);
          omega_WTC_surf_ex = omega_WTC-omega_WTC_b;
          if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(omega_WTC_surf_ex,"SURF.EX.: WTC BONDS");

          omega_sum += omega_WTC;
          omega_s_sum += omega_WTC_surf_ex;
       }

       if (Type_poly == WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){

          if (Lseg_densities){
             omega_WJDC=integrateInSpace_SumInComp(&integrand_WJDC_freen,Nel_hit2,x,Integration_profile);
             if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(omega_WJDC,"WJDC BONDS");

             omega_WJDC_b=integrateInSpace_SumInComp(&integrand_WJDC_freen_bulk,Nel_hit,x,Integration_profile);
             omega_WJDC_surf_ex = omega_WJDC-omega_WJDC_b;
             if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(omega_WJDC_surf_ex,"SURF.EX.: WJDC-DTERM");
          }
          else{
             omega_WJDC=integrateInSpace_SumInComp(&integrand_WJDCcomp_freen,Nel_hit2,x,Integration_profile);
             if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(omega_WJDC,"WJDC BONDS");

             omega_WJDC_b=integrateInSpace_SumInComp(&integrand_WJDCcomp_freen_bulk,Nel_hit,x,Integration_profile);
             omega_WJDC_surf_ex = omega_WJDC-omega_WJDC_b;
             if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(omega_WJDC_surf_ex,"SURF.EX.: WJDC-DTERM");
          }

          omega_sum += omega_WJDC;
          omega_s_sum += omega_WJDC_surf_ex;
       }

      } /* end of if(!first || !lfirst) */

      if (Proc==0 && Iwrite != NO_SCREEN){
        if(!first || !lfirst) {
        printf("\t----------------------------------------\n");
        print_to_screen(omega_sum,"TOTAL GRAND POTENTIAL");
        print_to_screen(omega_s_sum,"TOTAL SURF EX FREE ENERGY");
        printf("---------------------------------------------------------------\n");
        }
      }

	if (Proc==0 && Iwrite != NO_SCREEN){
                volume=1.0;
                for (idim=0;idim<Ndim; idim++) volume*=Size_x[0];
	        if (fp !=NULL && LBulk) print_to_file(fp,-omega_sum/volume,"pressure",first);
	        if (fp !=NULL && !LBulk) { 
                    print_to_file(fp,omega_sum,"omega",first);
		    if (Type_interface == UNIFORM_INTERFACE) print_to_file(fp,omega_s_sum,"omega_s",first);
                 }
	}
    }

                                    /* CMS FREE ENERGY */
    if (Type_poly == CMS || Type_poly==CMS_SCFT){
      if(!first || !lfirst) {
       omega_CMS=integrateInSpace_SumInComp(&integrand_CMS_freen,Nel_hit2,x,Integration_profile);
       if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(omega_CMS,"CMS FREE ENERGY");

       omega_CMS_b=integrateInSpace_SumInComp(&integrand_CMS_freen_bulk,Nel_hit2,x,Integration_profile);
       omega_CMS_surf_ex = omega_CMS-omega_CMS_b;
       if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(omega_CMS_surf_ex,"CMS FREE ENERGY DELTA");

       omega_sum += omega_CMS;
       omega_s_sum += omega_CMS_surf_ex;
      }

       if (Proc==0 && Iwrite != NO_SCREEN){
        if(!first || !lfirst) {
           printf("\t----------------------------------------\n");
           print_to_screen(omega_s_sum,"FREE ENERGY REL TO BULK");
           printf("---------------------------------------------------------------\n");
         }
       }
    
       if (Proc==0 && Iwrite != NO_SCREEN){
	 if (fp !=NULL) print_to_file(fp,omega_s_sum,"del_omega",first);
       }
    }

    energy = omega_s_sum;

    if (first && lfirst) first=FALSE;
 

    if (Integration_profile != NULL){
      print_freen_profile_1D(Integration_profile,"dft_freen_prof.dat");
      safe_free ((void *)&Integration_profile);
    }

  return(energy);
}
