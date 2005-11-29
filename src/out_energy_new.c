#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

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
       omega_CMS,omega_CMS_b,omega_CMS_surf_ex,
       omega_mu,omega_mu_b,omega_mu_surf_ex;
static int first=TRUE,loc_inode;

  if (Ndim==1 && Iwrite==VERBOSE){
    Integration_profile = (double *) array_alloc(1, Nnodes_per_proc, sizeof(double));
    for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
      Integration_profile[loc_inode]=0.0;
    }
  }
  else Integration_profile=NULL;

  if (Proc==0&&Iwrite != NO_SCREEN) printf("---------------------- FREE ENERGY ----------------------------\n");

    omega_sum=0.0;
    omega_s_sum=0.0;

    if (Type_poly==NONE || Type_poly ==WTC){
                                    /* IDEAL GAS CONTRIBUTIONS */
       integrateInSpace_SumInComp(&integrand_ideal_gas_freen,Nel_hit2,x,Integration_profile);
       omega_id=Temporary_sum;
       if (Proc==0 && Iwrite != NO_SCREEN){
             print_to_screen(omega_id,"IDEAL GAS");
        }

       integrateInSpace_SumInComp(&integrand_ideal_gas_freen_bulk,Nel_hit,x,Integration_profile);
       omega_id_b=Temporary_sum;
       omega_id_surf_ex = omega_id-omega_id_b;
       if (Proc==0 && Iwrite != NO_SCREEN){
              print_to_screen(omega_id_surf_ex,"SURF.EX.: IDEAL GAS");
       }
       omega_sum += omega_id_b;
       omega_s_sum += omega_id_surf_ex;
    
                                    /* CHEMICAL POTENTIAL CONTRIBUTIONS */
       integrateInSpace_SumInComp(&integrand_mu_freen,Nel_hit2,x,Integration_profile);
       omega_mu=Temporary_sum;
       if (Proc==0 && Iwrite != NO_SCREEN){
           print_to_screen(omega_mu,"CHEM.POTENTIAL");
       }

       integrateInSpace_SumInComp(&integrand_mu_freen_bulk,Nel_hit,x,Integration_profile);
       omega_mu_b=Temporary_sum;
       omega_mu_surf_ex = omega_mu-omega_mu_b;
       if (Proc==0 && Iwrite != NO_SCREEN){
            print_to_screen(omega_mu_surf_ex,"SURF.EX.: CHEM.POTENTIAL");
       }
       omega_sum += omega_mu_b;
       omega_s_sum += omega_mu_surf_ex;

                                    /* NEUTRAL EXTERNAL FIELD CONTRIBUTIONS */
       integrateInSpace_SumInComp(&integrand_vext_freen,Nel_hit2,x,Integration_profile);
       omega_vext=Temporary_sum;
       if (Proc==0 && Iwrite != NO_SCREEN){
             print_to_screen(omega_vext,"NEUTRAL EXT.FIELD");
       }

       omega_vext_surf_ex = omega_vext; /* note Vext=0 in the bulk */
       if (Proc==0 && Iwrite != NO_SCREEN){
             print_to_screen(omega_vext_surf_ex,"SURF.EX.: NEUTRAL EXT.FIELD");
       }
       omega_sum += omega_vext;
       omega_s_sum += omega_vext_surf_ex;

                                    /* HARD SPHERE CONTRIBUTIONS */
       if (Type_func != NONE){
          integrateInSpace(&integrand_hs_freen,0,Nel_hit,x,Integration_profile);
          omega_hs=Temporary_sum;
          if (Proc==0 && Iwrite != NO_SCREEN){
               print_to_screen(omega_hs,"HARD SPHERE TERM");
          }

          integrateInSpace(&integrand_hs_freen_bulk,0,Nel_hit,x,Integration_profile);
          omega_hs_b=Temporary_sum;
          omega_hs_surf_ex = omega_hs-omega_hs_b;
          if (Proc==0 && Iwrite != NO_SCREEN){
              print_to_screen(omega_hs_surf_ex,"SURF.EX.: HARD SPHERE TERM");
          }
          omega_sum += omega_hs;
          omega_s_sum += omega_hs_surf_ex;
      }

                                    /* ATTRACTION CONTRIBUTIONS */
      if (Type_attr != NONE){
         integrateInSpace_SumInComp(&integrand_att_freen,Nel_hit2,x,Integration_profile);
         omega_att = Temporary_sum;
         if (Proc==0 && Iwrite != NO_SCREEN){
               print_to_screen(omega_att,"ATTRACTIONS");
         }

         integrateInSpace_SumInComp(&integrand_att_freen_bulk,Nel_hit,x,Integration_profile);
         omega_att_b = Temporary_sum;
         omega_att_surf_ex = omega_att-omega_att_b;
         if (Proc==0 && Iwrite != NO_SCREEN){
               print_to_screen(omega_att_surf_ex,"SURF.EX.: ATTRACTIONS");
         }

         omega_sum += omega_att;
         omega_s_sum += omega_att_surf_ex;
      }

                                    /* COULOMB SYSTEMS */
      if (Type_coul != NONE){
                                    /* POINT CHARGE CONTRIBUTIONS */
         integrateInSpace_SumInComp(&integrand_elec_PB_freen,Nel_hit2,x,Integration_profile);
         omega_psirho = Temporary_sum;
         omega_psirho_surf_ex = omega_psirho; /* note elec. pot.=0 in the bulk */
         if (Proc==0 && Iwrite != NO_SCREEN){
             print_to_screen(omega_psirho,"PSI*RHO ELEC TERM");
         }

         omega_sum += omega_psirho;
         omega_s_sum += omega_psirho_surf_ex;

                                 /* CHARGED EXTERNAL FIELD CONTRIBUTIONS */
         if (Vext_coul != NULL){
         integrateInSpace_SumInComp(&integrand_vext_elec_freen,Nel_hit2,x,Integration_profile);
         omega_vext_elec=Temporary_sum;
         if (Proc==0 && Iwrite != NO_SCREEN){
              print_to_screen(omega_vext_elec,"CHARGED EXT.FIELD");
         }

         omega_vext_elec_surf_ex = omega_vext_elec; /* note Vext=0 in the bulk */
         if (Proc==0 && Iwrite != NO_SCREEN){
              print_to_screen(omega_vext_elec_surf_ex,"SURF.EX.: CHARGED EXT.FIELD");
         }
         omega_sum += omega_vext_elec;
         omega_s_sum += omega_vext_elec_surf_ex;
         }

         if (Type_coul == DELTAC){     /* MSA CORRECTIONS FOR ELECTROLYTES */
            integrateInSpace_SumInComp(&integrand_elec_MSAcorr_freen,Nel_hit2,x,Integration_profile);
            omega_MSA = Temporary_sum;
            if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(omega_MSA,"MSA CORRECTIONS");
  
            integrateInSpace_SumInComp(&integrand_elec_MSAcorr_freen_bulk,Nel_hit,x,Integration_profile);
            omega_MSA_b = Temporary_sum;
            omega_MSA_surf_ex = omega_MSA-omega_MSA_b;
            if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(omega_MSA_surf_ex,"SURF.EX.: MSA CORRECTIONS");
  
            omega_sum += omega_MSA;
            omega_s_sum += omega_MSA_surf_ex;
         }
      }

                                    /* WTC CONTRIBUTIONS */
       if (Type_poly == WTC){
          integrateInSpace_SumInComp(&integrand_WTC_freen,Nel_hit2,x,Integration_profile);
          omega_WTC = Temporary_sum;
          if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(omega_WTC,"WTC BONDS");

          integrateInSpace_SumInComp(&integrand_WTC_freen_bulk,Nel_hit,x,Integration_profile);
          omega_WTC_b = Temporary_sum;
          omega_WTC_surf_ex = omega_WTC-omega_WTC_b;
          if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(omega_WTC_surf_ex,"SURF.EX.: WTC BONDS");

          omega_sum += omega_WTC;
          omega_s_sum += omega_WTC_surf_ex;
       }
      if (Proc==0 && Iwrite != NO_SCREEN){
        printf("\t----------------------------------------\n");
        print_to_screen(omega_sum,"TOTAL GRAND POTENTIAL");
        print_to_screen(omega_s_sum,"TOTAL SURF EX FREE ENERGY");
        if (fp !=NULL) print_to_file(fp,omega_sum,"omega",first);
        if (fp !=NULL) print_to_file(fp,omega_s_sum,"omega_s",first);
        printf("---------------------------------------------------------------\n");
      }
    }

                                    /* CMS FREE ENERGY */
    if (Type_poly == CMS || Type_poly==CMS_GAUSSIAN || Type_poly==CMS_SCFT){
       integrateInSpace_SumInComp(&integrand_CMS_freen,Nel_hit2,x,Integration_profile);
       omega_CMS = Temporary_sum;
       if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(omega_CMS,"CMS FREE ENERGY");

       integrateInSpace_SumInComp(&integrand_CMS_freen_bulk,Nel_hit2,x,Integration_profile);
       omega_CMS_b = Temporary_sum;
       omega_CMS_surf_ex = omega_CMS-omega_CMS_b;
       if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(omega_CMS_surf_ex,"CMS FREE ENERGY DELTA");

       omega_sum += omega_CMS;
       omega_s_sum += omega_CMS_surf_ex;

       if (Proc==0 && Iwrite != NO_SCREEN){
           printf("\t----------------------------------------\n");
           print_to_screen(omega_s_sum,"FREE ENERGY REL TO BULK");
           if (fp !=NULL) print_to_file(fp,omega_s_sum,"del_omega",first);
           printf("---------------------------------------------------------------\n");
       }
    }

    if (first) first=FALSE;

    if (Integration_profile != NULL){
      print_freen_profile_1D(Integration_profile,"dft_freen_prof.dat");
      safe_free ((void *)&Integration_profile);
    }
  return;
}
