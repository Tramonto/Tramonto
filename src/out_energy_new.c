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
       omega_maxwell_stress,omega_surface_charge,
       omega_osmotic, omega_osmotic_b,omega_osmotic_surf_ex,
       omega_mu,omega_mu_b,omega_mu_surf_ex;
static int first=TRUE,loc_inode;
double energy;
int iunk;

       double L,L1,L2,sum,energy_RR,energy_RR_plus,energy_RR_minus,derivative_neg,lambda;
       int icomp;

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

       omega_id=omega_id_b=0.0;
       for (iunk=Phys2Unk_first[DENSITY];iunk<Phys2Unk_last[DENSITY];iunk++) {
          if (fabs(Charge_f[iunk-Phys2Unk_first[DENSITY]])<1.e-15){
              integrateInSpace(&integrand_ideal_gas_freen,iunk,Nel_hit2,x,Integration_profile);
              omega_id+=Temporary_sum;

              integrateInSpace(&integrand_ideal_gas_freen_bulk,iunk,Nel_hit,x,Integration_profile);
              omega_id_b+=Temporary_sum;
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

/*       if (Type_coul ==NONE){
       integrateInSpace_SumInComp(&integrand_ideal_gas_freen,Nel_hit2,x,Integration_profile);
       omega_id=Temporary_sum;
       if (Proc==0 && Iwrite != NO_SCREEN){
             print_to_screen(omega_id,"IDEAL GAS");
        }

       integrateInSpace_SumInComp(&integrand_ideal_gas_freen_bulk,Nel_hit,x,Integration_profile);
       omega_id_b=Temporary_sum;
       omega_id_surf_ex = omega_id-omega_id_b;
       if (Proc==0 && Iwrite != NO_SCREEN){
              print_to_screen(omega_id_b,"BULK TERM: IDEAL GAS");
              print_to_screen(omega_id_surf_ex,"SURF.EX.: IDEAL GAS");
       }
       omega_sum += omega_id;
       omega_s_sum += omega_id_surf_ex;
       }*/
    
                                    /* CHEMICAL POTENTIAL CONTRIBUTIONS */
/*       if (Type_coul==NONE){*/

       integrateInSpace_SumInComp(&integrand_mu_freen,Nel_hit2,x,Integration_profile);
       omega_mu=Temporary_sum;
       if (Proc==0 && Iwrite != NO_SCREEN){
           print_to_screen(omega_mu,"CHEM.POTENTIAL");
       }

       integrateInSpace_SumInComp(&integrand_mu_freen_bulk,Nel_hit,x,Integration_profile);
       omega_mu_b=Temporary_sum;
       omega_mu_surf_ex = omega_mu-omega_mu_b;
       if (Proc==0 && Iwrite != NO_SCREEN){
            print_to_screen(omega_mu_b,"BULK TERM: CHEM.POTENTIAL");
            print_to_screen(omega_mu_surf_ex,"SURF.EX.: CHEM.POTENTIAL");
       }
       omega_sum += omega_mu;
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

         /* Here we do some very simple analytical things to see if we can reproduce the sum rule of the most
            simple electrostatic systems */

            sum=0;
            for (icomp=0; icomp<Ncomp;icomp++) sum+=Charge_f[icomp]*Charge_f[icomp]*Rho_b[icomp];
            lambda = sqrt(Temp_elec/(4*PI*sum));
            printf("Debye wave length: lambda=%9.6f\n",lambda);


/* simple implementation of the electrostatic free energy in the Debye-Huckel limit according to Reiner and Radke.  Don't 
 * have the DH solution working yet - don't really want it anyway..... */
/*            L=Size_x[0]-2*WallParam[0];
            L1=L-0.001;
            L2=L+0.001;
            printf("cosh_L=%9.6f  cosh_L1=%9.6f  cosh_L2=%9.6f\n",cosh(L/(2*lambda)),cosh(L2/(2*lambda)),cosh(L1/(2*lambda)));
            printf("sinh_L=%9.6f  sinh_L1=%9.6f  sinh_L2=%9.6f\n",sinh(L/(2*lambda)),sinh(L2/(2*lambda)),sinh(L1/(2*lambda)));
            printf("coth_L=%9.6f  coth_L1=%9.6f  coth_L2=%9.6f\n",
                cosh(L/(2*lambda))/sinh(L/(2*lambda)),cosh(L2/(2*lambda))/sinh(L2/(2*lambda)),cosh(L1/(2*lambda))/sinh(L1/(2*lambda)));
    
            energy_RR=(PI*lambda*(4*Elec_param_w[0]*Elec_param_w[0])*(cosh(L/(2*lambda))/sinh(L/(2*lambda)) ))/(2*Temp_elec);
            energy_RR_plus=(PI*lambda*(4*Elec_param_w[0]*Elec_param_w[0])*(cosh(L2/(2*lambda))/sinh(L2/(2*lambda)) ))/(2*Temp_elec);
            energy_RR_minus=(PI*lambda*(4*Elec_param_w[0]*Elec_param_w[0])*(cosh(L1/(2*lambda))/sinh(L1/(2*lambda)) ))/(2*Temp_elec);
            derivative_neg= -(energy_RR_plus-energy_RR_minus)/(L2-L1);

             printf("L=%9.6f L1=%9.6f L2=%9.6f  energy_RR=%g energy_RR_plus=%g energy_RR_minus=%g derivative=%g\n",
                      L,L1,L2,energy_RR,energy_RR_plus,energy_RR_minus,derivative_neg);*/
  
         /* Reiner and Radke method for computing the free energy of a PB electrolyte near a charged surface */

                /* Maxwell Stress Term */
         integrateInSpace(&integrand_maxwell_stress_freen,0,Nel_hit,x,Integration_profile); 
         omega_maxwell_stress = Temporary_sum;
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
              integrateInSpace(&integrand_adsorption,iunk,Nel_hit2,x,Integration_profile);
              omega_osmotic-=Temporary_sum;

              integrateInSpace(&integrand_adsorption_bulk,iunk,Nel_hit,x,Integration_profile);
              omega_osmotic_b-=Temporary_sum;
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
         /* for now, to test this approach, just do a kludge.....need to implement surface integrals rigorously in dft_utils to get forces properly*/
/*         omega_surface_charge=1.58*x[Phys2Unk_first[POISSON]][20];*/
/*         omega_surface_charge=2*(-0.0903125)*x[Phys2Unk_first[POISSON]][10];*/

         printf("trying to call integrate in space for the surface charge integral\n");
         integrateOverSurface(&integrand_surface_charge,iunk,x,Integration_profile);
         printf("after call integrate in space for the surface charge integral\n");
         omega_surface_charge=Temporary_sum;

         if (Proc==0 && Iwrite != NO_SCREEN){
               print_to_screen(omega_surface_charge,"SURFACE CHARGE TERM");
         }
         omega_sum += omega_surface_charge;
         omega_s_sum += omega_surface_charge;

                                    /* POINT CHARGE CONTRIBUTIONS */
         /* term 1 based on Tang and Davis papers for electrostatics 
         integrateInSpace_SumInComp(&integrand_elec_PB_freen,Nel_hit2,x,Integration_profile);
         omega_psirho = Temporary_sum;
         omega_psirho_surf_ex = omega_psirho; * note elec. pot.=0 in the bulk *
         if (Proc==0 && Iwrite != NO_SCREEN){
             print_to_screen(omega_psirho,"PSI*RHO ELEC TERM");
         }

         omega_sum += omega_psirho;
         omega_s_sum += omega_psirho_surf_ex;
         */

                                 /* CHARGED EXTERNAL FIELD CONTRIBUTIONS */
         /* term 2 based on Tang and Davis papers for electrostatics 
         if (Vext_coul != NULL){
         integrateInSpace_SumInComp(&integrand_vext_elec_freen,Nel_hit2,x,Integration_profile);
         omega_vext_elec=Temporary_sum;
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
         */
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
    energy= omega_s_sum;

    if (first) first=FALSE;

    if (Integration_profile != NULL){
      print_freen_profile_1D(Integration_profile,"dft_freen_prof.dat");
      safe_free ((void *)&Integration_profile);
    }
  return(energy);
}
