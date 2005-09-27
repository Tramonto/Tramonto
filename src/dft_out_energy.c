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

#define BFD 0
#define FFD 1
#define CFD 2

double phispt_i(double *);
double phispt_bulk();
double int_stencil(double **,int, int,int);
void assemble_HS_free_energy(double **, double *, double *,double *);
double free_energy_charging_up(double **);
double energy_elec(double **,double *);
double charge_stress(double **,double *);
double calc_deriv_e(int,int,int,int *,double **,int);
double calc_deriv2(int,int,int,double **);
double calc_u_ideal(int, int *, double **, double *, double *);
double energy_elec_vext_vol(double **);

/****************************************************************************/
/****************************************************************************
 *calc_free_energy:  This routine calculates the contribution of processor  *
 *                   Nproc to the grand potential (omega) and the surface   *
 *                   free energy (omega + pV).                              */
/****************************************************************************/

double calc_free_energy(FILE *fp, double **x, double fac_area,
                              double fac_vol, int print_flag)
{
  int    icomp, ilist, loc_inode, inode,iunk,
         psiunk,i,inode_box,muunk;

  int    iwall,ijk[3],jwall;

  int    unk_xi2,unk_xi3,iseg,ibond,jseg,jcomp,unk_bond;
  double  wtc_term,y,wtc_b,wtc_sum[2];

  double rho_i,psi_i=0.0,omega_i,omega_i_b,area,stress_tensor,
         stress_tensor_p,energy_psi_q,int_psi_q,
         stress_in_wall,stress_in_wall_p,
         omega_sum[2],omega_s_sum[2], sum_phispt,sum_phispt_b,
         sum_phispt_b_old,psi_rho,ex_ads_t[2],ads_t[2],
         grand_free_energy[2]={0.0,0.0},surface_free_energy[2]={0.0,0.0},psi_rho_sum[2],
         vext_c_sum[2],hs_energy,ideal_sum[2],lj_sum[2],deltac_sum[2],vext_sum[2],
         ideal,ljterm,deltacterm,vext,ideal_b,lj_b,deltac_b,
         charging_free_energy,energy_charge,rr_omega_s[2],rr_omega_s2[2];
  double elec_energy_vol, elec_energy_vol_vext,prefac;
  double omega_b1,omega_b2,omega_b3;
  double o_id,o_v,o_att,o_mu,Uww_sum,fac_coulomb,vext_c;
  int ilink,jlink;
  static int first=TRUE;
  FILE   *fp10=NULL;

  for (i=0; i<2; i++){
     omega_sum[i] = 0.0;
     omega_s_sum[i] = 0.0;
     ideal_sum[i] = 0.0;
     lj_sum[i] = 0.0;
     deltac_sum[i] = 0.0;
     psi_rho_sum[i] = 0.0;
     vext_c_sum[i] = 0.0;
     vext_sum[i] = 0.0;
  }

  if (Proc==0 &&Iwrite!=NO_SCREEN) printf ("\n Energy Calculation summary before taking area/volume factors into account:\n");
  if (Ipot_ff_n != IDEAL_GAS){
     assemble_HS_free_energy(x, &sum_phispt, &sum_phispt_b, &sum_phispt_b_old);
     if (Proc ==0 && print_flag &&Iwrite!=NO_SCREEN) {
        printf("\t===== Hard Sphere Integrals ... integrating over domain volume =======\n");
        printf("\t Phispt: %g \n",sum_phispt);
        printf("\t Phispt_bulk: %g \n",sum_phispt_b);
        printf("\t===== Bulk Hard Sphere Integral ... integrating over fluid volume =======\n");
        printf("\t Phispt_bulk(fluid vol): %g \n",sum_phispt_b_old);
        printf("\t========================================================================\n");
      }
  }
  else{
     sum_phispt = sum_phispt_b = sum_phispt_b_old =0.0;
  }
  
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){

      inode_box = L2B_node[loc_inode];
      inode = L2G_node[loc_inode];
      node_to_ijk(inode,ijk);
   
      for (iunk=Phys2Unk_first[DENSITY];iunk<Phys2Unk_last[DENSITY];iunk++){ 
         icomp=Unk2Comp[iunk];
         for (i=0; i<Imax; i++){
            ilist = List[i];

	    if (Lsteady_state) muunk = Phys2Unk_first[DIFFUSION]+icomp;
	    if (Ipot_ff_c == 1) psiunk = Phys2Unk_first[POISSON];

            rho_i = x[iunk][inode_box];
            if (Ipot_ff_c == 1) psi_i = x[psiunk][inode_box];

            /* FIRST DETERMINE THE RHO_I TERMS */ 
            if (rho_i > Rho_b[icomp]*exp(-VEXT_MAX) && Vext[loc_inode][icomp] < VEXT_MAX){ 
               if (Lsteady_state){
                   omega_i = rho_i*( log(rho_i)-1.0 
                                 + Vext[loc_inode][icomp] 
                                 - x[muunk][inode_box]  );
                   ideal = rho_i*( log(rho_i)-1.0 - Betamu_RTF[icomp]  );
               }
               else{
                   omega_i = rho_i*( log(rho_i)-1.0 
                                 + Vext[loc_inode][icomp] 
                                 - Betamu[icomp]  );
                   ideal = rho_i*( log(rho_i)-1.0 - Betamu[icomp]  );
               }
               vext = Vext[loc_inode][icomp];

               if (Ipot_ff_n == LJ12_6) {  /* LJ fluid case */
                    ljterm = 0.5*rho_i*int_stencil(x,inode_box,icomp,U_ATTRACT);
                    omega_i += ljterm;
               }
               else ljterm=0.0;

              if (Type_poly_TC){
                    unk_xi2 = Phys2Unk_first[CAVITY_WTC];
                    unk_xi3 = Phys2Unk_first[CAVITY_WTC]+1;
                    iseg = iunk-Phys2Unk_first[DENSITY];
                    wtc_term=0.0;
                    for (ibond=0;ibond<Nbonds_SegAll[iunk];ibond++){
                         jseg = Bonds_SegAll[iunk][ibond];
                         jcomp = Unk2Comp[jseg+Phys2Unk_first[DENSITY]];
                         y = y_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],
                                   x[unk_xi2][inode_box],x[unk_xi3][inode_box]);
                         unk_bond = Poly_to_Unk_SegAll[iseg][ibond]+Phys2Unk_first[BOND_WTC];
                         wtc_term += 0.5*rho_i*(1.0-log(y*x[unk_bond][inode_box])) ;
                    }
                    omega_i += wtc_term;
              }
              else wtc_term=0.0;

   
               if (Ipot_ff_c == 1){ /* Coulomb fluid */
                  psi_rho = 0.5*rho_i*Charge_f[icomp]*psi_i;
                  if (Vol_charge_flag && Ndim==3){ 
                     vext_c = 0.5*rho_i*Vext_coul[loc_inode][icomp];
                  }
                  omega_i += (psi_rho+vext_c);

                  if (Sten_Type[THETA_CHARGE] == TRUE){
                     deltacterm = -0.5*rho_i*int_stencil(x,inode_box,icomp,THETA_CHARGE);
                     omega_i += deltacterm;
                   }
                   else deltacterm=0.0;
               }
               else{
                  psi_rho = 0.0;
                  vext_c=0.0;
               }
            }
            else {
               omega_i = 0.0;
               psi_rho = 0.0;
               vext_c=0.0;
               deltacterm=0.0;
               ljterm=0.0;
               ideal=0.0;
               vext=0.0;
            }
 
            /* NOW DETERMINE THE RHO_B TERMS */
            if (Lsteady_state) ideal_b = Rho_b_RTF[icomp]*(log(Rho_b_RTF[icomp])-1.0 - Betamu_RTF[icomp]);
            else   ideal_b = Rho_b[icomp]*(log(Rho_b[icomp])-1.0 - Betamu[icomp]);
            omega_i_b = ideal_b;
   
            if (Ipot_ff_n == LJ12_6) {  /* LJ fluid case */
               lj_b = 0.5*Rho_b[icomp]*Betamu_att[icomp];  /* note that Betamu_att=Betamu_att_RTF for steady state */
               omega_i_b += lj_b;
            }
            else lj_b=0.0;

            if (Type_poly_TC){
                wtc_b += 0.5*Rho_b[iseg]*Betamu_wtc[iseg];
            }
            else wtc_b=0.0;

            if (Ipot_ff_c == 1 && Sten_Type[THETA_CHARGE] == TRUE){
                if(Lsteady_state) deltac_b = -0.5*Rho_b_RTF[icomp]*Deltac_b[icomp];
                else deltac_b = -0.5*Rho_b[icomp]*Deltac_b[icomp];
                omega_i_b += deltac_b;
            }
            else deltac_b=0.0;

            omega_s_sum[i] += (omega_i*Nel_hit2[i][iunk][inode_box]
                               -omega_i_b*Nel_hit[i][iunk][inode_box])
                                   * Vol_el/((double)Nnodes_per_el_V);
            omega_sum[i]   +=  omega_i*Nel_hit2[i][iunk][inode_box]
                                   * Vol_el/((double)Nnodes_per_el_V);
            ideal_sum[i]  +=  (ideal*Nel_hit2[i][iunk][inode_box]-ideal_b*Nel_hit[i][iunk][inode_box])
                                   * Vol_el/((double)Nnodes_per_el_V);
            lj_sum[i]  +=  (ljterm*Nel_hit2[i][iunk][inode_box]-lj_b*Nel_hit[i][iunk][inode_box])
                                   * Vol_el/((double)Nnodes_per_el_V);
            vext_sum[i]  +=  vext*Nel_hit2[i][iunk][inode_box]
                                   * Vol_el/((double)Nnodes_per_el_V);
            if (Type_poly_TC)
               wtc_sum[i]  +=  wtc_term*Nel_hit2[i][iunk][inode_box]
                                   * Vol_el/((double)Nnodes_per_el_V);
            if (Type_coul > -1){
            psi_rho_sum[i]   +=  psi_rho*Nel_hit2[i][iunk][inode_box]
                                   * Vol_el/((double)Nnodes_per_el_V);
            vext_c_sum[i]  +=  vext_c*Nel_hit2[i][iunk][inode_box]
                                   * Vol_el/((double)Nnodes_per_el_V);
            deltac_sum[i]  +=  (deltacterm*Nel_hit2[i][iunk][inode_box]-deltac_b*Nel_hit[i][iunk][inode_box])
                                   * Vol_el/((double)Nnodes_per_el_V);
            }
        }      /* end of loop over lists */
     }         /* end of icomp loop */
  }            /* end of loc_inode loop */

  for (i=0; i<Imax; i++){
    grand_free_energy[i]   = gsum_double(omega_sum[i]);
    surface_free_energy[i] = gsum_double(omega_s_sum[i]);

    ideal_sum[i] = gsum_double(ideal_sum[i]);
    lj_sum[i] = gsum_double(lj_sum[i]);
    vext_sum[i] = gsum_double(vext_sum[i]);
    if (Type_poly_TC) wtc_sum[i]=gsum_double(wtc_sum[i]);
    if (Type_coul > -1){
      psi_rho_sum[i] = gsum_double(psi_rho_sum[i]);
      vext_c_sum[i] = gsum_double(vext_c_sum[i]);
      deltac_sum[i] = gsum_double(deltac_sum[i]);
    }
  }

  for (i=0; i<Imax; i++) {
     if (Proc==0 && print_flag &&Iwrite!=NO_SCREEN)
        printf("\n----------------------------------------------------------\n");

     grand_free_energy[i] += (sum_phispt);
     surface_free_energy[i] += (sum_phispt - sum_phispt_b);

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

     if (Lper_area && area >0.0) prefac = (fac_vol/fac_area)/area;
     else                        prefac = fac_vol;

     grand_free_energy[i] *= prefac; 
     surface_free_energy[i] *= prefac; 
     hs_energy = (sum_phispt-sum_phispt_b)*prefac; 
     ideal_sum[i] *= prefac; 
     vext_sum[i] *= prefac; 
     lj_sum[i] *= prefac; 
     if (Type_coul>-1){ 
        psi_rho_sum[i] *= prefac; 
        vext_c_sum[i] *= prefac; 
        deltac_sum[i] *= prefac; 
     }	         
     if (Type_poly_TC) wtc_sum[i] *= prefac; 
     
  }
  if (fp !=NULL && Proc==0 && print_flag){
    if (Iwrite!=NO_SCREEN)
     printf ("\n Free energy summary after volume and area factors \n");
    if (Lhard_surf) i=1;
    else            i=0;
    if (Iwrite!=NO_SCREEN){
       printf("\t===== Totals when Integrating over domain volume =====\n");
       printf("\t\t Grand Free Energy: %g \n",grand_free_energy[i]);
       printf("\t\t Surface Free Energy: %g \n",surface_free_energy[i]);
    }
    Energy = surface_free_energy[i]; /* return value for MC-DFT calculations */
    if (Iwrite!=NO_SCREEN){
       printf("\t====Contributions of each term to surface free energy =========\n");
       printf("\t\t Ideal Gas and Chem.pot term: %9.6f\n",ideal_sum[i]); 
       if (Type_func >=0) printf("\t\t Hard Sphere term: %9.6f\n",hs_energy); 
       if (Type_attr >=0) printf("\t\t Lennard-Jones term: %9.6f\n",lj_sum[i]); 
       if (Type_poly_TC) printf("\t\t WTC bond term: %9.6f\n",wtc_sum[i]);
       printf("\t\t External field (non-Coulomb) term: %9.6f\n",vext_sum[i]); 
       if (Type_coul > -1) {
         printf("\t\t Poisson-electrostatics term: %9.6f\n",psi_rho_sum[i]); 
         printf("\t\t External field (Coulomb) term: %9.6f\n",vext_c_sum[i]); 
         printf("\t\t Delta c correction term: %9.6f\n",deltac_sum[i]); 
       }
       if (Nwall==2 &&Lprint_pmf) printf("\t\t Direct wall-wall interaction: %9.6f\n",Uww[0][1]); 
       if (Nwall==2 &&Lprint_pmf) printf("\t\t Total surface free energy: %9.6f\n",
                                                surface_free_energy[i]+Uww[0][1]); 
    }
     if (Nwall > 2 && Lprint_pmf){ 
        for (iwall=0; iwall<Nwall-1;iwall++){
           for (jwall=iwall+1; jwall<Nwall;jwall++){
             Uww_sum += Uww[iwall][jwall];
             if (Iwrite != NO_SCREEN) 
                    printf("iwall=%d jwall=%d Uww=%9.6f\n",iwall,jwall,Uww[iwall][jwall]);
           }
        }
        if (Iwrite != NO_SCREEN) {
           printf("\t\t Direct wall-wall interaction: %9.6f\n",Uww_sum); 
           printf("\t\t Total surface free energy: %9.6f\n", surface_free_energy[i]+Uww_sum); 
        }
     }
     if (Iwrite != NO_SCREEN) {
     printf("\t========================================\n");
     if (Ipot_ff_n != IDEAL_GAS && Lhard_surf){
        printf("\t===== Totals when Integrating over fluid volume =====\n");
        printf("\t\t Grand Free Energy: %g \n",grand_free_energy[0]);
        printf("\t\t Surface Free Energy: %g \n",surface_free_energy[0]);
        printf("\t========================================\n");
     }
     }
     if (first){
        if (Lprint_pmf &&fp!=NULL){
           if (Nwall==2) fprintf(fp," Energy=%9.6f Surf_ex_energy=%9.6f  Uww[0][1]= %9.6f PMF=%9.6f ", 
                                 grand_free_energy[i],surface_free_energy[i],
                                 Uww[0][1],surface_free_energy[i]+Uww[0][1]);
           else if (Nwall>2) fprintf(fp," Eneregy=%9.6f Surf_ex_energy=%9.6f  Total_Uww=%9.6f Total PMF%9.6f  ", 
                                 grand_free_energy[i],surface_free_energy[i],
                                 Uww_sum,surface_free_energy[i]+Uww_sum);
           else fprintf(fp," Energy=%9.6f Surf_ex_energy=%9.6f  ", grand_free_energy[i],surface_free_energy[i]);
        }
        else    fprintf(fp," Energy=%9.6f Surf_ex_energy=%9.6f  ", grand_free_energy[i],surface_free_energy[i]);
        first=FALSE;
     }
     else{
        if (Lprint_pmf &&fp!=NULL){
           if (Nwall==2) fprintf(fp," %9.6f %9.6f   %9.6f %9.6f ", 
                                 grand_free_energy[i],surface_free_energy[i],
                                 Uww[0][1],surface_free_energy[i]+Uww[0][1]);
           else if (Nwall>2) fprintf(fp," %9.6f %9.6f   %9.6f %9.6f  ", 
                                 grand_free_energy[i],surface_free_energy[i],
                                 Uww_sum,surface_free_energy[i]+Uww_sum);
           else fprintf(fp," %9.6f %9.6f  ", grand_free_energy[i],surface_free_energy[i]);
        }
        else    fprintf(fp," %9.6f %9.6f  ", grand_free_energy[i],surface_free_energy[i]);
     }
     if (Iwrite !=NO_SCREEN) printf("----------------------------------------------------------\n");
  }
  return (surface_free_energy[Imax-1]);
}
/*******************************************************************************
assemble_HS_free_energy:  In this subroutine we calculate PHI on all of this   *
                          processor's internals and externals.  Then pass all  *
                          the values to processor 0 where they are collected   *
                          and summed.                                          */

void assemble_HS_free_energy(double **x, double *sum_phispt, double *sum_phispt_b,
                             double *sum_phispt_b_old)
{
  int    i,iel,ielement,iel_box,reflect_flag[NDIM_MAX],ijk[3];
  int    inode,inode_box,loc_inode,idim,nel_hit,nel_hit_b;
  double phi,fac,fac_b,sum_phi,sum_phi_b,sum_phi_b_old;
  double rho_bar[10];

/* integrate phi locally  */

  sum_phi       = 0.0;
  sum_phi_b     = 0.0;
  sum_phi_b_old = 0.0;
  for (idim=0;idim<Ndim;idim++) reflect_flag[idim]=FALSE;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){

     /* get box and global nodes */
      inode_box = L2B_node[loc_inode];
      inode = L2G_node[loc_inode];

      node_to_ijk(inode,ijk);

      /* Load rho-bar's into a single array*/
      for(i=0; i<Phys2Nunk[RHOBAR_ROSEN]; i++)
	      rho_bar[i] = x[i+Phys2Unk_first[RHOBAR_ROSEN]][inode_box];

      phi = phispt_i(rho_bar);

      nel_hit = Nnodes_per_el_V;
      for (iel=0; iel<Nnodes_per_el_V; iel++){
          ielement = node_to_elem(inode,iel,reflect_flag);
          /*if (ielement>=0)*/ iel_box = el_to_el_box(ielement);
          if (ielement != -2){
                if (ielement == -1 ||
                Wall_elems[Nlists_HW-1][iel_box] != -1) nel_hit--;
          }
      }     

      nel_hit_b = Nnodes_per_el_V;
      for (iel=0; iel<Nnodes_per_el_V; iel++){
          ielement = node_to_elem(inode,iel,reflect_flag);
          /*if (ielement>=0)*/ iel_box = el_to_el_box(ielement);

          if (ielement != -2){
                if (ielement == -1 ||
                Wall_elems[0][iel_box] != -1) nel_hit_b--;
          }

      }     

      for (idim=0; idim<Ndim; idim++){
         if ( (Type_bc[idim][0] == REFLECT || 
               Type_bc[idim][0] == IN_BULK || Type_bc[idim][0]==LAST_NODE) &&
                        ijk[idim] == 0) nel_hit /= 2;

         if ( (Type_bc[idim][1] == REFLECT || 
               Type_bc[idim][1] == IN_BULK || Type_bc[idim][1]==LAST_NODE) && 
                       ijk[idim] == Nodes_x[idim]-1) nel_hit /= 2;
      }
      fac   = Vol_el*((double)nel_hit)/((double)Nnodes_per_el_V);
      fac_b = Vol_el*((double)nel_hit_b)/((double)Nnodes_per_el_V);

      if (nel_hit > 0){
         sum_phi       += (phi*fac);
         sum_phi_b     += (phispt_bulk()*fac);
         sum_phi_b_old += (phispt_bulk()*fac_b);
      } 
  }    /* end of loop over local nodes */


  (*sum_phispt)   = gsum_double(sum_phi);
  (*sum_phispt_b) = gsum_double(sum_phi_b);
  (*sum_phispt_b_old)=gsum_double(sum_phi_b_old);

  return;
}
/*******************************************************************************
 energy_elec_vext_vol(x); compute the wall-fluid electrostatic  contributions
                            to the free energy. This is 1/2 times the integral
                            over the fixed charge distribution times the 
                            electrostatic potential. */

double energy_elec_vext_vol(double **x)
{
  int loc_inode, iunk, inode_box, inode, iel_box,ijk_box[3],jln,ielement,reflect_flag[3];
  double charge_at_node,energy;

  reflect_flag[0]=reflect_flag[1]=reflect_flag[2]=0;

 /* The only tricky issue here is that the volumetric charge distribution
    is assigned on an element by element basis while the electrostatic
    potentials are defined on the nodes of the calculation.  So, we
    compute the average charge at each node for this
    integral. */

    energy = 0.0;

    for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
       inode = L2G_node[loc_inode];
       inode_box = L2B_node[loc_inode];
   
      iunk = Phys2Unk_first[POISSON];
    
       charge_at_node = 0.0; 
       for (jln=0; jln< Nnodes_per_el_V; jln++) { 
          ielement = node_to_elem(inode,jln,reflect_flag);
          iel_box = el_to_el_box(ielement); 
          if (iel_box > 0){
             charge_at_node += Charge_vol_els[iel_box]/Nnodes_per_el_V;
          }
       }
      energy += 0.5*charge_at_node*x[iunk][inode_box]*Vol_el;

    }
    printf("PROC=%d: RETURNING energy=%9.6f \n",Proc,energy);
    return (energy*Vol_el);
}

/*******************************************************************************
 * energy_elec: For charged surfaces, find the electrostatic contribution   *
 *             to the force on each wall.                                  */

double energy_elec(double **x, double *sum3)
{
   int iunk,loc_inode,iwall,idim=0,ilist,
     iel_w,inode,surf_norm,ijk[3],inode_box;
   int wall_count[NWALL_MAX],blocked,nblock=0;
   double prefac,nodepos[3],geom_factor[3],deriv_x[3];
   double sum=0.0,wall_avg_psi[NWALL_MAX],
          sblock=0.0,sopen=0.0;
   double deriv_x_A[3],deriv_x_B[3];
   *sum3 = 0.;

   for (iwall=0; iwall<Nwall; iwall++){
        wall_avg_psi[iwall] = 0.;
        wall_count[iwall] = 0;
   }

   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){

       inode_box = L2B_node[loc_inode];
       inode = L2G_node[loc_inode];
       node_to_ijk(inode,ijk);

       for (ilist=0; ilist<Nlists_HW;ilist++) {
          iwall = Nodes_2_boundary_wall[ilist][inode_box];
          if (iwall !=-1 ) break;
       }
       /* iwall = Nodes_2_boundary_wall[Nlists_HW-1][inode_box];
          iwall = Nodes_2_boundary_wall[0][inode_box];*/

       if (iwall != -1){       /* sum up the boundary element contributions */

          node_to_position(inode,nodepos); 
/*          calc_geom_factor(iwall,nodepos,geom_factor);*/

          geom_factor[0]=geom_factor[1]=geom_factor[2]=1.0;

	  iunk = Phys2Unk_first[POISSON];

          for (idim=0; idim<Ndim; idim++)  deriv_x[idim] = 0.0;

            for (iel_w=0; iel_w<Nelems_S[ilist][loc_inode]; iel_w++){
              surf_norm = Surf_normal[ilist][loc_inode][iel_w];
              idim = abs(surf_norm) - 1;

              prefac = (double)(surf_norm/abs(surf_norm))*Area_surf_el[idim]/
                                                           (Nnodes_per_el_S);
 
              if (surf_norm < 0) {
                 deriv_x_A[idim] = calc_deriv_e(idim,inode_box,BFD,&blocked,x,ilist);
                 deriv_x_B[idim] = calc_deriv_e(idim,inode_box,FFD,&blocked,x,ilist);
              }
              else{               
                  deriv_x_A[idim] = calc_deriv_e(idim,inode_box,FFD,&blocked,x,ilist);
                  deriv_x_B[idim] = calc_deriv_e(idim,inode_box,BFD,&blocked,x,ilist);
              }

/*              printf("iel_w: %d  surf_norm: %d  idim: %d\n",iel_w,surf_norm,idim);*/

              sum += x[iunk][inode_box]*(deriv_x_A[idim]-deriv_x_B[idim])*prefac*geom_factor[idim];

              if (blocked){
                 sblock += deriv_x[idim]*prefac;
                 nblock++;
              }
              else
                 sopen += deriv_x[idim]*prefac;

              if (Type_bc_elec[WallType[iwall]]==2){
                       wall_avg_psi[iwall] += x[iunk][inode_box];
                       wall_count[iwall]++;
              }

          } /* end of surface element loop */ 
       }    /* end of test for boundary node */

       /* now calculate the volume term !! */

  /*   for (loc_node_el=0; loc_node_el< Nnodes_per_el_V ; loc_node_el++){

          iel = node_to_elem(inode, loc_node_el, reflect_flag);
          iel_box = el_to_el_box(iel);

*          if (iel >= 0 && Wall_elems[Nlists_HW-1][iel_box] != -1){*
            if (iel >= 0 && Wall_elemsF[0][iel_box] != -1){
*            if (iel >= 0 ){*

              for (idim=0; idim<Ndim; idim++){
                 deriv2_x[idim] = 0.0;
                 switch(Ndim)
                 {
                  case 3:
                    if ( 
                         (idim == 0 && (loc_node_el==0 || loc_node_el == 2
                                   || loc_node_el==4 || loc_node_el==6 )  
                                    && ijk[0] < Nodes_x[0]-1 )              ||
                         (idim == 1 && (loc_node_el==0 || loc_node_el == 1
                                   || loc_node_el==4 || loc_node_el==5 ) 
                                    && ijk[1] < Nodes_x[1]-1 )              ||
                         (idim == 2 && (loc_node_el==0 || loc_node_el == 1
                                   || loc_node_el==2 || loc_node_el==3 ) 
                                    && ijk[2] < Nodes_x[2]-1 )               
                         ){
                         deriv2_x[idim] = calc_deriv2(idim,inode_box,FFD,x);
                         sum2 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                    }
                    break;
                  case 2:
                    if ( 
                         (idim == 0 && (loc_node_el==0 || loc_node_el == 2 ) 
                                    && ijk[0] < Nodes_x[0]-1               ) ||
                         (idim == 1 && (loc_node_el==0 || loc_node_el == 1 ) 
                                    && ijk[1] < Nodes_x[1]-1               )
                         ){
                         deriv2_x[idim] = calc_deriv2(idim,inode_box,FFD,x);
                         sum2 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                    }
                    break;
                  case 1:
                    if ( loc_node_el==0 && ijk[0] < Nodes_x[0]-1 ){
                         deriv2_x[idim] = calc_deriv2(idim,inode_box,FFD,x);
                         sum2 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                    }
                    break;
                 }
              }   * end of idim loop * 
          }       * end of test for a wall element * 
       }          * end of loop over elements that include inode */


    }       /* end of loop over nodes on this processor */
    sum *= -Temp_elec/(8.0*PI);

/*     calculate charging free energy without derivative */
    for (iwall=0; iwall<Nwall; iwall++){
       if (Type_bc_elec[WallType[iwall]]==2) *sum3 = 0.5*wall_avg_psi[iwall]*
                    Area_surf_el[idim]*Elec_param_w[iwall]/ Nnodes_per_el_S;
           /*printf("wall_count %d    Nnodes_per_el_S %d   S_area_tot %6.3f\n",
             wall_count[iwall],Nnodes_per_el_S,S_area_tot[iwall]);*/
 /*  if (Proc==0 && nblock != 0)
      printf("Avg deriv: blocked %9.6f     open %9.6f\n",sblock/nblock,
               sopen/(wall_count[iwall]-nblock));*/
    }
/*if (Proc==0) printf("Surface integral: %9.6f from grad_psi.n    %9.6f from q\n",
          sum,*sum3);*/

/*  if (sum != 0.0 || sum2 != 0.0) 
    printf ("Proc: %d  surface term: %9.6f  volume term: %9.6f\n",
             Proc,sum,sum2*Temp_elec/(8.0*PI));*/

/*  if (Type_bc_elec[WallType[iwall]]==2) sum = *sum3;*/
    return (sum);
}
/*******************************************************************************
 * charge_stress: For charged surfaces, find the contribution of the stress *
 *             tensor to the local pressure.                                */

double charge_stress(double **x,double *sum2)
{
   int loc_inode,iwall,idim,
       inode,ijk[3],inode_box,iel_box;
   int loc_node_el,iel,reflect_flag[3]={FALSE,FALSE,FALSE};
   double deriv2_x[3];
   double sum3=0.0;
   *sum2=0.0;

   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){

       inode_box = L2B_node[loc_inode];
       inode = L2G_node[loc_inode];
       node_to_ijk(inode,ijk);

       iwall = Nodes_2_boundary_wall[Nlists_HW-1][inode_box];

       /* now calculate the volume term !! */

       for (loc_node_el=0; loc_node_el< Nnodes_per_el_V ; loc_node_el++){

          iel = node_to_elem(inode, loc_node_el, reflect_flag);
          iel_box = el_to_el_box(iel);

          if (iel >= 0){

              for (idim=0; idim<Ndim; idim++){
                 deriv2_x[idim] = 0.0;
                 switch(Ndim)
                 {
                  case 3:
                    if ( 
                         (idim == 0 && (loc_node_el==0 || loc_node_el == 2
                                   || loc_node_el==4 || loc_node_el==6 )  
                                    && ijk[0] < Nodes_x[0]-1 )              ||
                         (idim == 1 && (loc_node_el==0 || loc_node_el == 1
                                   || loc_node_el==4 || loc_node_el==5 ) 
                                    && ijk[1] < Nodes_x[1]-1 )              ||
                         (idim == 2 && (loc_node_el==0 || loc_node_el == 1
                                   || loc_node_el==2 || loc_node_el==3 ) 
                                    && ijk[2] < Nodes_x[2]-1 )               
                         ){
                         deriv2_x[idim] = calc_deriv2(idim,inode_box,FFD,x);
                         sum3 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                         if (Wall_elems[Nlists_HW-1][iel_box] != -1)
                           *sum2 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                    }
                    break;
                  case 2:
                    if ( 
                         (idim == 0 && (loc_node_el==0 || loc_node_el == 2 ) 
                                    && ijk[0] < Nodes_x[0]-1               ) ||
                         (idim == 1 && (loc_node_el==0 || loc_node_el == 1 ) 
                                    && ijk[1] < Nodes_x[1]-1               )
                         ){
                         deriv2_x[idim] = calc_deriv2(idim,inode_box,FFD,x);
                         sum3 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                         if (Wall_elems[Nlists_HW-1][iel_box] != -1)
                           *sum2 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                    }
                    break;
                  case 1:
                    if ( loc_node_el==0 && ijk[0] < Nodes_x[0]-1 ){
                         deriv2_x[idim] = calc_deriv2(idim,inode_box,FFD,x);
                         sum3 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                         if (Wall_elems[Nlists_HW-1][iel_box] != -1)
                           *sum2 += (deriv2_x[idim]*deriv2_x[idim]*Vol_el/
                                                 (0.5*Nnodes_per_el_V));
                    }
                    break;
                 }
              }   /* end of idim loop */

          }       /* end of test for a wall element */
       }          /* end of loop over elements that include inode */


    }       /* end of loop over nodes on this processor */

    *sum2 *= Temp_elec/(8.0*PI);
    sum3 *= Temp_elec/(8.0*PI);
/*  if (Proc==0) printf ("integral of (grad_phi)^2 over wall: %9.6f   total volume: %9.6f\n",*sum2,sum3);*/

    return (sum3);
}
/****************************************************************************
 * free_energy_charging_up:  In this routine we calculate the surface       *
 *                           integral that gives the free energy due to     *
 *                           charging up the surfaces                       */
double free_energy_charging_up(double **x)
{
  int loc_inode,idim,inode,iunk, inode_box;
  double psi_i,charging_int,charge_i=0.0;
  int i,iel_w,surf_norm,ilist;
  double psi_deriv[2][2],prefac,dphi_term;

  charging_int = 0.0;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
      inode = L2G_node[loc_inode];
      inode_box = L2B_node[loc_inode];

      if  (Nodes_2_boundary_wall[Nlists_HW-1][node_to_node_box(inode)] != -1){

         iunk = Phys2Unk_first[POISSON];

         ilist = Nlists_HW - 1;
         for (iel_w=0; iel_w<Nelems_S[ilist][loc_inode]; iel_w++){
              
              surf_norm = Surf_normal[ilist][loc_inode][iel_w];
              idim = abs(surf_norm) - 1;

              prefac = (double)(surf_norm/abs(surf_norm))*Area_surf_el[idim];

              if (prefac > 0.0) {  /* RHS*/
                   i = 1;
                   psi_deriv[i][0] = (
                               -    x[iunk][inode_box+2] 
                               + 4.*x[iunk][inode_box+1] 
                               - 3.*x[iunk][inode_box]         )/(2.*Esize_x[0]);
                   psi_deriv[i][1] = (x[iunk][inode_box] -
                                      x[iunk][inode_box-1])/Esize_x[0];
              }
              else {              /* LHS */
                    i = 0;
                    psi_deriv[i][0] = -(
                                  3.*x[iunk][inode_box] 
                                - 4.*x[iunk][inode_box-1]
                                +    x[iunk][inode_box-2])
                                    /(2.0*Esize_x[0]);
                    psi_deriv[i][1] = (x[iunk][inode_box+1] 
                                                      - x[iunk][inode_box])/Esize_x[0];
              }
              charge_i = (Temp_elec/(8.0*PI))*psi_deriv[i][0]; 


          } /* end of surface element loop */ 

/*         charge_i = 0.0;
         for (idim=0; idim<Ndim; idim++){
             charge_i -= Charge_w_sum_els[loc_inode][idim]*Area_surf_el[idim];
         }
 */
 

         psi_i = x[iunk][inode_box];
         charging_int -= (psi_i * charge_i);
/*         printf("loc_inode: %d  charge_i: %9.6f   psi_i: %9.6f  charging_int: %9.6f\n",
                 loc_inode,charge_i,psi_i,charging_int);   */

      }   /* end of test for boundary node */
  }       /* end of loc_inode loop */

/*  printf("PI=%9.6f\n",PI);   */
  dphi_term = (Temp_elec/(8.0*PI))*WallParam[0]
                  *psi_deriv[1][1]*psi_deriv[1][1];
  printf("dphi_term: %9.6f  charging_int: %9.6f \n",dphi_term,charging_int);
  charging_int += dphi_term;

  return charging_int;

  /* collate all the charging_int values from the different processors */
}
/****************************************************************************
int_stencil: Perform the integral sum(j)int rho_j(r')*weight[sten] */
 double int_stencil(double **x,int inode_box,int icomp,int sten_type)
{
  int isten,*offset,inode_sten,ijk_box[3],izone,idim;
  int j,jcomp,junk;
  double weight, sum;
  struct Stencil_Struct *current_sten;
  int **current_sten_offset, reflect_flag[NDIM_MAX];
  double *current_sten_weight;
  for (idim=0; idim<Ndim; idim++) reflect_flag[idim]=FALSE;

/*  izone = Nodes_to_zone_box[inode_box];*/
  izone = 0;

  sum = 0.0;
  node_box_to_ijk_box(inode_box,ijk_box);

  for (junk=Phys2Unk_first[DENSITY];junk<Phys2Unk_last[DENSITY];junk++){
     jcomp=Unk2Comp[junk];

     current_sten = &(Stencil[sten_type][izone][icomp+Ncomp*jcomp]);
     current_sten_offset = current_sten->Offset;
     current_sten_weight = current_sten->Weight;

     for (isten = 0; isten < current_sten->Length; isten++) {
        offset = current_sten_offset[isten];
        weight = current_sten_weight[isten];

         /* Find in the Stencil position on overall mesh */ 
        inode_sten =offset_to_node_box(ijk_box, offset, reflect_flag);

        if (inode_sten >= 0 && !Zero_density_TF[inode_sten][jcomp]) {
           if (Lhard_surf) {
               if (Nodes_2_boundary_wall[Nlists_HW-1][inode_sten]!=-1)
                  weight = HW_boundary_weight
                    (jcomp,Nlists_HW-1,current_sten->HW_Weight[isten], 
                                      inode_sten, reflect_flag);
           }

           if (inode_sten<Nnodes_box && inode_sten >=0){
               sum +=  weight*x[junk][inode_sten];
           }
        }
        else if (inode_sten<0){
             sum += weight*constant_boundary(junk,inode_sten);
        }

     }  /* end of loop over isten */ 
  }     /* end of loop over jcomp */ 
  return (sum);
}
/***************************************************************************
phispt_i: Calculate the hard sphere free energy contribution from 
           scaled particle theory at a given node i.                      */
double phispt_i(double *rho_bar)
{
  int idim;
  double rb0,rb1,rb2,rb3,rb2v[3],rb1v[3];
  double phi_s,phi_v,dot_12,dot_22;

  rb0 = rho_bar[3];
  rb1 = rho_bar[2];
  rb2 = rho_bar[1];
  rb3 = rho_bar[0];
  for (idim=0; idim<Ndim; idim++){
    rb1v[idim] = rho_bar[Nrho_bar_s+Ndim+idim];
    rb2v[idim] = rho_bar[Nrho_bar_s+idim];
  }

  if (rb3 < 1.0 && rb2 > 0.0){

     if (Type_func==0)
        phi_s = -rb0*log(1.0-rb3) + (rb1*rb2)/(1.0-rb3) +
                 (rb2*rb2*rb2)/(24.0*PI*(1.0-rb3)*(1.0-rb3));
     else
        phi_s = -rb0*log(1.0-rb3) + (rb1*rb2)/(1.0-rb3); 
   
     dot_12 = 0.0;
     dot_22 = 0.0;
     for (idim = 0; idim<Ndim; idim++){
         dot_12 += rb1v[idim]*rb2v[idim];
         dot_22 += rb2v[idim]*rb2v[idim];
     }
     if (Type_func==0)
        phi_v = - dot_12/(1.0-rb3) 
                - rb2*dot_22/(8.0*PI*(1.0-rb3)*(1.0-rb3));
     else
        phi_v = - dot_12/(1.0-rb3) 
                + POW_DOUBLE_INT(rb2-dot_22/rb2,3)/(24.0*PI*(1.0-rb3)*(1.0-rb3));
     return(phi_s + phi_v);
  }

  else return(0.0);

}
/***************************************************************************
phispt_bulk: Calculate the hard sphere free energy contribution from 
              scaled particle theory in a bulk fluid.                      */
double phispt_bulk()
{
  int icomp;
  double four_pi_R[NCOMP_MAX],four_pi_RSQ[NCOMP_MAX],pi_six_dCBD[NCOMP_MAX];
  double rb0=0.0,rb1=0.0,rb2=0.0,rb3=0.0,phi_s;
 
  /* 
   * first calculate the rho_bars at all the nodes 
   * rho_bar vectors are zero in a homogeneous bulk fluid 
   */
  for (icomp=0; icomp<Ncomp; icomp++){
     four_pi_R[icomp]   =  2.0 * PI *  Sigma_ff[icomp][icomp];
     four_pi_RSQ[icomp] =  PI * Sigma_ff[icomp][icomp]*Sigma_ff[icomp][icomp];
     pi_six_dCBD[icomp] =  (PI/6.0) * Sigma_ff[icomp][icomp]*Sigma_ff[icomp][icomp]
                                    * Sigma_ff[icomp][icomp];

     if (Lsteady_state){
        rb3 += pi_six_dCBD[icomp]*Rho_b_RTF[icomp];
        rb2 += four_pi_RSQ[icomp]*Rho_b_RTF[icomp];
        rb1 += Sigma_ff[icomp][icomp]*Rho_b_RTF[icomp]/2.0;
        rb0 += Rho_b_RTF[icomp];
     }
     else{
        rb3 += pi_six_dCBD[icomp]*Rho_b[icomp];
        rb2 += four_pi_RSQ[icomp]*Rho_b[icomp];
        rb1 += Sigma_ff[icomp][icomp]*Rho_b[icomp]/2.0;
        rb0 += Rho_b[icomp];
     }
  }

  /* given the rho_bars, calculate dphi_s */  
  phi_s = -rb0*log(1.0-rb3) + (rb1*rb2)/(1.0-rb3) +
           rb2*rb2*rb2/(24.0*PI*(1.0-rb3)*(1.0-rb3));

  return phi_s;
}
/*********************************************************************
 * calc_deriv_e : calculate a derivative of the electric potential!!   */

double calc_deriv_e(int idim,int inode0,int flag,int *blocked, double **x, 
                  int ilist)
{
   int inode1,inode2,offset1[3],offset2[3],
       iwall1,iwall2,
       jdim,ijk_box[3],reflect_flag[3],iunk;
   double deriv=0.0;

   node_box_to_ijk_box(inode0,ijk_box);

   for (jdim=0; jdim<Ndim; jdim++) {
       offset1[jdim] = 0;
       offset2[jdim] = 0;
   }

   switch(flag)
   {
      case CFD:
        offset1[idim] = -1;
        offset2[idim] =  1;
        break;

      case FFD:
        offset1[idim] =  1;
        offset2[idim] =  2;
        break;

      case BFD:
        offset1[idim] = -1;
        offset2[idim] = -2;
        break;
   }

   inode1 = offset_to_node_box(ijk_box,offset1,reflect_flag);
   inode2 = offset_to_node_box(ijk_box,offset2,reflect_flag);
   iunk = Phys2Unk_first[POISSON];

   iwall1 = Nodes_2_boundary_wall[ilist][inode1];
   iwall2 = Nodes_2_boundary_wall[ilist][inode2];

   if (iwall1 == -1 && iwall2 == -1) {
      *blocked = FALSE;
      switch(flag)
      {
         case CFD:
           deriv =    x[iunk][inode2] - x[iunk][inode1];break;
         case FFD:
           deriv = -3*x[iunk][inode0] + 4*x[iunk][inode1] - x[iunk][inode2];break;
         case BFD:
           deriv =  3*x[iunk][inode0] - 4*x[iunk][inode1] + x[iunk][inode2];break;
      }
      deriv /= (2.0*Esize_x[idim]);
   }
   else{
      *blocked = TRUE;
      switch(flag)
      {
         case FFD:
           deriv = x[iunk][inode1] - x[iunk][inode0]; break;
         case BFD:
           deriv = x[iunk][inode0] - x[iunk][inode1]; break;
      }
      deriv /= (Esize_x[idim]);
   }

   return deriv;
}
/***************************************************************************
 * calc_deriv2 : calculate a derivative of the electric potential (first order)!!   */

double calc_deriv2(int idim,int inode0,int flag,double **x)
{
   int inode1,offset[3], iunk,
       jdim,ijk_box[3],reflect_flag[3];
   double deriv=0.0;


   node_box_to_ijk_box(inode0,ijk_box);

   for (jdim=0; jdim<Ndim; jdim++) offset[jdim] = 0;

   switch(flag)
   {
      case FFD:
        offset[idim] =  1;
        break;

      case BFD:
        offset[idim] = -1;
        break;
   }

   inode1 = offset_to_node_box(ijk_box,offset,reflect_flag);
   iunk = Phys2Unk_first[POISSON];

   switch(flag)
   {
      case FFD:
        deriv = x[iunk][inode1] - x[iunk][inode0]; break;
      case BFD:
        deriv = x[iunk][inode0] - x[iunk][inode1];  break;
   }
   deriv /= (Esize_x[idim]);

   return deriv;
}
/***********************************************************************
 * calc_free_energy_polymer: calculates excess grand canonical free energy
                 see eq. 2.22 in Hooper, McCoy, Curro, J. Chem. Phys.,
                           112, p. 3090, (2000)                                 */

double calc_free_energy_polymer(FILE *fp,double **x,double fac_area,double fac_vol)
{
  double /*ads[NCOMP_MAX][2],*/area,free_energy,energy=0.0,boltz,ideal_nz;
  double ifluid,sfluid,ibulk,sbulk,vol;
  int loc_inode, icomp,jcomp, idim,iwall,pol_number;
  int reflect_flag[3],i_boltz, iunk;
  int izone,ijk_box[3],inode_box,i_box;
  double *freen_profile_1D;

  sfluid = 0.0;
  sbulk = 0.0;
  freen_profile_1D = (double *) array_alloc(1, Nnodes_per_proc, sizeof(double));
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
    freen_profile_1D[loc_inode]=0.0;
  }

  for (idim=0; idim<Ndim; idim++) reflect_flag[idim]=FALSE;

  /* loop to do double integral over c(r) */
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
      inode_box = L2B_node[loc_inode];
      for (icomp=0; icomp<Ncomp; icomp++){

	  iunk = Phys2Unk_first[DENSITY]+icomp;
          i_boltz = Phys2Unk_first[CMS_FIELD]+icomp;
          boltz = x[i_boltz][inode_box];
          node_box_to_ijk_box(inode_box, ijk_box);
          izone = 0;

          ideal_nz = calc_u_ideal(icomp,ijk_box,x,&ifluid,&ibulk);

          sfluid += 0.5*ifluid*x[iunk][inode_box]*Nel_hit2[0][iunk][inode_box]*Vol_el/((double)Nnodes_per_el_V);
          sbulk  += 0.5*ibulk*Nel_hit[0][iunk][inode_box]*Vol_el/((double)Nnodes_per_el_V);

         /* Collect the free energy as a function of position for pressure profiles ....
            Note that this funny code does the following...
                (1)find out if they type "icomp" is on a given polymer chain.  
                    If not, move on....if so, 
                (2) add the contribution of the density profile here */
          freen_profile_1D[loc_inode]+=0.5*(x[iunk][inode_box]*ifluid+ibulk);
          pol_number = 0;
          while (Nmer_t[pol_number][icomp]==0) pol_number++;  

          inode_box=L2B_node[loc_inode];
          i_box=loc_find(Phys2Unk_first[DENSITY]+icomp,inode_box,BOX);
          freen_profile_1D[loc_inode] -= (x[iunk][inode_box]-Rho_b[icomp])/Nmer[pol_number];

      }   /* end of icomp loop */
  }       /* end of loc_inode loop */

  /* write liquid-state theory (polymer) free energy */
  energy = sfluid+sbulk;
  free_energy = gsum_double(energy);

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

  if (Lper_area && area>0.0){
       free_energy /= area;
       free_energy *= fac_vol/fac_area;
  }
  else{ free_energy *= fac_vol;}

  /* note that the fac_vol and area adjustments were already 
     made for the excess adsorption terms in dft_out_ads.c */
  for (icomp=0; icomp<Ncomp; icomp++){
     if (icomp >= Ntype_mer){ free_energy -= Ads_ex[icomp][0];
/*        if (Proc==0) printf("adsorption contribution of icomp=%d is %9.6f\n",icomp,Ads_ex[icomp][0]);*/
     }
     else{
        pol_number = 0;
        while (Nmer_t[pol_number][icomp]==0) pol_number++;
        free_energy -= Ads_ex[icomp][0]/Nmer[pol_number];
/*        if (Proc==0) printf("adsorption contribution of icomp=%d is %9.6f  and Nmer is %d\n",icomp,
                     Ads_ex[icomp][0]/Nmer[pol_number],Nmer[pol_number]);*/
     }
  }

  print_freen_profile_1D(freen_profile_1D,"dft_freen_prof.dat");

  if (Proc == 0) {
     if (Nwall == 0 || Ndim == 1) {
        vol = Size_x[0];
        for (idim=1; idim < Ndim; idim++) vol *= Size_x[idim];
          printf("Delta(Free energy)/Volume (where Vol=%lf) = %lf\n",vol,free_energy/vol);
          printf("Delta(Free energy) = %lf\n",free_energy);
     }
     else
        printf("Delta(Free energy) = %lf \n",free_energy);
     if (fp !=NULL) fprintf(fp," %9.6f  %9.6f", free_energy, free_energy/vol);
  }
  return(free_energy);
}
/****************************************************************************/
/* calc_u_ideal:
             calculate the ideal field. Necessary only for nodes where the
             true field (walls) has already been incorporated into the
             boltz unknown                      */
double calc_u_ideal(int itype_mer, int *ijk_box, double **x, double *fluid, double *bulk)
{
  int   **sten_offset, *offset, isten;
  double *sten_weight,  weight, weight_bulk;
  struct Stencil_Struct *sten;

  double sign,ideal_nz;
  int jtype_mer, jlist, iunk;
  int reflect_flag[NDIM_MAX];
  int j_box, jnode_box;

  sign = 1.0;
  ideal_nz = 0.;
  *fluid = 0.;
  *bulk = 0.;
  for (jtype_mer=0; jtype_mer<Ncomp; jtype_mer++){
      if (Nlists_HW <= 2) jlist = 0;
      else                jlist = jtype_mer;

      sten = &(Stencil[POLYMER_CR][0][itype_mer+Ncomp*jtype_mer]);
      sten_offset = sten->Offset;
      sten_weight = sten->Weight;

      for (isten = 0; isten < sten->Length; isten++) {
        offset = sten_offset[isten];
        weight_bulk = sten_weight[isten];

         /* Find the Stencil point */
         jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);

         if (jnode_box >= 0 ) {
            if (Lhard_surf) {
                if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1)
                   weight = HW_boundary_weight
                    (jtype_mer,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
                else weight = weight_bulk;
            }
            else weight = weight_bulk;

	    iunk = Phys2Unk_first[DENSITY]+jtype_mer;

            *fluid +=  sign*weight*x[iunk][jnode_box];
            *bulk -= sign*weight_bulk*Rho_b[jtype_mer]*Rho_b[itype_mer];

            if (!Zero_density_TF[jnode_box][jtype_mer])  {
                  ideal_nz -=  sign*(weight*x[iunk][jnode_box]-
                                           weight*Rho_b[jtype_mer]);
            }
         }
         else if ( jnode_box == -2 ){     /*in wall*/
            *bulk -= sign*weight_bulk*Rho_b[jtype_mer]*Rho_b[itype_mer];
         /* *ideal_t +=  sign*weight_bulk*Rho_b[jtype_mer];*/
             /*  in wall does not contribute to ideal_nz     */
         }
         /* jnode_box == -1 = in bulk contributes nothing  */
      }
  }
  if ((Sten_Type[POLYMER_CR] == 2) && (ideal_nz < 0.5))
                                       ideal_nz = sqrt(1.-2.*ideal_nz);
  return ideal_nz;
}
/*****************************************************************************/

