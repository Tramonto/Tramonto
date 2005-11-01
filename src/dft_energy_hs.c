#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

/****************************************************************************/
double integrand_hs_freen(int iunk, int inode_box,double **x)
{
  int i;
  double integrand,rho_bar[10];

  for(i=0; i<Phys2Nunk[RHOBAR_ROSEN]; i++)
      rho_bar[i] = x[i+Phys2Unk_first[RHOBAR_ROSEN]][inode_box];

  integrand = phispt(rho_bar);
  return(integrand);
}
/*******************************************************************************/
double integrand_hs_freen_bulk(int iunk, int inode_box,double **x)
{
  double integrand;

  if (Lsteady_state) integrand = phispt(Rhobar_b_RTF);
  else               integrand = phispt(Rhobar_b);

  return(integrand);
}
/****************************************************************************/
/*phispt: Calculate the hard sphere free energy contribution from
           scaled particle theory at a given node i.                      */
double phispt(double *rho_bar)
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
     if (Type_func==ROSENFELD)
        phi_v = - dot_12/(1.0-rb3)
                - rb2*dot_22/(8.0*PI*(1.0-rb3)*(1.0-rb3));
     else if (Type_func==ROSENFELD2) 
        phi_v = - dot_12/(1.0-rb3)
                + POW_DOUBLE_INT(rb2-dot_22/rb2,3)/(24.0*PI*(1.0-rb3)*(1.0-rb3));
     else{
       printf("problem with type of HS / ROSENFELD FUNCTIONAL");
       exit(-1);
     }
     return(phi_s + phi_v);
  }
  else return(0.0);
}
/****************************************************************************/

