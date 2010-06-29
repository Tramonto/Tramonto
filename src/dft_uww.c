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
 *  FILE: dft_uww.c
 *
 *  This file contains routines to set up wall-wall interactions for
 *  cases where we want to include these energies in the total.  Note
 *  that these routines are in need of an overhaul.
 */

#include "dft_uww.h"

#define FLAG_LJ 0
#define FLAG_COULOMB 1

/******************************************************************************/

void setup_wall_wall_potentials() 

/* Set up any wall-wall potentials we may need for PMF calculations.
 *
 *    Authors:  Laura Frink
 */

{
   /* Local variable declarations */
   
   char *yo = "wall_wall_potentials";
   int iwall,jwall,itype_w,jtype_w;
   FILE *fp10, *fp11;
  
  /********************** BEGIN EXECUTION ************************************/

  if (Proc==0) printf("\n %s: Setting up Wall-Wall Potentials ... \n",yo);
  /*
   * Allocate and zero the arrays we will calculate here
   */

  if( (fp10  = fopen("dft_uww.dat","w")) == NULL) {
    printf("Can't open file dft_uww.dat\n");
    exit(1);
  }
  if( (fp11  = fopen("dft_uww_link.dat","w")) == NULL) {
    printf("Can't open file dft_uww.dat\n");
    exit(1);
  }

   Uww = (double **) array_alloc (2,Nwall,Nwall, sizeof(double));
   Uww_link = (double **) array_alloc (2,Nlink,Nlink, sizeof(double));

   for (iwall=0; iwall<Nwall; iwall++){
      for (jwall=0; jwall<Nwall; jwall++) {
          Uww[iwall][jwall]=0.0;
          Uww_link[Link[iwall]][Link[jwall]]=0.0;
       }
   }

  for (iwall=0; iwall<Nwall-1; iwall++){
     itype_w=WallType[iwall];
     for (jwall=iwall; jwall<Nwall; jwall++){
        jtype_w=WallType[jwall];

        /* compute wall-wall interactions for atoms */
        if (Ipot_ww_n[itype_w][jtype_w]==ATOM_CENTERS_WW){
            if (Type_uwwPot != PAIR_HARD) setup_atomic_ww(iwall,jwall,Type_uwwPot);
            if (Type_coul >= 0){ 
                 setup_atomic_ww(iwall,jwall,PAIR_COULOMB);
            }
         }

  } /* end loop over jwall */
  } /* end loop over iwall */

  for (iwall=0; iwall<Nwall-1; iwall++){
     for (jwall=iwall; jwall<Nwall; jwall++){
         if (Link[iwall] != Link[jwall])
         Uww_link[Link[iwall]][Link[jwall]] += Uww[iwall][jwall];
         Uww_link[Link[jwall]][Link[iwall]] += Uww[iwall][jwall];
         fprintf(fp10," %d  %d  %9.6f\n",iwall,jwall,Uww[iwall][jwall]);
     }
  }

  if (Nwall != Nlink){
     for (iwall=0; iwall<Nlink-1; iwall){
         for (jwall=iwall; jwall<Nlink; jwall){
            fprintf(fp11," %d  %d  %9.6f\n",iwall,jwall,Uww_link[iwall][jwall]);
         }
      }
  }

  fclose (fp10);
  fclose (fp11);
  return;
}
/******************************************************************************/
/* setup_atomic_ww: In this routine we base wall-wall potentials on the
                            centers of the atomic particles that make up the
                            surfaces.  */
void setup_atomic_ww(int iwall, int jwall,int type_uwwpot){


  int idim;
  double xi,xj,rsq,r;
  double param1,param2,param3,param4;

  pairPotparams_switch(type_uwwpot,WALL_WALL,iwall,jwall,&param1,&param2,&param3,&param4);

  rsq = 0.0;
  for (idim=0; idim<Ndim; idim++){
     xi = WallPos[idim][iwall];
     xj = WallPos[idim][jwall];
     rsq += (xi-xj)*(xi-xj);
  }
  r = sqrt (rsq);

  Uww[iwall][jwall] += pairPot_switch(r,param1,param2,param3,param4,Type_uwwPot);
  return;
}
/******************************************************************************/

