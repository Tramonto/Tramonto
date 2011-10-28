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

/* dft_plugin_user_continue.c:  This routine is a place holder where a user
   can implement their own problem specific continuation routines without 
   accessing archived files in Tramonto.  */

/* note that when implementing a new continuation type, a new case should
   be added to each of the routines below.   The new case must also be defined
   in dft_globals_const.h.  Refer to the routines in dft_switch_continuation.c
   and dft_plugins_archived_continue.c to see how variables can be implemented. */

/* here is a definition of parameters passed to these routines...
    int cont_type : This is the ID number associated with the new continuation type.  
                    User defined cases should start at number 200.

    int Loca_contID : This ID number indentifies which of the Loca functions a particular
                       call refers to.  Loca_contID=0 refers to Loca.contype1 (arc length),
                       Loca_contID=1 refers to Loca.conttype2 (binodal / spinodal).  

    double param :  This is the continuous variable that Loca will be adjusting during the run.

    FILE *fp : This file pointer refers to dft_output.dat and the output that you want in that file.

*/
                      
#include <stdio.h>
#include "dft_plugin_user_continue.h"
/*****************************************************************************/
double get_init_param_user_plugin(int cont_type,int Loca_contID)
{
  switch(cont_type){
   /*   user_param_type:*/ /*  ..... set the string for your continuation type - definine it in dft_globals_const.h */
         /* enter source code here to set the initial value of the parameter that will vary in your calculation - 
            let's assume it is a density, as defined by Loca_contID
           
            return Rho_b[Cont_ID[Loca_contID][0]];
            break;
          */

      default:
        printf("ERROR: Unknown Continuation parameter in user plugin %d\n",cont_type);
        exit(-1); 
        break;
  }
  return 0.0;
}
/*****************************************************************************/
void assign_param_user_plugin(int cont_type, int Loca_contID, double param,char *output_file1)
{
 double param_old;
 /* note that depending on which variables you are adjusting you may need to recalculate stencils, external fields, or 
    thermodynamics.  See code in dft_switch_continuation.c and dft_scale_variables.c for help. */

  switch(cont_type){

      /*user_param_type:*/ /* set the string for your new continuation type here - define it in dft_globals_const.h */

        /* enter source code here to change the value of the chose parameters.  
           As an example, assume that although you are varying Rho_b[0], you have another dependent density 
           that must also vary in the same ratio as the continuation density - 
           the user code would look something like this...
 
         param_old=Rho_b[Cont_ID[Loca_contID][0]];
         ratio=param/param_old
         Rho_b[Cont_ID[Loca_contID][0]] *= ratio;
         Rho_b[dependent_density_ID] *= ratio;
         break;
        */

      default:
        printf("ERROR_apt: Unknown Continuation parameter in user plugin %d\n",cont_type);
        exit(-1); 
        break;
  }

  /* on changing some parameters, it may be necessary to recompute stencils, thermodynamics, etc. 
     It is necessary to add the new continuation type to the routine adjust_dep_params() in 
     dft_switch_continuation.c.  That routine sets a series of logical to turn on various types
     of adjustments that may be needed for a given type of continuation */ 

  adjust_dep_params(cont_type,Loca_contID,param_old,param,output_file1);  
  return;
}
/*****************************************************************************/
void print_cont_type_user_plugin(int cont_type, FILE *fp,int Loca_contID)
{
  switch(cont_type){
      /*user_param_type:*/
          /* enter souce code to define labeling of your continuation parameter(s) in dft_output.dat 
          break;*/

      default:
        printf("ERROR_apt: Unknown Continuation parameter in user plugin %d\n",cont_type);
        exit(-1); 
        break;
  }
  return;
}
/*****************************************************************************/
void print_cont_variable_user_plugin(int cont_type,FILE *fp,int Loca_contID)
{
  switch(cont_type){
      /*user_param_type:*/
          /* enter souce code to print continuation parameter(s) in dft_output.dat 
          break;*/

      default:
        printf("ERROR_apt: Unknown Continuation parameter in user plugin %d\n",cont_type);
        exit(-1); 
        break;
  }
  return;
}
/*****************************************************************************/

