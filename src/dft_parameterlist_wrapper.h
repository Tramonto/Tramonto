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

#ifndef DFT_PARAMETERLIST_WRAPPER_H
#define DFT_PARAMETERLIST_WRAPPER_H
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

  /*****************************************************/
  /**                  dft_Parameterlist             **/
  /***************************************************/

  void * dft_parameterlist_create();

  void dft_parameterlist_destruct(void * parameterlistptr);

  void dft_parameterlist_set_int(void * parameterlistptr, char * name, int value);

  void dft_parameterlist_set_double(void * parameterlistptr, char * name, double value);

  void dft_parameterlist_set_int_array(void * parameterlistptr, char * name, int * values);

  void dft_parameterlist_set_double_array(void * parameterlistptr, char * name, double * values);

  void dft_parameterlist_set_all_aztec(void * parameterlistptr, int * options, double * params);

#ifdef __cplusplus
}
#endif

#endif /* DFT_PARAMETERLIST_WRAPPER_H */
