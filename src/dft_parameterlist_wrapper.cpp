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


#include "dft_parameterlist_wrapper.h"
#include "Teuchos_ParameterList.hpp"
#include "az_aztec_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

  /*****************************************************/
  /**                  dft_ParameterList             **/
  /***************************************************/

  void * dft_parameterlist_create() {
    Teuchos::ParameterList * parameterList_ = new Teuchos::ParameterList;
    return parameterList_;
  }

  void dft_parameterlist_destruct(void * parameterlistptr) {
    Teuchos::ParameterList * parameterList_ = (Teuchos::ParameterList *) parameterlistptr;
    delete parameterList_;
  }

  void dft_parameterlist_set_int(void * parameterlistptr, char * name, int value) {
    Teuchos::ParameterList * parameterList_ = (Teuchos::ParameterList *) parameterlistptr;
    parameterList_->set(name, value);
  }

  void dft_parameterlist_set_double(void * parameterlistptr, char * name, double value) {
    Teuchos::ParameterList * parameterList_ = (Teuchos::ParameterList *) parameterlistptr;
    parameterList_->set(name, value);
  }

  void dft_parameterlist_set_int_array(void * parameterlistptr, char * name, int * values) {
    Teuchos::ParameterList * parameterList_ = (Teuchos::ParameterList *) parameterlistptr;
    parameterList_->set(name, values);
  }

  void dft_parameterlist_set_double_array(void * parameterlistptr, char * name, double * values) {
    Teuchos::ParameterList * parameterList_ = (Teuchos::ParameterList *) parameterlistptr;
    parameterList_->set(name, values);
  }

  void dft_parameterlist_set_all_aztec(void * parameterlistptr, int * options, double * params) {
    Teuchos::ParameterList * parameterList_ = (Teuchos::ParameterList *) parameterlistptr;
    parameterList_->set("Solver", options[AZ_solver]);
    parameterList_->set("Scaling", options[AZ_scaling]);
    parameterList_->set("Precond", options[AZ_precond]);
    parameterList_->set("Conv", options[AZ_conv]);
    parameterList_->set("Output", options[AZ_output]);
    parameterList_->set("Pre_calc", options[AZ_pre_calc]);
    parameterList_->set("Max_iter", options[AZ_max_iter]);
    parameterList_->set("Poly_ord", options[AZ_poly_ord]);
    parameterList_->set("Overlap", options[AZ_overlap]);
    parameterList_->set("Type_overlap", options[AZ_type_overlap]);
    parameterList_->set("Kspace", options[AZ_kspace]);
    parameterList_->set("Orthog", options[AZ_orthog]);
    parameterList_->set("Aux_vec", options[AZ_aux_vec]);
    parameterList_->set("Reorder", options[AZ_reorder]);
    parameterList_->set("Keep_info", options[AZ_keep_info]);
    parameterList_->set("Subdomain_solve", options[AZ_subdomain_solve]);
    parameterList_->set("Graph_fill", options[AZ_graph_fill]);
    parameterList_->set("Init_guess", options[AZ_init_guess]);
    parameterList_->set("Keep_kvecs", options[AZ_keep_kvecs]);
    parameterList_->set("Apply_kvecs", options[AZ_apply_kvecs]);
    parameterList_->set("Orth_kvecs", options[AZ_orth_kvecs]);
    parameterList_->set("Ignore_scaling", options[AZ_ignore_scaling]);
    parameterList_->set("Check_update_size", options[AZ_check_update_size]);
    parameterList_->set("Extreme", options[AZ_extreme]);
    parameterList_->set("Diagnostics", options[AZ_diagnostics]);
    parameterList_->set("Tol", params[AZ_tol]);
    parameterList_->set("Drop", params[AZ_drop]);
    parameterList_->set("Ilut_fill", params[AZ_ilut_fill]);
    parameterList_->set("Omega", params[AZ_omega]);
    parameterList_->set("Rthresh", params[AZ_rthresh]);
    parameterList_->set("Athresh", params[AZ_athresh]);
    parameterList_->set("Update_reduction", params[AZ_update_reduction]);
    parameterList_->set("Ill_cond_thresh", params[AZ_ill_cond_thresh]); 
  }

#ifdef __cplusplus
}
#endif
