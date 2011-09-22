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


#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_TpetraBasicLinProbMgr.hpp"

#ifdef __cplusplus
extern "C" {
#endif

typedef Teuchos::Comm<int> COMM;

typedef dft_BasicLinProbMgr<double,int,int> BLPM;

  /*****************************************************/
  /**                  dft_BasicLinProbMgr            **/
  /*****************************************************/

  void * dft_basic_lin_prob_mgr_create(int numUnks, void * Parameterlist_list, MPI_Comm comm) {
    RCP<ParameterList> my_list = rcp( (ParameterList *) Parameterlist_list, false);
    const RCP<const COMM> my_comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
    BLPM * linprobmgr_ = new BLPM(numUnks, my_list, my_comm);
    return( (void *)linprobmgr_ );
  }

  void dft_linprobmgr_destruct(void * linprobmgr) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    delete linprobmgr_;
  }

  int dft_linprobmgr_setnodalrowmap(void * linprobmgr, int numgids, int * gids) {
    ArrayView<const int> gid_arr(gids, numgids);
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->setNodalRowMap(gid_arr);
    return( 0 );
  }

  int dft_linprobmgr_setnodalcolmap(void * linprobmgr, int numgids, int * gids) {
    ArrayView<const int> gid_arr(gids, numgids);
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->setNodalColMap(gid_arr);
    return( 0 );
  }

  int dft_linprobmgr_setcoarsenednodeslist(void * linprobmgr, int numgids, int * gids) {
    if(numgids == 0){
      return( 0 );
    }
    ArrayView<const int> gid_arr(gids, numgids);
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->setCoarsenedNodesList(gid_arr);
    return( 0 );
  }

  int dft_linprobmgr_finalizeblockstructure(void * linprobmgr) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->finalizeBlockStructure();
    return( 0 );
  }

  int dft_linprobmgr_initializeproblemvalues(void * linprobmgr) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->initializeProblemValues();
    return( 0 );
  }

  int dft_linprobmgr_insertrhsvalue (void * linprobmgr, int iunk, int inode, double value) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->insertRhsValue(iunk, inode, value);
    return( 0 );
  }

  int dft_linprobmgr_insertonematrixvalue (void * linprobmgr, int iunk, int ownednode,
					   int junk, int boxnode, double value) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->insertMatrixValue(iunk, ownednode, junk, boxnode, value);
    return( 0 );
  }

  int dft_linprobmgr_insertmultinodematrixvalues (void * linprobmgr, int iunk, int ownednode,
						  int junk, int * boxnodeindices, double *values, int numEntries) {
    ArrayView<const int> indices(boxnodeindices, numEntries);
    ArrayView<const double> vals(values, numEntries);

    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->insertMatrixValues(iunk, ownednode, junk, indices, vals);
    return( 0 );
  }

  int dft_linprobmgr_insertmultiphysicsmatrixvalues (void * linprobmgr, int iunk, int ownednode,
						    int *junkindices, int boxnode, double *values, int numEntries) {
    ArrayView<const int> indices(junkindices, numEntries);

    ArrayView<const double> vals(values, numEntries);
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->insertMatrixValues(iunk, ownednode, indices, boxnode, vals);
    return( 0 );
  }

  int dft_linprobmgr_finalizeproblemvalues(void * linprobmgr) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->finalizeProblemValues();
    return( 0 );
  }

  double dft_linprobmgr_getmatrixvalue (void * linprobmgr, int iunk, int ownednode,
					int junk, int boxnode) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    return( linprobmgr_->getMatrixValue(iunk, ownednode, junk, boxnode) );
  }

  int dft_linprobmgr_setrhs(void * linprobmgr, double** x) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    int numOwnedNodes = linprobmgr_->getNumOwnedNodes();
    int numUnknownsPerNode = linprobmgr_->getNumUnknownsPerNode();
    Array<ArrayView<const double> > my_x(numUnknownsPerNode);
    for(int i = 0; i < numUnknownsPerNode; i++){
      my_x[i] = ArrayView<const double>(x[i], numOwnedNodes);
    }
    ArrayView<ArrayView<const double> > ret_val(my_x());
    linprobmgr_->setRhs(ret_val);
    return( 0 );
  }

  int dft_linprobmgr_getlhs(void * linprobmgr, double** x) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    int numOwnedNodes = linprobmgr_->getNumOwnedNodes();
    int numUnknownsPerNode = linprobmgr_->getNumUnknownsPerNode();

    ArrayRCP<ArrayRCP<double> > data = linprobmgr_->getLhs();
    for(int j = 0; j < data.size(); j++){
      for(int k = 0; k < data[j].size(); k++){
	x[j][k] = data[j][k];
      }
    }
    return( 0 );
  }

  int dft_linprobmgr_getrhs(void * linprobmgr, double** x) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    int numOwnedNodes = linprobmgr_->getNumOwnedNodes();
    int numUnknownsPerNode = linprobmgr_->getNumUnknownsPerNode();

    ArrayRCP<ArrayRCP<double> > data = linprobmgr_->getRhs();
    for(int j = 0; j < data.size(); j++){
      for(int k = 0; k < data[j].size(); k++){
	x[j][k] = data[j][k];
      }
    }
    return( 0 );
  }

  int dft_linprobmgr_writeMatrix(void * linprobmgr, char * filename, char * matrixName, char * matrixDescription) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->writeMatrix(filename, matrixName, matrixDescription);
    return( 0 );
  }

  int dft_linprobmgr_setupsolver(void * linprobmgr) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->setupSolver();
    return( 0 );
  }

  int dft_linprobmgr_solve(void * linprobmgr) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->solve();
    return( 0 );
  }

  int dft_linprobmgr_applymatrix(void * linprobmgr, double**x, double** b) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    int numOwnedNodes = linprobmgr_->getNumOwnedNodes();
    int numUnknownsPerNode = linprobmgr_->getNumUnknownsPerNode();

    Array<ArrayView<const double> > x_views(numUnknownsPerNode);

    for(int i = 0; i < numUnknownsPerNode; i++){
      x_views[i] = Teuchos::arrayView(x[i], numOwnedNodes);
    }

    ArrayView<ArrayView<const double> > my_view(x_views);
    ArrayRCP<ArrayRCP<double> > b_data = linprobmgr_->applyMatrix(my_view);
    for(int i = 0; i < b_data.size(); i++){
      for(int j = 0; j < b_data[i].size(); j++){
	b[i][j] = b_data[i][j];
      }
    }
    return( 0 );
  }

  int dft_linprobmgr_importr2c(void * linprobmgr, double** x, double **b) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    int numOwnedNodes = linprobmgr_->getNumOwnedNodes();
    int numUnknownsPerNode = linprobmgr_->getNumUnknownsPerNode();

    ArrayRCP<ArrayRCP<const double> > my_x = Teuchos::arcp<ArrayRCP<const double> >(numUnknownsPerNode);
    for(int i = 0; i < numUnknownsPerNode; i++){
      my_x[i] = Teuchos::arcp<double>(x[i], 0, numOwnedNodes, false);
    }

    ArrayRCP<ArrayRCP<double> > my_b = linprobmgr_->importR2C(my_x);
    for(int i = 0; i < my_b.size(); i++){
      for(int j = 0; j < my_b[i].size(); j++){
	b[i][j] = my_b[i][j];
      }
    }
    return( 0 );
  }

  int dft_linprobmgr_importsingleunknownr2c(void * linprobmgr, double* x, double *b) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    int numOwnedNodes = linprobmgr_->getNumOwnedNodes();
    int numUnknownsPerNode = linprobmgr_->getNumUnknownsPerNode();

    ArrayRCP<double> my_x = Teuchos::arcp<double> (numOwnedNodes);

    my_x = Teuchos::arcp<double>(x, 0, numOwnedNodes, false);

    ArrayRCP<double> my_b = linprobmgr_->importR2C(my_x);

    for(int j = 0; j < my_b.size(); j++){
      b[j] = my_b[j];
    }
    return( 0 );
  }

  int dft_linprobmgr_importnodalr2c(void * linprobmgr, double* x, double *b) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    int numOwnedNodes = linprobmgr_->getNumOwnedNodes();
    int numUnknownsPerNode = linprobmgr_->getNumUnknownsPerNode();

    ArrayRCP<double> ret_x = Teuchos::arcp<double>(x, 0, numOwnedNodes, false);

    ArrayRCP<double> ret_val = linprobmgr_->importR2C(ret_x);

    for(int i = 0; i < ret_val.size(); i++){
      b[i] = ret_val[i];
    }
    return( 0 );
  }

#ifdef __cplusplus
}
#endif
