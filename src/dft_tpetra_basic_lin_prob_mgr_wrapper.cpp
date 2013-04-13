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

#if LINSOLVE_PREC == 0
// Solver precision
#define SCALAR                    float
// Matrix precision
#define MAT_SCALAR                SCALAR

#elif LINSOLVE_PREC == 1
// Solver precision
#define SCALAR                    double
// Matrix precision
#if MIXED_PREC == 1
#define MAT_SCALAR                float
#else
#define MAT_SCALAR                SCALAR
#endif

#elif LINSOLVE_PREC == 2
#include <qd/dd_real.h>
// Solver precision
#define SCALAR                    dd_real
// Matrix precision
#if MIXED_PREC == 1
#define MAT_SCALAR                double
#else
#define MAT_SCALAR                SCALAR
#endif

#elif LINSOLVE_PREC == 3
#include <qd/qd_real.h>
// Solver precision
#define SCALAR                    qd_real
// Matrix precision
#if MIXED_PREC == 1
#define MAT_SCALAR                dd_real
#else
#define MAT_SCALAR                SCALAR
#endif

#endif

typedef dft_BasicLinProbMgr<SCALAR,MAT_SCALAR,int,int> BLPM;

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
    ArrayView<const int> gid_arr((numgids == 0)? NULL : gids, numgids);
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->setNodalRowMap(gid_arr);
    return( 0 );
  }

  int dft_linprobmgr_setnodalcolmap(void * linprobmgr, int numgids, int * gids) {
    ArrayView<const int> gid_arr((numgids == 0)? NULL : gids, numgids);
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->setNodalColMap(gid_arr);
    return( 0 );
  }

  int dft_linprobmgr_setcoarsenednodeslist(void * linprobmgr, int numgids, int * gids) {
    ArrayView<const int> gid_arr((numgids == 0)? NULL : gids, numgids);
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
    SCALAR val = Teuchos::as<SCALAR>(value);
    linprobmgr_->insertRhsValue(iunk, inode, val);
    return( 0 );
  }

  int dft_linprobmgr_insertonematrixvalue (void * linprobmgr, int iunk, int ownednode,
					   int junk, int boxnode, double value) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    MAT_SCALAR val = Teuchos::as<MAT_SCALAR>(value);
    linprobmgr_->insertMatrixValue(iunk, ownednode, junk, boxnode, val);
    return( 0 );
  }

  int dft_linprobmgr_insertmultinodematrixvalues (void * linprobmgr, int iunk, int ownednode,
						  int junk, int * boxnodeindices, double *values, int numEntries) {
    ArrayView<const int> indices(boxnodeindices, numEntries);
    MAT_SCALAR *fvalues;
    fvalues = new MAT_SCALAR[numEntries];
    for (int i=0;i<numEntries;i++)
      fvalues[i] = Teuchos::as<MAT_SCALAR>(values[i]);
    ArrayView<const MAT_SCALAR> vals(fvalues, numEntries);
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->insertMatrixValues(iunk, ownednode, junk, indices, vals);
    delete [] fvalues;
    return( 0 );
  }

  int dft_linprobmgr_insertmultiphysicsmatrixvalues (void * linprobmgr, int iunk, int ownednode,
						    int *junkindices, int boxnode, double *values, int numEntries) {
    ArrayView<const int> indices(junkindices, numEntries);
    MAT_SCALAR *fvalues;
    fvalues = new MAT_SCALAR[numEntries];
    for (int i=0;i<numEntries;i++)
      fvalues[i] = Teuchos::as<MAT_SCALAR>(values[i]);
    ArrayView<const MAT_SCALAR> vals(fvalues, numEntries);
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    linprobmgr_->insertMatrixValues(iunk, ownednode, indices, boxnode, vals);
    delete [] fvalues;
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
    MAT_SCALAR retval1 = linprobmgr_->getMatrixValue(iunk, ownednode, junk, boxnode);
    double retval2 = Teuchos::as<double>(retval1);
    return retval2;
  }

  int dft_linprobmgr_setrhs(void * linprobmgr, double** x) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    int numOwnedNodes = linprobmgr_->getNumOwnedNodes();
    int numUnknownsPerNode = linprobmgr_->getNumUnknownsPerNode();

    SCALAR **fx;
    fx = new SCALAR*[numUnknownsPerNode];
    for (int i=0; i<numUnknownsPerNode;i++)
      fx[i] = new SCALAR[numOwnedNodes];
    for(int i=0;i<numUnknownsPerNode;i++)
      for(int j=0;j<numOwnedNodes;j++)
	fx[i][j] = Teuchos::as<SCALAR>(x[i][j]);

    Array<ArrayView<const SCALAR> > my_x(numUnknownsPerNode);
    for(int i = 0; i < numUnknownsPerNode; i++){
      my_x[i] = ArrayView<const SCALAR>(fx[i], numOwnedNodes);
    }
    ArrayView<ArrayView<const SCALAR> > ret_val(my_x());
    linprobmgr_->setRhs(ret_val);
    for (int i = 0; i < numUnknownsPerNode; ++i)
      delete [] fx[i];
    delete [] fx;

    return( 0 );
  }

  int dft_linprobmgr_getlhs(void * linprobmgr, double** x) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    int numUnknownsPerNode = linprobmgr_->getNumUnknownsPerNode();
    ArrayRCP<ArrayRCP<SCALAR> > data = linprobmgr_->getLhs();
    for(int j = 0; j < numUnknownsPerNode; j++){
      for(int k = 0; k < data[j].size(); k++){
	x[j][k] = Teuchos::as<double>(data[j][k]);
      }
    }
    return( 0 );
  }

  int dft_linprobmgr_getrhs(void * linprobmgr, double** x) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    int numUnknownsPerNode = linprobmgr_->getNumUnknownsPerNode();
    ArrayRCP<ArrayRCP<SCALAR> > data = linprobmgr_->getRhs();
    for(int j = 0; j < data.size(); j++){
      for(int k = 0; k < data[j].size(); k++){
	x[j][k] = Teuchos::as<double>(data[j][k]);
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
    int numBoxNodes = linprobmgr_->getNumBoxNodes();
    int numUnknownsPerNode = linprobmgr_->getNumUnknownsPerNode();

    SCALAR **fx;
    fx = new SCALAR*[numUnknownsPerNode];
    for (int i=0; i<numUnknownsPerNode;i++)
      fx[i] = new SCALAR[numBoxNodes];
    for(int i=0;i<numUnknownsPerNode;i++)
      for(int j=0;j<numBoxNodes;j++)
	fx[i][j] = Teuchos::as<SCALAR>(x[i][j]);

    Array<ArrayView<const SCALAR> > x_views(numUnknownsPerNode);

    for(int i = 0; i < numUnknownsPerNode; i++){
      x_views[i] = Teuchos::arrayView(fx[i], numBoxNodes);
    }

    ArrayView<ArrayView<const SCALAR> > my_view(x_views);
    ArrayRCP<ArrayRCP<SCALAR> > b_data = linprobmgr_->applyMatrix(my_view);
    for(int i = 0; i < b_data.size(); i++){
      for(int j = 0; j < b_data[i].size(); j++){
	b[i][j] = Teuchos::as<double>(b_data[i][j]);
      }
    }
    for (int i = 0; i < numUnknownsPerNode; ++i)
      delete [] fx[i];
    delete [] fx;

    return( 0 );
  }

  int dft_linprobmgr_importr2c(void * linprobmgr, double** x, double **b) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    int numOwnedNodes = linprobmgr_->getNumOwnedNodes();
    int numUnknownsPerNode = linprobmgr_->getNumUnknownsPerNode();

    SCALAR **fx;
    fx = new SCALAR*[numUnknownsPerNode];
    for (int i=0; i<numUnknownsPerNode;i++)
      fx[i] = new SCALAR[numOwnedNodes];
    for(int i=0;i<numUnknownsPerNode;i++)
      for(int j=0;j<numOwnedNodes;j++)
	fx[i][j] = Teuchos::as<SCALAR>(x[i][j]);

    ArrayRCP<ArrayRCP<const SCALAR> > my_x = Teuchos::arcp<ArrayRCP<const SCALAR> >(numUnknownsPerNode);
    for(int i = 0; i < numUnknownsPerNode; i++){
      my_x[i] = Teuchos::arcp<SCALAR>(fx[i], 0, numOwnedNodes, false);
    }

    ArrayRCP<ArrayRCP<SCALAR> > my_b = linprobmgr_->importR2C(my_x);
    for(int i = 0; i < my_b.size(); i++){
      for(int j = 0; j < my_b[i].size(); j++){
	b[i][j] = Teuchos::as<double>(my_b[i][j]);
      }
    }
    for (int i = 0; i < numUnknownsPerNode; ++i)
      delete [] fx[i];
    delete [] fx;

    return( 0 );
  }

  int dft_linprobmgr_importsingleunknownr2c(void * linprobmgr, double* x, double *b) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    int numOwnedNodes = linprobmgr_->getNumOwnedNodes();
    int numUnknownsPerNode = linprobmgr_->getNumUnknownsPerNode();

    SCALAR *fx;
    fx = new SCALAR[numOwnedNodes];
    for(int i=0;i<numOwnedNodes;i++)
      fx[i] = Teuchos::as<SCALAR>(x[i]);

    ArrayRCP<SCALAR> my_x = Teuchos::arcp<SCALAR> (numOwnedNodes);

    my_x = Teuchos::arcp<SCALAR>(fx, 0, numOwnedNodes, false);

    ArrayRCP<SCALAR> my_b = linprobmgr_->importR2C(my_x);
    for(int j = 0; j < my_b.size(); j++){
      b[j] = Teuchos::as<double>(my_b[j]);
    }
    delete [] fx;

    return( 0 );
  }

  int dft_linprobmgr_importnodalr2c(void * linprobmgr, double* x, double *b) {
    BLPM * linprobmgr_ = (BLPM *) linprobmgr;
    int numOwnedNodes = linprobmgr_->getNumOwnedNodes();
    int numUnknownsPerNode = linprobmgr_->getNumUnknownsPerNode();

    SCALAR *fx;
    fx = new SCALAR[numOwnedNodes];
    for(int i=0;i<numOwnedNodes;i++)
      fx[i] = Teuchos::as<SCALAR>(x[i]);

    ArrayRCP<SCALAR> ret_x = Teuchos::arcp<SCALAR>(fx, 0, numOwnedNodes, false);

    ArrayRCP<SCALAR> ret_val = linprobmgr_->importR2C(ret_x);
    for(int i = 0; i < ret_val.size(); i++){
      b[i] = Teuchos::as<double>(ret_val[i]);
    }
    delete [] fx;

    return( 0 );
  }

#ifdef __cplusplus
}
#endif
