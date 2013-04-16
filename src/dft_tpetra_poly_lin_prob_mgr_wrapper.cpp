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

#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_TpetraPolyLinProbMgr.hpp"

#ifdef __cplusplus
extern "C"
{
#endif

//template class dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal>;
//template class dft_PolyLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal>;

typedef Teuchos::Comm<int> COMM;

#if LINSOLVE_PREC == 0
// Solver precision
#define SCALAR                 float
// Matrix precision
#define MAT_SCALAR             SCALAR

#elif LINSOLVE_PREC == 1
// Solver precision
#define SCALAR                 double
// Matrix precision
#if MIXED_PREC == 1
#define MAT_SCALAR             float
#else
#define MAT_SCALAR             SCALAR
#endif

#elif LINSOLVE_PREC == 2
#include <qd/dd_real.h>
// Solver precision
#define SCALAR                 dd_real
// Matrix precision
#if MIXED_PREC == 1
#define MAT_SCALAR             double
#else
#define MAT_SCALAR             SCALAR
#endif

#elif LINSOLVE_PREC == 3
#include <qd/qd_real.h>
// Solver precision
#define SCALAR                 qd_real
// Matrix precision
#if MIXED_PREC == 1
#define MAT_SCALAR             dd_real
#else
#define MAT_SCALAR             SCALAR
#endif
#endif

#if NODE_TYPE == 0
#define NODE Kokkos::TPINode
#elif NODE_TYPE == 1
#define NODE Kokkos::TBBNode
#elif NODE_TYPE == 2
#define NODE Kokkos::OpenMPNode
#else
#define NODE Kokkos::SerialNode
#endif

#if PLATFORM_TYPE == 0
#define PLATFORM Tpetra::SerialPlatform<NODE>
#else
#define PLATFORM Tpetra::MpiPlatform<NODE>
#endif

#define THREAD_NUMBER 1

typedef dft_BasicLinProbMgr<SCALAR,MAT_SCALAR,int,int,NODE> BLPM;
typedef dft_PolyLinProbMgr<SCALAR,MAT_SCALAR,int,int,NODE> PLPM;

  /*****************************************************/
  /**                  dft_BasicLinProbMgr            **/
  /*****************************************************/

void *
dft_poly_lin_prob_mgr_create
(int numUnks, void * Parameterlist_list, MPI_Comm comm)
{
  RCP<ParameterList> my_list = rcp( (ParameterList *) Parameterlist_list, false);
  ParameterList nodeParams;
  nodeParams.set<int>("Num Threads", THREAD_NUMBER);
  RCP<NODE> my_node = rcp(new NODE(nodeParams));
  RCP<PLATFORM> platform = rcp(new PLATFORM(my_node));
  const RCP<const COMM> my_comm = platform->getComm();

  PLPM * tmp = new PLPM(numUnks, my_list, my_comm, my_node);
  BLPM * linprobmgr_ = dynamic_cast<BLPM *>(tmp);
  return((void *)linprobmgr_);
}

void *
dft_poly_lin_prob_mgr_create_debug
(int numUnks, void * Parameterlist_list, MPI_Comm comm)
{
  RCP<ParameterList> my_list = rcp( (ParameterList *) Parameterlist_list, false);
  ParameterList nodeParams;
  nodeParams.set<int>("Num Threads", THREAD_NUMBER);
  RCP<NODE> my_node = rcp(new NODE(nodeParams));
  RCP<PLATFORM> platform = rcp(new PLATFORM(my_node));
  const RCP<const COMM> my_comm = platform->getComm();

  PLPM * tmp = new PLPM(numUnks, my_list, my_comm, my_node, true);
  BLPM * linprobmgr_ = dynamic_cast<BLPM *>(tmp);
  return((void *)linprobmgr_);
}

void
dft_poly_lin_prob_mgr_destruct
(void * linprobmgr)
{
  BLPM * tmp = (BLPM *) linprobmgr;
  PLPM * linprobmgr_ = dynamic_cast<PLPM *>(tmp);
  delete linprobmgr_;
}

int
dft_poly_lin_prob_mgr_setgequationids
(void * linprobmgr, int numgids, int * gids)
{
  BLPM * tmp = (BLPM *) linprobmgr;
  PLPM * linprobmgr_ = dynamic_cast<PLPM *>(tmp);

  ArrayView<const int> gid_arr((numgids == 0)? NULL : gids, numgids);
  linprobmgr_->setGEquationIDs(gid_arr);
  return 0;
}

int
dft_poly_lin_prob_mgr_setginvequationids
(void * linprobmgr, int numgids, int * gids)
{
  BLPM * tmp = (BLPM *) linprobmgr;
  PLPM * linprobmgr_  = dynamic_cast<PLPM *>(tmp);

  ArrayView<const int> gid_arr((numgids == 0)? NULL : gids, numgids);
  linprobmgr_->setGInvEquationIDs(gid_arr);
  return 0;
}

int
dft_poly_lin_prob_mgr_setcmsequationids
(void * linprobmgr, int numgids, int * gids)
{
  BLPM * tmp = (BLPM *) linprobmgr;
  PLPM * linprobmgr_ = dynamic_cast<PLPM *>(tmp);

  ArrayView<const int> gid_arr((numgids == 0)? NULL : gids, numgids);
  linprobmgr_->setCmsEquationIDs(gid_arr);
  return 0;
}

int
dft_poly_lin_prob_mgr_setdensityequationids
(void * linprobmgr, int numgids, int * gids)
{
  BLPM * tmp = (BLPM *) linprobmgr;
  PLPM * linprobmgr_ = dynamic_cast<PLPM *>(tmp);

  ArrayView<const int> gid_arr((numgids == 0)? NULL : gids, numgids);
  linprobmgr_->setDensityEquationIDs(gid_arr);
  return 0;
}

int
dft_poly_lin_prob_mgr_setfieldondensityislinear
(void * linprobmgr, int isLinear)
{
  BLPM * tmp = (BLPM *) linprobmgr;
  PLPM * linprobmgr_ = dynamic_cast<PLPM *>(tmp);
  bool isLinear1 = (isLinear!=0);
  linprobmgr_->setFieldOnDensityIsLinear(isLinear1);
  return 0;
}

int
dft_poly_lin_prob_mgr_setpoissonequationids
(void * linprobmgr, int numgids, int * gids)
{
  BLPM * tmp = (BLPM *) linprobmgr;
  PLPM * linprobmgr_ = dynamic_cast<PLPM *>(tmp);

  ArrayView<const int> gid_arr((numgids == 0)? NULL : gids, numgids);
  linprobmgr_->setPoissonEquationIDs(gid_arr);
  return 0;
}

#ifdef __cplusplus
}
#endif
