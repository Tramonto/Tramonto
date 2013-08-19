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

#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "dft_TpetraHardSphereLinProbMgr.hpp"

#ifdef __cplusplus
extern "C"
{
#endif

typedef Teuchos::Comm<int> COMM;

TRAMONTO_TYPEDEF_HELPER(dft_BasicLinProbMgr,BLPM)
TRAMONTO_TYPEDEF_HELPER(dft_HardSphereLinProbMgr,HSLPM)

  /*****************************************************/
  /**                  dft_HardSphereLinProbMgr       **/
  /*****************************************************/
void *
dft_hardsphere_lin_prob_mgr_create
(int numUnks, void * Parameterlist_list, MPI_Comm comm)
{
  RCP<ParameterList> my_list = rcp( (ParameterList *) Parameterlist_list, false);
  ParameterList nodeParams;
  nodeParams.set<int>("Num Threads", NUM_THREADS);
  RCP<NODE> my_node = rcp(new NODE(nodeParams));
  RCP<PLATFORM> platform = rcp(new PLATFORM(my_node));
  const RCP<const COMM> my_comm = platform->getComm();

  HSLPM * tmp = new HSLPM(numUnks, my_list, my_comm, my_node);
  BLPM* linprobmgr_ = dynamic_cast<BLPM*>(tmp);
  return((void *)linprobmgr_);
}

void *
dft_hardsphere_lin_prob_mgr_create_debug
(int numUnks, void * Parameterlist_list, MPI_Comm comm)
{
  RCP<ParameterList> my_list = rcp((ParameterList *) Parameterlist_list, false);
  ParameterList nodeParams;
  nodeParams.set<int>("Num Threads", NUM_THREADS);
  RCP<NODE> my_node = rcp(new NODE(nodeParams));
  RCP<PLATFORM> platform = rcp(new PLATFORM(my_node));
  const RCP<const COMM> my_comm = platform->getComm();

  HSLPM * tmp = new HSLPM(numUnks, my_list, my_comm, my_node, true);
  BLPM * linprobmgr_ = dynamic_cast<BLPM *>(tmp);
  return((void *)linprobmgr_);
}

void
dft_hardsphere_lin_prob_mgr_destruct
(void * linprobmgr)
{
  BLPM * tmp = (BLPM *) linprobmgr;
  HSLPM * linprobmgr_ = dynamic_cast<HSLPM *>(tmp);
  delete linprobmgr_;
}

int
dft_hardsphere_lin_prob_mgr_setindnonlocalequationids
(void * linprobmgr, int numgids, int * gids)
{
  BLPM * tmp = (BLPM *) linprobmgr;
  HSLPM * linprobmgr_ = dynamic_cast<HSLPM *>(tmp);

  ArrayView<const int> gid_arr((numgids == 0)? NULL : gids, numgids);
  linprobmgr_->setIndNonLocalEquationIDs(gid_arr);
  return 0;
}

int
dft_hardsphere_lin_prob_mgr_setdepnonlocalequationids
(void * linprobmgr, int numgids, int * gids)
{
  BLPM * tmp = (BLPM *) linprobmgr;
  HSLPM * linprobmgr_ = dynamic_cast<HSLPM *>(tmp);

  ArrayView<const int> gid_arr((numgids == 0)? NULL : gids, numgids);
  linprobmgr_->setDepNonLocalEquationIDs(gid_arr);
  return 0;
}

int
dft_hardsphere_lin_prob_mgr_setdensityequationids
(void * linprobmgr, int numgids, int * gids)
{
  BLPM * tmp = (BLPM *) linprobmgr;
  HSLPM * linprobmgr_ = dynamic_cast<HSLPM *>(tmp);

  ArrayView<const int> gid_arr((numgids == 0)? NULL : gids, numgids);
  linprobmgr_->setDensityEquationIDs(gid_arr);
  return 0;
}

int
dft_hardsphere_lin_prob_mgr_seta22blockisdiagonal
(void * linprobmgr, int isa22diagonal)
{
  BLPM * tmp = (BLPM *) linprobmgr;
  HSLPM * linprobmgr_ =  dynamic_cast<HSLPM *>(tmp);
  bool isA22Diagonal = (isa22diagonal!=0);
  linprobmgr_->setA22BlockIsDiagonal(isA22Diagonal);
  return 0;
}

#ifdef __cplusplus
}
#endif
