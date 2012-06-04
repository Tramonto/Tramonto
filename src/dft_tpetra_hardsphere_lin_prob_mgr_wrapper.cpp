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

//template class dft_HardSphereLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal>;
//template class dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal>;

typedef Teuchos::Comm<int> COMM;

#if LINSOLVE_PREC == 0
// Use float
typedef dft_BasicLinProbMgr<float,int,int> BLPM;
typedef dft_HardSphereLinProbMgr<float,int,int> HSLPM;
#define WORKING_PREC float
#define WORKING_CAST( x ) float(x)
#define DOUBLE_CAST( x ) double(x)

#elif LINSOLVE_PREC == 1
// Use double
typedef dft_BasicLinProbMgr<double,int,int> BLPM;
typedef dft_HardSphereLinProbMgr<double,int,int> HSLPM;
#define WORKING_PREC double
#define WORKING_CAST( x ) (x)
#define DOUBLE_CAST( x ) (x)

#elif LINSOLVE_PREC == 2
// Use quad double
typedef dft_BasicLinProbMgr<qd_real,int,int> BLPM;
typedef dft_HardSphereLinProbMgr<qd_real,int,int> HSLPM;
#define WORKING_PREC qd_real
#include <qd/qd_real.h>
#define WORKING_CAST( x ) qd_real(x)
#define DOUBLE_CAST( x ) to_double(x)

#elif LINSOLVE_PREC == 3
// Use double double
typedef dft_BasicLinProbMgr<dd_real,int,int> BLPM;
typedef dft_HardSphereLinProbMgr<dd_real,int,int> HSLPM;
#define WORKING_PREC dd_real
#include <qd/dd_real.h>
#define WORKING_CAST( x ) dd_real(x)
#define DOUBLE_CAST( x ) to_double(x)
#endif


  /*****************************************************/
  /**                  dft_HardSphereLinProbMgr       **/
  /*****************************************************/
void *
dft_hardsphere_lin_prob_mgr_create
(int numUnks, void * Parameterlist_list, MPI_Comm comm)
{
  RCP<ParameterList> my_list = rcp( (ParameterList *) Parameterlist_list, false);
  RCP<const COMM> my_comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  HSLPM * tmp = new HSLPM(numUnks, my_list, my_comm);
  BLPM* linprobmgr_ = dynamic_cast<BLPM*>(tmp);
  return((void *)linprobmgr_);
}

void *
dft_hardsphere_lin_prob_mgr_create_debug
(int numUnks, void * Parameterlist_list, MPI_Comm comm)
{
  RCP<ParameterList> my_list = rcp((ParameterList *) Parameterlist_list, false);
  RCP<const COMM> my_comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  HSLPM * tmp = new HSLPM(numUnks, my_list, my_comm, true);
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

  ArrayView<const int> gid_arr(gids, numgids);
  linprobmgr_->setIndNonLocalEquationIDs(gid_arr);
  return 0;
}

int
dft_hardsphere_lin_prob_mgr_setdepnonlocalequationids
(void * linprobmgr, int numgids, int * gids)
{
  BLPM * tmp = (BLPM *) linprobmgr;
  HSLPM * linprobmgr_ = dynamic_cast<HSLPM *>(tmp);

  if(numgids == 0)
  {
    return 0;
  }

  ArrayView<const int> gid_arr(gids, numgids);
  linprobmgr_->setDepNonLocalEquationIDs(gid_arr);
  return 0;
}

int
dft_hardsphere_lin_prob_mgr_setdensityequationids
(void * linprobmgr, int numgids, int * gids)
{
  BLPM * tmp = (BLPM *) linprobmgr;
  HSLPM * linprobmgr_ = dynamic_cast<HSLPM *>(tmp);

  ArrayView<const int> gid_arr(gids, numgids);
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
