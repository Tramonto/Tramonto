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
 *  File:  dft_direct_solver_const.h
 *
 *  This file contains the Amesos solver constants used in dft_mesh
 *  and dft_BasicLinProbMgr.
 *
 */
/****************************************************************************/

/* Linear (direct) solver constants */
#define AM_lapack    15
#define AM_klu       16
#define AM_mumps     17
#define AM_umfpack   18
#define AM_superlu   19
#define AM_superludist 20
#define AM_pardiso   21
#define AM_taucs     22
//#define AM_none      0
