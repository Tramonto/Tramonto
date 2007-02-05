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

/* Data structures for parallel RCB

      rcb_dot has 4 required fields as shown below
      other fields can be added by user
      are just carried along as dots migrate to new processors

      examples:

  int       global;                global id # of dot
  int       local;                 local id # (memory loc) of dot before RCB
  int       proc;                  owner of this dot before RCB
*/

                             /* dot to balance on for RCB */ 
struct rcb_dot {	        /* dot = point in 3-space */
  double    x[3];		/* location of dot */
  double    weight;             /* weight of dot - if used must be > 0 */
};

struct rcb_tree {	     /* tree of RCB cuts */
  double    cut;        	/* position of cut */
  int       dim;	        /* dimension (012) of cut */
};

struct rcb_median {          /* RCB cut info */
  double    lototal, hitotal;   /* weight in each half of active partition */
  double    valuelo, valuehi;	/* position of dot(s) nearest to cut */
  double    wtlo, wthi;         /* total weight of dot(s) at that position */
  int       countlo, counthi;   /* # of dots at that position */
  int       proclo, prochi;	/* unique proc who owns a nearest dot */
};

struct rcb_box {       	     /* bounding box */
  double    lo[3], hi[3];	/* xyz lo/hi bounds */
};

