/*
//@HEADER
// ********************************************************************
// Copyright (2006) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000, there is a non-exclusive license for use of this
// work by or on behalf of the U.S. Government. Export of this program
// may require a license from the United States Government.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// ********************************************************************
//@HEADER
*/

/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/


/* function declarations for dynamic array allocation */

extern void *array_alloc (
	int numdim,
	...
);


extern void * array_alloc_1d ( size_t n1, size_t size );
extern void * array_alloc_2d ( size_t n1, size_t n2, size_t size );
extern void * array_alloc_3d ( size_t n1, size_t n2, size_t n3, size_t size );
extern void * array_alloc_4d ( size_t n1, size_t n2, size_t n3, size_t n4, size_t size );

extern void   safe_free ( void **ptr) ;

