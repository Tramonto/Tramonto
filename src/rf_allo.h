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

