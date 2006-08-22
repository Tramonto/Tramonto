
/*@HEADER
************************************************************************

************************************************************************
//@HEADER
*/

#ifndef TRAMONTO_CONFIGDEFS_H
#define TRAMONTO_CONFIGDEFS_H

#ifdef HAVE_CONFIG_H

/*
 * The macros PACKAGE, PACKAGE_NAME, etc, get defined for each Trilinos package
 * and need to be undef'd here to avoid warnings due to redefinition.
 * Kevin Long came up with this solution.
 */

#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

#include <Tramonto_config.h>

/*#else HAVE_CONFIG_H is not defined*/

#endif /*HAVE_CONFIG_H*/

#endif /* TRAMONTO_CONFIGDEFS_H */
