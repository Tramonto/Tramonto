dnl
dnl configure.ac
dnl
dnl This script sets the following variables involved
dnl in building Tramonto. The variables are set from
dnl the command line arguments the user passes to 
dnl the configure script, or automatically within this
dnl script.
dnl
dnl Variables:
dnl   CC       = the c compiler
dnl   CXX      = the c++ compiler
dnl   F77      = the fortran compiler
dnl   CFLAGS   = c compiler flags
dnl   CXXFLAGS = c++ compiler flags
dnl   LDFLAGS  = linker flags
dnl   FLIBS    = libraries to link fortran to c/c++
dnl   RANLIB   = ranlib command if available
dnl 
dnl Any of these variables can be set on the 
dnl command line when invoking configure with the 
dnl following formats:
dnl
dnl ./configure CC=/usr/bin/gcc
dnl
dnl or
dnl 
dnl ./configure --with-cc=/usr/bin/gcc
dnl
dnl
dnl Process this file with autoconf to produce the
dnl confiure script
dnl

dnl
dnl Initialize autoconf
dnl
AC_INIT(Tramonto, 1.0, ljfrink@sandia.gov, tramonto)
AC_PREREQ(2.57)
AC_CONFIG_SRCDIR([src/dft_main.c])

dnl
dnl see if a specifc C compiler was specified
dnl
AC_ARG_WITH( [cc],
               [AC_HELP_STRING([--with-cc=...], [use specified c compiler])], 
               [CC="${withval}", CC=""])

dnl
dnl see if a specifc C++ compiler was specified
dnl
AC_ARG_WITH( [cxx],
               [AC_HELP_STRING([--with-cxx=...], [use specified c++ compiler])], 
               [CXX="${withval}", CXX=""])

dnl
dnl see if a specifc F77 compiler was specified
dnl
AC_ARG_WITH( [f77],
               [AC_HELP_STRING([--with-f77=...], [use specified f77 compiler])], 
               [F77="${withval}", F77=""])

dnl
dnl see if any extra c compiler flags were specified
dnl
AC_ARG_WITH( [cflags],
               [AC_HELP_STRING([--with-cflags=...], [use cflags flags ... with c compiler])], 
               [CFLAGS="${CFLAGS} ${withval}", C_FLAGS=""])

dnl
dnl see if any extra c++ compiler flags were specified
dnl
AC_ARG_WITH( [cxxflags],
               [AC_HELP_STRING([--with-cxxflags=...], [use cxxflags flags ... with c++ compiler])],
               [CXXFLAGS="${CXXFLAGS} ${withval}"])

dnl
dnl see if any extra linker flags were specified
dnl
AC_ARG_WITH( [ldflags], 
               [AC_HELP_STRING([--with-ldflags=...], [use linker flags ...])], 
               [LIBS="${LIBS} ${withval}"])

dnl
dnl set up to use make
dnl
AC_PROG_MAKE_SET

dnl
dnl check for ranlib incase uses wants to build library
dnl
AC_PROG_RANLIB

dnl
dnl check for C++, C and Fotran compilers
dnl
AC_PROG_CXX( [mpic++ mpiCC c++ g++] )
AC_PROG_CC( [mpicc cc gcc] )
AC_PROG_F77( [mpif77 f77 g77] )

dnl
dnl check for fortran to C libraries
dnl
AC_F77_LIBRARY_LDFLAGS

dnl
dnl check for math library
dnl
AC_SEARCH_LIBS( [pow], 
                [dxml m],,
                [AC_MSG_FAILURE([Cannot link to math library functions.])] )

dnl
dnl check for blas library
dnl
AC_LANG_FORTRAN77
AC_SEARCH_LIBS( [idamax], 
                [blas dxml cxml csmath complib.sgimath aztec],,
                [AC_MSG_FAILURE([Cannot link to blas library functions.])] )

dnl
dnl check for lapack library
dnl
AC_SEARCH_LIBS( [DGETRF], 
                [lapack dxml cxml csmath complib.sgimath aztec],,
                [AC_MSG_FAILURE([Cannot link to lapack library functions.])] )

dnl
dnl check for AztecOO library
dnl
AC_LANG_CPLUSPLUS
AC_CHECK_LIB( [aztecoo], 
              [AZ_iterate],
              [LIBS="-laztecoo -ltriutils -lifpack -lepetraext -lepetra ${FLIBS} $LIBS"],
              [AC_MSG_FAILURE([Cannot link to AztecOO library.])], 
              [-ltriutils -lifpack -lepetraext -lepetra ${FLIBS}])

dnl
dnl If we're not using an mpi compiler, then see if we can find
dnl the mpi.h  mpi compilers automatically include the right
dnl path to find mpi.h
dnl
if test ${CC} != "mpicc"
then
  AC_CHECK_HEADER( [mpi.h],,
                   [AC_MSG_FAILURE([Cannot find mpi.h.  Please inlude its location using --with-cflags and --with-cxxflags])])
fi;


dnl
dnl piont to the make files
dnl
AC_CONFIG_FILES([Makefile src/Makefile])

dnl
dnl make the output
dnl
AC_OUTPUT
