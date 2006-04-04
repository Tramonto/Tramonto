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
dnl   FFLAGS   = F77 compiler flags (used for determining FLIBS)
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
AC_PREREQ(2.59)
AC_CONFIG_SRCDIR([src/dft_main.c])
AM_CONFIG_HEADER([src/Tramonto_config.h:src/Tramonto_config.h.in])

dnl Specify directory for auxillary build tools (e.g., install-sh,
dnl config.sub, config.guess) and M4 files.
AC_CONFIG_AUX_DIR(config)
dnl Allow users to specify their own "install" command.  If none is specified,
dnl the default is install-sh found in the config subdirectory.
                                                                                
AC_ARG_WITH(install,
 [AC_HELP_STRING([--with-install=INSTALL_PROGRAM],
 [Use the installation program INSTALL_PROGRAM rather the default that is provided.  For example --with-install="/path/install -p"])],
 [
   INSTALL=$withval
   INSTALL_PROGRAM=$withval
   INSTALL_SCRIPT=$withval
   INSTALL_DATA="$withval -m 644"
 ],)
                                                                                
dnl AM_MAINTAINER_MODE turns off maintainer-only makefile targets by
dnl default, and changes configure to understand a
dnl --enable-maintainer-mode option. --enable-maintainer-mode turns the
dnl maintainer-only targets back on. The maintainer-only makefile
dnl targets permit end users to clean automatically-generated files such
dnl as configure, which means they have to have autoconf and automake
dnl installed to repair the damage. AM_MAINTAINER_MODE makes it a bit
dnl harder for users to shoot themselves in the foot.
                                                                                
AM_MAINTAINER_MODE
                                                                                
dnl Define $build, $host, $target, etc
                                                                                
AC_CANONICAL_TARGET
                                                                                
dnl  Required version of automake.
AM_INIT_AUTOMAKE(1.8.2)

dnl
dnl check for ranlib incase uses wants to build library
dnl
AC_PROG_RANLIB

AC_ARG_WITH(mpi,
[AC_HELP_STRING([--with-mpi=MPIROOT],[use MPI root directory])],
[
  HAVE_PKG_MPI=yes
  MPI_DIR=${withval}
  AC_MSG_CHECKING(MPI directory)
  AC_MSG_RESULT([${MPI_DIR}])
]
)

AC_ARG_WITH(mpi-libs,
[AC_HELP_STRING([--with-mpi-libs="LIBS"],[MPI libraries @<:@"-lmpi"@:>@])],
[
  MPI_LIBS=${withval}
  AC_MSG_CHECKING(user-defined MPI libraries)
  AC_MSG_RESULT([${MPI_LIBS}])
]
)

AC_ARG_WITH(mpi-incdir,
[AC_HELP_STRING([--with-mpi-incdir=DIR],[MPI include directory @<:@MPIROOT/include@:>@  Do not use -I])],
[
  MPI_INC=${withval}
  AC_MSG_CHECKING(user-defined MPI includes)
  AC_MSG_RESULT([${MPI_INC}])
]
)

AC_ARG_WITH(mpi-libdir,
[AC_HELP_STRING([--with-mpi-libdir=DIR],[MPI library directory @<:@MPIROOT/lib@:>@  Do not use -L])],
[
  MPI_LIBDIR=${withval}
  AC_MSG_CHECKING(user-defined MPI library directory)
  AC_MSG_RESULT([${MPI_LIBDIR}])
]
)

dnl
dnl check for C++, C and Fotran compilers
dnl
AC_PROG_CXX( [mpicxx mpic++ mpiCC c++ g++ cxx CC] )
AC_PROG_CC( [mpicc cc gcc] )
AC_PROG_F77( [mpif77 f77 g77 f90 xlf90 f95] )

TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
TAC_ARG_WITH_FLAGS(ldflags, LDFLAGS)

dnl
dnl check for fortran to C libraries
dnl
AC_F77_LIBRARY_LDFLAGS

if test -n "${MPI_DIR}" && test -z "${MPI_INC}"; then
  MPI_INC="${MPI_DIR}/include"
fi

if test -n "${MPI_INC}"; then
  CPPFLAGS="${CPPFLAGS} -I${MPI_INC}"
fi

if test ${CC} != "mpicc"
then
  AC_CHECK_HEADER( [mpi.h],,
                   [AC_MSG_FAILURE([Cannot find mpi.h.  Please inlude its location using --with-cflags and --with-cxxflags])])
fi;

AC_LANG_CPLUSPLUS
AC_MSG_CHECKING(for mpi.h)
AC_TRY_CPP([#include "mpi.h"],
  [AC_MSG_RESULT(yes)],
  [
   AC_MSG_RESULT(no)
   echo "-----"
   echo "Cannot link simple MPI program."
   echo "Try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir"
   echo "to specify all the specific MPI compile options."
   echo "-----"
   AC_MSG_ERROR(MPI cannot link)
  ])

#Process any mpi directory information that was passed on the configure line
if test -n "${MPI_DIR}" && test -z "${MPI_LIBDIR}"; then
  MPI_LIBDIR="${MPI_DIR}/lib"
fi

if test -n "${MPI_LIBDIR}"; then
  LDFLAGS="${LDFLAGS} -L${MPI_LIBDIR}"
fi

if test -z "${MPI_LIBS}" && test -n "${MPI_LIBDIR}"; then
  MPI_LIBS="-lmpi"
fi


dnl ------------------------------------------------------------------------
dnl Checks if examples and GUI support should be built
dnl ------------------------------------------------------------------------

dnl If a test directory is ever added, add a check for tests to the
dnl Checks for libraries section below.                 
dnl TAC_ARG_ENABLE_FEATURE(tests, [Build Tramonto tests], TESTS, yes)
dnl AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_tests" != "Xno")
TAC_ARG_ENABLE_FEATURE(examples, [Build Tramonto examples], EXAMPLES, yes)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_examples" != "Xno")

#TAC_ARG_ENABLE_FEATURE(gui, [Build Tramonto GUI], GUI, yes)
#AM_CONDITIONAL(BUILD_GUI, test "X$ac_cv_use_gui" != "Xno")

# ------------------------------------------------------------------------
# Specify other directories
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS


dnl Option to turn off checks for third-party libraries

TAC_ARG_ENABLE_FEATURE(libcheck, [Check for some third-party libraries including BLAS and LAPACK.  (Cannot be disabled unless tests and examples are also disabled.)], LIBCHECK, yes)

dnl ------------------------------------------------------------------------
dnl Checks for libraries
dnl ------------------------------------------------------------------------
                                                                                
dnl If a test directory is ever added, add a check for tests here
dnl If examples and libcheck are disabled, we don't have to check
dnl for these libraries.
if test "X$ac_cv_use_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

dnl
dnl check for math library
dnl
AC_SEARCH_LIBS( [pow], 
                [dxml m],,
                [AC_MSG_FAILURE([Cannot link to math library functions.])] )

  ACX_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, specify a path using --with-ldflags=-I<DIR> (ex. --with-ldflags=-I/usr/path/lib) or a specific library using --with-blas=DIR/LIB (for example --with-blas=/usr/path/lib/libcxml.a)]])])
  ACX_LAPACK([], [AC_MSG_ERROR([[Cannot find lapack library, specify a path using --with-ldflags=-I<DIR> (ex. --with-dlflags=-I/usr/path/lib) or a specific library using --with-lapack=DIR/LIB (ex. --with-lapack=/usr/path/lib/libcxml.a)]])])

dnl
dnl check for AztecOO library
dnl
AC_LANG_CPLUSPLUS
AC_CHECK_LIB( [aztecoo], 
              [AZ_iterate],
              [LIBS="-laztecoo -ltriutils -lifpack -lepetraext -lepetra -lteuchos $LAPACK_LIBS  $BLAS_LIBS $LIBS $FLIBS"],
              [AC_MSG_FAILURE([Cannot link to AztecOO, Triutils, Ifpack, EpetraExt, Epetra, or Teuchos library.])], 
              [-ltriutils -lifpack -lepetraext -lepetra -lteuchos ${LAPACK_LIBS} ${BLAS_LIBS} ${LIBS} ${FLIBS}])

fi

dnl
dnl If we're not using an mpi compiler, then see if we can find
dnl the mpi.h  mpi compilers automatically include the right
dnl path to find mpi.h
dnl

dnl
dnl piont to the make files
dnl
AC_CONFIG_FILES([
		Makefile
		src/Makefile
		Examples/Makefile])

dnl
dnl make the output
dnl
AC_OUTPUT
