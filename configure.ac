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

dnl ------------------------------------------------------------------------
dnl Checks if examples and GUI support should be built
dnl ------------------------------------------------------------------------
                                                               
dnl If a test directory is ever added, add a check for tests to the
dnl Checks for libraries section below.                 
dnl TAC_ARG_ENABLE_FEATURE(tests, [Build Tramonto tests], TESTS, yes)
dnl AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_tests" != "Xno")                                                                                
TAC_ARG_ENABLE_FEATURE(examples, [Build Tramonto examples], EXAMPLES, yes)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_examples" != "Xno")

TAC_ARG_ENABLE_FEATURE(gui, [Build Tramonto GUI], GUI, yes)
AM_CONDITIONAL(BUILD_GUI, test "X$ac_cv_use_gui" != "Xno")

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

dnl
dnl check for blas library
dnl
AC_LANG_FORTRAN77
AC_SEARCH_LIBS( [idamax], 
                [blas dxml cxml csmath complib.sgimath],,
                [AC_MSG_FAILURE([Cannot link to blas library functions.])] )

dnl
dnl check for lapack library
dnl
AC_SEARCH_LIBS( [DGETRF], 
                [lapack dxml cxml csmath complib.sgimath],,
                [AC_MSG_FAILURE([Cannot link to lapack library functions.])] )

dnl
dnl check for AztecOO library
dnl
AC_LANG_CPLUSPLUS
AC_CHECK_LIB( [aztecoo], 
              [AZ_iterate],
              [LIBS="-laztecoo -ltriutils -lifpack -lepetraext -lepetra -lteuchos ${FLIBS} $LIBS"],
              [AC_MSG_FAILURE([Cannot link to AztecOO, Triutils, Ifpack, EpetraExt, Epetra, or Teuchos library.])], 
              [-ltriutils -lifpack -lepetraext -lepetra -lteuchos ${FLIBS}])

fi

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
AC_CONFIG_FILES([
		Makefile
		src/Makefile
		Examples/Makefile])

dnl
dnl make the output
dnl
AC_OUTPUT
