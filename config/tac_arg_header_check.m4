# This is a simplified version of the TRILINOS_HEADER macro
# used by NOX.
# We want to see if we can find <file name> in the available
# include directories.
# Author: Jim Willenbring jmwille@sandia.gov
# HEADER_INCLUDE(<file name>)
AC_DEFUN([HEADER_INCLUDE],
[
  done=no

AC_MSG_CHECKING(for [$1])

  AC_PREPROC_IFELSE([AC_LANG_SOURCE(
  [[
#undef PACKAGE
#undef PACKAGE_NAME
#undef PACKAGE_BUGREPORT
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#include "$1"
  ]])],[done=yes],[])

  if test $done = no; then
    echo "------"
    echo "Cannot preprocess the include file $1."
    echo "Try --with-incdir=\"-I<dir1> -I<dir2>\""
    echo "------"
    AC_MSG_ERROR([Cannot find $2])
  else
    AC_MSG_RESULT(yes)
  fi
])

