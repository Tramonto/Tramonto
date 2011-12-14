# This file will setup the the CMake cache for a Trilinos build that
# has all of the required Tramonto packages. To make use of this,
# configure Trilinos with the following CMake line:
#
# cmake [cmake-defines] -C Tramonto/config/TrilinosCache.cmake /path/to/Trilinos
#
# At a minimum, you probably want to add the flag:
# 
# -DCMAKE_INSTALL_PREFIX=/your/install/path
#
# To choose the installation prefix.

macro(setcache var value type)
  message("Setting ${var} to ${value} for the Tramonto build.")
  set(${var} ${value} CACHE ${type} "Value set for Tramonto")
endmacro()

setcache(TPL_ENABLE_MPI ON BOOL)
setcache(Trilinos_ENABLE_ALL_PACKAGES OFF BOOL)
setcache(Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES OFF BOOL)
setcache(Trilinos_ENABLE_Amesos ON BOOL)
setcache(Trilinos_ENABLE_AztecOO ON BOOL)
setcache(Trilinos_ENABLE_Triutils ON BOOL)
setcache(Trilinos_ENABLE_Ifpack ON BOOL)
setcache(Trilinos_ENABLE_EpetraExt ON BOOL)
setcache(Trilinos_ENABLE_Epetra ON BOOL)
setcache(NOX_ENABLE_LOCA ON BOOL)
setcache(Trilinos_ENABLE_ML ON BOOL)
setcache(Trilinos_ENABLE_NOX ON BOOL)
setcache(Trilinos_ENABLE_Teuchos ON BOOL)

