# @HEADER
#  ************************************************************************
#  
#                Tramonto:
#               Copyright () Sandia Corporation
#  
#  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
#  license for use of this work by or on behalf of the U.S. Government.
#  
#  This library is free software; you can redistribute it and or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation; either version 2.1 of the
#  License, or (at your option) any later version.
#   
#  This library is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#   
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
#  USA
#  Questions? Contact Laura Frink (ljfrink@sandia.gov) 
#  
#  ************************************************************************
# @HEADER

# NOTE: We also need to include the .h file as a "source" file, because this
# list is also used to package the distribution.

COMMON_H = \
	${srcdir}/Tramonto_ConfigDefs.h \
        ${srcdir}/dft_globals.h \
        ${srcdir}/dft_globals_const.h \
        ${srcdir}/dft_ldbal.h \
        ${srcdir}/rf_allo.h \
        ${srcdir}/loca_const.h \
        ${srcdir}/loca_util_const.h \
        ${srcdir}/az_aztec.h \
        ${srcdir}/az_aztec_defs.h \
        ${srcdir}/dft_BasicLinProbMgr.hpp \
        ${srcdir}/dft_PolyLinProbMgr.hpp \
        ${srcdir}/dft_HardSphereLinProbMgr.hpp \
        ${srcdir}/dft_PolyA11_Epetra_Operator.hpp \
        ${srcdir}/dft_PolyA22_Epetra_Operator.hpp \
        ${srcdir}/dft_A22Matrix_Epetra_Operator.hpp \
        ${srcdir}/dft_HardSphereA11_Epetra_Operator.hpp \
        ${srcdir}/dft_HardSphereA22_Epetra_Operator.hpp \
        ${srcdir}/dft_Schur_Epetra_Operator.hpp \
        ${srcdir}/dft_basic_lin_prob_mgr_wrapper.h \
        ${srcdir}/dft_poly_lin_prob_mgr_wrapper.h \
        ${srcdir}/dft_hardsphere_lin_prob_mgr_wrapper.h

COMMON_C = \
        ${srcdir}/dft_input.c \
        ${srcdir}/dft_main.c \
        ${srcdir}/dft_mesh.c \
        ${srcdir}/dft_mesh_surfaces.c \
        ${srcdir}/dft_mesh_lib.c \
        ${srcdir}/dft_newton.c \
        ${srcdir}/dft_ldbal.c \
        ${srcdir}/dft_fill_pde.c \
        ${srcdir}/dft_stencil.c \
        ${srcdir}/dft_thermo.c \
        ${srcdir}/dft_thermo_hs.c \
        ${srcdir}/dft_thermo_att.c \
        ${srcdir}/dft_thermo_cms.c \
        ${srcdir}/dft_thermo_elec_MSA.c \
        ${srcdir}/dft_thermo_id.c \
        ${srcdir}/dft_thermo_wtc.c \
        ${srcdir}/dft_out_main.c \
        ${srcdir}/dft_out_ads.c \
        ${srcdir}/out_energy_new.c \
        ${srcdir}/dft_energy_id.c \
        ${srcdir}/dft_energy_hs.c \
        ${srcdir}/dft_energy_att.c \
        ${srcdir}/dft_energy_wtc.c \
        ${srcdir}/dft_energy_elec.c \
        ${srcdir}/dft_energy_vext.c \
        ${srcdir}/dft_energy_mu.c \
        ${srcdir}/dft_energy_cms.c \
        ${srcdir}/dft_out_force.c \
        ${srcdir}/dft_out_flux.c \
        ${srcdir}/dft_out_profiles.c \
        ${srcdir}/dft_vext.c \
        ${srcdir}/dft_guess.c \
        ${srcdir}/dft_comm.c \
        ${srcdir}/dft_continuation.c \
        ${srcdir}/dft_utils.c \
        ${srcdir}/dft_fill_control.c \
        ${srcdir}/dft_fill_main.c \
        ${srcdir}/dft_fillp.c \
        ${srcdir}/dft_fill_rosen_rb.c \
        ${srcdir}/dft_fill_mf.c \
        ${srcdir}/dft_fill_wtc.c \
        ${srcdir}/dft_vext_coul.c \
        ${srcdir}/dft_uww.c \
        ${srcdir}/dft_images.c \
        ${srcdir}/dft_potentials.c \
        ${srcdir}/rf_allo.c \
        ${srcdir}/loca_lib.c \
        ${srcdir}/loca_bord.c \
        ${srcdir}/loca_util.c

COMMON_CXX = \
        ${srcdir}/dft_BasicLinProbMgr.cpp \
        ${srcdir}/dft_PolyLinProbMgr.cpp \
        ${srcdir}/dft_HardSphereLinProbMgr.cpp \
        ${srcdir}/dft_PolyA11_Epetra_Operator.cpp \
        ${srcdir}/dft_PolyA22_Epetra_Operator.cpp \
        ${srcdir}/dft_A22Matrix_Epetra_Operator.cpp \
        ${srcdir}/dft_HardSphereA11_Epetra_Operator.cpp \
        ${srcdir}/dft_HardSphereA22_Epetra_Operator.cpp \
        ${srcdir}/dft_Schur_Epetra_Operator.cpp \
        ${srcdir}/dft_basic_lin_prob_mgr_wrapper.cpp \
        ${srcdir}/dft_poly_lin_prob_mgr_wrapper.cpp \
        ${srcdir}/dft_hardsphere_lin_prob_mgr_wrapper.cpp

EXTRA_NOXLOCA_CXX = \
        ${srcdir}/NOXLOCA_Solver.cpp \
        ${srcdir}/NOXLOCA_Tramonto_Group.cpp \
        ${srcdir}/NOXLOCA_Tramonto_Vector.cpp

EXTRA_NOXLOCA_H = \
        ${srcdir}/NOXLOCA_Tramonto_Group.hpp \
        ${srcdir}/NOXLOCA_Tramonto_Vector.hpp

if HAVE_NOXLOCA

NOXLOCA_CXX = \
	$(EXTRA_NOXLOCA_CXX)

NOXLOCA_H = \
	$(EXTRA_NOXLOCA_H)

else

NOXLOCA_CXX =
NOXLOCA_H =

endif

# ------------------------------------------------------------------------
# Tramonto library specifications
# ------------------------------------------------------------------------

DFT_LIB = libdft.a

DFT_H = \
	$(COMMON_H) \
	$(NOXLOCA_H)

libdft_a_SOURCES = \
	$(COMMON_C) \
	$(COMMON_CXX) \
	$(NOXLOCA_CXX)

EXTRA_libdft_a_SOURCES = \
	$(NOXLOCA_H) \
	$(NOXLOCA_CXX)

#Trilinos has a solution for selecting alternate ar's if that is a problem

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
libdft_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# Tramonto program specifications
# ------------------------------------------------------------------------
                                                                                
bin_PROGRAMS = dft
                                                                                
dft_SOURCES = \
        $(srcdir)/dft_maindriver.c \
        $(COMMON_CXX) \
        $(COMMON_C)
                                                                                
dft_DEPENDENCIES = $(top_builddir)/src/libdft.a
                                                                                
dft_CXXFLAGS = -I$(srcdir)

dft_LDADD = $(top_builddir)/src/libdft.a $(LDFLAGS) $(LIBS) $(FLIBS)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(DFT_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(DFT_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Tramonto_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include

nodist_execinclude_HEADERS = Tramonto_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

