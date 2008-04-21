# @HEADER
# Tramonto: A molecular theory code for structured and uniform fluids
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
# @HEADER

# NOTE: We also need to include the .h file as a "source" file, because this
# list is also used to package the distribution.

COMMON_H = \
	${srcdir}/Tramonto_ConfigDefs.h \
        ${srcdir}/dft_ldbal.h \
        ${srcdir}/include_global/dft_globals.h \
	${srcdir}/dft_globals_const.h \
        ${srcdir}/rf_allo.h \
        ${srcdir}/loca_const.h \
        ${srcdir}/loca_util_const.h \
        ${srcdir}/dft_BasicLinProbMgr.hpp \
        ${srcdir}/dft_PolyLinProbMgr.hpp \
        ${srcdir}/dft_HardSphereLinProbMgr.hpp \
        ${srcdir}/dft_PolyA11_Epetra_Operator.hpp \
        ${srcdir}/dft_PolyA22_Epetra_Operator.hpp \
        ${srcdir}/dft_PolyA11_Coulomb_Epetra_Operator.hpp \
        ${srcdir}/dft_PolyA22_Coulomb_Epetra_Operator.hpp \
        ${srcdir}/dft_A22Matrix_Epetra_Operator.hpp \
        ${srcdir}/dft_HardSphereA11_Epetra_Operator.hpp \
        ${srcdir}/dft_HardSphereA22_Epetra_Operator.hpp \
        ${srcdir}/dft_Schur_Epetra_Operator.hpp \
        ${srcdir}/dft_basic_lin_prob_mgr_wrapper.h \
        ${srcdir}/dft_poly_lin_prob_mgr_wrapper.h \
        ${srcdir}/dft_hardsphere_lin_prob_mgr_wrapper.h \
        ${srcdir}/dft_parameterlist_wrapper.h \
        ${srcdir}/dft_direct_solver_const.h \
        ${srcdir}/dft_input.h \
        ${srcdir}/dft_main.h \
	${srcdir}/dft_maindriver.h \
        ${srcdir}/dft_mesh.h \
        ${srcdir}/dft_mesh_surfaces.h \
        ${srcdir}/dft_mesh_lib.h \
        ${srcdir}/dft_newton.h \
        ${srcdir}/dft_picard.h \
        ${srcdir}/dft_linsolver_setup.h \
        ${srcdir}/dft_ldbal.h \
        ${srcdir}/dft_fill_pde.h \
        ${srcdir}/dft_fill_pde_elec.h \
        ${srcdir}/dft_fill_pde_diffusion.h \
        ${srcdir}/dft_fill_coarse.h \
        ${srcdir}/dft_fill_EL.h \
        ${srcdir}/dft_fill_CMS.h \
        ${srcdir}/dft_fill_wjdc.h \
        ${srcdir}/dft_fill_CLSchain.h \
        ${srcdir}/dft_fill_CLSmf.h \
        ${srcdir}/dft_physics_FMT1.h \
        ${srcdir}/dft_physics_FMT2.h \
        ${srcdir}/dft_physics_FMT3.h \
        ${srcdir}/dft_defs_stencil.h \
        ${srcdir}/dft_defs_unknowns.h \
        ${srcdir}/dft_quadrature.h \
        ${srcdir}/dft_stencil.h \
        ${srcdir}/dft_stencil_delta_Bond.h \
        ${srcdir}/dft_stencil_delta_R.h \
        ${srcdir}/dft_stencil_theta_CrCMS.h \
        ${srcdir}/dft_stencil_theta_R.h \
        ${srcdir}/dft_stencil_theta_RPMmsa.h \
        ${srcdir}/dft_stencil_theta_Sigma.h \
        ${srcdir}/dft_stencil_theta_uattr.h \
        ${srcdir}/dft_thermo.h \
        ${srcdir}/dft_thermo_hs.h \
        ${srcdir}/dft_thermo_att.h \
        ${srcdir}/dft_thermo_cms.h \
        ${srcdir}/dft_thermo_elec_MSA.h \
        ${srcdir}/dft_thermo_id.h \
        ${srcdir}/dft_thermo_wtc.h \
        ${srcdir}/dft_thermo_wjdc.h \
        ${srcdir}/dft_out_main.h \
        ${srcdir}/dft_out_ads.h \
        ${srcdir}/dft_out_utils.h \
        ${srcdir}/dft_energy.h \
        ${srcdir}/dft_energy_id.h \
        ${srcdir}/dft_energy_hs.h \
        ${srcdir}/dft_energy_att.h \
        ${srcdir}/dft_energy_wtc.h \
        ${srcdir}/dft_energy_wjdc.h \
        ${srcdir}/dft_energy_elec.h \
        ${srcdir}/dft_energy_vext.h \
        ${srcdir}/dft_energy_mu.h \
        ${srcdir}/dft_energy_cms.h \
        ${srcdir}/dft_energy_CLSintsten.h \
        ${srcdir}/dft_out_force.h \
        ${srcdir}/dft_out_flux.h \
        ${srcdir}/dft_out_profiles.h \
        ${srcdir}/dft_vext.h \
        ${srcdir}/dft_guess.h \
        ${srcdir}/dft_guess_CHEMPOT.h \
        ${srcdir}/dft_guess_EL_DENSITY.h \
        ${srcdir}/dft_guess_HSFMT.h \
        ${srcdir}/dft_guess_MFATT.h \
        ${srcdir}/dft_guess_WTC.h \
        ${srcdir}/dft_guess_restart.h \
        ${srcdir}/dft_guess_CMS.h \
        ${srcdir}/dft_guess_WJDC.h \
        ${srcdir}/dft_guess_ELECPOT.h \
        ${srcdir}/dft_comm.h \
        ${srcdir}/dft_continuation.h \
        ${srcdir}/dft_utils.h \
        ${srcdir}/dft_fill_control.h \
        ${srcdir}/dft_fill_main.h \
        ${srcdir}/dft_fill_FMT.h \
        ${srcdir}/dft_fill_mf.h \
        ${srcdir}/dft_fill_wtc.h \
        ${srcdir}/dft_vext_coul.h \
        ${srcdir}/dft_uww.h \
        ${srcdir}/dft_images.h \
        ${srcdir}/dft_potential_params.h \
        ${srcdir}/dft_pairPot_LJ12_6.h \
        ${srcdir}/dft_pairPot_COULOMB.h \
        ${srcdir}/dft_pairPot_Yukawa.h \
        ${srcdir}/dft_switch_pairPot.h \
        ${srcdir}/dft_switch_vext1D.h \
        ${srcdir}/dft_switch_FMT.h \
        ${srcdir}/dft_switch_stencil.h \
        ${srcdir}/dft_vext_1D.h \
        ${srcdir}/dft_vext_integrated.h 

COMMON_C = \
        ${srcdir}/dft_input.c \
        ${srcdir}/dft_main.c \
        ${srcdir}/dft_mesh.c \
        ${srcdir}/dft_mesh_surfaces.c \
        ${srcdir}/dft_mesh_lib.c \
        ${srcdir}/dft_newton.c \
        ${srcdir}/dft_picard.c \
        ${srcdir}/dft_linsolver_setup.c \
        ${srcdir}/dft_ldbal.c \
        ${srcdir}/dft_fill_pde.c \
        ${srcdir}/dft_fill_pde_elec.c \
        ${srcdir}/dft_fill_pde_diffusion.c \
        ${srcdir}/dft_fill_coarse.c \
        ${srcdir}/dft_fill_EL.c \
        ${srcdir}/dft_fill_CMS.c \
        ${srcdir}/dft_fill_wjdc.c \
        ${srcdir}/dft_fill_CLSchain.c \
        ${srcdir}/dft_fill_CLSmf.c \
        ${srcdir}/dft_physics_FMT1.c \
        ${srcdir}/dft_physics_FMT2.c \
        ${srcdir}/dft_physics_FMT3.c \
        ${srcdir}/dft_defs_stencil.c \
        ${srcdir}/dft_defs_unknowns.c \
        ${srcdir}/dft_quadrature.c \
        ${srcdir}/dft_stencil.c \
        ${srcdir}/dft_stencil_delta_Bond.c \
        ${srcdir}/dft_stencil_delta_R.c \
        ${srcdir}/dft_stencil_theta_CrCMS.c \
        ${srcdir}/dft_stencil_theta_R.c \
        ${srcdir}/dft_stencil_theta_RPMmsa.c \
        ${srcdir}/dft_stencil_theta_Sigma.c \
        ${srcdir}/dft_stencil_theta_uattr.c \
        ${srcdir}/dft_thermo.c \
        ${srcdir}/dft_thermo_hs.c \
        ${srcdir}/dft_thermo_att.c \
        ${srcdir}/dft_thermo_cms.c \
        ${srcdir}/dft_thermo_elec_MSA.c \
        ${srcdir}/dft_thermo_id.c \
        ${srcdir}/dft_thermo_wtc.c \
        ${srcdir}/dft_thermo_wjdc.c \
        ${srcdir}/dft_out_main.c \
        ${srcdir}/dft_out_ads.c \
        ${srcdir}/dft_out_utils.c \
        ${srcdir}/dft_energy.c \
        ${srcdir}/dft_energy_id.c \
        ${srcdir}/dft_energy_hs.c \
        ${srcdir}/dft_energy_att.c \
        ${srcdir}/dft_energy_elec.c \
        ${srcdir}/dft_energy_wtc.c \
        ${srcdir}/dft_energy_wjdc.c \
        ${srcdir}/dft_energy_vext.c \
        ${srcdir}/dft_energy_mu.c \
        ${srcdir}/dft_energy_cms.c \
        ${srcdir}/dft_energy_CLSintsten.c \
        ${srcdir}/dft_out_force.c \
        ${srcdir}/dft_out_flux.c \
        ${srcdir}/dft_out_profiles.c \
        ${srcdir}/dft_vext.c \
        ${srcdir}/dft_guess.c \
        ${srcdir}/dft_guess_CHEMPOT.c \
        ${srcdir}/dft_guess_EL_DENSITY.c \
        ${srcdir}/dft_guess_HSFMT.c \
        ${srcdir}/dft_guess_MFATT.c \
        ${srcdir}/dft_guess_WTC.c \
        ${srcdir}/dft_guess_restart.c \
        ${srcdir}/dft_guess_CMS.c \
        ${srcdir}/dft_guess_WJDC.c \
        ${srcdir}/dft_guess_ELECPOT.c \
        ${srcdir}/dft_comm.c \
        ${srcdir}/dft_continuation.c \
        ${srcdir}/dft_utils.c \
        ${srcdir}/dft_fill_control.c \
        ${srcdir}/dft_fill_main.c \
        ${srcdir}/dft_fill_FMT.c \
        ${srcdir}/dft_fill_mf.c \
        ${srcdir}/dft_fill_wtc.c \
        ${srcdir}/dft_vext_coul.c \
        ${srcdir}/dft_uww.c \
        ${srcdir}/dft_images.c \
        ${srcdir}/dft_potential_params.c \
        ${srcdir}/dft_pairPot_LJ12_6.c \
        ${srcdir}/dft_pairPot_COULOMB.c \
        ${srcdir}/dft_pairPot_Yukawa.c \
        ${srcdir}/dft_switch_pairPot.c \
        ${srcdir}/dft_switch_vext1D.c \
        ${srcdir}/dft_switch_FMT.c \
        ${srcdir}/dft_switch_stencil.c \
        ${srcdir}/dft_vext_1D.c \
        ${srcdir}/dft_vext_integrated.c \
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
        ${srcdir}/dft_PolyA11_Coulomb_Epetra_Operator.cpp \
        ${srcdir}/dft_PolyA22_Coulomb_Epetra_Operator.cpp \
        ${srcdir}/dft_A22Matrix_Epetra_Operator.cpp \
        ${srcdir}/dft_HardSphereA11_Epetra_Operator.cpp \
        ${srcdir}/dft_HardSphereA22_Epetra_Operator.cpp \
        ${srcdir}/dft_Schur_Epetra_Operator.cpp \
        ${srcdir}/dft_basic_lin_prob_mgr_wrapper.cpp \
        ${srcdir}/dft_poly_lin_prob_mgr_wrapper.cpp \
        ${srcdir}/dft_hardsphere_lin_prob_mgr_wrapper.cpp \
        ${srcdir}/dft_parameterlist_wrapper.cpp

NOXLOCA_CXX = \
        ${srcdir}/NOXLOCA_Solver.cpp \
        ${srcdir}/NOXLOCA_Tramonto_Group.cpp \
        ${srcdir}/NOXLOCA_Tramonto_Vector.cpp \
        ${srcdir}/NOXLOCA_Tramonto_PTGroup.cpp \
        ${srcdir}/NOXLOCA_Tramonto_PTVector.cpp

NOXLOCA_H = \
        ${srcdir}/NOXLOCA_Tramonto_Group.hpp \
        ${srcdir}/NOXLOCA_Tramonto_Vector.hpp \
        ${srcdir}/NOXLOCA_Tramonto_PTGroup.hpp \
        ${srcdir}/NOXLOCA_Tramonto_PTVector.hpp

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

include $(top_builddir)/Makefile.export.tramonto

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(TRAMONTO_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(TRAMONTO_LIBS))
else
EXPORT_INCLUDES = $(TRAMONTO_INCLUDES)
EXPORT_LIBS = $(TRAMONTO_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)
#dft_CXXFLAGS = -I$(srcdir)

dft_LDADD = $(EXPORT_LIBS)

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

