
# @HEADER
# ********************************************************************
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
#  Questions? Contact Laura Frink (ljfrink@sandia.gov)
#  ************************************************************************
# @HEADER

## #######################################################################
## Options to automake (rarely used - don't worry about it)
## #######################################################################
AUTOMAKE_OPTIONS = foreign

## #######################################################################
## Aclocal command (rarely used - don't worry about it)
## #######################################################################
ACLOCAL	= aclocal -I config

#
#       I believe that by switching to AUX_DIR(../../config) one 
#       could get rid of these.
#
## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
EXTRA_DIST = config/strip_dup_incl_paths.pl config/strip_dup_libs.pl \
	test/utilities/report test/utilities/runbuild \
	test/utilities/runconfigure test/utilities/runmake \
	test/utilities/runtests test/harness/runharness \
	test/harness/config test/harness/invoke-configures/linux-beowulf \
	test/harness/invoke-configures/linux-hrothgar sampleConfigs/cygwin \
	sampleConfigs/cygwin2 sampleConfigs/irix64 sampleConfigs/liberty \
	sampleConfigs/linux sampleConfigs/linux2 sampleConfigs/macOsX \
	sampleConfigs/macOsX2 sampleConfigs/reddish sampleConfigs/sunos \
	sampleConfigs/tbird Examples/BC10_FMT1_WCA_ELEC_3D/dft_input.dat \
	COPYRIGHT LICENSE RELEASE_NOTES \
	Examples/BC10_FMT1_WCA_ELEC_3D/dft_surfaces.dat \
	Examples/BC10_FMT1_WCA_ELEC_3D/README \
	Examples/BC11_FMT1_WCA_ELEC_3D/dft_input.dat \
	Examples/BC11_FMT1_WCA_ELEC_3D/dft_surfaces.dat \
	Examples/BC11_FMT1_WCA_ELEC_3D/README \
	Examples/BC12_FMT1_WCA_ELEC_3D/dft_input.dat \
	Examples/BC12_FMT1_WCA_ELEC_3D/dft_surfaces.dat \
	Examples/BC12_FMT1_WCA_ELEC_3D/README \
	Examples/BC1_FMT1_WCA_ELEC_2D/dft_input.dat \
	Examples/BC1_FMT1_WCA_ELEC_2D/dft_surfaces.dat \
	Examples/BC1_FMT1_WCA_ELEC_2D/README \
	Examples/BC2_FMT1_WCA_ELEC_2D/dft_input.dat \
	Examples/BC2_FMT1_WCA_ELEC_2D/dft_surfaces.dat \
	Examples/BC2_FMT1_WCA_ELEC_2D/README \
	Examples/BC3_FMT1_WCA_ELEC_2D/dft_input.dat \
	Examples/BC3_FMT1_WCA_ELEC_2D/dft_surfaces.dat \
	Examples/BC3_FMT1_WCA_ELEC_2D/README \
	Examples/BC4_FMT1_WCA_ELEC_2D/dft_input.dat \
	Examples/BC4_FMT1_WCA_ELEC_2D/dft_surfaces.dat \
	Examples/BC4_FMT1_WCA_ELEC_2D/README \
	Examples/BC5_FMT1_WCA_ELEC_2D/dft_input.dat \
	Examples/BC5_FMT1_WCA_ELEC_2D/dft_surfaces.dat \
	Examples/BC5_FMT1_WCA_ELEC_2D/README \
	Examples/BC6_FMT1_WCA_ELEC_2D/dft_input.dat \
	Examples/BC6_FMT1_WCA_ELEC_2D/dft_surfaces.dat \
	Examples/BC6_FMT1_WCA_ELEC_2D/README \
	Examples/BC7_FMT1_WCA_ELEC_3D/dft_input.dat \
	Examples/BC7_FMT1_WCA_ELEC_3D/dft_surfaces.dat \
	Examples/BC7_FMT1_WCA_ELEC_3D/README \
	Examples/BC8_FMT1_WCA_ELEC_3D/dft_input.dat \
	Examples/BC8_FMT1_WCA_ELEC_3D/dft_surfaces.dat \
	Examples/BC8_FMT1_WCA_ELEC_3D/README \
	Examples/BC9_FMT1_WCA_ELEC_3D/dft_input.dat \
	Examples/BC9_FMT1_WCA_ELEC_3D/dft_surfaces.dat \
	Examples/BC9_FMT1_WCA_ELEC_3D/README \
	Examples/CONT1_FMT1_WCA_1D/dft_input.dat \
	Examples/CONT1_FMT1_WCA_1D/dft_surfaces.dat \
	Examples/CONT1_FMT1_WCA_1D/README \
	Examples/CONT2_FMT1_WCA_1D/dft_input.dat \
	Examples/CONT2_FMT1_WCA_1D/dft_surfaces.dat \
	Examples/CONT2_FMT1_WCA_1D/README \
	Examples/CONT3_FMT1_WCA_1D/dft_dens.dat \
	Examples/CONT3_FMT1_WCA_1D/dft_input.dat \
	Examples/CONT3_FMT1_WCA_1D/dft_surfaces.dat \
	Examples/CONT3_FMT1_WCA_1D/README \
	Examples/CONT4_FMT1_WCA_1D/dft_dens2.dat \
	Examples/CONT4_FMT1_WCA_1D/dft_dens.dat \
	Examples/CONT4_FMT1_WCA_1D/dft_input.dat \
	Examples/CONT4_FMT1_WCA_1D/dft_surfaces.dat \
	Examples/CONT4_FMT1_WCA_1D/README \
	Examples/HS1_FMT1_1D/dft_input.dat \
	Examples/HS1_FMT1_1D/dft_surfaces.dat \
	Examples/HS1_FMT1_1D/README \
	Examples/HS2_FMT1_1D/dft_input.dat \
	Examples/HS2_FMT1_1D/dft_surfaces.dat \
	Examples/HS2_FMT1_1D/README \
	Examples/HS3_FMT1_1D/dft_input.dat \
	Examples/HS3_FMT1_1D/dft_surfaces.dat \
	Examples/HS3_FMT1_1D/README \
	Examples/HS4_FMT1_1D/dft_input.dat \
	Examples/HS4_FMT1_1D/dft_surfaces.dat \
	Examples/HS4_FMT1_1D/README \
	Examples/HS5_FMT1_1D/dft_input.dat \
	Examples/HS5_FMT1_1D/dft_surfaces.dat \
	Examples/HS5_FMT1_1D/README \
	Examples/HS6_FMT1_1D/dft_input.dat \
	Examples/HS6_FMT1_1D/dft_surfaces.dat \
	Examples/HS6_FMT1_1D/README \
	Examples/HS7_FMT1_1D/dft_input.dat \
	Examples/HS7_FMT1_1D/dft_surfaces.dat \
	Examples/HS7_FMT1_1D/README \
	Examples/LJ1_FMT1_WCA_1D/dft_input.dat \
	Examples/LJ1_FMT1_WCA_1D/dft_surfaces.dat \
	Examples/LJ1_FMT1_WCA_1D/README \
	Examples/LJ2_FMT1_WCA_1D/dft_input.dat \
	Examples/LJ2_FMT1_WCA_1D/dft_surfaces.dat \
	Examples/LJ2_FMT1_WCA_1D/README \
	Examples/LJ3_FMT1_WCA_1D/dft_input.dat \
	Examples/LJ3_FMT1_WCA_1D/dft_surfaces.dat \
	Examples/LJ3_FMT1_WCA_1D/README \
	Examples/LJ4_FMT1_WCA_1D/dft_input.dat \
	Examples/LJ4_FMT1_WCA_1D/dft_surfaces.dat \
	Examples/LJ4_FMT1_WCA_1D/README \
	Examples/MIX1_FMT1_1D/dft_input.dat \
	Examples/MIX1_FMT1_1D/dft_surfaces.dat \
	Examples/MIX1_FMT1_1D/README \
	Examples/MIX2_FMT1_WCA_1D/dft_input.dat \
	Examples/MIX2_FMT1_WCA_1D/dft_surfaces.dat \
	Examples/MIX2_FMT1_WCA_1D/README \
	Examples/MIX3_ELEC_1D/dft_input.dat \
	Examples/MIX3_ELEC_1D/dft_surfaces.dat \
	Examples/MIX3_ELEC_1D/README \
	Examples/MIX4_FMT1_ELEC_1D/dft_input.dat \
	Examples/MIX4_FMT1_ELEC_1D/dft_surfaces.dat \
	Examples/MIX4_FMT1_ELEC_1D/README \
	Examples/MIX5_FMT1_ELEC_1D/dft_input.dat \
	Examples/MIX5_FMT1_ELEC_1D/dft_surfaces.dat \
	Examples/MIX5_FMT1_ELEC_1D/README \
	Examples/MIX6_FMT1_ELEC_1D/dft_input.dat \
	Examples/MIX6_FMT1_ELEC_1D/dft_surfaces.dat \
	Examples/MIX6_FMT1_ELEC_1D/README \
	Examples/MIX7_FMT1_WCA_ELEC_1D/dft_input.dat \
	Examples/MIX7_FMT1_WCA_ELEC_1D/dft_surfaces.dat \
	Examples/MIX7_FMT1_WCA_ELEC_1D/README \
	Examples/MIX8_FMT1_WCA_ELEC_1D/dft_input.dat \
	Examples/MIX8_FMT1_WCA_ELEC_1D/dft_surfaces.dat \
	Examples/MIX8_FMT1_WCA_ELEC_1D/README \
	Examples/MIX9_FMT1_WCA_ELEC_1D/dft_input.dat \
	Examples/MIX9_FMT1_WCA_ELEC_1D/dft_surfaces.dat \
	Examples/MIX9_FMT1_WCA_ELEC_1D/README \
	Examples/POLY1_CMS_1D/cr8.2.8_0.711eff1 \
	Examples/POLY1_CMS_1D/dft_input.dat \
	Examples/POLY1_CMS_1D/dft_surfaces.dat \
	Examples/POLY1_CMS_1D/lin_18 \
	Examples/POLY1_CMS_1D/README \
	Examples/POLY2_CMS_1D/cr8.2.8_0.711eff1 \
	Examples/POLY2_CMS_1D/dft_input.dat \
	Examples/POLY2_CMS_1D/dft_surfaces.dat \
	Examples/POLY2_CMS_1D/lin_18 \
	Examples/POLY2_CMS_1D/README \
	Examples/POLY3_CMS_1D/cr8.2.8_0.375_1_0.375 \
	Examples/POLY3_CMS_1D/dft_dens.dat \
	Examples/POLY3_CMS_1D/dft_dens.datg \
	Examples/POLY3_CMS_1D/dft_input.dat \
	Examples/POLY3_CMS_1D/dft_surfaces.dat \
	Examples/POLY3_CMS_1D/lin_18_1 \
	Examples/POLY3_CMS_1D/README \
	Examples/POLY4_CMS_1D/cr8.2.8_0.375_1_0.375 \
	Examples/POLY4_CMS_1D/dft_dens.dat \
	Examples/POLY4_CMS_1D/dft_dens.datg \
	Examples/POLY4_CMS_1D/dft_input.dat \
	Examples/POLY4_CMS_1D/dft_surfaces.dat \
	Examples/POLY4_CMS_1D/lin_18_1 \
	Examples/POLY4_CMS_1D/README \
	Examples/POLY5_WTC_FMT1_1D/dft_input.dat \
	Examples/POLY5_WTC_FMT1_1D/dft_surfaces.dat \
	Examples/POLY5_WTC_FMT1_1D/lin_3mer \
	Examples/POLY5_WTC_FMT1_1D/README \
	Examples/POLY6_WTC_FMT1_1D/dens_1Dpoly6_chapmanisaft.dat \
	Examples/POLY6_WTC_FMT1_1D/dft_input.dat \
	Examples/POLY6_WTC_FMT1_1D/dft_surfaces.dat \
	Examples/POLY6_WTC_FMT1_1D/lin_5mer \
	Examples/POLY6_WTC_FMT1_1D/README \
	Examples/POLY7_WTC_FMT1_1D/dft_input.dat \
	Examples/POLY7_WTC_FMT1_1D/dft_surfaces.dat \
	Examples/POLY7_WTC_FMT1_1D/lin_5mer \
	Examples/POLY7_WTC_FMT1_1D/README \
	Examples/POLY8_WTC_FMT1_ELEC_1D/dft_input.dat \
	Examples/POLY8_WTC_FMT1_ELEC_1D/dft_surfaces.dat \
	Examples/POLY8_WTC_FMT1_ELEC_1D/lin_5mer \
	Examples/POLY8_WTC_FMT1_ELEC_1D/README \
	Examples/POLY9_WTC_FMT1_1D/dft_input.dat \
	Examples/POLY9_WTC_FMT1_1D/dft_surfaces.dat \
	Examples/POLY9_WTC_FMT1_1D/lin_8_lin_4 \
	Examples/POLY9_WTC_FMT1_1D/README

#
#       I believe that by switching to AUX_DIR(../../config) one 
#       could get rid of these.
#
## #######################################################################
## Tools in the auxillary directory 
## #######################################################################
AUX_DIST = config/install-sh config/missing config/mkinstalldirs 
#
#  Again, I hope that AUX_DIR(../../config) eliminates these
#  config/install-sh config/missing config/mkinstalldirs 

## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 autom4te.cache/* \
	configure config.status config.log \
	src/common/config-h.in src/common/stamp-h.in \
	$(AUX_DIST) 

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################
#  BUILD_TESTS, BUILD_EXAMPLES, and BUILD_GUI are defined in configure.ac.

#if BUILD_TESTS
#TEST_SUBDIR=TEST
#endif

if BUILD_EXAMPLES
EXAMPLE_SUBDIR=Examples
endif

#if BUILD_GUI
#GUI_SUBDIR=TramontoGUI
#endif

SUBDIRS = src $(EXAMPLE_SUBDIR)

