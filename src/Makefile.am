## this file is used when automake is run in the parent directory. 
## There should be no need to process it by hand.
# NEVER edit the Makefile.in produced by automake or the Makefile produced
# by configure, as it will get overwritten when automake or configure
# are run -- edit Makefile.am and re-run "automake --foreign"

#Define the path to the source code
srcdir=@srcdir@
VPATH=@srcdir@

# the subdirectories that need to be made first, with '.' at the end to make
# this directory last
SUBDIRS= .

# The result of makes are programs  and libraries.  In this directory,
# we produce a program, dft, to be installed in "bindir" when make install
# is run
bin_PROGRAMS=dft

lib_LIBRARIES=libdft.a

# NOTE: We also need to include the .h file as a "source" file, because this
# list is also used to package the distribution.
common_SOURCES= \
	dft_input.c \
	dft_main.c \
	dft_mesh.c \
	dft_mesh_surfaces.c \
	dft_mesh_lib.c \
	dft_newton.c \
	dft_ldbal.c \
	dft_fill_msr.c \
	dft_fill_pde.c \
	dft_stencil.c \
	dft_thermo.c \
	dft_out_main.c \
	dft_out_ads.c \
	dft_out_energy.c \
	dft_out_force.c \
	dft_out_flux.c \
	dft_out_profiles.c \
	dft_vext.c \
	dft_guess.c \
	dft_comm.c \
	dft_continuation.c \
	dft_fill_shared.c \
	dft_fill_main.c \
	dft_fillp.c \
	dft_fill_rosen.c \
	dft_fill_rosen_rb.c \
	dft_fill_mf.c \
	dft_vext_coul.c \
	dft_uww.c \
	dft_images.c \
	dft_potentials.c \
	rf_allo.c \
	loca_con_lib.c \
	loca_con_bord.c \
	loca_util.c \
	dft_globals.h \
	dft_globals_const.h \
	dft_ldbal.h
	rf_allo.h \
	loca_con_const.h \ 
	loca_util_const.h \
	az_aztec.h \
	az_aztec_defs.h \
	$(srcdir)/include/dft_fill_msr.h \
	$(srcdir)/include/dft_fill_shared.h  \
	$(srcdir)/include/dft_main.h \
	$(srcdir)/include/dft_out_ads.h \
	$(srcdir)/include/dft_out_profiles.h \
	$(srcdir)/include/dft_thermo.h \
	$(srcdir)/include/dft_comm.h \
	$(srcdir)/include/dft_fill_pde.h \     
	$(srcdir)/include/dft_guess.h \
	$(srcdir)/include/dft_mesh.h \
	$(srcdir)/include/dft_out_energy.h \
	$(srcdir)/include/dft_postavs.h \
	$(srcdir)/include/dft_uww.h \
	$(srcdir)/include/dft_continuation.h \
	$(srcdir)/include/dft_fillp.h \
	$(srcdir)/include/dft_images.h \
	$(srcdir)/include/dft_mesh_lib.h \
	$(srcdir)/include/dft_out_flux.h	\
	$(srcdir)/include/dft_potentials.h  \
	$(srcdir)/include/dft_vext.h \
	$(srcdir)/include/dft_fill_main.h \
	$(srcdir)/include/dft_fill_rosen.h \
	$(srcdir)/include/dft_input.h \
	$(srcdir)/include/dft_mesh_surfaces.h \
	$(srcdir)/include/dft_out_force.h \
	$(srcdir)/include/dft_quadratue.h	 \
	$(srcdir)/include/loca_con_bord.h \
	$(srcdir)/include/dft_fill_mf.h \
	$(srcdir)/include/dft_fill_rosen_rb.h \
	$(srcdir)/include/dft_ldbal.h \
	$(srcdir)/include/dft_newton.h \
	$(srcdir)/include/dft_out_main.h \
	$(srcdir)/include/dft_stencil.h \
	$(srcdir)/include/loca_util.h

dft_SOURCES= dft_maindriver.c $(common_SOURCES) 
	
libdft_a_SOURCES=$(common_SOURCES)
EXTRA_DIST = $(srcdir)/az_aztec_defs.h $(srcdir)/az_aztec.h $(srcdir)/loca_con_const.h $(srcdir)/loca_util_const.h $(srcdir)/rf_con_const.h $(srcdir)/rf_allo.h $(srcdir)/include
#
# let the configure script figure out these parameters
#

CC = @CC@
CXX = @CXX@
MPICXX = @MPICXX@
CFLAGS = @CFLAGS@ @TRAMONTO_CFLAGS@ -I.
DFT_FILL_PDE_CFLAGS = @DFT_FILL_PDE_CFLAGS@
TRAMONTO_LDFLAGS = -L../../lib/@BUILD_MACHINE@ @TRAMONTO_LDFLAGS@ @LDFLAGS@ @LIBS@ @FTOCLIBS@
EXE = dft_@BUILD_MACHINE@

# Targets
dft$(EXEEXT): $(dft_OBJECTS) 
	${MPICXX} ${CFLAGS} -o ${EXE} ${dft_OBJECTS} ${TRAMONTO_LDFLAGS} 
	cp ${EXE} ../.
	
libdft.a: ${libdft_a_OBJECTS}
	ar cru libdft.a ${libdft_a_OBJECTS}
	${RANLIB} libdft.a

# This file gets special treatment on pc-linux platforms
# on all other platforms it's compiled like any other 
# c file
dft_fill_pde.o: dft_fill_pde.c
	${CC} ${DFT_FILL_PDE_CFLAGS} -c -o $@ $<
