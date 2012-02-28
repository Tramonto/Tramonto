include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#include_directories(${MPI_INCLUDE_PATH})
set(common_h
	Tramonto_ConfigDefs.h 
        dft_ldbal.h 
		dft_globals_const.h 
        include_global/dft_globals.h 
        rf_allo.h 
        loca_const.h 
        loca_util_const.h 
        dft_BasicLinProbMgr.hpp 
        dft_PolyLinProbMgr.hpp 
        dft_HardSphereLinProbMgr.hpp 
        dft_PolyA11_Epetra_Operator.hpp 
        dft_PolyA22_Epetra_Operator.hpp 
        dft_PolyA11_Coulomb_Epetra_Operator.hpp 
        dft_PolyA22_Coulomb_Epetra_Operator.hpp 
        dft_A22Matrix_Epetra_Operator.hpp 
        dft_HardSphereA11_Epetra_Operator.hpp 
        dft_HardSphereA22_Epetra_Operator.hpp 
        dft_Schur_Epetra_Operator.hpp 
        dft_basic_lin_prob_mgr_wrapper.h 
        dft_poly_lin_prob_mgr_wrapper.h 
        dft_hardsphere_lin_prob_mgr_wrapper.h 
        dft_parameterlist_wrapper.h 
        dft_direct_solver_const.h 
        dft_input.h 
        dft_setup_params.h
        dft_input_broadcast.h
        dft_chain_architecture.h 
        dft_main.h 
	dft_maindriver.h 
        dft_mesh.h 
        dft_surfaces.h
        dft_surfaces_planar.h
        dft_surfaces_block.h
        dft_surfaces_sphere.h
        dft_surfaces_cylinder2D.h
        dft_surfaces_cylinder3D.h
        dft_surfaces_1elem.h
        dft_surfaces_atoms.h
        dft_surfaces_cylindricalPore2D.h
        dft_surfaces_sphericalCavity3D.h
        dft_surfaces_slitPore2D.h
        dft_surfaces_cylindricalPore3D.h
        dft_surfacesAdjust_angleCutout.h
        dft_surfacesAdjust_roughness.h
        dft_surfacesAdjust_periodic.h
        dft_surfacesAdjust_linear.h
        dft_mesh_lib.h 
        dft_newton.h 
        dft_picard.h 
        dft_linsolver_setup.h 
        dft_ldbal.h 
        dft_fill_pde.h 
        dft_fill_pde_elec.h 
        dft_fill_pde_diffusion.h 
        dft_fill_coarse.h 
        dft_fill_EL.h 
        dft_fill_CMS.h 
        dft_fill_wjdc.h 
        dft_fill_CLSchain.h 
        dft_fill_CLSmf.h 
        dft_physics_FMT1.h 
        dft_physics_FMT2.h 
        dft_physics_FMT3.h 
        dft_physics_FMT4.h 
        dft_defs_stencil.h 
        dft_defs_unknowns.h 
        dft_quadrature.h 
        dft_stencil.h 
        dft_stencil_delta_Bond.h 
        dft_stencil_delta_R.h 
        dft_stencil_theta_CrCMS.h 
        dft_stencil_theta_R.h 
        dft_stencil_theta_RPMmsa.h 
        dft_stencil_theta_Sigma.h 
        dft_stencil_theta_uattr.h 
        dft_stencil_theta_GENmsa.h
        dft_thermo.h 
        dft_thermo_hs.h 
        dft_thermo_att.h 
        dft_thermo_cms.h 
        dft_thermo_elec_MSA.h 
        dft_thermo_elec_MSAgeneral.h 
        dft_thermo_id.h 
        dft_thermo_wtc.h 
        dft_thermo_wjdc.h 
        dft_out_main.h 
        dft_out_ads.h 
        dft_out_utils.h 
        dft_energy.h 
        dft_energy_id.h 
        dft_energy_hs.h 
        dft_energy_att.h 
        dft_energy_wtc.h 
        dft_energy_wjdc.h 
        dft_energy_elec.h 
        dft_energy_vext.h 
        dft_energy_mu.h 
        dft_energy_cms.h 
        dft_energy_CLSintsten.h 
        dft_out_force.h 
        dft_out_flux.h 
        dft_out_profiles.h 
        dft_vext.h 
        dft_guess.h 
        dft_guess_CHEMPOT.h 
        dft_guess_EL_DENSITY.h 
        dft_guess_HSFMT.h 
        dft_guess_MFATT.h 
        dft_guess_WTC.h 
        dft_guess_restart.h 
        dft_guess_CMS.h 
		dft_guess_SCF.h 
        dft_guess_WJDC.h 
        dft_guess_ELECPOT.h 
        dft_comm.h 
        dft_continuation.h 
        dft_switch_continuation.h 
        dft_plugin_archived_continue.h 
        dft_plugin_user_continue.h 
        dft_scale_variables.h 
        dft_utils.h 
        dft_fill_control.h 
        dft_fill_main.h 
        dft_fill_FMT.h 
        dft_fill_mf.h 
		dft_fill_SCF.h 
        dft_fill_wtc.h 
        dft_vext_coul.h 
        dft_uww.h 
        dft_images.h 
        dft_potential_params.h 
        dft_pairPot.h 
        dft_pairPot_LJ12_6.h 
        dft_pairPot_LJ12_6_SIGTORCUT.h 
        dft_pairPot_COULOMB.h 
        dft_pairPot_Yukawa.h 
        dft_pairPot_LJ12_6YukawaSum.h
        dft_pairPot_r12YukawaSum.h
        dft_pairPot_r18YukawaSum.h
        dft_pairPot_nvarYukawaSum.h
        dft_pairPot_EXP.h 
		dft_pairPot_SW.h 
        dft_switch_pairPot.h 
        dft_switch_vext1D.h 
        dft_switch_FMT.h 
        dft_switch_stencil.h 
        dft_vext_1D.h 
        dft_vext_integrated.h )

set(common_c
        dft_input.c 
        dft_setup_params.c
        dft_input_broadcast.c
        dft_chain_architecture.c 
        dft_main.c 
        dft_mesh.c 
        dft_surfaces.c
        dft_surfaces_planar.c
        dft_surfaces_block.c
        dft_surfaces_sphere.c
        dft_surfaces_cylinder2D.c
        dft_surfaces_cylinder3D.c
        dft_surfaces_1elem.c
        dft_surfaces_atoms.c
        dft_surfaces_cylindricalPore2D.c
        dft_surfaces_sphericalCavity3D.c
        dft_surfaces_slitPore2D.c
        dft_surfaces_cylindricalPore3D.c
        dft_surfacesAdjust_angleCutout.c
        dft_surfacesAdjust_roughness.c
        dft_surfacesAdjust_periodic.c
        dft_surfacesAdjust_linear.c
        dft_mesh_lib.c 
        dft_newton.c 
        dft_picard.c 
        dft_linsolver_setup.c 
        dft_ldbal.c 
        dft_fill_pde.c 
        dft_fill_pde_elec.c 
        dft_fill_pde_diffusion.c 
        dft_fill_coarse.c 
        dft_fill_EL.c 
        dft_fill_CMS.c 
        dft_fill_wjdc.c 
        dft_fill_CLSchain.c 
        dft_fill_CLSmf.c 
        dft_physics_FMT1.c 
        dft_physics_FMT2.c 
        dft_physics_FMT3.c 
        dft_physics_FMT4.c 
        dft_defs_stencil.c 
        dft_defs_unknowns.c 
        dft_quadrature.c 
        dft_stencil.c 
        dft_stencil_delta_Bond.c 
        dft_stencil_delta_R.c 
        dft_stencil_theta_CrCMS.c 
        dft_stencil_theta_R.c 
        dft_stencil_theta_RPMmsa.c 
        dft_stencil_theta_Sigma.c 
        dft_stencil_theta_uattr.c 
        dft_stencil_theta_GENmsa.c
        dft_thermo.c 
        dft_thermo_hs.c 
        dft_thermo_att.c 
        dft_thermo_cms.c 
        dft_thermo_elec_MSA.c 
        dft_thermo_elec_MSAgeneral.c 
        dft_thermo_id.c 
        dft_thermo_wtc.c 
        dft_thermo_wjdc.c 
        dft_out_main.c 
        dft_out_ads.c 
        dft_out_utils.c 
        dft_energy.c 
        dft_energy_id.c 
        dft_energy_hs.c 
        dft_energy_att.c 
        dft_energy_elec.c 
        dft_energy_wtc.c 
        dft_energy_wjdc.c 
        dft_energy_vext.c 
        dft_energy_mu.c 
        dft_energy_cms.c 
        dft_energy_CLSintsten.c 
        dft_out_force.c 
        dft_out_flux.c 
        dft_out_profiles.c 
        dft_vext.c 
        dft_guess.c 
        dft_guess_CHEMPOT.c 
        dft_guess_EL_DENSITY.c 
        dft_guess_HSFMT.c 
        dft_guess_MFATT.c 
        dft_guess_WTC.c 
        dft_guess_restart.c 
        dft_guess_CMS.c 
		dft_guess_SCF.c 
        dft_guess_WJDC.c 
        dft_guess_ELECPOT.c 
        dft_comm.c 
        dft_continuation.c 
        dft_switch_continuation.c 
        dft_plugin_archived_continue.c 
        dft_plugin_user_continue.c 
        dft_scale_variables.c 
        dft_utils.c 
        dft_fill_control.c 
        dft_fill_main.c 
        dft_fill_FMT.c 
        dft_fill_mf.c 
		dft_fill_SCF.c 
        dft_fill_wtc.c 
        dft_vext_coul.c 
        dft_uww.c 
        dft_images.c 
        dft_potential_params.c 
        dft_pairPot.c 
        dft_pairPot_LJ12_6.c 
        dft_pairPot_LJ12_6_SIGTORCUT.c 
        dft_pairPot_COULOMB.c 
        dft_pairPot_Yukawa.c 
        dft_pairPot_LJ12_6YukawaSum.c
        dft_pairPot_r12YukawaSum.c
        dft_pairPot_r18YukawaSum.c
        dft_pairPot_nvarYukawaSum.c
        dft_pairPot_EXP.c 
		dft_pairPot_SW.c 
        dft_switch_pairPot.c 
        dft_switch_vext1D.c 
        dft_switch_FMT.c 
        dft_switch_stencil.c 
        dft_vext_1D.c 
        dft_vext_integrated.c 
        rf_allo.c 
        loca_lib.c 
        loca_bord.c 
        loca_util.c )

set(common_cxx
        dft_BasicLinProbMgr.cpp 
        dft_PolyLinProbMgr.cpp 
        dft_HardSphereLinProbMgr.cpp 
        dft_PolyA11_Epetra_Operator.cpp 
        dft_PolyA22_Epetra_Operator.cpp 
        dft_PolyA11_Coulomb_Epetra_Operator.cpp 
        dft_PolyA22_Coulomb_Epetra_Operator.cpp 
        dft_A22Matrix_Epetra_Operator.cpp 
        dft_HardSphereA11_Epetra_Operator.cpp 
        dft_HardSphereA22_Epetra_Operator.cpp 
        dft_Schur_Epetra_Operator.cpp 
        dft_basic_lin_prob_mgr_wrapper.cpp 
        dft_poly_lin_prob_mgr_wrapper.cpp 
        dft_hardsphere_lin_prob_mgr_wrapper.cpp 
        dft_parameterlist_wrapper.cpp 
        )
set(noxloca_cxx
        NOXLOCA_Solver.cpp 
        NOXLOCA_Tramonto_Group.cpp 
        NOXLOCA_Tramonto_Vector.cpp 
        NOXLOCA_Tramonto_PTGroup.cpp 
        NOXLOCA_Tramonto_PTVector.cpp )

set(noxloca_headers
        NOXLOCA_Tramonto_Group.hpp 
        NOXLOCA_Tramonto_Vector.hpp 
        NOXLOCA_Tramonto_PTGroup.hpp 
        NOXLOCA_Tramonto_PTVector.hpp)

#Add Trilinos information to the include and link lines
include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
link_directories(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} )
# /Library/Frameworks/QtCore.framework /Library/Frameworks/QtGui.framework)

#list of necessary Trilinos packages:
#ml zoltan amesos loca nox aztecoo ifpack epetraext epetra teuchos muelu optika 
set(Tramonto_Trilinos_LIBRARIES ${ML_LIBRARIES} ${Amesos_LIBRARIES}
${NOX_LIBRARIES} ${AztecOO_LIBRARIES} ${Ifpack_LIBRARIES} ${MueLu_LIBRARIES}
${EpetraExt_LIBRARIES} ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES} )
#${Optika_LIBRARIES} /Library/Frameworks/QtCore.framework/QtCore /Library/Frameworks/QtGui.framework/QtGui)

#remove duplicates from the link line
LIST(REVERSE Tramonto_Trilinos_LIBRARIES )
LIST(REMOVE_DUPLICATES Tramonto_Trilinos_LIBRARIES )
LIST(REVERSE Tramonto_Trilinos_LIBRARIES )

#verify all necessary Trilinos packages are available
LIST(FIND Trilinos_PACKAGE_LIST ML ML_FOUND)
IF(ML_FOUND GREATER -1)
  MESSAGE("Trilinos package ML found.")
ELSE()
  MESSAGE(FATAL_ERROR "Could not find Trilinos package ML!")
ENDIF()

LIST(FIND Trilinos_PACKAGE_LIST Amesos AMESOS_FOUND)
IF(AMESOS_FOUND GREATER -1)
  MESSAGE("Trilinos package Amesos found.")
ELSE()
  MESSAGE(FATAL_ERROR "Could not find Trilinos package Amesos!")
ENDIF()

LIST(FIND Trilinos_PACKAGE_LIST NOX NOX_FOUND)
IF(NOX_FOUND GREATER -1)
  MESSAGE("Trilinos package NOX found.")
ELSE()
  MESSAGE(FATAL_ERROR "Could not find Trilinos package NOX!")
ENDIF()

LIST(FIND Trilinos_PACKAGE_LIST AztecOO AZTECOO_FOUND)
IF(AZTECOO_FOUND GREATER -1)
  MESSAGE("Trilinos package AztecOO found.")
ELSE()
  MESSAGE(FATAL_ERROR "Could not find Trilinos package AztecOO!")
ENDIF()

LIST(FIND Trilinos_PACKAGE_LIST Ifpack IFPACK_FOUND)
IF(IFPACK_FOUND GREATER -1)
  MESSAGE("Trilinos package Ifpack found.")
ELSE()
  MESSAGE(FATAL_ERROR "Could not find Trilinos package Ifpack!")
ENDIF()

LIST(FIND Trilinos_PACKAGE_LIST EpetraExt EPETRAEXT_FOUND)
IF(EPETRAEXT_FOUND GREATER -1)
  MESSAGE("Trilinos package EpetraExt found.")
ELSE()
  MESSAGE(FATAL_ERROR "Could not find Trilinos package EpetraExt!")
ENDIF()

LIST(FIND Trilinos_PACKAGE_LIST Epetra EPETRA_FOUND)
IF(EPETRA_FOUND GREATER -1)
  MESSAGE("Trilinos package Epetra found.")
ELSE()
  MESSAGE(FATAL_ERROR "Could not find Trilinos package Epetra!")
ENDIF()

LIST(FIND Trilinos_PACKAGE_LIST Teuchos TEUCHOS_FOUND)
IF(TEUCHOS_FOUND GREATER -1)
  MESSAGE("Trilinos package Teuchos found.")
ELSE()
  MESSAGE(FATAL_ERROR "Could not find Trilinos package Teuchos!")
ENDIF()

if (BUILD_GUI)
   message("BUILD_GUI defined. Compiling against Optika code.")
  LIST(FIND Trilinos_PACKAGE_LIST Optika OPTIKA_FOUND)
  IF(OPTIKA_FOUND GREATER -1)
    MESSAGE("Trilinos package Optika found.")
  ELSE()
    MESSAGE(FATAL_ERROR "Could not find Trilinos package Optika!")
  ENDIF()
  list(APPEND common_h GUI/dft_GUI.hpp)
  list(APPEND common_h GUI/dft_GUI.h)
  list(APPEND common_cxx GUI/dft_GUI_mesh.cpp)
  list(APPEND common_cxx GUI/dft_GUI_output.cpp)
  list(APPEND common_cxx GUI/dft_GUI_potentialsWF.cpp)
  list(APPEND common_cxx GUI/dft_GUI_numericalMethods.cpp)
  list(APPEND common_cxx GUI/dft_GUI_functionals.cpp)
  list(APPEND common_cxx GUI/dft_GUI_vextType.cpp)
  list(APPEND common_cxx GUI/dft_GUI_toTramonto.cpp)
  list(APPEND common_cxx GUI/dft_GUI_main.cpp)
  list(APPEND common_cxx GUI/dft_GUI_surfaces.cpp)
  list(APPEND common_cxx GUI/dft_GUI_potentialsFF.cpp)
  list(APPEND common_cxx GUI/dft_GUI_statePoint.cpp)
  list(APPEND common_cxx GUI/dft_GUI_Polymer.cpp)
  list(APPEND common_cxx GUI/dft_GUI_continuation.cpp)

  list(APPEND Tramonto_Trilinos_LIBRARIES ${Optika_LIBRARIES} )
  ADD_DEFINITIONS(-DBUILD_GUI=1)
endif(BUILD_GUI)

# Handle selective compilation of Tpetra source
if (USE_TPETRA)
  message("USE_TPETRA defined. Compiling against Tpetra code.")
#  list(APPEND common_h Tpetra_OperatorApplyInverse.hpp)

  list(REMOVE_ITEM common_h dft_BasicLinProbMgr.hpp)
  list(REMOVE_ITEM common_cxx dft_BasicLinProbMgr.cpp)
  list(APPEND common_h dft_TpetraBasicLinProbMgr.hpp)
  list(APPEND common_cxx dft_TpetraBasicLinProbMgr.cpp)

  list(REMOVE_ITEM common_h dft_PolyLinProbMgr.hpp)
  list(REMOVE_ITEM common_cxx dft_PolyLinProbMgr.cpp)
  list(APPEND common_h dft_TpetraPolyLinProbMgr.hpp)
  list(APPEND common_cxx dft_TpetraPolyLinProbMgr.cpp)

  list(REMOVE_ITEM common_h dft_HardSphereLinProbMgr.hpp)
  list(REMOVE_ITEM common_cxx dft_HardSphereLinProbMgr.cpp)
#  list(APPEND common_h dft_TpetraHardSphereLinProbMgr.hpp)
#  list(APPEND common_cxx dft_TpetraHardSphereLinProbMgr.cpp)

  list(REMOVE_ITEM common_cxx  dft_basic_lin_prob_mgr_wrapper.cpp)
  list(APPEND common_cxx dft_tpetra_basic_lin_prob_mgr_wrapper.cpp)

  list(REMOVE_ITEM common_cxx  dft_poly_lin_prob_mgr_wrapper.cpp)
  list(APPEND common_cxx dft_tpetra_poly_lin_prob_mgr_wrapper.cpp)

  list(REMOVE_ITEM common_cxx  dft_hardsphere_lin_prob_mgr_wrapper.cpp)
  list(APPEND common_cxx dft_tpetra_hardsphere_lin_prob_mgr_wrapper.cpp)

#  list(APPEND common_h dft_A22Matrix_Tpetra_Operator.hpp)
#  list(APPEND common_cxx dft_A22Matrix_Tpetra_Operator.cpp)

  list(APPEND common_h dft_PolyA11_Tpetra_Operator.hpp)
  list(APPEND common_cxx dft_PolyA11_Tpetra_Operator.cpp)
  list(APPEND common_h dft_PolyA11_Coulomb_Tpetra_Operator.hpp)
  list(APPEND common_cxx dft_PolyA11_Coulomb_Tpetra_Operator.cpp)

  list(APPEND common_h dft_PolyA22_Tpetra_Operator.hpp)
  list(APPEND common_cxx dft_PolyA22_Tpetra_Operator.cpp)
  list(APPEND common_h dft_PolyA22_Coulomb_Tpetra_Operator.hpp)
  list(APPEND common_cxx dft_PolyA22_Coulomb_Tpetra_Operator.cpp)

  list(APPEND common_h dft_Schur_Tpetra_Operator.hpp)
  list(APPEND common_cxx dft_Schur_Tpetra_Operator.cpp)

#  list(APPEND common_h dft_HardSphereA11_Tpetra_Operator.hpp)
#  list(APPEND common_cxx dft_HardSphereA11_Tpetra_Operator.cpp)
#  list(APPEND common_h dft_HardSphereA22_Tpetra_Operator.hpp)
#  list(APPEND common_cxx dft_HardSphereA22_Tpetra_Operator.cpp)

   list(APPEND Tramonto_Trilinos_LIBRARIES ${Tpetra_LIBRARIES} ${Ifpack2_LIBRARIES} )

   LIST(FIND Trilinos_PACKAGE_LIST Tpetra TPETRA_FOUND)
   IF(TPETRA_FOUND GREATER -1)
     MESSAGE("Trilinos package Tpetra found.")
   ELSE()
     MESSAGE(FATAL_ERROR "Could not find Trilinos package Tpetra!")
   ENDIF()

   LIST(FIND Trilinos_PACKAGE_LIST Ifpack2 IFPACK2_FOUND)
   IF(IFPACK2_FOUND GREATER -1)
     MESSAGE("Trilinos package Ifpack2 found.")
   ELSE()
     MESSAGE(FATAL_ERROR "Could not find Trilinos package Ifpack2!")
   ENDIF()

#  Select precision of linear solver
   if (LINSOLVE_PREC_FLOAT)   
     message("LINSOLVE_PREC_FLOAT defined. Using float precision for linear solver code.")
     ADD_DEFINITIONS(-DLINSOLVE_PREC=0)
   elseif (LINSOLVE_PREC_DOUBLE)   
     message("LINSOLVE_PREC_DOUBLE defined. Using double precision for linear solver code.")
     ADD_DEFINITIONS(-DLINSOLVE_PREC=1)
   elseif (LINSOLVE_PREC_QUAD_DOUBLE)   
     message("LINSOLVE_PREC_QUAD_DOUBLE defined. Using quad double precision for linear solver code.")
     ADD_DEFINITIONS(-DLINSOLVE_PREC=2)
   elseif (LINSOLVE_PREC_DOUBLE_DOUBLE)   
     message("LINSOLVE_PREC_DOUBLE_DOUBLE defined. Using double double precision for linear solver code.")
     ADD_DEFINITIONS(-DLINSOLVE_PREC=3)
   else () 
     message("No precision specified. Using double precision for linear solver code.")
     ADD_DEFINITIONS(-DLINSOLVE_PREC=1)   
   endif (LINSOLVE_PREC_FLOAT)

endif(USE_TPETRA)

# add dft library
add_library(dft ${common_h} ${noxloca_h} ${common_c} ${common_cxx} ${noxloca_cxx} )

#add dft_exec executable
#FIX Epetra_TPL_LIBRARIES should be changed to Trilinos_TPL_LIBRARIES when
#necessary fix is applied to Trilinos
add_executable(dft_exec dft_maindriver.c)
target_link_libraries(dft_exec dft ${Tramonto_Trilinos_LIBRARIES}
${Trilinos_TPL_LIBRARIES} )

INSTALL(TARGETS dft DESTINATION lib)
INSTALL(TARGETS dft_exec DESTINATION bin)
INCLUDE(CPack)

