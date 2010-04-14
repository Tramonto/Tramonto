include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TRILINOS_PATH}/include)
include_directories(${MPI_INCLUDE_PATH})
link_directories(${TRILINOS_PATH}/lib)
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
        dft_chain_architecture.h 
        dft_main.h 
	dft_maindriver.h 
        dft_mesh.h 
        dft_mesh_surfaces.h 
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
        dft_pairPot_LJ12_6.h 
        dft_pairPot_LJ12_6_SIGTORCUT.h 
        dft_pairPot_COULOMB.h 
        dft_pairPot_Yukawa.h 
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
        dft_chain_architecture.c 
        dft_main.c 
        dft_mesh.c 
        dft_mesh_surfaces.c 
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
        dft_pairPot_LJ12_6.c 
        dft_pairPot_LJ12_6_SIGTORCUT.c 
        dft_pairPot_COULOMB.c 
        dft_pairPot_Yukawa.c 
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
        dft_parameterlist_wrapper.cpp )
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

set(TRILINOS_LIBRARIES  ml amesos loca nox aztecoo ifpack epetraext epetra teuchos )

add_library(dft ${common_h} ${noxloca_h} ${common_c} ${common_cxx} ${noxloca_cxx} )
target_link_libraries(dft ${TRILINOS_LIBRARIES} ${LAPACK_LIBRARIES}
${BLAS_LIBRARIES} ${MPI_LIBRARIES} )

add_executable(dft_exec dft_maindriver.c)
target_link_libraries(dft_exec dft ${TRILINOS_LIBRARIES} ${LAPACK_LIBRARIES}
${BLAS_LIBRARIES} ${MPI_LIBRARIES})
