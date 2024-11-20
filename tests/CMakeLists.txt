# ==========================
# A collection of unit tests
# ==========================

include("${CMAKE_SOURCE_DIR}/cmake/macros/macro_setup_test.cmake")
configure_file(globals.cpp.in globals.cpp @ONLY)

# Enable testing functionality
enable_testing()

# torch tests
if (${Torch_FOUND})
  #setup_test(test_eos)
  #setup_test(test_weno)
  #setup_test(test_reconstruct)
  #setup_test(test_riemann)
  #setup_test(test_hydro)
  setup_test(test_sedimentation)
endif()

# athena tests
setup_test(test_mesh)

if (${NVAPOR} EQUAL 3 AND ${NCLOUD} EQUAL 5)
  setup_test(test_jup_thermo)
endif()

if (${NVAPOR} EQUAL 1)
  setup_test(test_earth_thermo)
endif()

setup_test(test_glog)
setup_test(test_absorber)
setup_test(test_read_stellar_flux)
setup_test(test_read_cia_ff)
setup_test(test_read_rayleigh)
setup_test(test_read_yaml)
setup_test(test_helios_ck_read)
setup_test(test_cs_velocity_rotation)
#setup_test(test_rad_unit)

if (PVFMM_OPTION STREQUAL "ENABLE_PVFMM")
  setup_test(test_nbody)
endif()

if (${NVAPOR} EQUAL 0)
  if (NOT HYDROSTATIC OR NOT DEFINED HYDROSTATIC)
    #setup_test(test_convective_adjustment)
  endif()
endif()

if (${NVAPOR} EQUAL 1 AND ${NCLOUD} EQUAL 1 AND ${NPHASE_LEGACY} EQUAL 2)
  #setup_test(test_moist_adiabat)
endif()

if (${NVAPOR} EQUAL 2)
  if (${NCLOUD} EQUAL 4)
    setup_test(test_impl)
    setup_test(test_air_parcel)
    setup_test(test_radiation)
    setup_test(test_microwave_opacity)
    setup_test(test_microphysics)
  elseif (${NCLOUD} EQUAL 5)
    setup_test(test_ammonium_hydrosulfide)
  endif()
endif()

# copy yaml files
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/tests/)
endforeach()
