# define default parameters

list(APPEND CMAKE_PREFIX_PATH "/usr/include/mpich-x84_64" "/usr/lib64/mpich/")

# populate the default values
set_if_empty(NVAPOR 0)

set_if_empty(NCLOUD 0)

set_if_empty(NPHASE_LEGACY 3)

set_if_empty(NCHEMISTRY 0)

set_if_empty(NTRACER 0)

set_if_empty(NTURBULENCE 0)

set_if_empty(NSTATIC 0)

set_if_empty(NINT_PARTICLE_DATA 0)

set_if_empty(NREAL_PARTICLE_DATA 0)

# set_if_empty(RSOLVER lmars)

# set_if_empty(NUMBER_GHOST_CELLS 3)

set_if_empty(EOS "ideal_moist")

if(EOS STREQUAL "shallow_xy")
  set(NVAPOR 0)
  set(NON_BAROTROPIC_EOS 0)
  set_if_empty(RSOLVER "roe_shallow_xy")
elseif(EOS STREQUAL "shallow_yz")
  set(NVAPOR 0)
  set(NON_BAROTROPIC_EOS 0)
  set_if_empty(RSOLVER "roe_shallow_yz")
elseif(EOS STREQUAL "ideal_moist")
  set(NON_BAROTROPIC_EOS 1)
else()
  message(FATAL_ERROR "Unknown EquationOfState (EOS) : ${EOS}")
endif()
set(EQUATION_OF_STATE ${EOS})

if(NOT AFFINE OR NOT DEFINED AFFINE)
  set(AFFINE_OPTION "NOT_AFFINE")
else()
  set(AFFINE_OPTION "AFFINE")
  set(COORDINATE_SYSTEM "affine_coordinate")
endif()

if(NOT CUBED_SPHERE OR NOT DEFINED CUBED_SPHERE)
  set(CUBED_SPHERE_OPTION "NOT_CUBED_SPHERE")
else()
  set(CUBED_SPHERE_OPTION "CUBED_SPHERE")
  set(COORDINATE_SYSTEM "gnomonic_equiangle")
endif()

if(NOT NETCDF OR NOT DEFINED NETCDF)
  set(NETCDF_OPTION "NO_NETCDFOUTPUT")
else()
  set(NETCDF_OPTION "NETCDFOUTPUT")
  find_package(NetCDF REQUIRED)
endif()

if(NOT PNETCDF OR NOT DEFINED PNETCDF)
  set(PNETCDF_OPTION "NO_PNETCDFOUTPUT")
else()
  set(PNETCDF_OPTION "PNETCDFOUTPUT")
  find_package(PNetCDF REQUIRED)
endif()

if(NOT FITS OR NOT DEFINED FITS)
  set(FITS_OPTION "NO_FITSOUTPUT")
else()
  set(FITS_OPTION "FITSOUTPUT")
  find_package(cfitsio REQUIRED)
endif()

if(NOT HYDROSTATIC OR NOT DEFINED HYDROSTATIC)
  set(HYDROSTATIC_OPTION "NOT_HYDROSTATIC")
else()
  set(HYDROSTATIC_OPTION "HYDROSTATIC")
endif()

if(NOT DISORT OR NOT DEFINED DISORT)
  set(DISORT_OPTION "NOT_RT_DISORT")
else()
  set(DISORT_OPTION "RT_DISORT")
endif()

if(NOT RRTMG_SW OR NOT DEFINED RRTMG_SW)
  set(RRTMG_SW_OPTION "NOT_RT_RRTMG_SW")
else()
  set(RRTMG_SW_OPTION "RT_RRTMG_SW")
endif()

if(NOT RRTMG_LW OR NOT DEFINED RRTMG_LW)
  set(RRTMG_LW_OPTION "NOT_RT_RRTMG_LW")
else()
  set(RRTMG_LW_OPTION "RT_RRTMG_LW")
endif()

if(NOT PVFMM OR NOT DEFINED PVFMM)
  set(PVFMM_OPTION "DISABLE_PVFMM")
else()
  set(PVFMM_OPTION "ENABLE_PVFMM")
  set(MPI ON)
endif()

if(NOT MPI OR NOT DEFINED MPI)
  set(MPI_OPTION "NOT_MPI_PARALLEL")
else()
  set(MPI_OPTION "MPI_PARALLEL")
  find_package(
    MPI
    COMPONENTS C
    REQUIRED)
endif()

if(NOT GLOG OR NOT DEFINED GLOG)
  set(GLOG_OPTION "DISABLE_GLOG")
else()
  set(GLOG_OPTION "ENABLE_GLOG")
  find_package(glog REQUIRED)
  set(GLOG_LIBRARIES glog::glog)
endif()

if(NOT PYTHON_BINDINGS OR NOT DEFINED PYTHON_BINDINGS)
  set(PYTHON_BINDINGS_OPTION "NO_PYTHON_BINDINGS")
else()
  if(NOT DEFINED PYTHON_VERSION)
    find_package(Python3 3.8 REQUIRED COMPONENTS Interpreter Development.Module)
  else()
    set(Python3_FIND_STRATEGY VERSION)
    find_package(Python3 ${PYTHON_VERSION} EXACT REQUIRED
                 COMPONENTS Interpreter Development.Module)
  endif()
  set(PYTHON_BINDINGS_OPTION "PYTHON_BINDINGS")
endif()
