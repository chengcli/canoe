# define default parameters

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

list(APPEND CMAKE_PREFIX_PATH "/usr/include/mpich-x84_64" "/usr/lib64/mpich/")

# populate the default values
set_if_empty(NVAPOR 0)

set_if_empty(NCLOUD 0)

set_if_empty(NCHEMISTRY 0)

set_if_empty(NTRACER 0)

set_if_empty(NSTATIC 0)

set_if_empty(EQUATION_OF_STATE "ideal_moist")

set_if_empty(TASKLIST TimeIntegratorTaskList)

set_if_empty(PLANET "UNKNOWN")

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

if(NOT MPI OR NOT DEFINED MPI)
  set(MPI_OPTION "NOT_MPI_PARALLEL")
else()
  set(MPI_OPTION "MPI_PARALLEL")
  find_package(MPI REQUIRED)
endif()
