# configure file for hydrogen-water world

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

# canoe configure set(CUBED_SPHERE ON) set(NVAPOR 2) set(NCLOUD 4)
set(NVAPOR 1)
set(NCLOUD 2)
set(NPHASE_LEGACY 3)
set(PNETCDF ON)
set(MPI ON)
set(DISORT ON)
set(PYTHON_BINDINGS ON)
# set_if_empty(RSOLVER hllc_transform)
set(RSOLVER lmars)
