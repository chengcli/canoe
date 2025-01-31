# configuration for titan hydrodynamcis

# athena variables
set(NUMBER_GHOST_CELLS 3)
# set(EQUATION_OF_STATE adiabatic)
set(EQUATION_OF_STATE ideal_moist)
set(NON_BAROTROPIC_EOS 1)
set(RSOLVER lmars)

# canoe variables
set(MPI ON)
set(PNETCDF ON)
