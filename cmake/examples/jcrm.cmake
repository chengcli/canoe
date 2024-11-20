# configure file for test jvortex

# athena variables
set(NUMBER_GHOST_CELLS 3)

# canoe configure
set(NVAPOR 1)
set(NCLOUD 1)
set(NPRECIP 1)
set(NETCDF ON)
set(PNETCDF ON)
set(MPI ON)
set(EQUATION_OF_STATE ideal_moist)
set(RSOLVER lmars)
# set(RSOLVER hllc_transform)
