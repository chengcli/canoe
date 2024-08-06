# configure file for bryan test case

# athena variables
set(NUMBER_GHOST_CELLS 3)

# canoe configure
set(NVAPOR 1)
set(NCLOUD 1)
set(NPHASE_LEGACY 2)
set(NETCDF OFF)
set(PNETCDF ON)
set(MPI ON)
set(EQUATION_OF_STATE ideal_moist)
# set(RSOLVER lmars)
set(RSOLVER hllc_transform)
