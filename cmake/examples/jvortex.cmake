# configure file for test jvortex

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

# canoe configure
set(NVAPOR 2)
set(NCLOUD 4)
set(NPHASE_LEGACY 3)
set(NETCDF ON)
set(PNETCDF ON)
set(MPI ON)
set(TASKLIST ImplicitHydroTasks)
set(RSOLVER lmars)
