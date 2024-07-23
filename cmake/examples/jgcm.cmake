# configure file for test cubedsphere jcrm

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

# canoe configure
set(CUBED_SPHERE ON)
set(NVAPOR 2)
set(NCLOUD 4)
set(NPHASE_LEGACY 3)
set(NETCDF ON)
set(PNETCDF OFF)
set(MPI ON)
set(TASKLIST ImplicitHydroTasks)
set(RSOLVER hllc_transform)
