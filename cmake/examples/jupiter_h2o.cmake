# configure file for test jupiter crm

# athena variables
set(NUMBER_GHOST_CELLS 3)

# canoe configure
set(NVAPOR 1)
set(NCLOUD 2)
set(NPHASE_LEGACY 3)
set(PNETCDF ON)
set(MPI ON)
set(EQUATION_OF_STATE ideal_moist)
set(TASKLIST ImplicitHydroTasks)
# set_if_empty(RSOLVER hllc_transform)
set(RSOLVER lmars)
