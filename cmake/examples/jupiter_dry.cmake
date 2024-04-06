# configure file for jupiter_dry

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

# canoe configure
set(PNETCDF ON)
set(DISORT ON)
set(PYTHON_BINDINGS ON)
set(MPI ON)
set(TASKLIST ImplicitHydroTasks)
set(RSOLVER lmars)
