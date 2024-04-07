# configuration for Jupiter Juno forward and inversion calculation

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 2)
set_if_empty(NVAPOR 2)

# canoe variables
set_if_empty(NCLOUD 4)
set_if_empty(NTRACER 2)

# canoe configure
set(HYDROSTATIC ON)
set(NETCDF ON)
set(RadianceSolver lambert)
set(PYTHON_BINDINGS ON)
