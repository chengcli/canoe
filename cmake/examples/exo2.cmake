# configuration for shallow water model on cubed sphere

# athena variables
set(NUMBER_GHOST_CELLS 3)
set(NTRACER 3)
set(EQUATION_OF_STATE shallow_yz)
set(NON_BAROTROPIC_EOS 0)
set(RSOLVER roe_shallow_yz)

# canoe variables
set(CUBED_SPHERE ON)
set(HYDROSTATIC ON)
set(NETCDF ON)
set(MPI ON)
set(PNETCDF ON)
