# Instructions on running 2D shallow water (Global Cubed Sphere) model
## Building the code
First, create a build directory to host this build. For example, `mkdir build`, and then `cd build`.
In the build directory, run
```
cmake .. -DTASK=exo2
```
And then run the make command. Run
```
make
```
Or if parallel execution is desired, use
```
make -j8
```
Instead. If no error shows up, then the build is complete. The desired executables built from the problem generators will show up with the name `pgen_name.release`. Input files will also be copied over. Execute the code using
```
mpirun -np [Number of processors] ./pgen_name.release -i pgen_name.inp
```

## Writing problem generator
In writing the problem generator, what you need to do is mostly similar to what other files in canoe needs to be done. In using the gnomonic equiangle coordinate, we have provided utilities to help identify where the current grid is. We suggest using `examples\2023-Chen-exo3\W92.cpp` as the example. Detailed suggestions are here:
First, the exocubed headers are needed:
```
// exo3
#include <exo3/cubed_sphere.hpp>
#include <exo3/cubed_sphere_utility.hpp>
```
Then in problem generator, get the Cubed Sphere pointer using
```
auto pexo3 = pimpl->pexo3;
```
Note that in forcing, we need to redirect from the meshblock
```
auto pexo3 = pmb->pimpl->pexo3;
```
Then at each point, you can get the latitude and longitude in radians
```
Real lat, lon;
pexo3->GetLatLon(&lat, &lon, k, j, i);
```
Note that the vectors also need to be transformed. The common practice is to transform the vectors to lat-lon grid, work in the lat-lon grid, and then transform the results back to the gnomonic equiangle grid. For example,
```
Real U, V;
pexo3->GetUV(&U, &V, w(IVY, k, j, i), w(IVZ, k, j, i), k, j, i);
Real ll_acc_U = f * V;
Real ll_acc_V = -f * U;
Real acc1, acc2, acc3;
pexo3->GetVyVz(&acc2, &acc3, ll_acc_U, ll_acc_V, k, j, i);
```
An important side note is that, the primitive variables are in the contravariant basis, while the conservative variables are in covariant basis. So the calculations made from the primitives must be transformed once more to obtain the change in the covariant basis.
```
pexo3->ContravariantVectorToCovariant(j, k, acc2, acc3, &acc2, &acc3);
```
A complete implementation of the Coriolis forcing is available in `examples\2023-Chen-exo3\W92.cpp` as a reference.
To join the problem generator written to the list of compilation, put the pgen file to `examples\2023-Chen-exo3\` and modify the `examples\2023-Chen-exo3\CMakeLists.txt`. Add one line after the block within the if statement
```
if(${EQUATION_OF_STATE} STREQUAL "shallow_yz")
```
It would become apparent what to do the next. Add
```
setup_problem(your_pgen_name)
```
After the current examples.

## Writing input file
Writing the input file is simple. Using the same example, we take a look at the `examples\2023-Chen-exo3\W92.inp` file:
```
<mesh>
nx1        = 3            # Number of zones in X1-direction
x1min      = 6371000.           # minimum value of X1
x1max      = 6371001.           # maximum value of X1
ix1_bc     = reflecting   # inner-X1 boundary flag
ox1_bc     = reflecting   # outer-X1 boundary flag

nx2        = 96            # Number of zones in X2-direction
x2min      = 0.           # minimum value of X2
x2max      = 1.E3         # maximum value of X2
ix2_bc     = periodic   # inner-X2 boundary flag
ox2_bc     = periodic   # outer-X2 boundary flag

nx3        = 144           # Number of zones in X3-direction
x3min      = 0.           # minimum value of X3
x3max      = 1.E3         # maximum value of X3
ix3_bc     = periodic   # inner-X3 boundary flag
ox3_bc     = periodic   # outer-X3 boundary flag
```

In the mesh, x2min, x2max, and x3min, x3max are not important. But make sure that the boundary conditions for x2 and x3 are set to be `periodic`.
The way that we allocate the grids in athena makes the panels be placed 2 by 3 in the mesh. So if you want to do a simulation with N by N grids in each panel, with $`N\times N\times 6`$ cells in total, specify `nx2` to be $`2N`$ and `nx3` to be $`3N`$.
It is required that each panel is divided into equal number of meshblocks along `x2` and `x3` directions. This is because the `x2` direction in one panel potentially becomes the `x3` direction in another. Additionally, the number of meshblock must be equal to the number of processors used. This requirement comes from the MPI communication used in the implementation.
