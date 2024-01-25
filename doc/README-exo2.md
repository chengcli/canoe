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


## Writing input file