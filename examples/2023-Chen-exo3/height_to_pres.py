import netCDF4 as nc
import numpy as np
from scipy.interpolate import interp1d
from tqdm import tqdm

# Open the original NetCDF file
with nc.Dataset("cart_hotjupiter-main.nc", "r") as src:
    # Read the dimensions
    t_dim = src.dimensions["time"]
    x2_dim = src.dimensions["lat"]
    x3_dim = src.dimensions["lon"]
    x1 = src.variables["x1"][:]  # Original vertical coordinates
    press_var = src.variables["press"][:]  # Original pressure variable

    # Assume 'new_press_levels' is provided as a 1D numpy array of the new pressure levels
    new_press_levels = np.linspace(1e5, 0.01e5, 100)

    # Create a new NetCDF file
    with nc.Dataset("pres_hotjupiter.nc", "w") as dst:
        # Copy dimensions from the source to destination, except for 'x1' which is replaced by 'press'
        dst.createDimension("time", len(t_dim))
        dst.createDimension("lat", len(x2_dim))
        dst.createDimension("lon", len(x3_dim))
        dst.createDimension("press", len(new_press_levels))

        # Copy all variables from the source to destination, except those that depend on 'x1'
        for name, variable in src.variables.items():
            if "x1" not in variable.dimensions:
                dst.createVariable(name, variable.datatype, variable.dimensions)
                dst.variables[name][:] = src.variables[name][:]

        # Define the new pressure coordinate variable
        new_press_var = dst.createVariable("press", new_press_levels.dtype, ("press",))
        new_press_var[:] = new_press_levels

        # Interpolate the variables that depend on 'x1'
        for name, variable in src.variables.items():
            if (
                name == "time"
                or name == "lat"
                or name == "lon"
                or name == "press"
                or name == "x1"
            ):
                continue
            print("Interpolating variable: " + name)
            if "x1" in variable.dimensions:
                # Create the new variable in the destination file
                new_dimensions = tuple(
                    "press" if dim == "x1" else dim for dim in variable.dimensions
                )
                interp_var = dst.createVariable(name, variable.datatype, new_dimensions)

                # Loop over the additional dimensions
                for t in tqdm(range(len(t_dim)), desc="Processing time steps"):
                    for x2 in range(len(x2_dim)):
                        for x3 in range(len(x3_dim)):
                            # Extract the slice of pressure values for the current point
                            press_slice = press_var[t, :, x2, x3]
                            # Extract the original data slice
                            original_data = variable[t, :, x2, x3]
                            # Create an interpolation function based on the original pressure and data
                            f = interp1d(
                                press_slice,
                                original_data,
                                kind="linear",
                                bounds_error=False,
                                fill_value="extrapolate",
                            )
                            # Interpolate to the new press levels
                            interp_data = f(new_press_levels)
                            # Insert the interpolated data into the variable
                            interp_var[t, :, x2, x3] = interp_data
