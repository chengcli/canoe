import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata
from tqdm import tqdm
import os

# Path to your single combined .nc file
filepath = "hotjupiter-main.nc"

# Variables to process
variables_to_process = ["temp", "theta", "rho", "press", "vlat", "vlon", "vel1"]


# Function to interpolate and extrapolate data
def interpolate_data(lat, lon, data, grid_lats, grid_lons, height):
    new_data = np.empty((height, grid_lats.shape[0], grid_lats.shape[1]))

    for z in range(height):
        points = np.column_stack((lat[z, :].ravel(), lon[z, :].ravel()))
        values = data[z, :].ravel()
        grid_z = griddata(points, values, (grid_lats, grid_lons), method="linear")

        nan_indices = np.isnan(grid_z)
        if np.any(nan_indices):
            grid_z[nan_indices] = griddata(
                points,
                values,
                (grid_lats[nan_indices], grid_lons[nan_indices]),
                method="nearest",
            )

        new_data[z, :, :] = grid_z

    return new_data


def process_timestep(t):
    with Dataset(filepath, mode="r") as nc:
        lat = nc.variables["lat"][t, :, :] * 180 / np.pi
        lon = nc.variables["lon"][t, :, :] * 180 / np.pi
        x1 = nc.variables["x1"][:]
        time_value = nc.variables["time"][t]

        data_arrays = {
            var: nc.variables[var][t, :, :, :] for var in variables_to_process
        }

    # Interpolation to regular lat-lon grid
    lat_min, lat_max = -90, 90
    lon_min, lon_max = 0, 360
    K = 50
    new_lats = np.linspace(lat_min, lat_max, K)
    new_lons = np.linspace(lon_min, lon_max, K)
    grid_lons, grid_lats = np.meshgrid(new_lons, new_lats)

    new_data_arrays = {}
    for var, data in data_arrays.items():
        new_data_arrays[var] = interpolate_data(
            lat, lon, data, grid_lats, grid_lons, len(x1)
        )

    return new_data_arrays


# Extract number of time steps in the file
with Dataset(filepath, "r") as nc:
    num_time_steps = len(nc.variables["time"][:])

# Process each time step and store results in a list
results = []
for t in tqdm(range(num_time_steps), desc="Processing time steps"):
    results.append(process_timestep(t))

with Dataset(filepath, "r") as nc:
    num_time_steps = len(nc.variables["time"][:])
    x1 = nc.variables["x1"][:]

K = 50  # You can adjust this value as needed

# Save the processed data to a new NetCDF file
output_filepath = "polar_" + os.path.basename(filepath)
with Dataset(output_filepath, mode="w") as new_nc:
    # Create the dimensions
    new_nc.createDimension("time", None)  # Unlimited dimension (usually time)
    new_nc.createDimension("x1", len(x1))
    new_nc.createDimension("lat", K)
    new_nc.createDimension("lon", K)

    # Create and assign the time variable
    time_var = new_nc.createVariable("time", np.float64, ("time",))
    with Dataset(filepath, "r") as nc:
        time_var[:] = nc.variables["time"][:]

    # Create and assign the x1 variable
    new_nc.createVariable("x1", np.float32, ("x1",))
    new_nc.variables["x1"][:] = x1

    lat_var = new_nc.createVariable("lat", np.float32, ("lat",))
    lat_var[:] = np.linspace(-90, 90, K)

    lon_var = new_nc.createVariable("lon", np.float32, ("lon",))
    lon_var[:] = np.linspace(0, 360, K)

    # Create and assign the processed data
    for var in variables_to_process:
        new_var = new_nc.createVariable(var, np.float32, ("time", "x1", "lat", "lon"))
        for t in range(num_time_steps):
            new_var[t, :, :, :] = results[t][var]
        # Create and assign the latitude and longitude variables


print(f"Processing completed. Results saved to {output_filepath}.")
