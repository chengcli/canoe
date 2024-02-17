import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import interp1d
import os
from tqdm import tqdm

# Set the filepath to your single combined .nc file
filepath = "pres_hotjupiter.nc"  # Change this to your file path

# Set the name of the output file where the results will be saved
output_file = "averages_and_products.nc"

# Variables to process
variables_to_process = ["temp", "vlat", "vlon", "vel1"]

# Initialize accumulators for summing data across all time steps
data_sums = {}
upvp_sum = None
upwp_sum = None
uv_variance_sum = None

timestep_count = 0

P0 = 1e5
R = 3779
g = 8.0

pressure_levels = np.linspace(1e5, 0.01e5, 100)  # Replace with actual pressure levels

# Extract number of time steps in the file
with Dataset(filepath, "r") as nc:
    num_time_steps = len(nc.variables["time"][:])


# Process each time step
start_t = 100
Rp = 1e8
for t in tqdm(range(start_t, num_time_steps), desc="Processing time steps"):
    with Dataset(filepath, mode="r") as nc:
        if t == start_t:
            x1 = nc.variables["press"][:]
            latitudes = nc.variables["lat"][:]
        # If this is the first time step, initialize data_sums and other accumulators
        if timestep_count == 0:
            for var in variables_to_process:
                # Initialize the sum arrays
                data_sums[var] = np.zeros(
                    (nc.variables[var].shape[1], nc.variables[var].shape[2])
                )

            upvp_sum = np.zeros(
                (nc.variables["vlat"].shape[1], nc.variables["vlat"].shape[2])
            )
            upwp_sum = np.zeros(
                (nc.variables["vlat"].shape[1], nc.variables["vlat"].shape[2])
            )
            uv_variance_sum = np.zeros(
                (nc.variables["vlat"].shape[1], nc.variables["vlat"].shape[2])
            )

        # Update sums
        mean_rho = np.mean(nc.variables["rho"][t], axis=2)
        for var in variables_to_process:
            data = nc.variables[var][t]  # We take the first index of the time dimension
            zonal_mean = np.mean(
                data, axis=2
            )  # Zonal mean over longitude, reducing dimension
            data_sums[var] += zonal_mean  # Summing up the zonal mean data

            if var in ["vlat", "vlon", "vel1"]:
                prime = data - np.expand_dims(
                    zonal_mean, axis=2
                )  # Subtracting zonal mean, keeping dimensions consistent

                if var == "vlat":
                    u_prime = prime
                else:
                    if var == "vlon":
                        v_prime = prime
                    else:
                        w_prime = prime
        # Calculate and sum u'v' and (u'^2 + v'^2) / 2 for each time step
        upvp = np.mean(u_prime * v_prime, axis=2)  # Zonal mean of the product
        upwp = np.mean(u_prime * w_prime, axis=2)  # Zonal mean of the product
        uv_variance = np.mean(
            (u_prime**2 + v_prime**2 + w_prime**2) / 2, axis=2
        )  # Zonal mean of the variance

        upvp_sum += upvp
        upwp_sum += upwp
        uv_variance_sum += uv_variance

        timestep_count += 1  # Update the timestep count

# After processing all files, calculate the averages
averages = {
    var: data_sum / timestep_count for var, data_sum in data_sums.items()
}  # Averaging over all files
upvp_avg = upvp_sum / timestep_count  # Averaging over all files
upwp_avg = upwp_sum / timestep_count  # Averaging over all files
uv_variance_avg = uv_variance_sum / timestep_count  # Averaging over all files

# Save the results to a new .nc file
with Dataset(output_file, mode="w") as new_nc:
    # Create dimensions
    new_nc.createDimension("pressure", len(pressure_levels))
    new_nc.createDimension("lat", averages["temp"].shape[1])

    # Create pressure variable
    pressure_var = new_nc.createVariable("pressure", np.float32, ("pressure",))
    pressure_var[:] = pressure_levels
    pressure_var.units = "Pa"

    # Create latitude variable
    lat_var = new_nc.createVariable("lat", np.float32, ("lat",))
    lat_var[:] = np.linspace(-90, 90, averages["temp"].shape[1])
    lat_var.units = "degrees_north"

    # Create variables for the zonal means and other variables
    for var, data in averages.items():
        new_var = new_nc.createVariable(var + "_avg", np.float32, ("pressure", "lat"))
        new_var[:] = data

    # Create variables for 'upvp' and 'uv_variance'
    upvp_var = new_nc.createVariable("upvp", np.float32, ("pressure", "lat"))
    upvp_var[:] = upvp_avg

    upwp_var = new_nc.createVariable("upwp", np.float32, ("pressure", "lat"))
    upwp_var[:] = upwp_avg

    uv_variance_var = new_nc.createVariable(
        "uv_variance", np.float32, ("pressure", "lat")
    )
    uv_variance_var[:] = uv_variance_avg

    # Optionally, add descriptions, units, or other metadata as attributes to the variables

# Print out a message indicating the script has finished
print(f"Data processing completed. Results saved to {output_file}.")
