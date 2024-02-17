import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import interp1d
import os
from tqdm import tqdm

# Set the filepath to your single combined .nc file
filepath = "cart_hotjupiter-main.nc"  # Change this to your file path

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
            x1 = nc.variables["x1"][:]
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

temp_avg = averages["temp"]  # Temperature average over all files


# Calculate the pressure at each height level for each latitude using the hydrostatic equation
# We assume that the first level is the surface (P = P0) and that 'x1_levels' is an array of height levels in meters
x1_levels = x1
pressures = np.zeros((temp_avg.shape[0], temp_avg.shape[1]))  # [height, lat]
pressures[0, :] = P0

for j in range(temp_avg.shape[1]):  # Loop over all latitudes
    for i in range(1, temp_avg.shape[0]):  # Loop over all height levels
        dz = x1_levels[i] - x1_levels[i - 1]  # Difference in height between two levels
        T_avg = 0.5 * (
            temp_avg[i, j] + temp_avg[i - 1, j]
        )  # Average temperature between two levels
        pressures[i, j] = pressures[i - 1, j] * np.exp(
            -g * dz / (R * T_avg)
        )  # Hydrostatic equation

# Interpolate the variables to the new pressure levels for each latitude
averages_interpolated = {}
for var, avg in averages.items():
    averages_interpolated[var] = np.zeros(
        (len(pressure_levels), temp_avg.shape[1])
    )  # [pressure, lat]
    for j in range(temp_avg.shape[1]):
        f = interp1d(pressures[:, j], avg[:, j], fill_value="extrapolate")
        averages_interpolated[var][:, j] = f(pressure_levels)

# Also interpolate 'upvp' and 'uv_variance'
upvp_interpolated = np.zeros(
    (len(pressure_levels), temp_avg.shape[1])
)  # [pressure, lat]
upwp_interpolated = np.zeros(
    (len(pressure_levels), temp_avg.shape[1])
)  # [pressure, lat]
uv_variance_interpolated = np.zeros(
    (len(pressure_levels), temp_avg.shape[1])
)  # [pressure, lat]

for j in range(temp_avg.shape[1]):
    upvp_function = interp1d(pressures[:, j], upvp_avg[:, j], fill_value="extrapolate")
    upwp_function = interp1d(pressures[:, j], upwp_avg[:, j], fill_value="extrapolate")
    uv_variance_function = interp1d(
        pressures[:, j], uv_variance_avg[:, j], fill_value="extrapolate"
    )
    upvp_interpolated[:, j] = upvp_function(pressure_levels)
    upwp_interpolated[:, j] = upwp_function(pressure_levels)
    uv_variance_interpolated[:, j] = uv_variance_function(pressure_levels)

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
    for var, data in averages_interpolated.items():
        new_var = new_nc.createVariable(var + "_avg", np.float32, ("pressure", "lat"))
        new_var[:] = data

    # Create variables for 'upvp' and 'uv_variance'
    upvp_var = new_nc.createVariable("upvp", np.float32, ("pressure", "lat"))
    upvp_var[:] = upvp_interpolated

    upwp_var = new_nc.createVariable("upwp", np.float32, ("pressure", "lat"))
    upwp_var[:] = upwp_interpolated

    uv_variance_var = new_nc.createVariable(
        "uv_variance", np.float32, ("pressure", "lat")
    )
    uv_variance_var[:] = uv_variance_interpolated

    # Optionally, add descriptions, units, or other metadata as attributes to the variables

# Print out a message indicating the script has finished
print(f"Data processing completed. Results saved to {output_file}.")
