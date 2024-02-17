import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

# File and variable setup
nc_file_path = "averages_and_products.nc"  # Specify the path to your .nc file
variables_to_plot = [
    "temp_avg",
    "vlat_avg",
    "vlon_avg",
    "upvp",
    "upwp",
    "uv_variance",
]  # Specify the variable names

# Plot customizations
colorbar_ranges = {
    "temp_avg": (1200, 1700),
    "vlat_avg": (-800, 1200),
    "vlon_avg": (-50, 50),
    "upvp": (-20 * 1000, 20 * 1000),
    "upwp": (-250, 250),
    "uv_variance": (0, 180 * 1000),
}

colormaps = {
    "temp_avg": "jet",
    "vlat_avg": "jet",
    "vlon_avg": "seismic",
    "upvp": "seismic",
    "upwp": "seismic",
    "uv_variance": "jet",
}

aspect_ratio = (6, 4)  # Width, Height
x_label = "Latitude"
y_label = "Pressure (Pa)"

# Specific colorbar ticks for each variable
colorbar_ticks = {
    "temp_avg": np.linspace(1200, 1700, 6),  # Example specific ticks
    "vlat_avg": np.linspace(-800, 1200, 11),
    "vlon_avg": np.linspace(-50, 50, 11),
    "upvp": np.linspace(-20 * 1000, 20 * 1000, 9),
    "upwp": np.linspace(-250, 250, 6),
    "uv_variance": np.linspace(0, 180 * 1000, 7),
}

# Open the NetCDF file
with Dataset(nc_file_path, mode="r") as nc:
    # Get the pressure and latitude values
    pressure = nc.variables["pressure"][:]
    latitude = nc.variables["lat"][:]

    # Create a meshgrid for the plots (needed for contourf)
    X, Y = np.meshgrid(latitude, pressure / pressure[0])

    # Plot each variable
    for var_name in variables_to_plot:
        # Extract the data for the current variable
        data = nc.variables[var_name][:]

        # Create a new figure for the current variable
        plt.figure(figsize=aspect_ratio)

        # Create the contour plot
        contour = plt.contourf(
            X,
            Y,
            data,
            levels=np.linspace(
                colorbar_ranges[var_name][0], colorbar_ranges[var_name][1], 100
            ),
            cmap=colormaps[var_name],
            extend="both",
        )

        # Customize axes and labels
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(var_name)

        # Create the colorbar with specified ticks
        cbar = plt.colorbar(contour, label=var_name, ticks=colorbar_ticks[var_name])
        cbar.ax.set_yticklabels(
            [f"{tick:.0f}" for tick in colorbar_ticks[var_name]]
        )  # Set tick labels with desired format

        plt.gca().invert_yaxis()  # Invert the y-axis to have the surface (lower pressure) at the top
        plt.grid(True)  # Turn on the grid

        # Save the figure
        plt.savefig(f"{var_name}_plot.png", bbox_inches="tight")
        plt.close()
