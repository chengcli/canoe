# Description: Calculate the streamfunction from the velocity field, and save it to a new .nc file together with the original variables.

import netCDF4 as nc
import numpy as np
from scipy.interpolate import interp1d
from tqdm import tqdm


infile = "polar_hotjupiter-main.nc"
tgtfile = "polar_hotjupiter-psi.nc"
Rp = 1e8

# Read the dimensions
with nc.Dataset(infile, "r") as src:
    t_dim = src.dimensions["time"]
    lat = src.variables["lat"][:]
    lon = src.dimensions["lon"]
    dlon = src.variables["lon"][1] - src.variables["lon"][0]  # assume uniform grid
    x1 = src.variables["x1"][:]
    psi = np.zeros((len(t_dim), len(x1), len(lat), len(lon)))
    rho = src.variables["rho"][:]
    vel1 = src.variables["vel1"][:]
    vlon = src.variables["vlon"][:]
    for i in range(1, len(x1)):
        for j in range(len(lat)):
            psi[:, i, j, :] = psi[:, i - 1, j, :] + vlon[:, i, j, :] * (
                x1[i] - x1[i - 1]
            ) * rho[:, i, j, :] * dlon * np.pi / 180.0 * Rp * np.cos(
                lat[j] * np.pi / 180.0
            )
    for i in range(1, len(lat)):
        psi[:, :, i, :] = psi[:, :, i - 1, :] - vel1[:, :, i, :] * (
            lat[i] - lat[i - 1]
        ) * np.pi / 180.0 * Rp * rho[:, :, i, :] * dlon * np.pi / 180.0 * Rp * np.cos(
            lat[i] * np.pi / 180.0
        )

with nc.Dataset(infile) as src, nc.Dataset(tgtfile, "w") as dst:
    # copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)
    # copy dimensions
    for name, dimension in src.dimensions.items():
        dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None)
        )
    # copy all file data except for the excluded
    for name, variable in src.variables.items():
        x = dst.createVariable(name, variable.datatype, variable.dimensions)
        dst[name][:] = src[name][:]
        # copy variable attributes all at once via dictionary
        dst[name].setncatts(src[name].__dict__)

    psi_var = dst.createVariable("psi", psi.dtype, ("time", "x1", "lat", "lon"))
    psi_var[:] = psi
    psi_var.units = "kg/s"
