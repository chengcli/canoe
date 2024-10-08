#! /usr/bin/env python3
import sys

sys.path.append("../python")
sys.path.append(".")

import canoe
from canoe import load_configure
from canoe.harp import radiation_band
from netCDF4 import Dataset
from scipy.interpolate import interp1d
from multiprocessing import Pool, cpu_count
from typing import List
import numpy as np


def get_gauss_legendre(n: int):
    """
    Get the Gauss-Legendre quadrature points and weights

    Parameters
    ----------
    n : int
        The number of quadrature points

    Returns
    -------
    p : np.ndarray
        The quadrature points
    w : np.ndarray
        The quadrature weights
    """

    p, w = np.polynomial.legendre.leggauss(n)  # range is -1 to 1
    # now change it from 0 to 1
    return (p + 1.0) / 2.0, w / 2.0


class CorrelatedKtable:
    """
    A base class to generate a correlated k-table

    Attributes
    ----------
    name : str
        The name of the absorbing species
    """

    def __init__(self, species: List[str]):
        """
        Initialize the class

        Parameters
        ----------
        name : str
            The name of the absorbing species
        """
        self.species = species
        self.name = "-".join(species)

    def load_opacity(self, fname: str) -> None:
        """
        Load the opacity data from a file.
        This is a virtual function that should be implemented by the derived class.

        Parameters
        ----------
        fname : str
            The name of the file
        """
        pass

    def make_ck_coeff(
        self,
        kcoeff: np.ndarray,
        ilayer: int,
        nbins: int,
        npoints: int,
        log_opacity: bool = True,
    ) -> np.ndarray:
        """
        Make the correlated k-table axis and sort the k-coefficients

        Parameters
        ----------
        kcoeff : np.ndarray
            The absorption coefficients with shape (nwave, nlayer)
        ilayer : int
            The layer index to use for the bin_divides
        nbins : int
            The number of bins to divide the spectral band
        npoints : int
            The number of points in each bin
        log_opacity : bool
            If True, the absorption coefficients are in log scale

        Returns
        -------
        ckcoeff : np.ndarray
            The correlated k-coefficients with shape (nbins * npoints, nlayer)
        """
        bin_divides = np.zeros(nbins + 1)
        bin_divides[0] = 0.0
        bin_divides[-1] = 1.0

        # use ilayer to determine the bin_divides
        if log_opacity:
            lnkmax = kcoeff[:, ilayer].max()
            lnkmin = kcoeff[:, ilayer].min()
        else:
            lnkmax = np.log(kcoeff[:, ilayer].max())
            lnkmin = np.log(kcoeff[:, ilayer].min())

        nwave = kcoeff.shape[0]
        nlayer = kcoeff.shape[1]

        for i in range(nlayer):
            kcoeff[:, i] = np.sort(kcoeff[:, i])

        for i in range(1, nbins):
            lnk = lnkmin + (lnkmax - lnkmin) * i / nbins
            if log_opacity:
                bin_divides[i] = np.searchsorted(kcoeff[:, ilayer], lnk) / nwave
            else:
                bin_divides[i] = np.searchsorted(kcoeff[:, ilayer], exp(lnk)) / nwave

        # print('bin_divides:', bin_divides)

        gaxis = np.zeros(nbins * npoints)
        weights = np.zeros(nbins * npoints)
        ckcoeff = np.zeros((nbins * npoints, nlayer))

        gg, ww = get_gauss_legendre(npoints)
        for i in range(nbins):
            gaxis[i * npoints : (i + 1) * npoints] = (
                gg * (bin_divides[i + 1] - bin_divides[i]) + bin_divides[i]
            )
            weights[i * npoints : (i + 1) * npoints] = ww * (
                bin_divides[i + 1] - bin_divides[i]
            )

        for j in range(nlayer):
            kcoeff_func = interp1d(np.arange(nwave), kcoeff[:, j])
            ckcoeff[:, j] = kcoeff_func(gaxis * (nwave - 1))

        self.gaxis = gaxis
        self.weights = weights

        return ckcoeff

    def write_opacity(self, fname: str):
        """
        Write the correlated k-table to a file

        Parameters
        ----------
        fname : str
            The name of the file
        """
        ncfile = Dataset(fname, "w")
        ncfile.createDimension("gaxis", len(self.gaxis))
        dim = ncfile.createVariable("gaxis", "f8", ("gaxis",))
        dim[:] = self.gaxis
        dim.long_name = "gaussian quadrature points"
        dim.units = "1"

        ncfile.createDimension("Wavenumber", len(self.gaxis))
        dim = ncfile.createVariable("Wavenumber", "f8", ("Wavenumber",))
        dim[:] = self.wave
        dim.long_name = "gaussian quadrature equivalent wavenumber"
        dim.units = "1/cm"

        ncfile.createDimension("weights", len(self.weights))
        dim = ncfile.createVariable("weights", "f8", ("weights",))
        dim[:] = self.weights
        dim.long_name = "gaussian quadrature weights"
        dim.units = "1"

        ncfile.createDimension("Pressure", len(self.pres))
        dim = ncfile.createVariable("Pressure", "f8", ("Pressure",))
        dim[:] = self.pres
        dim.long_name = "reference pressure"
        dim.units = "pa"

        ncfile.createDimension("TempGrid", len(self.temp_grid))
        dim = ncfile.createVariable("TempGrid", "f8", ("TempGrid",))
        dim[:] = self.temp_grid
        dim.long_name = "temperature anomaly grid"
        dim.units = "K"

        var = ncfile.createVariable("Temperature", "f8", ("Pressure"))
        var[:] = self.temp
        var.long_name = "reference temperature"
        var.units = "K"

        for name in self.species:
            var = ncfile.createVariable(
                name, "f8", ("Wavenumber", "Pressure", "TempGrid")
            )
            var[:] = self.ckcoeff[name]
            var.long_name = "correlated k-coefficients"
            var.units = self.kunits[name]

        ncfile.close()
        print("Correlated k-table written to", fname)


class HitranCorrelatedKtable(CorrelatedKtable):
    """
    Derived class to generate a correlated k-table from HITRAN line-by-line opacity
    """

    def load_opacity(self, fname: str):
        """
        Load the opacity data from a file. Overrides the base class method.

        Parameters
        ----------
        fname : str
            The name of the file
        """
        data = Dataset(fname, "r")
        self.kcoeff = {}
        self.kunits = {}

        for name in self.species:
            self.kcoeff[name] = data.variables[name][:]
            self.kunits[name] = data.variables[name].units
        self.pres = data.variables["Pressure"][:]
        self.temp = data.variables["Temperature"][:]
        self.temp_grid = data.variables["TempGrid"][:]

    def make_cktable(self, wmin: float, wmax: float, nbins: int = 1, npoints: int = 50):
        """
        Make the correlated k-table. This function will call make_ck_axis for each temperature
        grid point.

        Parameters
        ----------
        nbins : int
            The number of bins to divide the spectral band
        npoints : int
            The number of points in each bin
        """
        nlayer = len(self.pres)
        ntemp = len(self.temp_grid)

        self.ckcoeff = {}
        for name in self.species:
            self.ckcoeff[name] = np.zeros((nbins * npoints, nlayer, ntemp))
            for i in range(ntemp):
                self.ckcoeff[name][:, :, i] = self.make_ck_coeff(
                    self.kcoeff[name][:, :, i], nlayer // 2, nbins, npoints
                )
        self.wave = wmin + self.gaxis * (wmax - wmin)


def run_cktable_one_band(bname: str):
    band = radiation_band(bname, config)

    species = list(map(str, config[bname]["opacity"]))
    wmin, wmax = band.get_range()
    ab_ck = HitranCorrelatedKtable(species)
    ab_ck.load_opacity(opacity_input + "-" + bname + ".nc")
    ab_ck.make_cktable(wmin, wmax)
    ab_ck.write_opacity(opacity_output + "-" + bname + ".nc")


if __name__ == "__main__":
    canoe.start()

    opacity_input = "giants-lbl"
    opacity_output = "giants-ck"
    opacity_config = "giants-op.yaml"
    max_threads = cpu_count()

    config = load_configure(opacity_config)

    bnames = list(map(str, config["bands"]))
    for bname in bnames:
        run_cktable_one_band(bname)
    # pool = Pool(max_threads)
    # pool.map(run_cktable_one_band, bnames)
