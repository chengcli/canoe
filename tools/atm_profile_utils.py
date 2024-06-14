#! /usr/bin/env python3
from netCDF4 import Dataset
from typing import Tuple, Optional, Dict, List
from dataclasses import dataclass
import argparse, datetime, os, sys
import numpy as np


@dataclass
class __Constants:
    """
    A class to hold physical constants

    This class is not meant to be used outside of this module.

    Attributes:
    -----------
    kRgas: float
        Ideal gas constant in J/(mol K)
    kUnitConversions: Dict[str, float]
        A dictionary of unit conversion factors
    """

    kRgas: float = 8.31446261815324

    # Define conversion factors
    kUnitConversions = {
        "km": 1e3,
        "m": 1,
        "cm": 1e-2,
        "bar": 1e5,
        "Pa": 1,
        "hPa": 1e2,
        "mbar": 1e2,
        "pa": 1,
        "hpa": 1e2,
        "ppmv": 1e-6,
        "ppbv": 1e-9,
        "ppm": 1e-6,
        "ppb": 1e-9,
        "mol/mol": 1,
        "1": 1,
    }


def get_command_string() -> str:
    # Get the current script file path
    script_path = os.path.basename(__file__)

    # Get all command-line arguments
    args = sys.argv[1:]

    # Form the command string
    command_string = f"python3 {script_path} {' '.join(args)}"

    return command_string


def parse_column_name(col_name: str) -> Tuple[str, Optional[str]]:
    """
    Parse the column name and return the name and unit.

    Parameters:
    -----------
    col_name: str
        The name of the column.

    Returns:
    --------
    name, unit: Tuple[str, Optional[str]]
        The name and unit of the column.
    """

    if "[" in col_name and "]" in col_name:
        name = col_name[: col_name.index("[")]
        unit = col_name[col_name.index("[") + 1 : col_name.index("]")]
    else:
        name = col_name
        unit = None
    return name, unit


def read_atm_profile_nc(
    filename: str, dry_air_name: Optional[str], vapors: List[str], wghts: List[float]
) -> Dict[str, np.ndarray]:
    """
    Read an atmospheric profile from a netCDF file.

    Parameters:
    -----------
    filename: str
        The name of the netCDF file to read.
    dry_air_name: Optional[str]
        The name of the dry air species.
    vapors: List[str]
        A list of vapor names to read from the file.
    wghts: List[float]
        A list of molecular weights for all species (except for the dry air) in g/mol.

    Returns:
    --------
    atm: Dict[str, np.ndarray]
        A dictionary containing the atmospheric profile
    """

    if len(vapors) != len(wghts):
        raise ValueError(
            "The number of species names and molecular weights must be the same."
        )

    data = Dataset(filename, "r")

    # Initialize atmosphere dictionary
    atm = {}

    # index
    atm["IDX"] = np.arange(data.dimensions["x1"].size) + 1

    # altitude in m
    atm["HGT"] = data["x1"][:]

    # pressure in Pa
    atm["PRE"] = data["press"][0, :, 0, 0]

    # temperature in K
    atm["TEM"] = data["temp"][0, :, 0, 0]

    # mean molecular weight in kg/mol
    mmw = data["rho"][0, :, 0, 0] * __Constants.kRgas * atm["TEM"] / atm["PRE"]

    # mole fraction of dry air
    dry_mf = np.ones(atm["HGT"].shape)

    # mole fractions
    if len(vapors) == 1:
        atm[vapors[0]] = data["vapor"][0, :, 0, 0] / wght[0] * mmw
        dry_mf -= atm[vapors[0]]
    else:
        for i, spec in enumerate(vapors):
            atm[spec] = data["vapor" + str(i + 1)][0, :, 0, 0] / wght[i] * mmw
            dry_mf -= atm[spec]

    # dry air species
    if dry_air_name is not None:
        atm[dry_air_name] = dry_mf

    return atm


def read_atm_profile_txt(filename: str) -> Dict[str, np.ndarray]:
    """
    Read an atmospheric profile from a text file.

    Parameters:
    -----------
    filename: str
        The name of the text file to read.

    Returns:
    --------
    data: Dict[str, np.ndarray]
        A dictionary containing the atmospheric profile
    """

    # Initialize atmosphere dictionary
    data = {}

    # Read the file
    with open(filename, "r") as f:
        lines = f.readlines()

    # Process the header
    header = None
    for j, line in enumerate(lines):
        if line.startswith("#"):
            continue
        if header is None:
            header = line.strip().split()
            for col in header:
                name, unit = parse_column_name(col)
                if (unit is not None) and (unit not in __Constants.kUnitConversions):
                    raise ValueError(f"Unknown unit: {unit}")
                else:
                    data[name] = []
        else:
            break

    header = []
    # Process the data
    for line in lines[j:]:
        if line.startswith("#"):
            continue
        values = line.strip().split()

        if len(header) == 0:
            header = values.copy()
            continue

        for i, value in enumerate(values):
            name, unit = parse_column_name(header[i])

            # default unit for HGT is km
            if name == "IDX":
                unit = "1"
            elif name == "HGT" and unit is None:
                unit = "km"
            elif name == "PRE" and unit is None:
                unit = "mbar"
            elif name == "TEM" and unit is None:
                unit = "K"
            elif unit is None:
                unit = "ppmv"

            if name in data:
                data[name].append(
                    float(value) * __Constants.kUnitConversions.get(unit, 1)
                )
            else:
                data[name] = [float(value) * __Constants.kUnitConversions.get(unit, 1)]

    if "IDX" not in data:
        raise ValueError("Missing required column: IDX")

    if "HGT" not in data:
        raise ValueError("Missing required column: HGT")

    if "PRE" not in data:
        raise ValueError("Missing required column: PRE")

    if "TEM" not in data:
        raise ValueError("Missing required column: TEM")

    # Convert lists to numpy arrays
    for key in data:
        if key == "IDX":
            data[key] = np.array(data[key], dtype=int)
        else:
            data[key] = np.array(data[key])

    return data


def write_atm_profile(
    data: Dict[str, np.ndarray],
    filename: str,
    units: List[str] = ["m", "pa"],
    comment: Optional[str] = None,
) -> None:
    """
    Write an atmospheric profile to a text file.

    Parameters:
    -----------
    data: Dict[str, np.ndarray]
        A dictionary containing the atmospheric profile
    filename: str
        The name of the text file to write.
    units: List[str]
        A list of units for altitude, pressure, and mole fractions.
    comment: Optional[str]
        A comment string to write to the file.

    Returns:
    --------
    None
    """
    with open(filename, "w") as f:
        # write comments
        f.write(f"# File generated by atm_profile_utils.py\n")
        f.write(f"# Date: {datetime.datetime.now()}\n")
        f.write(f"# Width of first column: 7 characters\n")
        f.write(f"# Width of the other columns: 12 characters\n")
        if comment:
            f.write(f"# {comment}\n\n")

        # Write header
        headers = []
        for key in data.keys():
            if key == "IDX":
                headers.append(f"{key}".ljust(6))
            elif key == "HGT":
                headers.append(f"{key}[{units[0]}]".ljust(11))
            elif key == "PRE":
                headers.append(f"{key}[{units[1]}]".ljust(11))
            elif key == "TEM":
                headers.append(f"{key}[K]".ljust(11))
            else:
                headers.append(f"{key}[{units[2]}]".ljust(11))
        f.write(" ".join(headers) + "\n")

        # Write data
        num_rows = len(data["IDX"])
        for i in range(num_rows):
            row = []
            for key in data.keys():
                if key == "HGT":
                    value = data[key][i] / __Constants.kUnitConversions.get(units[0], 1)
                elif key == "PRE":
                    value = data[key][i] / __Constants.kUnitConversions.get(units[1], 1)
                elif key in ["IDX", "TEM"]:
                    value = data[key][i]
                else:
                    value = data[key][i] / __Constants.kUnitConversions.get(units[2], 1)

                if key == "IDX":
                    row.append(f"{value}".ljust(6))
                else:
                    formatted_value = f"{value:.6g}"
                    if len(formatted_value) > 10:
                        formatted_value = f"{value:.4e}"
                    row.append(formatted_value.ljust(11))
            f.write(" ".join(row) + "\n")
    print(f"Atmosphere profile written to {filename}")


def read_atm_profile(
    filename: str,
    dry_air_name: Optional[str] = None,
    vapors: List[str] = [],
    wghts: List[float] = [],
) -> Dict[str, np.ndarray]:
    """
    Read an atmospheric profile from a file.

    This function checks the file extension and calls the appropriate function to read the file.
    The supported file formats are:
    - .atm: text file, calls read_atm_profile_txt
    - .txt: text file, calls read_atm_profile_txt
    - .nc: netCDF file, calls read_atm_profile_nc

    Parameters:
    -----------
    filename: str
        The name of the file to read.
    species: Optional[List[str]]
        A list of species names to read from the file.
    wght: Optional[List[float]]
        A list of molecular weights for the species in g/mol.
    num_vapors: Optional[int]
        Total number of condensible species in the atmosphere
    dry_molfrac: Optional[List[float]]
        A list of the mole fractions of the dry species, in the dry atmosphere

    Returns:
    --------
    data: Dict[str, np.ndarray]
        A dictionary containing the atmospheric profile
    """
    # split the filename into the root and the extension
    extension = os.path.splitext(filename)[1]
    if extension == ".atm" or extension == ".txt":
        data = read_atm_profile_txt(filename)
    elif extension == ".nc":
        data = read_atm_profile_nc(filename, dry_air_name, vapors, wghts)
    return data


def get_species_names(atm: Dict[str, np.ndarray]) -> List[str]:
    """
    Get the species names from an atmospheric profile.

    Parameters:
    -----------
    atm: Dict[str, np.ndarray]
        A dictionary containing the atmospheric profile

    Returns:
    --------
    species: List[str]
        A list of species names
    """
    species = []
    for key in atm.keys():
        if key not in ["IDX", "HGT", "PRE", "TEM"]:
            species.append(key)
    return species


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Read and process an atmospheric profile file. \
                                     Supported file formats are .atm, .txt, and .nc. \
                                     If the file is an .nc file, the vapor names and molecular weights \
                                     must be provided. Use the -v and -m options to \
                                     provide the vapor names. \
                                     The length of the vapor names and molecular weights must be the same. \
                                     The orderings of the vapor names and molecular \
                                     weights matter. The vapor names must be in the same order as they appear \
                                     in the netCDF file described by vapor1, vapor2, etc. \
                                     The molecular weights must be in the same order as \
                                     the vapor names. \
                                     Use the -o option to write the output to a text file."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="The name of the atmospheric profile file to process.",
    )
    parser.add_argument(
        "-o", "--output", required=False, help="The name of the output file."
    )
    parser.add_argument(
        "-u",
        "--units",
        required=False,
        default="m,pa,1",
        help="A list of units for altitude, pressure, and mole fractions, seperated by comma.",
    )
    parser.add_argument(
        "-v",
        "--vapor",
        required=False,
        help="A list of vapor names to read from the file, sperated by comma.",
    )
    parser.add_argument(
        "-w",
        "--wght",
        required=False,
        help="A list of molecular weights in g/mol for the vapors, seperated by comma.",
    )
    parser.add_argument(
        "-d", "--dry", required=False, help="name of the dry air species"
    )
    args = vars(parser.parse_args())

    if args["dry"] is not None:
        dry = args["dry"]
        print("Dry air species:", dry)
    else:
        dry = None

    if args["vapor"] is not None:
        vapor = args["vapor"].split(",")
        print("Vapors:", vapor)
    else:
        vapor = None

    if args["wght"] is not None:
        # g/mol to kg/mol
        wght = [float(m) * 1.0e-3 for m in args["wght"].split(",")]
        print("Molecular weights:", [f"{w:.4g}" for w in wght])
    else:
        wght = None

    units = args["units"].split(",")
    if len(units) != 3:
        raise ValueError("The number of unit strings must be 3, e.g., km,pa,1")

    data = read_atm_profile(args["input"], dry_air_name=dry, vapors=vapor, wghts=wght)

    if args["output"]:
        comment = "Command: " + get_command_string()
        write_atm_profile(data, args["output"], units=units, comment=comment)
