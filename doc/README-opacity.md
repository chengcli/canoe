# Instruction on generating opacity

## Files required
1. atmospheric profile, e.g. `jupiter.atm` 
2. opacity definition file, e.g. `jupiter.yaml`

## Atmospheric profile
The atmospheric profile is a text file with named columns. The particular order of the columns is not important,
except for the first column, which is assumed to be an index of the vertical levels.
The vertical level index starts from 1 at the bottom of the atmosphere and progresses upwards. 
The rest of the columns are altitude, pressure, temperature, and gas mole fractions.
The following names are required to be present in the file:

- 'IDX': index of the layer
- 'ALT': altitude (km)
- 'PRE': pressure (mbar)
- 'TEM': temperature (K)

If no units are provided, the code assumes `km` for altitude, `mbar` for pressure, and `K` for temperature.
It is possible to use different units by specifying them after the name of the column
(without spaces), and enclosed in square brackets, e.g. 'ALT[km]'. The code will convert the units automatically to SI units.
The recognized units are: `km`, `m`, `cm`, `bar`, `Pa`, `hPa`, `mbar`, `pa`, `hpa`, `K`

Following the required columns, the file can contain any number of additional columns with gas mole fractions.
The names of the columns should correspond to the names of the gases in the opacity definition file.
The units of the mole fractions are by default in `ppmv` (parts per million by volume), but it is possible to specify
different units by enclosing them in square brackets, e.g. 'H2[ppbv]' or `H2[1]`, where `1` is the unitless mole fraction.
The recognized units are: `ppmv`, `ppbv`, `ppm`, `ppb`, `mol/mol`, `1`

It is required that the atmospheric profile file contains all the gases that are defined in the opacity definition file.
The gases that are not present in the atmospheric profile are assumed to be radiatively inactive.
The mean molecular weight of the atmosphere is calculated based on the pressure and
altitude levels in the atmospheric profile.

Here is an example of an atmospheric profile file for Jupiter:
```text
# Sample atmospheric profile for Jupiter
IDX ALT[km] PRE[mbar] TEM[K] H2[ppmv] He[ppmv] CH4[ppmv]
1 -315.8 6289.6 294.2 861993 136001 1809
2 -310.2 5501.5 282.4 861993 136001 1809
```

You may add comments to the file by starting the line with `#`.
Any text after the `#` symbol wil be ignored.

We provide a python script `atm_profile_utils.py` that parses the atmospheric profile file and returns the data in a dictionary.
The dictionary contains the following keys:
- 'IDX': index of the layer
- 'ALT': altitude (m)
- 'PRE': pressure (Pa)
- 'TEM': temperature (K)
and additional keys for the names of the gases present in the file.
The entries are numpy arrays with the column data in SI units.
The script checks the units of the columns and converts them to SI units if necessary.
It will raise an error if the units are not recognized or the file format is incorrect.

The script can be used as follows:
```python
from atm_profile_utils import read_atm_profile
data = read_atm_profile('jupiter.atm')
```

Or you can use the script as a command line tool:
```bash
python3 atm_profile_utils.py -i [input_file] -o [output_file] -u [units]
```
where `[input_file]` is the name of the input atmospheric profile file, `[output_file]` is the name of the output file,
The format of the input file is determined by the extension of the file name.
Currently, the script supports `.atm`, `.txt` and `.nc` file formats.

## Generate line-by-line opacity using RFM
