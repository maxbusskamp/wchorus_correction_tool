# WCHORUS Correction Tool

This tool can be used to generate corrected wurst shape files for the use in the CHORUS sequence. It generates shapefiles for the use with Simpson and Bruker Topspin.

It is part of my doctoral thesis.

It is recommended to run the script on a Linux machine. Windows is not fully tested, but should work if SIMPSON was correctly installed.

## Dependencies:
- SIMPSON [Download](https://inano.au.dk/about/research-centers-and-projects/nmr/software/simpson)
- Python >3.7
- PySimpleGUI
- Anaconda Distribution, or:
  - os
  - subprocess
  - matplotlib
  - numpy
  - scipy

## Installation:
Install the Anaconda Distribution, or any Python installation of your choosing. Check if the packages name in the dependencies are available in your Python shell.
Install Simpson by unpacking the archive downloaded from the *Danish Center for Ultrahigh Field NMR Spectroscopy* website and running the `install.sh` as sudo.

## Usage
The WCHORUS Correction Tool can be started by running `python wchorus_correction_tool.py` or, ddependding on your system `python3 wchorus_correction_tool.py`. The tool can also be accessed in other folders, by supplying the absolute path.

Afterwards the type of pulse and experiment to simulate are chosen and the simulation parameters can be specified.
The results of the simulations are shown and saved to shape files for SIMPSON and TopSpin. Logs files containing all simulation parameters and the produced figure are also saved. All files are encoded with a unique string, to allow the simulation of multiple parameters in quick succession.

# License and Copyright Info

The WCHORUS Correction Tool is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
