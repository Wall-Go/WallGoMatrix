![WallGo](wallgo.svg)

# WallGoMatrix

**WallGoMatrix** is a Mathematica package for computing 2-to-2 scattering matrix
elements for arbitrary quantum field theories. It is part of the WallGo project
for computing the bubble wall speed for cosmological phase transitions.

## Status

**Version** v0.1.0 (alpha)

## About this project

**WallGo** is an open source code for the computation of the bubble wall
velocity and bubble wall width in first-order cosmological phase transitions.
The main WallGo Python package determines the wall velocity and width by
solving the scalar field(s) equation of motion, the Boltzmann equations and
energy-momentum conservation for the fluid velocity and temperature.

**WallGo** is accompanied by two subsidiary software packages:
- **WallGoMatrix** computes the relevant matrix elements for the
out-of-equilibrium particles, and is written in Mathematica.
It builds on existing Mathematica packages DRalgo and GroupMath.
- **WallGoCollision** performs the higher-dimensional integrals to obtain the
 collision terms in the Boltzmann equations, and is written in C++. It also
 has Python bindings so that it can be called directly from Python, but
 still benefits from the speedup from compiled C++ code.

Users can implement their own models. For WallGoMatrix, the model building
routines are taken from DRalgo, and involve constructing coupling tensors.
The WallGoMatrix routines then contract these coupling tensors with kinematic
factors appropriately to yield the 2-to-2 scattering matrix elements.

## Installation

**WallGoMatrix** can be installed as a Mathematica package by downloading the
package

![Download](download.svg)

The WallGoMatrix directory should be put inside either the base or the user
Mathematica Applications directory. The locations of these directories can be found
by inspecting **$BaseDirectory** or **$UserBaseDirectory** within Mathematica.
Often the paths are the following ones

**Linux**
- /usr/share/Mathematica/Applications
- ~/.Mathematica/Applications

**macOS**
- ~/Library/Mathematica/Applications/GroupMath

**Windows**
- C:\ProgramData\Mathematica\Applications 
- C:\Users\(computer name)\AppData\Roaming\Mathematica\Applications 

Once the WallGoMatrix directory has been put inside the Mathematica applications
directory, the package can be loaded from within Mathematica using

    <<WallGoMatrix`

To see how WallGoMatrix is used in practice, we recommend taking a look at the
[examples](https://github.com/Wall-Go/WallGoMatrix/tree/main/examples).

### Requirements

WallGoMatrix is written in Wolfram Mathematica, and depends on the Mathematica package GroupMath. It has been tested on the following versions.

- [Mathematica](https://www.wolfram.com/mathematica/) versions 12.x, 13.x and 14.0
    - [GroupMath](https://renatofonseca.net/groupmath) version 1.1.2

WallGoMatrix builds on [DRalgo](https://github.com/DR-algo/DRalgo) version 1.2, but the required elements are included directly in the WallGoMatrix package, so separate installation of DRalgo is not necessary.


## License

Copyright (c) 2024 Andreas Ekstedt, Oliver Gould, Joonas Hirvonen,
Benoit Laurent, Lauri Niemi, Philipp Schicho, and Jorinde van de Vis.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
