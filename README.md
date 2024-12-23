<img src="https://raw.githubusercontent.com/Wall-Go/WallGo/refs/heads/main/docs/source/figures/wallgo.svg" alt="WallGoLogo" width="200"/>

# WallGoMatrix

**WallGoMatrix** is a Mathematica package for computing 2-to-2 scattering matrix
elements for arbitrary quantum field theories. It is part of the WallGo project
for computing the bubble wall speed for cosmological phase transitions.

## Status

[![Version](https://img.shields.io/github/v/tag/Wall-Go/WallGoMatrix?label=Version)](https://github.com/Wall-Go/WallGoMatrix/releases/latest/)

## About this project

[**WallGo**](https://github.com/Wall-Go) is an open source code for the
computation of the bubble wall velocity and bubble wall width in first-order
cosmological phase transitions. There is a main Python package together with
two subsidiary software packages. The physical and mathematical details are explained in [the associated paper](https://arxiv.org/abs/2411.04970).

- [**WallGo**](https://github.com/Wall-Go/WallGo) is the name of the main Python
package. This determines the wall velocity and width by
solving the scalar field(s) equation of motion, the Boltzmann equations and
energy-momentum conservation for the fluid velocity and temperature.
- [**WallGoMatrix**](https://github.com/Wall-Go/WallGoMatrix) computes the relevant matrix elements for the
out-of-equilibrium particles, and is written in Mathematica.
It builds on existing Mathematica packages [DRalgo](https://github.com/DR-algo/DRalgo) and [GroupMath](https://renatofonseca.net/groupmath).
- [**WallGoCollision**](https://github.com/Wall-Go/WallGoCollision) performs the higher-dimensional integrals to obtain the
 collision terms in the Boltzmann equations, and is written in C++. It also
 has Python bindings so that it can be called directly from Python, but
 still benefits from the speedup from compiled C++ code.

Users can implement their own models.
For WallGoMatrix, the model building
routines are taken from DRalgo and involve constructing coupling tensors.
The WallGoMatrix routines then contract these coupling tensors with
the appropriate kinematic factors to yield the 2-to-2 scattering matrix elements.

## Installation

### Paclet Installation

**WallGoMatrix** can be installed as a Wolfram Paclet by running either of the following two commands in Mathematica
#### [Wolfram repository](https://resources.wolframcloud.com/PacletRepository/resources/WallGo/WallGoMatrix/)
```mathematica
PacletInstall["WallGo/WallGoMatrix"];
```
#### GitHub repository
```mathematica
PacletInstall["https://github.com/Wall-Go/WallGoMatrix/releases/latest/download/WallGoMatrix.paclet"];
```
Note that the dependencies of **WallGoMatrix** must also be installed, for which see the [Requirements section](#requirements) below.

### Manual Installation

**WallGoMatrix** can alternatively be installed manually by downloading the
zip file from the following link:

[![Download zip](https://custom-icon-badges.demolab.com/badge/-Download-blue?style=for-the-badge&logo=download&logoColor=white "Download zip")](https://github.com/Wall-Go/WallGoMatrix/releases/latest/download/WallGoMatrix.zip)

After unpacking the zip file, place the **WallGoMatrix** directory inside the **Applications** folder within either the base or user-specific **Mathematica Applications** directory. These directories store Mathematica packages and can be located by evaluating the variables **`$BaseDirectory`** and **`$UserBaseDirectory`** within a Mathematica session.
To find these directories, you can run the following command in Mathematica:

```mathematica
Print["Base Directory: ", FileNameJoin[{$BaseDirectory, "Applications"}]]
Print["User Base Directory: ", FileNameJoin[{$UserBaseDirectory, "Applications"}]]
```

For versions of Mathematica before 14.1, the paths are commonly the
following ones

**Linux**
- /usr/share/Mathematica/Applications
- ~/.Mathematica/Applications

**macOS**
- ~/Library/Mathematica/Applications

**Windows**
- C:\ProgramData\Mathematica\Applications 
- C:\Users\(computer name)\AppData\Roaming\Mathematica\Applications

Note that from Mathematica 14.1 onwards,
these directories are renamed to contain Wolfram instead of Mathematica.
For example, the path for macOS becomes `~/Library/Wolfram/Applications`.

### Requirements

WallGoMatrix is written in the Wolfram Mathematica language, and depends on the Mathematica package GroupMath.
It has been tested on the following versions.

- [Mathematica](https://www.wolfram.com/mathematica/) versions 12.x, 13.x and 14.x
    - [GroupMath](https://renatofonseca.net/groupmath) version 1.1.2

GroupMath can be either manually obtained from the link above or by setting the following flag before loading WallGoMatrix in Mathematica
```mathematica
WallGo`WallGoMatrix`$InstallGroupMath=True;
```
WallGoMatrix builds on [DRalgo](https://github.com/DR-algo/DRalgo) version 1.2, but the required elements are included directly in the WallGoMatrix package, so separate installation of DRalgo is not necessary.

## Running

Once the WallGoMatrix directory has been installed, the package can be loaded from within Mathematica using
```mathematica
<<WallGo`WallGoMatrix`
```

To see how WallGoMatrix is used in practice, we recommend taking a look at the
[examples](https://github.com/Wall-Go/WallGoMatrix/tree/main/examples).

### Running the examples

Within **WallGo**, **WallGoMatrix** can be executed using 
[Wolframscript](https://www.wolfram.com/wolframscript/).
Wolframscript provides the core computational capabilities of Wolfram Mathematica and allows Wolfram Language scripts to be run without needing a full Mathematica installation.
To run the example files, you can use the following command:

    $ wolframscript -file examples/qcd.m 

One requirement for the above command is an active `WolframKernel`.


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
