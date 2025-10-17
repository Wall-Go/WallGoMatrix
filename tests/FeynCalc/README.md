# Matrix Element Computation with `FeynRules`, `FeynArts`, and `FeynCalc`

This directory provides a lightweight workflow to cross-check matrix elements computed by the WallGoMatrix Elements package using Mathematica-based tools:

- `FeynRules`: model and Feynman rule generation
- `FeynArts`: diagram generation and amplitudes
- `FeynCalc`: algebra, simplification, and squared matrix elements

## Prerequisites

- [Mathematica](https://www.wolfram.com/mathematica/) (use `14.2` or earlier; `FeynRules` `fr-2.3.49` is not compatible with newer versions)
- [FeynRules](https://feynrules.irmp.ucl.ac.be/)
- [FeynArts](https://feynarts.de/)
- [FeynCalc](https://feyncalc.github.io/)

Refer to each package’s documentation for installation and setup.

## Directory layout

As an example, see the Yukawa model under `tests/FeynCalc/yukawa`:
- `yukawa.fr`: FeynRules model file
- `feynRules.m`: generates FeynArts files (`.mod`, `.gen`) and parameters (`.pars`)
- `yukawa.matrixElements.m`: builds diagrams (FeynArts) and computes squared MEs (FeynCalc)

## Quick start

1. Open `feynRules.m` in Mathematica and evaluate it to generate:
   - `<model>.mod`, `<model>.gen` (FeynArts model files, e.g., `u1-higgs-yukawa.mod`)
   - `<model>.pars` (model parameters)
2. Restart the Mathematica kernel before proceeding (avoids stale symbol/context issues).
3. Convert the FeynRules output to the conventions expected by this repo using the helper script (run it on `.mod` or `.gen` files, not on directories):
   - macOS Terminal (from the repository root):
     ```
     sh tests/FeynCalc/replaceFeynCalc.sh tests/FeynCalc/yukawa/u1-higgs-yukawa.mod
     ```
     - Optionally also run on the corresponding `.gen`:
       ```
       sh tests/FeynCalc/replaceFeynCalc.sh tests/FeynCalc/yukawa/u1-higgs-yukawa.gen
       ```
   - Notes:
     - The script updates files in place; consider backing up originals.
     - See the script header for usage details and optional arguments.
4. Open and run `yukawa.matrixElements.m` in Mathematica to generate diagrams (FeynArts) and compute/simplify squared matrix elements (FeynCalc).

## Helper: replaceFeynCalc.sh

Path: `tests/FeynCalc/replaceFeynCalc.sh`

Purpose: Post-processes FeynRules `.mod`/`.gen` exports so they load cleanly in our FeynArts/FeynCalc workflow (file naming and symbol/convention adjustments expected by the scripts in this directory).

Typical usage:
- From repo root:
  ```
  sh tests/FeynCalc/replaceFeynCalc.sh tests/FeynCalc/<model-dir>/<model>.mod

  sh tests/FeynCalc/replaceFeynCalc.sh tests/FeynCalc/<model-dir>/<model>.gen
  ```
- From within the model directory:
  ```
  sh ../../replaceFeynCalc.sh <model>.mod
  ```

## Compatibility

FeynRules `fr-2.3.49` is not compatible with Mathematica versions newer than `14.2`. Use Mathematica `14.2` or earlier.

For further details, consult the documentation of the respective packages.