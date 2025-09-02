# Matrix Element Computation with `FeynRules`, `FeynArts`, and `FeynCalc`

This document describes an alternative approach for computing matrix elements, used to cross-check results from the `WallGoMatrix` Elements package.

## Prerequisites

The calculations require the following `Mathematica` packages:

- [FeynRules](https://feynrules.irmp.ucl.ac.be/)
- [FeynArts](https://feynarts.de/)
- [FeynCalc](https://feyncalc.github.io/)

Refer to each package's documentation for installation and usage details.

## Workflow Overview

The computation proceeds in several steps. As an example, consider the Yukawa model in the `yukawa` directory:

1. **Model Definition**  
    Create a plain text `FeynRules` model file (e.g., `yukawa.fr`).

2. **Feynman Rule Generation**  
    Run `feynRules.m` in `Mathematica` to generate Feynman rules and export files for `FeynArts` (`yukawa.mod`, `yukawa.gen`). Model parameters are saved in `yukawa.pars`.

3. **Matrix Element Calculation**  
    Execute `yukawa.matrixElements.m` to construct diagrams and amplitudes with `FeynArts`, then compute and simplify squared matrix elements using `FeynCalc`.

**Note:**  
After running `feynRules.m`, restart the `Mathematica` kernel before proceeding to matrix element calculations.

## Compatibility

As of now, `FeynRules` `fr-2.3.49` is not compatible with `Mathematica` versions newer than `14.2`. Use `Mathematica` `14.2` or earlier for reliable results.

---
For further details, consult the documentation of the respective packages.