======================================================
Matrix elements with FeynRules/Arts/Calc
======================================================

Here we present an alternative computation of the matrix elements, which we have
used as a cross check of the WallGo MatrixElements package.

Requirements
============
The matrix element calculations presented here use some packages written in
Mathematica. For more information on these packages, please see their
documentation.

- `FeynRules <https://feynrules.irmp.ucl.ac.be/>`_
- `FeynArts <https://feynarts.de/>`_
- `FeynCalc <https://feyncalc.github.io/>`_

Running the computation
=======================
The calculation proceeds in steps as follows. Here we give the example of a
Yukawa model in the directory `yukawa`.

1. Write a plain text FeynRules model file. The example here is called `yukawa.fr`.
2. Run `feynRules.m` in Mathematica, to compute the Feynman rules, and output for
FeynArts. These are files called `yukawa.mod` and `yukawa.gen`, as well as `yukawa.pars`
which holds the parameters of the model.
3. Run the `matrixElements.m` file to produce the required matrix elements. This file
consists of first constructing the diagrams and amplitudes using FeynArts, and then
computing and simplifying the matrix elements squared, using FeynCalc. Note that
you will need to quit the Mathematica kernel after running `feynRules.m`.