======================================================
Yukawa model matrix elements with FeynRules/Arts/Calc
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
The calculation proceeds in steps as follows.

1. Write a plain text FeynRules model file. The example here is called `yukawa.fr`.
2. Run `feynRules.nb` in Mathematica, to compute the Feynman rules, and output for
FeynArts. These are files called `yukawa.mod` and `yukawa.gen`, as well as `yukawa.pars`
which holds the parameters of the model.
3. Run the `matrixElements.nb` file to produce the required matrix elements. This
consists of first constructing the diagrams and amplitudes using FeynArts, and then
computing and simplifying the matrix elements squared, using FeynCalc.