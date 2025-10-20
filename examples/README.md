# Examples
This folder contains examples of model implementations in `WallGoMatrix`.
Instructions can be found in the [manual](https://arxiv.org/abs/2411.04970).

# Loading WallGoMatrix
The line 
    ```Get["WallGo`WallGoMatrix`"]```
in the first cell of each example loads the installed version of WallGoMatrix. See the [manual](https://arxiv.org/abs/2411.04970) or the [WallGoMatrix repo](https://github.com/Wall-Go/WallGoMatrix) for installation instructions.
If you want to run with the version of WallGoMatrix in the local folder (e.g. to test changes in the source code), use instead:
    ```Get["../Kernel/WallGoMatrix.m"]```

# Models implemented in WallGo
Some of the model files in this folder are related to models implemented in the [WallGo repository](https://github.com/Wall-Go/WallGo):
- `2hdm.m` gives the matrix elements for the Two Higgs Doublet Model
- `2scalars.m` gives the matrix elements for the pure scalar sector of the Two Higgs Doublet Model
- `SMFull.m` gives the matrix elements for the StandardModel
- `SMSimplified.m` gives the matrix elements for the StandardModel (only quarks + gauge bosons)
- `idm.m` gives the matrix elements for InertDoubletModel
- `qcd.m` gives the matrix elements for ManySinglets and SingletStandardModel_Z2
- `xsm.m` gives the matrix elements for SingletStandardModel_Z2 beyond the QCD sector
- `weak.m` gives the matrix elements for QCD and W-bosons
- `yukawa.m` gives the matrix elements for the Yukawa Model

The WallGo repository contains copies of these example files.