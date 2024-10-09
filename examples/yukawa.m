(* ::Package:: *)

Quit[];


SetDirectory[NotebookDirectory[]];
(*Put this if you want to create multiple model-files with the same kernel*)
$GroupMathMultipleModels=True;
$LoadGroupMath=True;
<<../src/WallGoMatrix.m


(* ::Chapter:: *)
(*Yukawa Model*)


(* ::Section:: *)
(*Model*)


Group={"U1"};
RepAdjoint={0};
RepScalar={{{0},"R"}};
CouplingName={g1};


RepFermion={{{0},"L"},{{0},"R"}};


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


(* \[Sigma] \[Phi] *)
InputInv={{1},{True}};
LinearTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VLinear=0;(* turned off *)\[Sigma] LinearTerm
\[Lambda]1=GradTadpole[VLinear];


(* 1/2m^2\[Phi]^2 *)
InputInv={{1,1},{True,True}}; 
MassTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VMass=msq/2 MassTerm
\[Mu]ij=GradMass[VMass]//SparseArray;


(* 1/6\[Gamma]\[Phi]^3 *)
InputInv={{1,1,1},{True,True,True}};
CubicTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VCubic=gamma/6 CubicTerm
\[Lambda]3=GradCubic[VCubic];


(* 1/24\[Lambda]\[Phi]^4 *)
InputInv={{1,1,1,1},{True,True,True,True}};
QuarticTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VQuartic=lam/24 MassTerm^2
\[Lambda]4=GradQuartic[VQuartic];


(* m(Subscript[\[Psi], R]Subscript[\[Psi], L]+Subscript[\[Psi]^+, L]Subscript[\[Psi], R]^+)*)
InputInv={{2,1},{True,True}}; (*Subscript[\[Psi], R]^+Subscript[\[Psi], L]*)
MassTerm1=CreateInvariantFermion[Group,RepFermion,InputInv][[1]]
InputInv={{1,2},{False,False}};  (*Subscript[\[Psi]^+, L]Subscript[\[Psi], R]*)
MassTerm2=CreateInvariantFermion[Group,RepFermion,InputInv][[1]]


\[Mu]IJ=mpsi*GradMassFermion[MassTerm1];
\[Mu]IJC=mpsi*GradMassFermion[MassTerm2];


(* y \[Phi](Subscript[\[Psi], R]Subscript[\[Psi], L]+Subscript[\[Psi]^+, L]Subscript[\[Psi], R]^+)*)
InputInv={{1,2,1},{True,True,True}};
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
InputInv={{1,1,2},{True,False,False}}; 
YukawaDoublet2=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;


Ysff= y*GradYukawa[YukawaDoublet1];
YsffC=y*GradYukawa[YukawaDoublet2];


ImportModel[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Section:: *)
(*A model with 1 Dirac fermions and 1 scalar bosons*)


(* ::Subsection:: *)
(*SymmetryBreaking*)


gvss//MatrixForm


WallGoMatrix`SimplifySparse


vev={v};
(*SymmetryBreaking[vev,VevDependentCouplings->True]*) (*uncomment if you want vev-dependent couplings*)
SymmetryBreaking[vev]


(* ::Subsection:: *)
(*UserInput*)


(*
In DRalgo fermions are Weyl.
So to create one Dirac we need
one left-handed and
one right-handed fermion
*)


(*Below
rep 1-2 are fermions,
(*rep 3 is a scalar*)
*)
(* scalar *)
RepScalar=CreateParticle[{1},"S"];

(* left-handed fermion *)
RepFermionL=CreateParticle[{1},"F"];

(* right-handed fermion *)
RepFermionR=CreateParticle[{2},"F"];

(*Vector bosons*)
RepZ=CreateParticle[{1},"V"];


CreateParticle[{1,2},"F"]


(*
These particles do not necessarily have to be out of equilibrium
the remainin particle content is set as light
*)
ParticleList={RepScalar,RepFermionL,RepFermionR};


(*Defining various masses and couplings*)


VectorMass=Table[mv2,{i,1,RepZ[[1]]//Length}];
FermionMass=Table[mf2,{i,1,Length[gvff[[1]]]}];
ScalarMass=Table[ms2,{i,1,Length[gvss[[1]]]}];
ParticleMasses={VectorMass,FermionMass,ScalarMass};
(*
up to the user to make sure that the same order is given in the python code
*)
UserMasses={ms2,mf2,mf2};
UserCouplings=Variables@Normal@{Ysff,gvss,gvff,gvvv,\[Lambda]4,\[Lambda]3}//DeleteDuplicates;


(*
	output of matrix elements
*)
SetDirectory[NotebookDirectory[]];
OutputFile="output/matrixElements.yukawa";
ParticleName={"Phi","PsiL","PsiR"};
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	UserMasses,
	UserCouplings,
	ParticleName,
	ParticleMasses,
	{TruncateAtLeadingLog->True,Format->{"json","txt"}}];


MatrixElements