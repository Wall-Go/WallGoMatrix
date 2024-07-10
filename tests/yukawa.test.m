(* ::Package:: *)

Quit[];


(* Check Mathematica version *)
If[$VersionNumber < 13.3,
  Print["The Mathematica testing framework requires Mathematica version ", requiredVersion," or higher. You are using version ", currentVersion, "."];
  Abort[]
];

SetDirectory[NotebookDirectory[]];
(*Put this if you want to create multiple model-files with the same kernel*)
$GroupMathMultipleModels=True;
$LoadGroupMath=True;
<<../DRalgo/DRalgo.m
<<../src/matrixElements.m


(* ::Chapter:: *)
(*Yukawa Model*)


(* ::Section::Closed:: *)
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
VCubic=\[Gamma]/6 CubicTerm
\[Lambda]3=GradCubic[VCubic];


(* 1/24\[Lambda]\[Phi]^4 *)
InputInv={{1,1,1,1},{True,True,True,True}};
QuarticTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VQuartic=\[Lambda]/24 MassTerm^2
\[Lambda]4=GradQuartic[VQuartic];


(* m(Subscript[\[Psi], R]Subscript[\[Psi], L]+Subscript[\[Psi]^+, L]Subscript[\[Psi], R]^+)*)
InputInv={{2,1},{True,True}}; (*Subscript[\[Psi], R]^+Subscript[\[Psi], L]*)
MassTerm1=CreateInvariantFermion[Group,RepFermion,InputInv][[1]]
InputInv={{1,2},{False,False}};  (*Subscript[\[Psi]^+, L]Subscript[\[Psi], R]*)
MassTerm2=CreateInvariantFermion[Group,RepFermion,InputInv][[1]]


\[Mu]IJ=m\[Psi]*GradMassFermion[MassTerm1];
\[Mu]IJC=m\[Psi]*GradMassFermion[MassTerm2];


(* y \[Phi](Subscript[\[Psi], R]Subscript[\[Psi], L]+Subscript[\[Psi]^+, L]Subscript[\[Psi], R]^+)*)
 InputInv={{1,2,1},{True,True,True}};
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
InputInv={{1,1,2},{True,False,False}}; 
YukawaDoublet2=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;


Ysff= y*GradYukawa[YukawaDoublet1];
YsffC=y*GradYukawa[YukawaDoublet2];


(* ::Section:: *)
(*User Input*)


(* ::Subsection:: *)
(*SymmetryBreaking*)


vev={v};
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
RepScalar=CreateOutOfEq[{1},"S"];

(* left-handed fermion *)
RepFermionL=CreateOutOfEq[{1},"F"];

(* right-handed fermion *)
RepFermionR=CreateOutOfEq[{2},"F"];

(*Vector bosons*)
RepZ=CreateOutOfEq[{1},"V"];


(*
These particles do not necessarily have to be out of equilibrium
the remainin particle content is set as light
*)
ParticleList={RepScalar,RepFermionL,RepFermionR};


(*Defining various masses and couplings*)


VectorMass=Table[mv,{i,1,RepZ[[1]]//Length}];
FermionMass=Table[mf,{i,1,Length[gvff[[1]]]}];
ScalarMass=Table[ms,{i,1,Length[gvss[[1]]]}];
ParticleMasses={VectorMass,FermionMass,ScalarMass};
(*
up to the user to make sure that the same order is given in the python code
*)
UserMasses={ms,mf,mf};
UserCouplings={CouplingName,\[Lambda],\[Gamma],y}//Flatten;


(*
	output of matrix elements
*)
OutputFile="matrixElements.yukawa";
SetDirectory[NotebookDirectory[]];
ParticleName={"Phi","PsiL","PsiR"};
RepOptional={};
MatrixElements=ExportMatrixElements[OutputFile,ParticleList,UserMasses,UserCouplings,ParticleName,ParticleMasses,RepOptional];


MatrixElements


(* ::Section:: *)
(*Tests*)


(* ::Subsection:: *)
(*Test hard*)


testList={};


(* scalar-scalar scattering*)
TestCreate[M[0,0,0,0]/.MatrixElements/.msq[i_]->0,
	(*(c[1] + c[2]^2 (1/s + 1/t + 1/u))^2*) (* all channels *)
	(c[2]^2 (1/t + 1/u))^2 (* just IR sensitive terms *)
];


(* scalar to fermions *)
AppendTo[testList,
TestCreate[M[0,0,1,1]+M[0,0,2,2]+M[0,0,1,2]/.MatrixElements/.msq[i_]->0,
	(*c[3]^4 (2 t/u + 2 u/t - 4) + c[2]^2 c[3]^2 (2/s)*)(* all channels *)
	 c[3]^4 (2 t/u + 2 u/t) (* just IR sensitive terms *)
]];


(* fermions to scalar *)
AppendTo[testList,
TestCreate[M[1,1,0,0]+M[1,2,0,0]+M[2,1,0,0]+M[2,2,0,0]/.MatrixElements/.msq[i_]->0,
	(*c[3]^4 (4 t/u + 4 u/t - 8) + c[2]^2 c[3]^2 (4/s)*)(* all channels *)
	 c[3]^4 (4 t/u + 4 u/t)(* just IR sensitive terms *)
]];


(* scalar-fermion scattering *)
AppendTo[testList,
TestCreate[M[0,1,0,1]+M[0,2,0,2]/.MatrixElements/.msq[i_]->0,
	(*c[3]^4 (4 s/u + 4 u/s - 8) + c[2]^2 c[3]^2 (4/t)*)(* all channels *)
	 c[3]^4 (4 s/u) + c[2]^2 c[3]^2 (4/t)(* just IR sensitive terms *)
]];


(* fermion-scalar scattering *)
AppendTo[testList,
TestCreate[M[1,0,0,1]+M[2,0,0,2]/.MatrixElements/.msq[i_]->0,
	(*c[3]^4 (4 s/t + 4 t/s - 8) + c[2]^2 c[3]^2 (4/u)*)(* all channels *)
	 c[3]^4 (4 s/t + 4 t/s) + c[2]^2 c[3]^2 (4/u)(* just IR sensitive terms *)
]];


(* fermion-fermion scattering*)
AppendTo[testList,
TestCreate[M[1,2,1,2]+M[2,1,1,2]/.MatrixElements/.msq[i_]->0,
	(*c[3]^4 (24)*)(* all channels *)
	0 (* just IR sensitive terms *)
]];


report=TestReport[testList]
report["ResultsDataset"]





