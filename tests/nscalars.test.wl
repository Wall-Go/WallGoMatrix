(* ::Package:: *)

(*Quit[];*)


If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
(*Put this if you want to create multiple model-files with the same kernel*)
WallGo`WallGoMatrix`$GroupMathMultipleModels=True;
WallGo`WallGoMatrix`$LoadGroupMath=True;
Check[
    Get["../Kernel/WallGoMatrix.m"],
    Message[Get::noopen, "WallGo`WallGoMatrix` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*N scalar Model*)


(* ::Section:: *)
(*Model*)


(* number of scalars *)
n=3;


Group={"U1"};
RepAdjoint={0};
RepScalar=Table[{{0},"R"},n];
RepFermion={};
CouplingName={g1};


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


(* 1/2 msq \[Phi]^2 *)
(*VMass=Sum[1/2!msq[Sort[{i,j}]/.List->Sequence]CreateInvariant[Group,RepScalar,{{i,j},{True,True,True}}][[1]],{i,n},{j,n}];
\[Mu]ij=GradMass[VMass[[1]]]//Simplify//SparseArray;*)


(* 1/6 g \[Phi]^3 *)
VCubic=Sum[1/3! g[Sort[{i,j,k}]/.List->Sequence]CreateInvariant[Group,RepScalar,{{i,j,k},{True,True,True}}][[1]],{i,n},{j,n},{k,n}];
\[Lambda]3=GradCubic[VCubic];


(* 1/24 lam \[Phi]^3 *)
VQuartic=Sum[1/4! lam[Sort[{i,j,k,l}]/.List->Sequence]CreateInvariant[Group,RepScalar,{{i,j,k,l},{True,True,True,True}}][[1]],{i,n},{j,n},{k,n},{l,n}];
\[Lambda]4=GradQuartic[VQuartic];


ImportModel[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Section:: *)
(*User Input*)


(* ::Subsection:: *)
(*UserInput*)


(*
In DRalgo fermions are Weyl.
So to create one Dirac we need
one left-handed and
one right-handed fermion
*)


(*Below
rep 1 is a scalar
rep 2-3 are fermions,
*)
(* scalar *)
RepScalars=Table[
	CreateParticle[{i}, "S", ms, StringJoin[{"Phi",ToString[i]}]],
	{i,n}
]
(*RepScalars=CreateParticle[{1,2},"S",ms,"Phi"]*)

(* fermion *)
RepFermion={};

(*Vector bosons*)
RepZ=CreateParticle[{1}, "V", mv, "LightParticle"]


(*
These particles do not necessarily have to be out of equilibrium
RepZ can later be assumed to be a light particle
*)
ParticleList=RepScalars;
LightParticleList={RepZ};


(*Defining various masses and couplings*)
UserMasses={ms,ms,mv};


(*
	output of matrix elements
*)
OutputFile="nscalar.test";
SetDirectory[NotebookDirectory[]];
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{
		TruncateAtLeadingLog->False,
		Format->{"json","txt"}
	}
];


MatrixElements


(* ::Chapter:: *)
(*Tests*)


(* ::Section:: *)
(*Importing results from WallGo*)


If[Not[ValueQ[MatrixElements]],
{particleNames,parameters,MatrixElements}=ImportMatrixElements["nscalars.test.json"]
];


(* ::Section:: *)
(*Importing results from FeynCalc*)


file=FileNameJoin[{NotebookDirectory[],"nscalars.feyncalc.test.json"}];
{particleNamesFeyn,parametersFeyn,MatrixElementsFeyn}=ImportMatrixElements[file];


(* ::Section:: *)
(*Comparison tests*)


insertCouplings={};


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})


UserMasses={ms,mv};
fixConvention[arg_]:=symmetriseTU[arg/.Thread[UserMasses->0]/.{s->(-t-u)}/.insertCouplings]//Expand//Simplify//Expand


removeMissing[arg_]:=arg/.M[__]->0/.Missing["KeyAbsent", _]->0


(* ::Subsection:: *)
(*Test hard*)


testList={};


(* the whole thing scattering*)
test["WallGo"]=Sum[M[a,b,c,d]/.MatrixElements//fixConvention//removeMissing,{a,0,n-1},{b,0,n-1},{c,0,n-1},{d,0,n-1}]
test["Feyn"]=Sum[M[a,b,c,d]/.MatrixElementsFeyn//fixConvention//removeMissing,{a,0,n-1},{b,0,n-1},{c,0,n-1},{d,0,n-1}];
Simplify[test["WallGo"]-test["Feyn"]]
AppendTo[testList,
	TestCreate[test["WallGo"],test["Feyn"]]
];


report=TestReport[testList]
report["ResultsDataset"]





