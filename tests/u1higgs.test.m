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
<<../src/WallGoMatrix.m


(* ::Chapter:: *)
(*Abelian Higgs*)


(* ::Section:: *)
(*Model*)


Group={"U1"};
CouplingName={g1};
RepAdjoint={0};
Higgs={{1},"C"}; (* Y\[Phi] = 1 *)
RepScalar={Higgs};


RepFermion={};


(* ::Text:: *)
(*The first element is the vector self-interaction matrix:*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


InputInv={{1,1},{True,False}}; (*This specifies that we want a \[Phi]^+\[Phi] term*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;


VMass=msq*MassTerm1[[1]];(*This is the \[Phi]^+\[Phi] term written in component form*)


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


QuarticTerm1=MassTerm1[[1]]^2; (*Because MassTerm1=\[Phi]^+\[Phi], we can write (\[Phi]^+\[Phi])^2=MassTerm1^2*)


VQuartic=\[Lambda]*QuarticTerm1;


\[Lambda]4=GradQuartic[VQuartic];


ImportModel[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Section:: *)
(*User Input*)


(* ::Subsection:: *)
(*SymmetryBreaking*)


(*vev={v,0};
SymmetryBreaking[vev]*)


(* ::Subsection:: *)
(*UserInput*)


(* vector *)
RepVector=CreateParticle[{1},"V"];

(* fermion *)
RepFermion={};

(* scalar *)
RepScalar=CreateParticle[{1},"S"];



(*
These particles do not necessarily have to be out of equilibrium
the remainin particle content is set as light
*)
ParticleList={RepScalar,RepVector};


(*Defining various masses and couplings*)


VectorMass=Table[mv,{i,1,Length[RepVector[[1]]]}];
FermionMass={};
ScalarMass=Table[ms,{i,1,Length[RepScalar[[1]]]}];
ParticleMasses={VectorMass,FermionMass,ScalarMass};
(*
up to the user to make sure that the same order is given in the python code
*)
UserMasses={ms,ms,mv};
UserCouplings={CouplingName,\[Lambda]}//Flatten;


(*
	output of matrix elements
*)
SetDirectory[NotebookDirectory[]];
OutputFile="matrixElements.ah";
ParticleName={"Phi","Vector"};
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	UserMasses,
	UserCouplings,
	ParticleName,
	ParticleMasses,
	{TruncateAtLeadingLog->False,Format->{"json","txt"}}];


MatrixElements


(* ::Section:: *)
(*Tests*)


file=FileNameJoin[{NotebookDirectory[],"u1higgs.test.json"}];
{particleNames,parameters,FeynMatrixElements}=ImportMatrixElements[file];


insertCouplings={g1->e,\[Lambda]->lam};


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})


fixConvention[arg_]:=symmetriseTU[arg/.Thread[UserMasses->0]/.{s->(-t-u)}/.insertCouplings]//Expand//Simplify//Expand


removeMissing[arg_]:=arg/.M[__]->0/.Missing["KeyAbsent", _]->0


(* ::Subsection:: *)
(*Test hard*)


testList={};


(* scalar-scalar scattering*)
AppendTo[testList,
TestCreate[
	M[0,0,0,0]/.MatrixElements//fixConvention//removeMissing//Simplify,
	1/2 Sum[M[a,b,c,d],{a,0,1},{b,0,1},{c,0,1},{d,0,1}]/.FeynMatrixElements//fixConvention//removeMissing//Simplify(* explicit 1/2 is due to average over leg 1 *)
]];


(* scalar to vector *)
AppendTo[testList,
TestCreate[
	M[0,0,1,1]/.MatrixElements//fixConvention//removeMissing,
	1/2 Sum[M[a,b,2,2],{a,0,1},{b,0,1}]/.FeynMatrixElements//fixConvention//removeMissing(* explicit 1/2 is due to average over leg 1 *)
]];


(* vector to scalar *)
AppendTo[testList,
TestCreate[
	M[1,1,0,0]/.MatrixElements//fixConvention//removeMissing,
	1/2 Sum[M[2,2,a,b],{a,0,1},{b,0,1}]/.FeynMatrixElements//fixConvention//removeMissing (* explicit 1/2 is due to average over leg 1 *)
]];


(* scalar-vector scattering *)
AppendTo[testList,
TestCreate[
	M[0,1,1,0]+M[0,1,0,1]/.MatrixElements//fixConvention//removeMissing,
	1/2 Sum[M[a,2,2,b]+M[a,2,b,2],{a,0,1},{b,0,1}]/.FeynMatrixElements//fixConvention//removeMissing(* explicit 1/2 is due to average over leg 1 *)
]];


(* vector-scalar scattering *)
AppendTo[testList,
TestCreate[
	M[1,0,1,0]+M[1,0,0,1]/.MatrixElements//fixConvention//removeMissing,
	1/2 Sum[M[2,a,2,b]+M[2,a,b,2],{a,0,1},{b,0,1}]/.FeynMatrixElements//fixConvention//removeMissing (* explicit 1/2 is due to average over leg 1 *)
]];


(* vector-vector scattering*)
AppendTo[testList,
TestCreate[
	M[1,1,1,1]/.MatrixElements//fixConvention//removeMissing,
	1/2 M[2,2,2,2]/.FeynMatrixElements//fixConvention//removeMissing (* explicit 1/2 is due to average over leg 1 *)
]];


report=TestReport[testList]
report["ResultsDataset"]




