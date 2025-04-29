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
(*Abelian-Higgs-Yukawa Model*)


(* ::Section:: *)
(*Model*)


Group={"SU2","U1"}; (* I've added the U1 field because without it the Higgs is only given 2 degrees of freedom *)
CouplingName={g,gu1};
RepAdjoint={{2},0};
Higgs1={{{1},0},"C"}; (* fundamental *)
RepScalar={Higgs1};


Rep1={{{1},0},"L"};
Rep2={{{0},0},"R"};
RepFermion={Rep1,Rep2};


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


InputInv={{1,1},{True,False}}; (*This specifies that we want a \[Phi]^+\[Phi] term*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;


VMass=msq*MassTerm1[[1]];(*This is the \[Phi]^+\[Phi] term written in component form*)


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


QuarticTerm1=MassTerm1[[1]]^2; (*Because MassTerm1=\[Phi]^+\[Phi], we can write (\[Phi]^+\[Phi])^2=MassTerm1^2*)


VQuartic=\[Lambda]*QuarticTerm1;


\[Lambda]4=GradQuartic[VQuartic];


(* y \[Phi]^*(Subscript[\[Psi], L]Subscript[\[Xi], R]^++Subscript[\[Psi], R]Subscript[\[Xi], L]^+)*)
InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv]//Simplify;


Ysff=-GradYukawa[yt*YukawaDoublet[[1]]];
YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];


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


vev={0,0,0,v};
SymmetryBreaking[vev]


(*Below
rep 1-2 are fermions,
(*rep 3 is a scalar*)
*)
(* scalar *)
RepPhi=CreateParticle[{{1,1},{1,2}},"S",ms,"Phi"]

(* left-handed fermion *)
RepPsi=CreateParticle[{{1,1},{2,1}},"F",mf,"Psi"]
RepXi=CreateParticle[{{1,2}},"F",mf,"Xi"]

(*Vector bosons*)
RepA=CreateParticle[{1},"V",mv,"VectorSU2"]
RepB=CreateParticle[{2},"V",mv,"VectorU1"]


(*
These particles do not necessarily have to be out of equilibrium
*)
ParticleList={RepPhi,RepPsi,RepXi,RepA,RepB};
LightParticleList={};


(*Defining various masses and couplings*)
UserMasses={mv,mf,ms};


(*
	output of matrix elements
*)
OutputFile="output/matrixElements.su2_higgs_yukawa";
SetDirectory[NotebookDirectory[]];
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{
		TruncateAtLeadingLog->False,
		NormalizeWithDOF->False,
		Format->{"json","txt"}}];


MatrixElements;


(* ::Section:: *)
(*Tests*)


(* ::Section:: *)
(*Importing results from WallGo*)


{particles,parameters,MatrixElements}=ImportMatrixElements["output/matrixElements.su2_higgs_yukawa.json"];


(* ::Section:: *)
(*Importing results from FeynCalc*)


{particlesFeyn,parametersFeyn,MatrixElementsFeyn}=ImportMatrixElements["sunhiggs.test.json"];


(* ::Section:: *)
(*Comparison tests*)


insertCouplings={Global`g->g,\[Lambda]->lam,SUNN->2,gu1->0};


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})


fixConvention[arg_]:=symmetriseTU[arg/.Thread[UserMasses->0]/.{s->(-t-u)}/.insertCouplings]//Expand//Simplify//Expand


removeMissing[arg_]:=arg/.M[__]->0/.Missing["KeyAbsent", _]->0


testsRulesWallGo[arg_]:=arg/.Flatten[particles]/.MatrixElements//fixConvention//removeMissing;
testsRulesFeynCalc[arg_]:=arg/.Flatten[particlesFeyn]/.MatrixElementsFeyn//fixConvention//removeMissing
testWallGo[particlesA_,particlesB_,particlesC_,particlesD_]:=Sum[
	M[a,b,c,d],{a,particlesA},{b,particlesB},{c,particlesC},{d,particlesD}
]//testsRulesWallGo
testFeynCalc[particlesA_,particlesB_,particlesC_,particlesD_]:=Sum[
	M[a,b,c,d],{a,particlesA},{b,particlesB},{c,particlesC},{d,particlesD}
]//testsRulesFeynCalc


(* ::Subsection:: *)
(*Test hard*)


particles
particlesFeyn


testList={};


(* scalar-scalar scattering*)
test["WallGo"][1]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"Phi"},
	{"Phi"}
]
test["FeynCalc"][1]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Phi","Phibar"}
]
AppendTo[testList,TestCreate[test["WallGo"][1],test["FeynCalc"][1],TestID->"WallGo vs FeynCalc: SS->SS"]];


(* scalar to vector *)
test["WallGo"][2]=testWallGo[
	{"Psi"},
	{"Psi"},
	{"VectorSU2"},
	{"VectorSU2"}
]
test["FeynCalc"][2]=testFeynCalc[
	{"Psi","Psibar"},
	{"Psi","Psibar"},
	{"A"},
	{"A"}
]
AppendTo[testList,TestCreate[test["WallGo"][2],test["FeynCalc"][2],TestID->"WallGo vs FeynCalc: SS->VV"]];


report=TestReport[testList]
report["ResultsDataset"]



