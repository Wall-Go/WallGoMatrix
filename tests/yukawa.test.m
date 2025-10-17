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


(* ::Title:: *)
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


ImportModel[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


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
rep 1 is a scalar
rep 2-3 are fermions,
*)
(* scalar *)
RepScalar=CreateParticle[{1},"S",ms,"Phi"];

(* left-handed fermion *)
RepFermionL=CreateParticle[{1},"F",mf,"PsiL"];

(* right-handed fermion *)
RepFermionR=CreateParticle[{2},"F",mf,"PsiR"];

(*Vector bosons*)
RepZ=CreateParticle[{1},"V",mv,"LightParticle"];


(*
These particles do not necessarily have to be out of equilibrium
RepZ can later be assumed to be a light particle
*)
ParticleList={RepScalar,RepFermionL,RepFermionR};
LightParticleList={RepZ};


(*Defining various masses and couplings*)
UserMasses={ms,mf};


(*
	output of matrix elements
*)
OutputFile="output/yukawa.test";
SetDirectory[NotebookDirectory[]];
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{
		TruncateAtLeadingLog->False,
		Format->{"json","txt"}}];


MatrixElements


(* ::Chapter:: *)
(*Tests*)


(* ::Section:: *)
(*Importing results from WallGo*)


{particleNames,parameters,MatrixElements}=ImportMatrixElements["output/yukawa.test.json"];


(* ::Section:: *)
(*Importing results from FeynCalc*)


file=FileNameJoin[{NotebookDirectory[],"testFiles/yukawa.feyncalc.test.json"}];
{particleNamesFeyn,parametersFeyn,MatrixElementsFeyn}=ImportMatrixElements[file];


(* ::Section:: *)
(*Comparison tests*)


(* ::Subsection::Closed:: *)
(*Translate input*)


insertCouplings={\[Lambda]->lam,\[Gamma]->(g+lam v),y->y};
customCouplings={ms2->mPhi^2};


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})


fixConvention[arg_]:=symmetriseTU[
	arg/.customCouplings/.Thread[UserMasses->0]/.insertCouplings/.v->0/.{s->(-t-u)}
	]//Simplify//Expand


removeMissing[arg_]:=arg/.M[__]->0/.Missing["KeyAbsent", _]->0


testsRulesWallGo[arg_]:=arg/.Flatten[particleNames]/.MatrixElements//fixConvention//removeMissing;
testsRulesFeynCalc[arg_]:=arg/.Flatten[particleNamesFeyn]/.MatrixElementsFeyn//fixConvention//removeMissing;
testWallGo[particlesA_,particlesB_,particlesC_,particlesD_]:=Sum[
	M[a,b,c,d],{a,particlesA},{b,particlesB},{c,particlesC},{d,particlesD}
]//testsRulesWallGo
testFeynCalc[particlesA_,particlesB_,particlesC_,particlesD_]:=Sum[
	M[a,b,c,d],{a,particlesA},{b,particlesB},{c,particlesC},{d,particlesD}
]//testsRulesFeynCalc


(* ::Subsection::Closed:: *)
(*Initialize tests*)


particleNames
particleNamesFeyn


testList={};


(* ::Subsubsection::Closed:: *)
(*SStoSS*)


(* scalar-scalar scattering*)
process="SS->SS"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"Phi"},
	{"Phi"}
]//Simplify
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Phi","Phibar"}
]//Simplify
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*SStoFF*)


(* scalar to fermions *)
process="SS->FF"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"PsiL","PsiR"},
	{"PsiL","PsiR"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Psi","Psibar"},
	{"Psi","Psibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*FFtoSS*)


(* fermions to scalar *)
process="FF->SS"
test["WallGo"][process]=testWallGo[
	{"PsiL","PsiR"},
	{"PsiL","PsiR"},
	{"Phi"},
	{"Phi"}
]
(* explicit 1/2 is due to average over leg 1 *)
test["FeynCalc"][process]=1/2*testFeynCalc[
	{"Psi","Psibar"},
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"Phi","Phibar"}
]//Expand
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*SFtoSF*)


(* scalar-fermion scattering *)
process="SF->SF"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"PsiL","PsiR"},
	{"Phi"},
	{"PsiL","PsiR"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"Psi","Psibar"}
]//Expand
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*FStoFS*)


(* fermion-scalar scattering *)
process="FS->FS"
test["WallGo"][process]=testWallGo[
	{"PsiL","PsiR"},
	{"Phi"},
	{"PsiL","PsiR"},
	{"Phi"}
]
(* explicit 1/2 is due to average over leg 1 *)
test["FeynCalc"][process]=1/2*testFeynCalc[
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"Psi","Psibar"},
	{"Phi","Phibar"}
]//Expand
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*FFtoFF*)


(* fermion-fermion scattering*)
process="F1F1->F1F1"
test["WallGo"][process]=testWallGo[
	{"PsiL","PsiR"},
	{"PsiL","PsiR"},
	{"PsiL","PsiR"},
	{"PsiL","PsiR"}
]
(* explicit 1/2 is due to average over leg 1 *)
test["FeynCalc"][process]=1/2*testFeynCalc[
	{"Psi","Psibar"},
	{"Psi","Psibar"},
	{"Psi","Psibar"},
	{"Psi","Psibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];
		
(s1*test["WallGo"][process]-s2*test["FeynCalc"][process]//Simplify)/.{s1-s2->0}//fixConvention//Simplify


(* ::Subsection:: *)
(*Test report*)


report=TestReport[testList]
report["ResultsDataset"]



