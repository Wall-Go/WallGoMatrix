(* ::Package:: *)

Quit[];


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


Group={"SU5","U1"}; (* I've added the U1 field because without it the Higgs is only given 2 degrees of freedom *)
CouplingName={g,gU1};
RepAdjoint={{1,0,0,1},0};
Higgs1={{{1,0,0,0},0},"C"}; (* fundamental *)
RepScalar={Higgs1};


su5Reps = RepsUpToDimN[SU5,25];
Grid[Prepend[{#,RepName[SU5,#]}&/@ su5Reps,{"Dynkin coefficients","Name"}],
Frame->All,FrameStyle->LightGray]


Rep1={{{1,0,0,0},0},"L"};
Rep2={{{0,0,0,0},0},"R"};
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
InputInv={{1,1,2},{True,False,False}};
Yukawa1=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
(*InputInv={{1,2,3},{True,False,False}};
Yukawa2=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;*)


Ysff=-y*GradYukawa[Yukawa1+Yukawa2];
YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{y>0}]];


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


vev={0*Range[2*5-1],v}//Flatten
SymmetryBreaking[vev]


(*Below
rep 1-2 are fermions,
(*rep 3 is a scalar*)
*)
(* scalar *)
RepPhi=CreateParticle[{{1,1},{1,2}},"S",ms2,"Phi"]

(* left-handed fermion *)
RepPsiL=CreateParticle[{{1,1}},"F",mf2,"PsiL"]
RepPsiR=CreateParticle[{{1,2}},"F",mf2,"PsiR"]
RepXi=CreateParticle[{2},"F",mf2,"Xi"]


(*Vector bosons*)
RepA=CreateParticle[{1},"V",mv2,"VectorSU2"]
RepB=CreateParticle[{2},"V",mv2,"VectorU1"]


(*
These particles do not necessarily have to be out of equilibrium
*)
ParticleList={RepPhi,RepPsiL,RepPsiR,RepXi,RepA,RepB};
LightParticleList={};


(*Defining various masses and couplings*)
UserMasses={mv2,mf2,ms2};


(*
	output of matrix elements
*)
OutputFile="output/matrixElements.su3_higgs_yukawa";
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


{particles,parameters,MatrixElements}=ImportMatrixElements["output/matrixElements.su3_higgs_yukawa.json"];


(* ::Section:: *)
(*Importing results from FeynCalc*)


{particlesFeyn,parametersFeyn,MatrixElementsFeyn}=ImportMatrixElements["sun-higgs-yukawa.test.json"];


(* ::Section:: *)
(*Comparison tests*)


(* ::Subsection:: *)
(*Translate input*)


insertCouplings={Global`g->g,\[Lambda]->lam,SUNN->5,gu1->0,mChi->0,mPhi->0};
customCouplings={ms2->mPhi^2};


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})


fixConvention[arg_]:=symmetriseTU[
	arg/.customCouplings/.{s->(-t-u)}/.insertCouplings/.Thread[UserMasses->0]
	]//(*Expand//*)Simplify//Expand


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


(* ::Subsubsection:: *)
(*VVtoVV*)


test["WallGo"][1]


process="VV->VV"
test["WallGo"][process]=testWallGo[
	{"VectorSU2"},
	{"VectorSU2"},
	{"VectorSU2"},
	{"VectorSU2"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"A"},
	{"A"},
	{"A"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* vector-vector scattering versus AMY *)
process="VV->VV"
test["WallGo"][process]=testWallGo[
	{"VectorSU2"},
	{"VectorSU2"},
	{"VectorSU2"},
	{"VectorSU2"}
]
test["AMY"][process]=fixConvention[16 dA CA^2 g^4 (3-(s u)/t^2-(s t)/u^2-(t u)/s^2)/.{dA->5^2-1,CA->5}]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["AMY"][process],
		TestID->"WallGo vs AMY:"<>process]];


(* ::Subsubsection::Closed:: *)
(*FFtoFF*)


process="F1F1->F1F1"
test["WallGo"][process]=testWallGo[
	{"PsiL","PsiR"},
	{"PsiL","PsiR"},
	{"PsiL","PsiR"},
	{"PsiL","PsiR"}
]
test["FeynCalc"][process]=testFeynCalc[
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


(* F2F2->F2F2 *)
process="F2F2->F2F2"
test["WallGo"][process]=testWallGo[
	{"Xi"},
	{"Xi"},
	{"Xi"},
	{"Xi"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Chi","Chibar"},
	{"Chi","Chibar"},
	{"Chi","Chibar"},
	{"Chi","Chibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F1F2->F1F2"
test["WallGo"][process]=testWallGo[
	{"PsiL","PsiR"},
	{"Xi"},
	{"PsiL","PsiR"},
	{"Xi"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Chi","Chibar"},
	{"Psi","Psibar"},
	{"Chi","Chibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F1F1->F2F2"
test["WallGo"][process]=testWallGo[
	{"PsiL","PsiR"},
	{"PsiL","PsiR"},
	{"Xi"},
	{"Xi"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Psi","Psibar"},
	{"Chi","Chibar"},
	{"Chi","Chibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*FFtoVV*)


process="F1F1->VV"
test["WallGo"][process]=testWallGo[
	{"PsiL","PsiR"},
	{"PsiL","PsiR"},
	{"VectorSU2"},
	{"VectorSU2"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Psi","Psibar"},
	{"A"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*FVtoFV*)


process="F1V->F1V"
test["WallGo"][process]=testWallGo[
	{"PsiL","PsiR"},
	{"VectorSU2"},
	{"PsiL","PsiR"},
	{"VectorSU2"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"A"},
	{"Psi","Psibar"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*SStoSS*)


process="SS->SS"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"Phi"},
	{"Phi"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Phi","Phibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection:: *)
(*SStoFF*)


process="SS->F1F1"
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


2*3*8/40*4/3


(* ::Subsubsection:: *)
(*FStoFS*)


process="F1S->F1S"
test["WallGo"][process]=testWallGo[
	{"PsiL","PsiR"},
	{"Phi"},
	{"PsiL","PsiR"},
	{"Phi"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"Psi","Psibar"},
	{"Phi","Phibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*FStoFV*)


process="F1S->F1V"
test["WallGo"][process]=testWallGo[
	{"PsiL","PsiR"},
	{"Phi"},
	{"PsiL","PsiR"},
	{"VectorSU2"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"Psi","Psibar"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*FFtoSV*)


process="F1F1->SV"
test["WallGo"][process]=testWallGo[
	{"PsiL","PsiR"},
	{"PsiL","PsiR"},
	{"Phi"},
	{"VectorSU2"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*SStoVV*)


process="SS->VV"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"VectorSU2"},
	{"VectorSU2"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"A"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*SVtoSV*)


process="SV->SV"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"VectorSU2"},
	{"Phi"},
	{"VectorSU2"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"A"},
	{"Phi","Phibar"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*SStoSV*)


process="SS->SV"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"Phi"},
	{"VectorSU2"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection:: *)
(*Test report*)


report=TestReport[testList]
report["ResultsDataset"]




