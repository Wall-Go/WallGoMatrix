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
(*SU(3)+Higgs-Yukawa Model*)


(* ::Section:: *)
(*Model*)


Group={"SU3","U1"}; (* I've added the U1 field because without it the Higgs is only given 2 degrees of freedom *)
CouplingName={g,gU1};
RepAdjoint={{1,1},0};
Higgs1={{{1,0},0},"C"}; (* fundamental *)
RepScalar={Higgs1};


su3Reps = RepsUpToDimN[SU3,3];
Grid[Prepend[{#,RepName[SU3,#]}&/@ su3Reps,{"Dynkin coefficients","Name"}],
Frame->All,FrameStyle->LightGray]


Rep1={{{1,0},0},"L"};
Rep2={{{0,0},0},"R"};
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


vev={0,0,0,0,0,v};
SymmetryBreaking[vev]


(*Below
rep 1-2 are fermions,
(*rep 3 is a scalar*)
*)
(* scalar *)
RepPhi=CreateParticle[{{1,1},{1,2}},"S",ms2,"Phi"]

(* left-handed fermion *)
RepPsiL=CreateParticle[{{1,1}},"F",mf2,"PsiL"]
RepXiL=CreateParticle[{{1,2}},"F",mf2,"XiL"]
RepPsiR=CreateParticle[{2},"F",mf2,"PsiR"]


(*Below
rep 1-2 are fermions,
(*rep 3 is a scalar*)
*)
(* scalar *)
RepPhi=CreateParticle[{{1,1},{1,2}},"S",ms2,"Phi"]

(* left-handed fermion *)
RepPsiL=CreateParticle[{{1,1}},"F",mf2,"PsiL"]
RepXiL=CreateParticle[{{1,2}},"F",mf2,"XiL"]
RepPsiR=CreateParticle[{2},"F",mf2,"PsiR"]

(*Vector bosons*)
RepA=CreateParticle[{1},"V",mv2,"Vector1"]
RepB=CreateParticle[{2},"V",mv2,"VectorU1"]


(*
These particles do not necessarily have to be out of equilibrium
*)
ParticleList={RepPhi,RepPsiL,RepPsiR,RepXiL,RepA,RepB};
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


(* ::Chapter:: *)
(*Tests*)


(* ::Section:: *)
(*Importing results from WallGo*)


{particles,parameters,MatrixElements}=ImportMatrixElements["output/matrixElements.su3_higgs_yukawa.json"];


(* ::Section:: *)
(*Importing results from FeynCalc*)


{particlesFeyn,parametersFeyn,MatrixElementsFeyn}=ImportMatrixElements["testFiles/sun-higgs-yukawa.test.json"];


testsRulesWallGo[arg_]:=arg/.Flatten[particles]/.MatrixElements//fixConvention//removeMissing;
testsRulesFeynCalc[arg_]:=arg/.Flatten[particlesFeyn]/.MatrixElementsFeyn//fixConvention//removeMissing
testWallGo[particlesA_,particlesB_,particlesC_,particlesD_]:=Sum[
	M[a,b,c,d],{a,particlesA},{b,particlesB},{c,particlesC},{d,particlesD}
]//testsRulesWallGo
testFeynCalc2[particlesA_,particlesB_,particlesC_,particlesD_]:=Sum[
	M[a,b,c,d],{a,particlesA},{b,particlesB},{c,particlesC},{d,particlesD}
]


(* ::Section:: *)
(*Comparison tests*)


(* ::Subsection::Closed:: *)
(*Translate input*)


insertCouplings={g->g,\[Lambda]->lam,SUNN->3,gu1->0,mChi->0,mPhi->0,mPsi->0};
customCouplings={ms2->mPhi^2};


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})


fixConvention[arg_]:=symmetriseTU[
	arg/.customCouplings/.Thread[UserMasses->0]/.{s->(-t-u)}/.insertCouplings
	]//Expand//Simplify//Expand


removeMissing[arg_]:=arg/.M[__]->0/.Missing["KeyAbsent", _]->0


testsRulesWallGo[arg_]:=arg/.Flatten[particles]/.MatrixElements//fixConvention//removeMissing;
testsRulesFeynCalc[arg_]:=arg/.Flatten[particlesFeyn]/.MatrixElementsFeyn//fixConvention//removeMissing
testWallGo[particlesA_,particlesB_,particlesC_,particlesD_]:=Sum[
	M[a,b,c,d],{a,particlesA},{b,particlesB},{c,particlesC},{d,particlesD}
]//testsRulesWallGo
testFeynCalc[particlesA_,particlesB_,particlesC_,particlesD_]:=Sum[
	M[a,b,c,d],{a,particlesA},{b,particlesB},{c,particlesC},{d,particlesD}
]//testsRulesFeynCalc


(* ::Subsection::Closed:: *)
(*Initialize tests*)


particles
particlesFeyn


testList={};


(* ::Subsubsection::Closed:: *)
(*VVtoVV*)


test["WallGo"][1]


process="VV->VV"
test["WallGo"][process]=testWallGo[
	{"Vector1"},
	{"Vector1"},
	{"Vector1"},
	{"Vector1"}
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
	{"Vector1"},
	{"Vector1"},
	{"Vector1"},
	{"Vector1"}
]
test["AMY"][process]=fixConvention[16 dA CA^2 g^4 (3-(s u)/t^2-(s t)/u^2-(t u)/s^2)/.{dA->CA^2-1}/.{CA->3}]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["AMY"][process],
		TestID->"WallGo vs AMY:"<>process]];


(* ::Subsubsection::Closed:: *)
(*FFtoFF*)


process="F1F1->F1F1"
test["WallGo"][process]=testWallGo[
	{"PsiL","XiL"},
	{"PsiL","XiL"},
	{"PsiL","XiL"},
	{"PsiL","XiL"}
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
	{"PsiR"},
	{"PsiR"},
	{"PsiR"},
	{"PsiR"}
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
	{"PsiL","XiL"},
	{"PsiR"},
	{"PsiL","XiL"},
	{"PsiR"}
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
	{"PsiL","XiL"},
	{"PsiL","XiL"},
	{"PsiR"},
	{"PsiR"}
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
	{"PsiL","XiL"},
	{"PsiL","XiL"},
	{"Vector1"},
	{"Vector1"}
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


process="F2F2->VV"
test["WallGo"][process]=testWallGo[
	{"PsiR"},
	{"PsiR"},
	{"Vector1"},
	{"Vector1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Chi","Chibar"},
	{"Chi","Chibar"},
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
	{"PsiL","XiL"},
	{"Vector1"},
	{"PsiL","XiL"},
	{"Vector1"}
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


(* ::Subsubsection::Closed:: *)
(*SStoFF*)


process="SS->F1F1"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"PsiL","XiL"},
	{"PsiL","XiL"}
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
(*FStoFS*)


process="F1S->F1S"
test["WallGo"][process]=testWallGo[
	{"PsiL","XiL"},
	{"Phi"},
	{"PsiL","XiL"},
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
	{"PsiL","XiL"},
	{"Phi"},
	{"PsiL","XiL"},
	{"Vector1"}
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


process="F1S->F2V"
test["WallGo"][process]=testWallGo[
	{"PsiL","XiL"},
	{"Phi"},
	{"PsiR"},
	{"Vector1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"Chi","Chibar"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F1S->F2V"
test["WallGo"][process]=testWallGo[
	{"PsiL","XiL"},
	{"Phi"},
	{"PsiR"},
	{"Vector1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"Chi","Chibar"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F2S->F1V"
test["WallGo"][process]=testWallGo[
	{"PsiR"},
	{"Phi"},
	{"PsiL","XiL"},
	{"Vector1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Chi","Chibar"},
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
	{"PsiL","XiL"},
	{"PsiL","XiL"},
	{"Phi"},
	{"Vector1"}
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


process="F1F2->SV"
test["WallGo"][process]=testWallGo[
	{"PsiL","XiL"},
	{"PsiR"},
	{"Phi"},
	{"Vector1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Chi","Chibar"},
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
	{"Vector1"},
	{"Vector1"}
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
	{"Vector1"},
	{"Phi"},
	{"Vector1"}
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
	{"Vector1"}
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


(* ::Subsection:: *)
(*Test report*)


report=TestReport[testList]
report["ResultsDataset"]


(* ::Section:: *)
(*Full test*)


(* ::Subsection::Closed:: *)
(*Initialize tests*)


totalWallGo=Sum[M[a,b,c,d],{a,0,4},{b,0,4},{c,0,4},{d,0,4}]/.MatrixElements/.Thread[UserMasses->0]//removeMissing//fixConvention;
totalFeyn=Sum[M[a,b,c,d],{a,0,6},{b,0,6},{c,0,6},{d,0,6}]/.MatrixElementsFeyn//removeMissing//fixConvention;


Collect[totalFeyn,{g,y,lam}];


Collect[s1*totalWallGo-s2*totalFeyn,{g,y,lam},Simplify[fixConvention[#]]&]/.{s1-s2->0}


testList={};


(* everything *)
AppendTo[testList,
TestCreate[
	totalWallGo,
	totalFeyn
]];


(* ::Subsection:: *)
(*Test report*)


report=TestReport[testList]
report["ResultsDataset"]



