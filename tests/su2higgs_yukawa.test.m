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


Group={"SU2","U1"}; (* I've added the U1 field because without it the Higgs is only given 2 degrees of freedom *)
CouplingName={g,gU1};
RepAdjoint={{2},0};
Higgs1={{{1},0},"C"}; (* fundamental *)
RepScalar={Higgs1};


su2Reps = RepsUpToDimN[SU2,2];
Grid[Prepend[{#,RepName[SU2,#]}&/@ su2Reps,{"Dynkin coefficients","Name"}],
Frame->All,FrameStyle->LightGray]


Rep1={{{1},0},"L"};
Rep2={{{1},0},"R"};
Rep3={{{0},0},"L"};
Rep4={{{0},0},"R"};
RepFermion={Rep1,Rep2,Rep3,Rep4};


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
InputInv={{1,1,4},{False,False,True}};
Yukawa1=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
InputInv={{1,2,3},{False,False,True}};
Yukawa2=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
InputInv={{1,4,1},{True,False,True}};
Yukawa1HC=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
InputInv={{1,3,2},{True,False,True}};
Yukawa2HC=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;


Ysff=-y*GradYukawa[Yukawa1+Yukawa2];
YsffC=-y*GradYukawa[Yukawa1HC+Yukawa2HC];
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


vev={0,0,0,v};
SymmetryBreaking[vev]


(*Below
rep 1-2 are fermions,
(*rep 3 is a scalar*)
*)
(* scalar *)
RepPhi=CreateParticle[{{1,1},{1,2}},"S",ms2,"Phi"]

(* left-handed fermion *)
RepPsi=CreateParticle[{1,2},"F",mf2,"Psi"]
RepXi=CreateParticle[{3,4},"F",mf2,"Xi"]

(*Vector bosons*)
RepA=CreateParticle[{1},"V",mv2,"VectorSU2"]
RepB=CreateParticle[{2},"V",mv2,"VectorU1"]


(*
These particles do not necessarily have to be out of equilibrium
*)
ParticleList={RepPhi,RepPsi,RepXi,RepA,RepB};
LightParticleList={};


(*Defining various masses and couplings*)
UserMasses={mv2,mf2,ms2};


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
		Verbose->True,
		Format->{"json","txt"}}];


MatrixElements;


M[0,0,1,1]/.MatrixElements
fixConvention[%]


(* ::Section:: *)
(*Tests*)


(* ::Section:: *)
(*Importing results from WallGo*)


{particles,parameters,MatrixElements}=ImportMatrixElements["output/matrixElements.su2_higgs_yukawa.json"];


(* ::Section:: *)
(*Importing results from FeynCalc*)


{particlesFeyn,parametersFeyn,MatrixElementsFeyn}=ImportMatrixElements["sun-higgs-yukawa.test.json"];


(* ::Section:: *)
(*Comparison tests*)


(* ::Subsection:: *)
(*Translate input*)


insertCouplings={Global`g->g,\[Lambda]->lam,SUNN->2,gu1->0,mChi->0,mPhi->0};
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


(* ::Subsection:: *)
(*Test hard*)


particles
particlesFeyn


testList={};


(* ::Subsubsection::Closed:: *)
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
test["AMY"][process]=fixConvention[16 dA CA^2 g^4 (3-(s u)/t^2-(s t)/u^2-(t u)/s^2)/.{dA->3,CA->2}]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["AMY"][process],
		TestID->"WallGo vs AMY:"<>process]];


(* ::Subsubsection::Closed:: *)
(*FFtoFF*)


process="F1F1->F1F1"
test["WallGo"][process]=testWallGo[
	{"Psi"},
	{"Psi"},
	{"Psi"},
	{"Psi"}
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
	{"Psi"},
	{"Xi"},
	{"Psi"},
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
	{"Psi"},
	{"Psi"},
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
	{"Psi"},
	{"Psi"},
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
	{"Psi"},
	{"VectorSU2"},
	{"Psi"},
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


(* ::Subsubsection::Closed:: *)
(*SStoFF*)


process="SS->F1F1"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"Psi"},
	{"Psi"}
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
	{"Psi"},
	{"Phi"},
	{"Psi"},
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
	{"Psi"},
	{"Phi"},
	{"Psi"},
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
	{"Psi"},
	{"Psi"},
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




