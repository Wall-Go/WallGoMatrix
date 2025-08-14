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
(*Full Standard Model*)


(* ::Section:: *)
(*Model*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
HiggsDoublet={{{0,0},{1},1/2},"C"};
RepScalar={HiggsDoublet};
CouplingName={gs,gw,gY};


Rep1={{{1,0},{1},1/6},"L"};
Rep2={{{1,0},{0},2/3},"R"};
Rep3={{{1,0},{0},-1/3},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3};


RepFermion3Gen={RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


(* ::Text:: *)
(*The first element is the vector self - interaction matrix :*)


InputInv={{1,1},{True,False}};
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;


VMass=m2*MassTerm1;


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


QuarticTerm1=MassTerm1^2;


VQuartic=lam1H*QuarticTerm1;


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


Ysff=-GradYukawa[yt*YukawaDoublet[[1]]];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];


ImportModel[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Section:: *)
(*SM quarks + gauge bosons*)


(* ::Subsection:: *)
(*SymmetryBreaking*)


vev={0,v,0,0};
SymmetryBreaking[vev]


(* ::Subsection:: *)
(*UserInput*)


(*Third generation of fermions*)
ReptL=CreateParticle[{{1,1}},"F",mq2,"TopL"];
RepbL=CreateParticle[{{1,2}},"F",mq2,"BotL"];
ReptR=CreateParticle[{{2,1}},"F",mq2,"TopR"];
RepbR=CreateParticle[{3},"F",mq2,"BotR"];


(*Vector bosons*)
RepGluon=CreateParticle[{1},"V",mg2,"Gluon"]; (*Gluons*)
RepW=CreateParticle[{{2,1}},"V",mW2,"W"]; (*SU2 gauge bosons*)
RepB=CreateParticle[{3},"V",mB2,"B"]; (*U1 gauge boson*)


(*Scalars bosons*)
RepHiggs=CreateParticle[{{1,2}},"S",mH2,"Higgs"];
RepGoldstone=CreateParticle[{{1,1}},"S",mG2,"Goldstone"];


ParticleList={
	ReptL,RepbL,ReptR,RepbR,
	RepGluon,RepW,RepB,RepHiggs,RepGoldstone
	};
(*
Light particles are never incoming particles 
*)
LightParticleList={};


(*
	output of matrix elements
*)
OutputFile="matrixElements.smYukawa.test";

MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{
		TruncateAtLeadingLog->False,
		NormalizeWithDOF->False,
		Verbose->True,
		Format->{"json","txt"}}];


(* ::Chapter:: *)
(*Tests*)


(* ::Section:: *)
(*Importing results from WallGo*)


{particles,parameters,MatrixElements}=ImportMatrixElements["matrixElements.smYukawa.test.json"];


(* ::Section:: *)
(*Comparison tests*)


insertCouplings={lam1H->lam,gY->0 (* TEMPORARILY TURNING OFF U(1)*),gw->gW,gs->gs};


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})


UserMasses={mq2,ml2,mg2,mW2,mB2,mH2,mG2};
fixConvention[arg_]:=symmetriseTU[arg/.Thread[UserMasses->0]/.{s->(-t-u)}/.insertCouplings/.v->0]//Expand//Simplify//Expand


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
particlesFeyn={"W"->0,"Wbar"->1,"Z"->2,"gamma"->3,"H"->4,"G0"->5,"Gp"->6,"Gpbar"->7}


testList={};


(* ::Subsubsection:: *)
(*Electroweak bosons sector*)


(*Importing results from FeynCalc*)
{particlesFeyn,parametersFeyn,MatrixElementsFeyn}=ImportMatrixElements["SMQCD.EWbosons.test.json"];


process="{Higgs},{Higgs}->{Higgs},{Higgs}"
test["WallGo"][process]=testWallGo[
	{"Higgs"},
	{"Higgs"},
	{"Higgs"},
	{"Higgs"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"H"},
	{"H"},
	{"H"},
	{"H"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="{Higgs},{Higgs}->{W,B},{W,B}"
test["WallGo"][process]=testWallGo[
	{"Higgs"},
	{"Higgs"},
	{"W","B"},
	{"W","B"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"H"},
	{"H"},
	{"W","Wbar","Z","photon"},
	{"W","Wbar","Z","photon"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="{H},{W,B}->{H},{W,B}"
test["WallGo"][process]=testWallGo[
	{"Higgs"},
	{"W","B"},
	{"Higgs"},
	{"W","B"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"H"},
	{"W","Wbar","Z","photon"},
	{"H"},
	{"W","Wbar","Z","photon"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="{W,B},{W,B}->{W,B},{W,B}"
test["WallGo"][process]=testWallGo[
	{"W","B"},
	{"W","B"},
	{"W","B"},
	{"W","B"}
]//Simplify//Expand
test["FeynCalc"][process]=testFeynCalc[
	{"W","Wbar","Z","photon"},
	{"W","Wbar","Z","photon"},
	{"W","Wbar","Z","photon"},
	{"W","Wbar","Z","photon"}
]//Simplify//Expand
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];
		
(s1*test["WallGo"][process]-s2*test["FeynCalc"][process]//Simplify)/.{s1-s2->0}


process="{Higgs},{Higgs}->{W,Z,B,g},{W,Z,B,g}"
test["WallGo"][process]=testWallGo[
	{"Higgs"},
	{"Higgs"},
	{"Gluon","W","B"},
	{"Gluon","W","B"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"H"},
	{"H"},
	{"W","Wbar","Z","g"},
	{"W","Wbar","Z","g"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="{Higgs},{W,Z,B,g}->{Higgs},{W,Z,B,g}"
test["WallGo"][process]=testWallGo[
	{"Higgs"},
	{"Gluon","W","B"},
	{"Higgs"},
	{"Gluon","W","B"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"H"},
	{"W","Wbar","Z","g"},
	{"H"},
	{"W","Wbar","Z","g"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="{W,Z},{W,Z}->{W,Z},{W,Z}"
test["WallGo"][process]=testWallGo[
	{"W"},
	{"W"},
	{"W"},
	{"W"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"W","Wbar","Z"},
	{"W","Wbar","Z"},
	{"W","Wbar","Z"},
	{"W","Wbar","Z"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];
		
(s1*test["WallGo"][process]-s2*test["FeynCalc"][process]//Simplify)/.{s1-s2->0}


(* ::Subsubsection:: *)
(*Higgs tau sector*)


(*Importing results from FeynCalc*)
{particlesFeyn,parametersFeyn,MatrixElementsFeyn}=ImportMatrixElements["SMQCD.HiggsTau.test.json"];


process="{H,H}->{tau,tau}"
test["WallGo"][process]=testWallGo[
	{"Higgs"},
	{"Higgs"},
	{"TauL","TauR"},
	{"TauL","TauR"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"H"},
	{"H"},
	{"tau","taubar"},
	{"tau","taubar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection:: *)
(*Higgs top sector*)


(*Importing results from FeynCalc*)
{particlesFeyn,parametersFeyn,MatrixElementsFeyn}=ImportMatrixElements["SMQCD.HiggsTop.test.json"];


process="{H,H}->{t,t}"
test["WallGo"][process]=testWallGo[
	{"Higgs"},
	{"Higgs"},
	{"TopL","BotL","TopR"},
	{"TopL","BotL","TopR"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"H"},
	{"H"},
	{"t","tbar"},
	{"t","tbar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];
(s1*test["WallGo"][process]-s2*test["FeynCalc"][process]//Simplify)/.{s1-s2->0}


(* ::Subsubsection:: *)
(*Top bottom QCD sector *)


(*Importing results from FeynCalc*)
{particlesFeyn,parametersFeyn,MatrixElementsFeyn}=ImportMatrixElements["SMQCD.tbg.test.json"];


process="{t,b}->{t,b}"
test["WallGo"][process]=testWallGo[
	{"TopL","TopR","BotL","BotR"},
	{"TopL","TopR","BotL","BotR"},
	{"TopL","TopR","BotL","BotR"},
	{"TopL","TopR","BotL","BotR"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"t","tbar","b","bbar"},
	{"t","tbar","b","bbar"},
	{"t","tbar","b","bbar"},
	{"t","tbar","b","bbar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];

(* doesn't cancel exactly, and the difference involves the Yukawa coupling *)
(s1*test["WallGo"][process]-s2*test["FeynCalc"][process]//Simplify)/.{s1-s2->0}//fixConvention//Simplify


(* ::Subsubsection:: *)
(*QCD sector*)


(*Importing results from FeynCalc*)
{particlesFeyn,parametersFeyn,MatrixElementsFeyn}=ImportMatrixElements["SMQCD.tbg.test.json"];


process="{g},{g}->{g},{g}"
test["WallGo"][process]=testWallGo[
	{"Gluon"},
	{"Gluon"},
	{"Gluon"},
	{"Gluon"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"g"},
	{"g"},
	{"g"},
	{"g"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];
(s1*test["WallGo"][process]-s2*test["FeynCalc"][process]//Simplify)/.{s1-s2->0}


(* vector-vector scattering versus AMY *)
process="AMY{g},{g}->{g},{g}"
test["WallGo"][process]=testWallGo[
	{"Gluon"},
	{"Gluon"},
	{"Gluon"},
	{"Gluon"}
]
test["AMY"][process]=fixConvention[16 dA CA^2 g^4 (3-(s u)/t^2-(s t)/u^2-(t u)/s^2)/.{dA->8,CA->3,g->gs}]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["AMY"][process],
		TestID->"WallGo vs AMY:"<>process]];


report=TestReport[testList]
report["ResultsDataset"]



