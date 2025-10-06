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
(*IDM*)


(*See 2211.13142 for implementation details*)


(* ::Section:: *)
(*Model*)


Group={"SU3","SU2"};
RepAdjoint={{1,1},{2}};
HiggsDoublet1={{{0,0},{1}},"C"};
HiggsDoublet2={{{0,0},{1}},"C"};
RepScalar={HiggsDoublet1,HiggsDoublet2};
CouplingName={g3,gw};


Rep1={{{1,0},{1}},"L"};
Rep2={{{1,0},{0}},"R"};
Rep3={{{1,0},{0}},"R"};
Rep4={{{0,0},{1}},"L"};
Rep5={{{0,0},{0}},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};


(* ::Text:: *)
(*The input for the gauge interactions to DRalgo are then given by*)


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The first element is the vector self-interaction matrix:*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


InputInv={{1,1},{True,False}};
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;
InputInv={{2,2},{True,False}};
MassTerm2=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;
InputInv={{1,2},{True,False}};
MassTerm3=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;
InputInv={{2,1},{True,False}};
MassTerm4=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;


VMass=(
	+m1*MassTerm1
	+m2*MassTerm2
	);


\[Mu]ij=GradMass[VMass[[1]]]//Simplify//SparseArray;


QuarticTerm1=MassTerm1[[1]]^2;
QuarticTerm2=MassTerm2[[1]]^2;
QuarticTerm3=MassTerm1[[1]]*MassTerm2[[1]];
QuarticTerm4=MassTerm3[[1]]*MassTerm4[[1]];
QuarticTerm5=(MassTerm3[[1]]^2+MassTerm4[[1]]^2)//Simplify;


VQuartic=(
	+lam1H/2*QuarticTerm1
	+lam2H/2*QuarticTerm2
	+lam3H*QuarticTerm3
	+lam4H*QuarticTerm4
	+lam5H/2*QuarticTerm5
	);


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


Ysff=-GradYukawa[yt1*YukawaDoublet1[[1]]];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt1>0}]];


ImportModel[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Section:: *)
(*MatrixElements*)


(*
In DRalgo fermions are Weyl.
So to create one Dirac we need
one left-handed and
one right-handed fermoon
*)


(* ::Subsection:: *)
(*TopL, TopR*)


vev={0,v,0,0,0,0,0,0};
(*using this multiple times, should not change the outcome -> needs fixing *)
SymmetryBreaking[vev,VevDependentCouplings->True] (*uncomment if you want vev-dependent couplings*)
(*SymmetryBreaking[vev]*)


(*Third generation of fermions*)
ReptL=CreateParticle[{{1,1}},"F",mq2,"TopL"]
RepbL=CreateParticle[{{1,2}},"F",mq2,"BotL"]
ReptR=CreateParticle[{{2,1}},"F",mq2,"TopR"];
RepbR=CreateParticle[{{3,1}},"F",mq2,"BotR"];


(*Vector bosons*)
RepGluon=CreateParticle[{1},"V",mg2,"Gluon"]; (*Gluons*)
RepW=CreateParticle[{{2,1}},"V",mW2,"W"]; (*SU2 gauge bosons*)


(*Scalars bosons*)
RepHiggsh=CreateParticle[{{1,2}},"S",mh2,"Higgs"]; (*Higgs*)
RepGoldstoneGpR={{1},"S",mG2,"GoldstoneGpR"}; (*real charged Goldstone*)
RepGoldstoneGpI={{3},"S",mG2,"GoldstoneGpI"}; (*imag charged Golstone*)
RepGoldstoneGp0={{4},"S",mG2,"GoldstoneG0"}; (*neutral Goldstone*)
RepHiggsH=CreateParticle[{{2,2}},"S",mH2,"H"]; (*CP-even inert scalar*)
RepGoldstoneA=CreateParticle[{{2,3}},"S",mA2,"A"]; (*CP-odd inert and charged scalars *)
RepGoldstoneHpR={{5},"S",mHp,"RepGoldstoneHpR"}; (*real charged inert scalar*)
RepGoldstoneHpI={{7},"S",mHp,"RepGoldstoneHpI"}; (*imag charged inert scalar*)


(*Light fermions*)
LightQuarks=CreateParticle[{6,7,8,11,12,13},"F",mq2,"LightQuarks"];
LightLeptons=CreateParticle[{4,5,9,10,14,15},"F",ml2,"LightLeptons"];


ParticleList={
	ReptL,RepbL,ReptR,RepbR,
	RepGluon,RepW,
	RepHiggsh,RepGoldstoneGp0,RepGoldstoneGpR,RepGoldstoneGpI,
	RepHiggsH,RepGoldstoneA,RepGoldstoneHpR,RepGoldstoneHpI
	};
LightParticleList={
	LightQuarks,
	LightLeptons
	};


(* ::Subsection::Closed:: *)
(*Test W out of equilibrium only*)


ParticleListTesting={
	RepW
	};
LightParticleListTesting={
	ReptL,RepbL,ReptR,RepbR,
	RepGluon,
	RepHiggsh,RepGoldstoneGp0,RepGoldstoneGpR,RepGoldstoneGpI,
	RepHiggsH,RepGoldstoneA,RepGoldstoneHpR,RepGoldstoneHpI,
	LightQuarks,
	LightLeptons
	};


(*Allow only W on the first entry for speed*)
OutputFileTesting="output/matrixElementsTesting.idm";

MatrixElements=ExportMatrixElements[
	OutputFileTesting,
	ParticleListTesting,
	LightParticleListTesting,
	{
		TruncateAtLeadingLog->False,
		Replacements->{lam4H->0,lam5H->0},
		Format->{"json","txt"},
		NormalizeWithDOF->False,
		Verbose->True
	}
];


M[0,14,0,14]/.MatrixElements


(* ::Subsection:: *)
(*Test t out of equilibrium only*)


ParticleListTesting={
	ReptL
	};
LightParticleListTesting={
	RepbL,ReptR,RepbR,
	RepGluon,RepW,
	RepHiggsh,RepGoldstoneGp0,RepGoldstoneGpR,RepGoldstoneGpI,
	RepHiggsH,RepGoldstoneA,RepGoldstoneHpR,RepGoldstoneHpI,
	LightQuarks,
	LightLeptons
	};


(*Allow only W on the first entry for speed*)
OutputFileTesting="output/matrixElementsTesting.idm";

MatrixElements=ExportMatrixElements[
	OutputFileTesting,
	ParticleListTesting,
	LightParticleListTesting,
	{
		TruncateAtLeadingLog->False,
		TagLeadingLog->True,
		Replacements->{lam4H->0,lam5H->0},
		Format->{"json","txt"},
		NormalizeWithDOF->False,
		Verbose->True
	}
];


(* Due to the mixing between up and down-type quarks via W bosons we get power divergences *)
M[0,1,0,1]/.MatrixElements


testtensor=Table[gvff[[;;,Particle1,Particle2]],
		{Particle1,{{1,3,5},{2,4,6},{1,3,5},{2,4,6}}},
		{Particle2,{{1,3,5},{2,4,6},{1,3,5},{2,4,6}}}];


(* off diagonal entries can be seen here *)
gvff[[9]]//MatrixForm


(* or here *)
testtensor[[1,2]]//MatrixForm


(* ::Subsection::Closed:: *)
(*Full Test*)


(*
	output of matrix elements
*)
OutputFile="output/matrixElements.idm";

MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{
		TruncateAtLeadingLog->True,
		TagLeadingLog->True,
		Replacements->{lam4H->0,lam5H->0},
		Format->{"json","txt"},
		NormalizeWithDOF->False,
		Verbose->True
	}
];


M[5,14,5,14]/.MatrixElements


(*Show power divergences at variaous places that not only involves cubic scalar coupling*)
powDivEntries=Select[MatrixElements, ! FreeQ[#[[2]], powDiv] &]
Length[powDivEntries]


powDivEntries[[1]]


(* ::Chapter:: *)
(*Tests*)


(* ::Section:: *)
(*Importing results from WallGo*)


{particles,parameters,MatrixElements}=ImportMatrixElements["output/matrixElements.idm.json"];


(* ::Section:: *)
(*Comparison tests*)


(* ::Subsection:: *)
(*Translate input*)


groupFactors={ };
UserMasses={mq2,ml2,mg2,mw2};
customParameters={ms2->mPhi^2,mf2->mPsi^2,g->g,\[Lambda]->lam,logDiv->1,powDiv->1};
setMassesToZero={Thread[UserMasses->0],Thread[{mChi,mPhi,mPsi}->0],Thread[{mA2,mW2,mG2,mh2}->0]}//Flatten;
comparisonReplacements={
	groupFactors,
	customParameters,
	setMassesToZero
	}//Flatten;


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})

fixConvention[arg_]:=symmetriseTU[
	arg/.comparisonReplacements/.{s->(-t-u)}
	]//Expand//Simplify//Expand

removeMissing[arg_]:=arg/.M[__]->0/.Missing["KeyAbsent", _]->0

(*truncateAtLeadingLog[arg_]:=Module[{res,tterms,uterms,crossterms,constterms},
res[1]=fixConvention[arg];
(*res[1]=arg;*)
tterms=Normal[Series[res[1],{t,0,-1}]];
uterms=Normal[Series[res[1],{u,0,-1}]];
crossterms=(1/(t u))SeriesCoefficient[res[1],{t,0,-1},{u,0,-1}];
res[2]=tterms+uterms-crossterms;
constterms=SeriesCoefficient[res[2],{t,0,0},{u,0,0}];
res[3]=res[2]-constterms
];

dropTUCrossTerm[arg_]:=Module[{res},
(* the reference clearly drops 1/(t u) terms in the scalar matrix elements, so we do too *)
res[1]=fixConvention[arg];
res[2]=Normal[Series[res[1],{t,0,-2}]]+Normal[Series[res[1],{u,0,-2}]]
]*)


generateWallGo[particlesA_,particlesB_,particlesC_,particlesD_]:=Sum[
	M[a,b,c,d],{a,particlesA},{b,particlesB},{c,particlesC},{d,particlesD}
]/.Flatten[particles]/.MatrixElements/.customParameters//removeMissing;

testWallGo[particlesA_,particlesB_,particlesC_,particlesD_]:=
	generateWallGo[particlesA,particlesB,particlesC,particlesD]//fixConvention


result=(s1 (128s^2)/(t s)+s2 (128t^2
 )/(s t)+s3 (128s^2
 )/(u s)+s4 (128t^2
 )/(u t)+s5 (128u^2
 )/(s u)+s6 (128u^2 )/(t u)
 )
 result(*//fixConvention*)//truncateAtLeadingLog
 %(*/.Thread[{s1,s2,s3,s4,s5,s6}->1]*)//symmetriseTU


(* ::Subsection:: *)
(*Initialize tests*)


particles


testList={};


(* ::Subsubsection:: *)
(*FFtoVV*)


process="tt->gg"
test["WallGo"][process]=testWallGo[
	{"TopL"},
	{"TopL"},
	{"Gluon"},
	{"Gluon"}
](*//truncateAtLeadingLog*)
test["Reference"][process]=(
		128/3*g3^4(u/(t-mq2^2)+t/(u-mq2^2))
	)/.setMassesToZero//fixConvention(*//truncateAtLeadingLog*)
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


process="tq->tq"
test["WallGo"][process]=testWallGo[
	{"TopL"},
	{"LightQuarks","BotL","BotR"},
	{"TopL"},
	{"LightQuarks","BotL","BotR"}
]/.{gw->0,yt1->0}(*//truncateAtLeadingLog*)
test["Reference"][process]=(
		+160*g3^4*(s^2+u^2)/(t-mq2^2)^2
	)/.setMassesToZero//fixConvention(*//truncateAtLeadingLog*)
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


process="Wg->qq"
test["WallGo"][process]=testWallGo[
	{"W"},
	{"Gluon"},
	{"LightQuarks","BotL","BotR","TopL","TopR"},
	{"LightQuarks","BotL","BotR","TopL","TopR"}
](*//truncateAtLeadingLog*)
(* have to put factor 4 to agree *)
test["Reference"][process]=4*(
		-72*g3^2*gw^2*s*t/(t-mq2)^2
	)/.setMassesToZero//fixConvention(*//truncateAtLeadingLog*)
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


process="WW->ff"
test["WallGo"][process]=testWallGo[
	{"W"},
	{"W"},
	{"LightLeptons","LightQuarks","BotL","TopL"},
	{"LightLeptons","LightQuarks","BotL","TopL"}
](*//truncateAtLeadingLog*)
(* have to put factor 4 to agree *)
test["Reference"][process]=4*(
		-27/2*gw^4*(3*s/(t-mq2)+s/(t-mq2))
	)/.setMassesToZero//fixConvention(*//truncateAtLeadingLog*)
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


(* ::Subsubsection::Closed:: *)
(*FVtoFV*)


process="tg->tg"
test["WallGo"][process]=testWallGo[
	{"TopL"},
	{"Gluon"},
	{"TopL"},
	{"Gluon"}
](*//truncateAtLeadingLog*)
test["Reference"][process]=(
		-128/3*g3^4*(s*u/(u-mg2^2)^2)
		+96*g3^4*(s^2+u^2)/(t-mq2^2)^2
	)/.setMassesToZero//fixConvention(*//truncateAtLeadingLog*)
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


process="Wq->qg"
test["WallGo"][process]=testWallGo[
	{"W"},
	{"LightQuarks","BotL","BotR","TopL","TopR"},
	{"LightQuarks","BotL","BotR","TopL","TopR"},
	{"Gluon"}
]//truncateAtLeadingLog
test["Reference"][process]=2*(
		-72*g3^2*gw^2*s/(t-mq2)
	)/.setMassesToZero//fixConvention//truncateAtLeadingLog
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


process="fW->fW"
test["WallGo"][process]=testWallGo[
	{"W"},
	{"LightLeptons","LightQuarks","BotL","BotR","TopL","TopR"},
	{"W"},
	{"LightLeptons","LightQuarks","BotL","BotR","TopL","TopR"}
]/.setMassesToZero//fixConvention//truncateAtLeadingLog
(* have to replace 360 -> 280 in first term *)
test["Reference"][process]=(
		+s1*288*gw^4*(s^2+u^2)/(t-mw2)^2
		-s2*27/2*gw^4*(3*s/(u-mq2)+s/(u-ml2))
	)/.{s1->1,s2->2}/.setMassesToZero//fixConvention//truncateAtLeadingLog
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


(* ::Subsubsection::Closed:: *)
(*FFtoSV*)


process="tt->gh"
test["WallGo"][process]=testWallGo[
	{"TopL"},
	{"TopR"},
	{"Gluon"},
	{"Higgs"}
](*//truncateAtLeadingLog*)
test["Reference"][process]=(
		8*yt1^2*g3^2(u/(t-mq2)+t/(u-mq2))
	)/.setMassesToZero//fixConvention(*//truncateAtLeadingLog*)
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


process="tt->g G0"
test["WallGo"][process]=testWallGo[
	{"TopL"},
	{"TopR"},
	{"Gluon"},
	{"GoldstoneG0"}
](*//truncateAtLeadingLog*)
test["Reference"][process]=(
		8*yt1^2*g3^2(u/(t-mq2)+t/(u-mq2))
	)/.setMassesToZero//fixConvention(*//truncateAtLeadingLog*)
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


(* ::Subsubsection::Closed:: *)
(*FVtoFS*)


process="tg->th"
test["WallGo"][process]=testWallGo[
	{"TopL"},
	{"Gluon"},
	{"TopR"},
	{"Higgs"}
](*//truncateAtLeadingLog*)
test["Reference"][process]=(
		-8*yt1^2*g3^2(s/(t-mq2))
	)/.setMassesToZero//fixConvention(*//truncateAtLeadingLog*)
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


process="tg->tG0"
test["WallGo"][process]=testWallGo[
	{"TopL"},
	{"Gluon"},
	{"TopR"},
	{"GoldstoneG0"}
]//truncateAtLeadingLog
test["Reference"][process]=(
		-8*yt1^2*g3^2(s/(t-mq2))
	)/.setMassesToZero//fixConvention//truncateAtLeadingLog
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


process="tg->tG+"
test["WallGo"][process]=testWallGo[
	{"TopL","TopR"},
	{"Gluon"},
	{"BotL"},
	{"GoldstoneGpR"}
]//truncateAtLeadingLog
test["Reference"][process]=(
		-8*yt1^2*g3^2(s/(t-mq2))
	)/.setMassesToZero//fixConvention//truncateAtLeadingLog
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


(* ::Subsubsection::Closed:: *)
(*FStoFV*)


process="tG-->bg"
test["WallGo"][process]=testWallGo[
	{"TopL","TopR"},
	{"GoldstoneGpI"},
	{"BotL"},
	{"Gluon"}
](*//truncateAtLeadingLog*)
test["Reference"][process]=(
		-8*yt1^2*g3^2(s/(t-mq2))
	)/.setMassesToZero//fixConvention(*//truncateAtLeadingLog*)
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


(* ::Subsubsection:: *)
(*FFtoSS*)


process="tb->h,G+"
(* here we have only u/t but not t/u contribution but there should be both such
that for similar particles between 12 and 34, the result is correct *)
test["WallGo"][process]=testWallGo[
	{"TopL","TopR"},
	{"BotL","BotR"},
	{"Higgs"},
	{"GoldstoneGpR","GoldstoneGpI"}
](*//Coefficient[#,logDiv]&//Simplify//Expand*)
test["Reference"][process]=1/2*(
		+6 u/t yt1^4
		(*8*yt1^2*g3^2(u/(t-mq2)+t/(u-mq2)*)
	)/.setMassesToZero//fixConvention(*//truncateAtLeadingLog*)
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


(*this is the Feyncalc result and it shows that only the yt^4 contribution is LL *)
resFC=(3*u*(gw^2*t + 2*s*yt^2)^2)/(2*s^2*t)//Expand
%//fixConvention(*//truncateAtLeadingLog*)


(*there is no more difference in putting factors before or after LL trunction *)
s1 resWG- s2 resFC/.{yt->yt1,lam1H->0}//fixConvention//truncateAtLeadingLog
%/.{s1->2,s2->1}
s1 resWG- s2 resFC/.{s1->2,s2->1}/.{yt->yt1,lam1H->0}//fixConvention//truncateAtLeadingLog


(* ::Subsubsection:: *)
(*SStoSS*)


process="HH->AA"
test["WallGo"][process]=testWallGo[
	{"Higgs"},
	{"Higgs"},
	{"A"},
	{"A"}
](*//truncateAtLeadingLog//dropTUCrossTerm*)
(* the Reference result drops the 1/(t u) cross term here *)
test["Reference"][process]=2*(
		lam3H^4*v^4/2*(1/(t-mA2)^2+1/(u-mA2)^2)
	)/.setMassesToZero//fixConvention(*//truncateAtLeadingLog*)
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


process="AH->HA"
test["WallGo"][process]=testWallGo[
	{"A"},
	{"Higgs"},
	{"Higgs"},
	{"A"}
]/.{lam1H->0}(*//truncateAtLeadingLog//dropTUCrossTerm*)
(* the Reference result drops the 1/(t u) cross term here *)
test["Reference"][process]=2*(
		lam3H^4*v^4/2*(1/(t-mA2)^2)
	)/.setMassesToZero//fixConvention(*//truncateAtLeadingLog*)
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


(* ::Subsection:: *)
(*Test report*)


report=TestReport[testList]
report["ResultsDataset"]




(* ::Section:: *)
(*Test if the propagator mass is correct*)


fermions = {{"TopL"},{"BotL"},{"TopR"},{"BotR"},{"LightQuarks"},{"LightLeptons"}};
vectors = {{"Gluon"},{"W"}};


(* ::Subsection:: *)
(*FF -> VV*)


testList = {};
For[l = 1, l<= Length[vectors], l++,
	For[k = 1, k <= Length[vectors],k++,
		For[j = 1, j<= Length[fermions],j++,
			For[i = 1, i<= Length[fermions],i++,
				matrixElement = generateWallGo[fermions[[i]],fermions[[j]],vectors[[k]],vectors[[l]]];
				(*If the external vectors are identical, there can be an s-channel gauge boson exchange, or t- or u-channel quark exchange*)
				If[vectors[[k]]=={"Gluon"}&&vectors[[l]]=={"Gluon"},
					AppendTo[testList,matrixElement -(SeriesCoefficient[matrixElement, {1/(-mg2+s)^2,0,1}]/(-mg2+s)^2+SeriesCoefficient[matrixElement, {1/(-mq2+t)^2,0,1}]/(-mq2+t)^2+SeriesCoefficient[matrixElement, {1/(-mq2+u)^2,0,1}]/(-mq2+u)^2)//FullSimplify ]
				];
				(*If the external vectors are different, there can only be t- or u-channel quark exchange*)
				If[(vectors[[k]]=={"Gluon"}&&vectors[[l]]=={"W"}) || (vectors[[k]]=={"W"} && vectors[[l]]=={"Gluon"}),
					AppendTo[testList,matrixElement -(SeriesCoefficient[matrixElement, {1/(-mq2+t)^2,0,1}]/(-mq2+t)^2+SeriesCoefficient[matrixElement, {1/(-mq2+u)^2,0,1}]/(-mq2+u)^2)//FullSimplify ]
				];
				(*If the external vectors are identical, there can be an s-channel gauge boson exchange, or t- or u-channel quark exchange*)
				If[vectors[[k]]=={"W"}&&vectors[[l]]=={"W"},
					AppendTo[testList,matrixElement -(SeriesCoefficient[matrixElement, {1/(-mW2+s)^2,0,1}]/(-mW2+s)^2+SeriesCoefficient[matrixElement, {1/(-mq2+t)^2,0,1}]/(-mq2+t)^2+SeriesCoefficient[matrixElement, {1/(-mq2+u)^2,0,1}]/(-mq2+u)^2)//FullSimplify ]
				];


				]
		]
	]
];
testmasses=TestCreate[Total[testList - Table[0,{i,Length[testList]}]]==0];
TestEvaluate[testmasses]
