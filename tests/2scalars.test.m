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
(*2scalar test*)


(*Pure scalar sector of 2211.13142*)


(* ::Section:: *)
(*Model*)


Group={"SU2"};
RepAdjoint={{2}};


HiggsDoublet1={{{1}},"C"};
HiggsDoublet2={{{1}},"C"};
RepScalar={HiggsDoublet1,HiggsDoublet2};
CouplingName={gw};


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
(*using this multiple times, should not changre the outcome -> needs fixing *)
SymmetryBreaking[vev,VevDependentCouplings->True] (*uncomment if you want vev-dependent couplings*)
(*SymmetryBreaking[vev]*)


(*Vector bosons*)
RepW=CreateParticle[{{1,1}},"V",mW2,"W"]; (*SU2 gauge bosons*)


(*Scalars bosons*)
RepHiggsh=CreateParticle[{{1,2}},"S",mh2,"Higgs"]; (*Higgs*)
RepGoldstoneGpR={{1},"S",mG2,"GoldstoneGPR"}; (*real charged Goldstone*)
RepGoldstoneGpI={{3},"S",mG2,"RepGoldstoneGpI"}; (*imag charged Golstone*)
RepGoldstoneGp0={{4},"S",mG2,"RepGoldstoneGp0"}; (*neutral Goldstone*)
RepHiggsH=CreateParticle[{{2,2}},"S",mH2,"H"]; (*CP-even inert scalar*)
RepGoldstoneA=CreateParticle[{{2,3}},"S",mA2,"A"]; (*CP-odd inert scalar*)
RepGoldstoneHpR={{5},"S",mHp,"RepGoldstoneHpR"}; (*real charged inert scalar*)
RepGoldstoneHpI={{7},"S",mHp,"RepGoldstoneHpI"}; (*imag charged inert scalar*)


(*Defining various masses and couplings*)
UserMasses={mw2,mG2,mh2,mH2,mA2,mHp};


ParticleList={
	RepHiggsh,RepGoldstoneGp0,RepGoldstoneGpR,RepGoldstoneGpI,
	RepHiggsH,RepGoldstoneA,RepGoldstoneHpR,RepGoldstoneHpI};
(*
Light particles are never incoming particles 
*)
LightParticleList={RepW};


(*
	output of matrix elements
*)
OutputFile="output/matrixElements.2scalars";

MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{
		TruncateAtLeadingLog->True,
		Replacements->{gw->0,lam4H->0,lam5H->0},
		Format->{"json","txt"},
		NormalizeWithDOF->True
	}
];


(* ::Chapter:: *)
(*Tests*)


(* ::Section:: *)
(*Importing results from WallGo*)


{particles,parameters,MatrixElements}=ImportMatrixElements["output/matrixElements.2scalars.json"];


(* ::Section:: *)
(*Comparison tests*)


(* ::Subsection::Closed:: *)
(*Translate input*)


groupFactors={ };
customParameters={ms2->mPhi^2,mf2->mPsi^2,g->g,\[Lambda]->lam};
setMassesToZero={Thread[UserMasses->0],Thread[{mChi,mPhi,mPsi}->0]}//Flatten;
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


generateWallGo[particlesA_,particlesB_,particlesC_,particlesD_]:=Sum[
	M[a,b,c,d],{a,particlesA},{b,particlesB},{c,particlesC},{d,particlesD}
]/.Flatten[particles]/.MatrixElements/.customParameters//removeMissing;

testWallGo[particlesA_,particlesB_,particlesC_,particlesD_]:=
	generateWallGo[particlesA,particlesB,particlesC,particlesD]//fixConvention


(* ::Subsection::Closed:: *)
(*Initialize tests*)


particles


testList={};


(* ::Subsubsection::Closed:: *)
(*SStoSS*)


process="HH->AA"
test["WallGo"][process]=testWallGo[
	{"Higgs"},
	{"Higgs"},
	{"A"},
	{"A"}
]
test["Reference"][process]=2*(
		lam3H^4*v^4/2*(1/(t-mA2)^2+1/(u-mA2)^2)
	)/.setMassesToZero//fixConvention
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
]/.{lam1H->0,v^2->0}
test["Reference"][process]=2*(
		lam3H^4*v^4/2*(1/(t-mA2)^2)
	)/.setMassesToZero//fixConvention
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["Reference"][process],
		TestID->"WallGo vs Reference: "<>process]];


(* ::Subsection:: *)
(*Test report*)


report=TestReport[testList]
report["ResultsDataset"]



