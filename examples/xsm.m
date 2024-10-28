(* ::Package:: *)

(*Quit[];*)


If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
(*Put this if you want to create multiple model-files with the same kernel*)
$GroupMathMultipleModels=True;
$LoadGroupMath=True;
Check[
    Get["WallGoMatrix`"],
    Message[Get::noopen, "WallGoMatrix` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*SM+scalar Singlet*)


(*see 2102.11145 [hep-ph]*)


(* ::Section::Closed:: *)
(*Model*)


HypY={Yl,Ye,Yq,Yu,Yd,Y\[Phi],Y\[Eta]};
repY=Thread[HypY->{-1,-2,1/3,4/3,-(2/3),1,2}];


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
scalar1={{{0,0},{1},Y\[Phi]/2},"C"};
scalar2={{{0,0},{0},0},"R"};
RepScalar={scalar1,scalar2}/.repY;
CouplingName={g3,gw,g1};


Rep1={{{1,0},{1},Yq/2},"L"};
Rep2={{{1,0},{0},Yu/2},"R"};
Rep3={{{1,0},{0},Yd/2},"R"};
Rep4={{{0,0},{1},Yl/2},"L"};
Rep5={{{0,0},{0},Ye/2},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5}/.repY;


(* ::Text:: *)
(*The input for the gauge interactions to DRalgo are then given by*)


RepFermion1Gen={Rep1,Rep2,Rep3,Rep1,Rep2,Rep3,Rep1,Rep2,Rep3,Rep4,Rep5}/.repY;
RepFermion3Gen={RepFermion1Gen}//Flatten[#,1]&;
(*RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;*)


(* ::Text:: *)
(*The first element is the vector self-interaction matrix:*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar]/.repY;


InputInv={{1,1},{True,False}};
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;
InputInv={{2,2},{True,True}};
MassTerm2=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;


VMass=(
	+m1*MassTerm1
	+ms/2*MassTerm2
	);


\[Mu]ij=GradMass[VMass[[1]]]//Simplify//SparseArray;


QuarticTerm1=MassTerm1[[1]]^2;
QuarticTerm2=MassTerm2[[1]]^2;
QuarticTerm3=MassTerm1[[1]]*MassTerm2[[1]];


VQuartic=(
	+lam1H*QuarticTerm1
	+lams/4*QuarticTerm2
	+lamm/2*QuarticTerm3
	);


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{True,False,True}};
CubicTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2,2},{True,True,True}};
CubicTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VCubic=(
	+mum/2*CubicTerm1
	+mu3/3*CubicTerm2
	);


\[Lambda]3=GradCubic[VCubic];


InputInv={{2},{True}};
TadpoleTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VTadpole=\[Mu]1*TadpoleTerm1;


\[Lambda]1=GradTadpole[VTadpole];


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


Ysff=-yt1*GradYukawa[YukawaDoublet1[[1]]];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt1>0}]];


ImportModel[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Section:: *)
(*MatrixElements*)


(*
In DRalgo fermions are Weyl.
So to create one Dirac we need
one left-handed and
one right-handed fermion
*)


(* ::Subsection:: *)
(*TopL, TopR*)


vev={0,v,0,0,0};
SymmetryBreaking[vev]


(*left-handed top-quark*)
ReptL=CreateParticle[{{1,1}},"F","TopL"];

(*right-handed top-quark*)
ReptR=CreateParticle[{2},"F","TopR"];

(*left-handed bottom-quark*)
RepbL=CreateParticle[{{1,2}},"F","BotL"];

(*light quarks*)
RepLight=CreateParticle[Range[3,2*4+3],"F","Light"];

(*Vector bosons*)
RepGluon=CreateParticle[{1},"V","Gluon"];
RepW=CreateParticle[{{2,1}},"V","W"];
RepB=CreateParticle[{{3,1}},"V","B"];

(*Higgs*)
RepH = CreateParticle[{1},"S","H"];
RepS = CreateParticle[{2},"S","S"];


VectorMass=Join[
	Table[mg2,{i,1,RepGluon[[1]]//Length}],
	Table[mW2,{i,1,RepW[[1]]//Length}],
	Table[mB2,{i,1,RepB[[1]]//Length}]];
FermionMass=Table[mq2,{i,1,Length[gvff[[1]]]}];
ScalarMass=Table[ms2,{i,1,Length[gvss[[1]]]}];
ParticleMasses={VectorMass,FermionMass,ScalarMass};


(*
These particles do not have out-of-eq contributions
*)
ParticleList={ReptL,ReptR,RepGluon,RepW,RepB,RepH,RepS};


(*
	output of matrix elements
*)
OutputFile="output/matrixElements.xsm";
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	ParticleMasses,
	{
		Replacements->{gw->0,g1->0},
		Format->{"json","txt"}}];


(* ::Subsection:: *)
(*QCD -like*)


vev={0,v,0,0,0};
SymmetryBreaking[vev]


(*left-handed top-quark*)
ReptL=CreateParticle[{{1,1}},"F","TopL"];

(*right-handed top-quark*)
ReptR=CreateParticle[{2},"F","TopR"];

(*join topL and topR into one rep*)
Rept=CreateParticle[{{1,1},2},"F","Top"];

(*left-handed bottom-quark*)
RepbL=CreateParticle[{{1,2}},"F","BotL"];

(*light quarks*)
RepLight=CreateParticle[Range[3,2*4+3],"F","Light"];

(*Vector bosons*)
RepGluon=CreateParticle[{1},"V","Gluon"];
RepW=CreateParticle[{2},"V","W"];
RepB=CreateParticle[{3},"V","B"];

(*Higgs*)
RepH = CreateParticle[{1},"S","H"];
RepS = CreateParticle[{2},"S","S"];


VectorMass=Join[
	Table[mg2,{i,1,RepGluon[[1]]//Length}],
	Table[mW2,{i,1,RepW[[1]]//Length}],
	Table[mB2,{i,1,RepB[[1]]//Length}]];
FermionMass=Table[mq2,{i,1,Length[gvff[[1]]]}];
ScalarMass=Table[ms2,{i,1,Length[gvss[[1]]]}];
ParticleMasses={VectorMass,FermionMass,ScalarMass};


(*
These particles do not have out-of-eq contributions
*)
ParticleList={Rept,RepGluon,RepW,RepB,RepH,RepS};


(*
	output of matrix elements
*)
OutputFile="output/matrixElements.xsm.qcd";
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	ParticleMasses,
	{
		Replacements->{gw->0,g1->0},
		Format->{"json","txt"}}];


MatrixElements//Expand
