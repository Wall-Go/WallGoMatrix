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
    Get["../Kernel/WallGoMatrix.m"],
    Message[Get::noopen, "WallGoMatrix` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*2HDM*)


(*See 2211.13142 for implementation details*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3","SU2"};
RepAdjoint={{1,1},{2},0};
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
(*using this multiple times, should not changre the outcome -> needs fixing *)
SymmetryBreaking[vev,VevDependentCouplings->True] (*uncomment if you want vev-dependent couplings*)
(*SymmetryBreaking[vev]*)


(*Third generation of fermions*)
ReptL=CreateParticle[{{1,1}},"F",mq2,"TopL"];
RepbL=CreateParticle[{{1,2}},"F",mq2,"BotL"];
ReptR=CreateParticle[{{2,1}},"F",mq2,"TopR"];


(*Vector bosons*)
RepGluon=CreateParticle[{1},"V",mg2,"Gluon"]; (*Gluons*)
RepW=CreateParticle[{{2,1}},"V",mW2,"W"]; (*SU2 gauge bosons*)


(*Scalars bosons*)
RepHiggsh=CreateParticle[{{1,2}},"S",mh2,"Higgs"]; (*Higgs*)
RepGoldstoneGpR={{1},"S",mG2,"GoldstoneGpR"}; (*real charged Goldstone*)
RepGoldstoneGpI={{3},"S",mG2,"GoldstoneGpI"}; (*imag charged Golstone*)
RepGoldstoneGp0={{4},"S",mG2,"GoldstoneG0"}; (*neutral Goldstone*)
RepHiggsH=CreateParticle[{{2,2}},"S",mH2,"H"]; (*CP-even inert scalar*)
RepGoldstoneA=CreateParticle[{{2,3},{2,1}},"S",mA2,"A"]; (*CP-odd inert and charged scalars *)


(*Light fermions*)
LightFermions=CreateParticle[{3,4,5,6,7,8,9,10,11,12,13,14,15},"F",mq2,"LightFermions"];


ParticleList={
	ReptL,RepbL,ReptR,
	RepGluon,RepW,
	RepHiggsh,RepGoldstoneGp0,RepGoldstoneGpR,RepGoldstoneGpI,
	RepHiggsH,RepGoldstoneA};
LightParticleList={LightFermions};


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
		Replacements->{lam4H->0,lam5H->0},
		Format->{"json","txt"},
		NormalizeWithDOF->False}];

