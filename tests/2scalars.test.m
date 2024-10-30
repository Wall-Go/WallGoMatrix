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
		NormalizeWithDOF->False}];


