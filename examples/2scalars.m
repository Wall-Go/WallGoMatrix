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
RepW=CreateParticle[{{1,1}},"V"]; (*SU2 gauge bosons*)


(*Scalars bosons*)
RepHiggsh=CreateParticle[{{1,2}},"S"]; (*Higgs*)
RepGoldstoneGpR={{1},"S"}; (*real charged Goldstone*)
RepGoldstoneGpI={{3},"S"}; (*imag charged Golstone*)
RepGoldstoneGp0={{4},"S"}; (*neutral Goldstone*)
RepHiggsH=CreateParticle[{{2,2}},"S"]; (*CP-even inert scalar*)
RepGoldstoneA=CreateParticle[{{2,3}},"S"]; (*CP-odd inert scalar*)
RepGoldstoneHpR={{5},"S"}; (*real charged inert scalar*)
RepGoldstoneHpI={{7},"S"}; (*imag charged inert scalar*)


(*Defining various masses and couplings*)
VectorMass=Join[
	Table[mW2,{i,1,RepW[[1]]//Length}]
	];
FermionMass={};
ScalarMass={mG2,mh2,mG2,mG2,mHp,mH2,mHp,mA2};
ParticleMasses={VectorMass,FermionMass,ScalarMass};

UserCouplings=Variables@Normal@{Ysff,gvss,gvff,gvvv,\[Lambda]4,\[Lambda]3}//DeleteDuplicates;


ParticleList={
	RepHiggsh,RepGoldstoneGp0,RepGoldstoneGpR,RepGoldstoneGpI,
	RepHiggsH,RepGoldstoneA,RepGoldstoneHpR,RepGoldstoneHpI};
ParticleName={
	"Higgs","GoldstoneG0","GoldstoneGpR","GoldstoneGpI",
	"H","A","GoldstoneHpR","GolstoneHpI"};


(*
	output of matrix elements
*)
OutputFile="output/matrixElements.2scalars";

MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	UserCouplings,
	ParticleName,
	ParticleMasses,
	{
		TruncateAtLeadingLog->True,
		Replacements->{gw->0,lam4H->0,lam5H->0},
		Format->{"json","txt"},
		NormalizeWithDOF->False}];

