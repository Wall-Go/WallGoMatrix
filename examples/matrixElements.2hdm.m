(* ::Package:: *)

Quit[];


SetDirectory[NotebookDirectory[]];
(*Put this if you want to create multiple model-files with the same kernel*)
$GroupMathMultipleModels=True;
$LoadGroupMath=True;
<<../src/WallGoMatrix.m


(* ::Chapter:: *)
(*2HDM*)


(*See 2211.13142 for implementation details*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
HiggsDoublet1={{{0,0},{1},1/2},"C"};
HiggsDoublet2={{{0,0},{1},1/2},"C"};
RepScalar={HiggsDoublet1,HiggsDoublet2};
CouplingName={g3,gw,g1};


Rep1={{{1,0},{1},1/6},"L"};
Rep2={{{1,0},{0},2/3},"R"};
Rep3={{{1,0},{0},-1/3},"R"};
Rep4={{{0,0},{1},-1/2},"L"};
Rep5={{{0,0},{0},-1},"R"};
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
SymmetryBreaking[vev,VevDependentCouplings->True] (*uncomment if you want vev-dependent couplings*)
(*SymmetryBreaking[vev]*)


(*Third generation of fermions*)
ReptL=CreateParticle[{{1,1}},"F"];
RepbL=CreateParticle[{{1,2}},"F"];
ReptR=CreateParticle[{{2,1}},"F"];


(*Vector bosons*)
RepGluon=CreateParticle[{1},"V"]; (*Gluons*)
RepW=CreateParticle[{{2,1}},"V"]; (*SU2 gauge bosons*)
RepB=CreateParticle[{3},"V"]; (*U1 gauge boson*)


(*Scalars bosons*)
RepHiggs1=CreateParticle[{{1,2}},"S"]; (*Higgs*)
RepGoldstone1=CreateParticle[{{1,1}},"S"]; (*Goldstone*)
RepHiggs2=CreateParticle[{{2,1}},"S"]; (*CP-even inert scalar*)
RepGoldstone21=CreateParticle[{{2,2}},"S"]; (*CP-odd inert scalar*)
RepGoldstone22=CreateParticle[{{2,3}},"S"]; (*charged inert scalar*)


(*left+right top-quark*)
Rept=CreateParticle[{{1,1},{2,1}},"F"];

(*left-handed bottom-quark*)
Repb=CreateParticle[{{1,2}},"F"];

(*scalar reps*)
Reph=CreateParticle[{{1,2}},"S"];
RepG=CreateParticle[{{1,1}},"S"];

RepH=CreateParticle[{{2,2}},"S"];
RepA=CreateParticle[{{2,3}},"S"];
RepHpm=CreateParticle[{{2,1}},"S"];

(*Vector bosons*)
RepGluon=CreateParticle[{1},"V"];
RepW=CreateParticle[{{2,1}},"V"];
RepZ=CreateParticle[{{3,1}},"V"];


(*Defining various masses and couplings*)
VectorMass=Join[
	Table[mg2,{i,1,RepGluon[[1]]//Length}],
	Table[mw2,{i,1,RepW[[1]]//Length}],{mb2}]; (*mb2 is the mass of the U(1) gauge field*)
FermionMass=Table[mq2,{i,1,Length[gvff[[1]]]}];
ScalarMass={mG2,mh2,mG2,mG2,mH2,mA2,mHp,mHp};
ParticleMasses={VectorMass,FermionMass,ScalarMass};

UserMasses={mq2,mg2,mw2,mb2,mG2,mh2,mH2,mA2,mHp};
UserCouplings=Variables@Normal@{Ysff,gvss,gvff,gvvv,\[Lambda]4,\[Lambda]3}//DeleteDuplicates;


ParticleList={
	ReptL,RepbL,ReptR,
	RepGluon,RepW,RepB,
	RepHiggs1,RepGoldstone1,
	RepHiggs2,RepGoldstone21,RepGoldstone22};
ParticleName={
	"TopL","BotL","TopR",
	"Gluon","W","B",
	"Higgs","Goldstone",
	"Higgs2","Goldstone21","Goldstone22"};


(*
	output of matrix elements
*)
OutputFile="matrixElements.2hdm";
SetDirectory[NotebookDirectory[]];

MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	UserMasses,
	UserCouplings,
	ParticleName,
	ParticleMasses,
	{
		TruncateAtLeadingLog->True,
		Replacements->{lam4H->0,lam5H->0},
		Format->{"json","txt"},
		NormalizeWithDOF->False}];
