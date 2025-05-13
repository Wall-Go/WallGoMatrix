(* ::Package:: *)

(* ::Section:: *)
(*Matrix elements  in SU(2) Higgs + Yukawa theory*)


Quit[];


(* ::Subsection:: *)
(*Loading packages*)


$LoadAddOns={"FeynArts"};
Get["FeynCalc`"];
$FAVerbose=0;
FCCheckVersion[9,3,0];
SetOptions[FourVector,FeynCalcInternal->False];
$KeepLogDivergentScalelessIntegrals=True;


(*FCEnableTraditionalFormOutput[]*)


FCDisableTraditionalFormOutput[]


(* ::Subsection:: *)
(*Choosing model*)


(* choose model *)
model = "sun-higgs-yukawa";


(* power counting *)
insertx={
	\[Sigma]->x^2 \[Sigma],
	msq->x^2 msq,
	mPhiv->x mPhiv,
	g->x^2 g,
	lam->x^2 lam,
	mPsi->x mPsi,
	mPsiv->x mPsiv,
	mChi->x mChi,
	mChiv->x mChiv,
	y->x y};
powerx[arg_,n_]:=SeriesCoefficient[arg/.insertx,{x,0,n},Assumptions->x>0];


(* this is where the FeynArts files made by FeynRules can be found, ending in .gen and .mod *)
modLoc=FileNameJoin[{NotebookDirectory[],model}];
genericLoc=FileNameJoin[{NotebookDirectory[],model}];
InitializeModel[modLoc,GenericModel->genericLoc];


Get[modLoc<>".pars"]
M$ExtParams//MatrixForm
parameters=Map[#[[1]]&,M$ExtParams]


particleTypes=DeleteCases[F$Classes,U[1]|-U[1]]
makeProcesses[firstParticle_]:=DeleteDuplicates[Flatten[Table[{firstParticle,a}->{b,c},{a,particleTypes},{b,particleTypes},{c,particleTypes}],2]];
makeProcesses[S[1]];
%[[1;;Min[10,Length[%]]]]//MatrixForm


particleNames=Map[ToString[TheLabel[#]]<>If[Head[#]===Times,"bar",""]&,particleTypes]


(*mapParticleToInteger:={S[1]->0,F[1]->1,-F[1]->2};*)
mapParticleToInteger:=Thread[particleTypes->(Range[Length[particleTypes]]-1)];
mapParticleToInteger
makeMName[process_]:=M[process[[1,1]],process[[1,2]],process[[2,1]],process[[2,2]]]/.mapParticleToInteger


(* relevant processes for Boltzmann collision integrals with given particle in position p1 *)
makeProcesses[firstParticle_]:=Flatten[
	Table[{firstParticle,a}->bc,
		{a,{S[1],-S[1],V[1]}},
		{bc,{{S[1],-S[1]},{S[1],V[1]},{-S[1],V[1]},{V[1],V[1]}}}
		],1]
makeProcesses[S[1]]//MatrixForm
processesByHand={
	{S[1],S[1]}->{S[1],S[1]},
	{S[1],S[1]}->{V[1],V[1]},
	{V[1],V[1]}->{V[1],V[1]}
	};
