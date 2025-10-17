(* ::Package:: *)

(* ::Section:: *)
(*Matrix elements  in SU(M) + SU(N) + Higgs + Yukawa theory*)


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


(* ::Subsection::Closed:: *)
(*Choosing model*)


(* choose model *)
model = "sum-sun-higgs-yukawa_replaced";


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


particleTypes


particleTypes=DeleteCases[F$Classes,U[1]|-U[1]]
makeProcesses[firstParticle_]:=DeleteDuplicates[Flatten[
	Table[{firstParticle,a}->{b,c},
	{a,particleTypes},
	{b,particleTypes},
	{c,particleTypes}],2]];
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
makeProcesses[firstParticle_]:=Flatten[
	Table[{firstParticle,a}->bc,
		{a,{S[1],-S[1],F[1],-F[1],F[2],-F[2],V[1]}},
		{bc,{{S[1],-S[1]},{S[1],V[1]},{-S[1],V[1]},{V[1],V[1]},{F[1],-F[2]}}}
		],1]
makeProcesses[S[1]]//MatrixForm
makeProcesses[F[1]]//MatrixForm
makeProcesses[V[1]]//MatrixForm
processesByHand={
	{S[1],S[1]}->{S[1],S[1]},
	{S[1],S[1]}->{V[1],V[1]},
	{V[1],V[1]}->{V[1],V[1]}
	};


(* ::Section::Closed:: *)
(*Automating steps*)


tops=CreateTopologies[0,2->2];
massFn[particle_]:=0(*Which[particle===F[1]||particle===-F[1],mPsi,particle===S[1]||particle===-S[1],0,particle===V[1]||particle===-V[1],0]*);


(* converting between particle index types *)
momenta={Momentum[p1,4],Momentum[p2,4],Momentum[p3,4],Momentum[p4,4]};
particleIndices={{1,1},{1,2},{2,1},{2,2}};
integerToParticleIndex[integer_]:=particleIndices[[integer]];
particleIndexToInteger[particleIndex_]:=Position[particleIndices,particleIndex][[1]][[1]];


(* polarization sums *)
polsums[x_,vec_,aux_,spinfac_]:=
	x//Collect2[#,Pair[_,Momentum[Polarization[vec,__]]]]&//Isolate[#,{Polarization[vec,__]}]&//DoPolarizationSums[#,vec,aux,ExtraFactor->spinfac]&//FixedPoint[ReleaseHold,#]&


interstingParticles={S[1],-S[1],V[1],V[2],F[1],-F[1],F[2],-F[2]};
makeProcesses[firstParticle_]:=Flatten[Table[{firstParticle,a}->{b,c},
	{a,interstingParticles},
	{b,interstingParticles},
	{c,interstingParticles}]
	,2];


ContractMultipleSUN[expression_]:=Module[{exprI=expression},
	exprI=exprI/.{
		SUNFDelta[
			SUNFIndex[Index[Fundamental1, s1_]],
			SUNFIndex[Index[Fundamental1, s2_]]
			] -> 
		SUN1[SUNFDelta[
			SUNFIndex[Index[Fundamental1, s1]],
			SUNFIndex[Index[Fundamental1, s2]]
			]],
		SUNFDelta[
			SUNFIndex[Index[Fundamental2, s1_]],
			SUNFIndex[Index[Fundamental2, s2_]]
			] -> 
		SUN2[SUNFDelta[
			SUNFIndex[Index[Fundamental2, s1]],
			SUNFIndex[Index[Fundamental2, s2]]
			]],
		SUNF1[x__]->SUN1[SUNF[x]],
		SUNTF1[x__]->SUN1[SUNTF[x]],
		SUNF2[x__]->SUN2[SUNF[x]],
		SUNTF2[x__]->SUN2[SUNTF[x]]
	}//SUNSimplify;
	exprI
]


ClearAll[makeAmplitude]
makeAmplitude[process_,channels_:All]:=Module[
	{
		diags,ampFA,ampFC,ampSq, ampMsq,ampMsqSummed,masses,transverse,i
	},
(* assumes 2->2 topology *)
(* create diagrams for process *)
diags=InsertFields[tops[[channels]],process,InsertionLevel->{Classes},Model->modLoc,GenericModel->modLoc];

(* create FeynArts amplitude from diagrams *)
ampFA[1]=CreateFeynAmp[diags,PreFactor->1,
	Truncated->False,
	GaugeRules->{
		GaugeXi[S[___]]->1,
		GaugeXi[V[___]]->1,
		FAGaugeXi[S[___]]->1,
		FAGaugeXi[V[___]]->1
	}];

(* introducing FeynAmpDenominator, as required by FeynCalc (not sure why this is missing in the first place) *)
ampFA[2]=ampFA[1]/.FAPropagatorDenominator[x_,y_]:>FAFeynAmpDenominator[PropagatorDenominator[x,y]];

(* converting to FeynCalc form *)
transverse={};
For[i=1,i<=4,i++,
	If[process[[particleIndices[[i]]/.List->Sequence]]==V[1],
	AppendTo[transverse,momenta[[i]]]
]];
ampFC[1]=FCFAConvert[ampFA[2],
	IncomingMomenta->{Momentum[p1,4],Momentum[p2,4]},
	OutgoingMomenta->{Momentum[p3,4],Momentum[p4,4]},
	LorentzIndexNames->{\[Mu],\[Nu],\[Rho],\[Sigma]},
	(*SUNIndexNames->{a1,a2,a3,a4,a5,a6},
	SUNFIndexNames->{f1,f2,f3,f4,f5,f6},*)
	ChangeDimension->4,
	List->False,
	SMP->False,
	Contract->False,
	DropSumOver->True,
	UndoChiralSplittings->True,
	TransversePolarizationVectors->transverse
];
ampFC[2]=ampFC[1](*/.PropagatorDenominator[x_,y_]\[RuleDelayed]FeynAmpDenominator[PropagatorDenominator[x,y]]*);
ampFC[3]=ampFC[2]/.Conjugate[PolarizationVector][_,a_,b_]->Conjugate[PolarizationVector[a,b]]/.PolarizationVector[_,a_,b_]->PolarizationVector[a,b];
ampFC[4]=ampFC[3](*/.Spinor[a__]Spinor[b__]->Spinor[a] . Spinor[b]*);
ampFC[5]=ampFC[4](*/.{SUNF[a_,b_,c_,d_] -> SUNF[a,b,Index[Adjoint,5]] SUNF[Index[Adjoint,5],c,d]} *);
(*/.{SUNT[SUNIndex[Index[Adjoint, a_]]] . SUNT[SUNIndex[Index[Fundamental, \[Beta]_]]] . SUNT[SUNIndex[Index[Fundamental, \[Gamma]_]]]->SUNTF[{SUNIndex[adj[a]]},SUNFIndex[fun[\[Beta]]],SUNFIndex[fun[\[Gamma]]]]} ;*)
(*Print[ampFC[5]//Contract//Expand];*)

ampFC[5]=ContractMultipleSUN[ampFC[5]];
(*ampFC[5]=ampFC[5]/.{
	SUNFDelta[
		SUNFIndex[Index[Fundamental2, s1_]],
		SUNFIndex[Index[Fundamental2, s2_]]
		] -> 
	SUN2[SUNFDelta[
		SUNFIndex[Index[Fundamental2, s1]],
		SUNFIndex[Index[Fundamental2, s2]]
		]],
	SUNFDelta[
		SUNFIndex[Index[Fundamental1, s1_]],
		SUNFIndex[Index[Fundamental1, s2_]]
		] -> 
	SUN1[SUNFDelta[
		SUNFIndex[Index[Fundamental1, s1]],
		SUNFIndex[Index[Fundamental1, s2]]
		]],
	SUNF1[x__]->SUN1[SUNF[x]],
	SUNTF1[x__]->SUN1[SUNTF[x]],
	SUNF2[x__]->SUN2[SUNF[x]],
	SUNTF2[x__]->SUN2[SUNTF[x]]}//SUNSimplify;*)
	
Return[ampFC[5]//Contract]
]


(* ::Subsubsection:: *)
(*Amplitude squared*)


ClearAll[makeAmplitudeSquared]
makeAmplitudeSquared[process_,channels_:All,dropUnphysicalPolarizations_:True]:=
Module[
	{
		ampM,ampSq, ampMsq,masses,stu,i,j,indexPair
	},
(* assumes 2->2 topology *)
masses={0,0,0,0};
stu=Plus@@(masses^2);
ampM[1]=makeAmplitude[process,channels];

(*/.PropagatorDenominator[a_,b_]->PropagatorDenominator[a,0]*);
(* making amplitude squared *)
ampSq[1]=(*1/2*)(* 1/2 is a symmetry factor *)ComplexConjugate[ampM[1]]ampM[1];
ampSq[2]=FeynAmpDenominatorExplicit[ampSq[1]];
ampSq[3]=ampSq[2];
FCClearScalarProducts[];

(* MANDELSTAM VARIABLES *)
(*FCSetScalarProducts[{SPD[momenta[[1]]],SPD[momenta[[2]]],SPD[momenta[[3]]],SPD[momenta[[4]]]},{0,0,0,0}];*)
SetMandelstam[s,t,u,p1,p2,-p3,-p4,(masses/.List->Sequence)];
(*For[i=1,i<=4,i++,
(* setting on-shell momenta *)
ScalarProduct[momenta[[i]],momenta[[i]]]:=masses[[i]];
];*)
(* SIMPLIFYING SPINOR INDICES *)
ampSq[4]=FermionSpinSum[ampSq[3]]/.DiracTrace[aaa__]->DiracTrace[aaa,Mandelstam->{s,t,u,stu}];
ampSq[5]=DiracSimplify[ampSq[4]];

(* SIMPLIFYING VECTOR INDICES *)
For[i=1,i<=4,i++,
	For[j=1,j<=4,j++,
		ampSq[i+5]=ampSq[i+4];
		(* summing over polarizations *)
		If[j==i,
			If[
				Not[FreeQ[process[[integerToParticleIndex[i]/.List->Sequence]],V[_]]],
				ampSq[i+5]=ampSq[i+4]//polsums[#,momenta[[i]],0,1]&;
				Break[];
			];
			Continue[];
		];
		
		If[
			Not[FreeQ[process[[integerToParticleIndex[i]/.List->Sequence]],V[_]]]&&
			Not[FreeQ[process[[integerToParticleIndex[j]/.List->Sequence]],V[_]]],
			
			ampSq[i+5]=ampSq[i+4]//
				If[dropUnphysicalPolarizations,
					polsums[#,momenta[[i]],momenta[[j]],1]&,
					polsums[#,momenta[[i]],0,1]&
				];
			Break[];
		];
	];
];
(* SIMPLIFYING COLOUR INDICES *)
ampSq[10]=ampSq[9];
(*ampSq[11]=SUNSimplify[ampSq[10],SUNNToCACF->False](*/.SUNN->2*);*)
ampSq[11]=SUNSimplify[ampSq[10]/.{SUN1[x__]->x},SUNNToCACF->False]/.{SUNN->SUNN1,CA->CA1,CF->CF1};
ampSq[11]=SUNSimplify[ampSq[11]/.{SUN2[x__]->x},SUNNToCACF->False]/.{SUNN->SUNN2,CA->CA2,CF->CF2};
(*ampSq[11]=ContractMultipleSUN[ampSq[10]];*)
(*ampSq[11]=ampSq[10];*)
(*Print[ampSq[11]];*)
(* FINAL SIMPLIFICATIONS *)
ampMsq[1]=PropagatorDenominatorExplicit[ampSq[11]];
(* simplify *)
ampMsq[3]=Expand[TrickMandelstam[ampMsq[1],{s,t,u,Plus@@(masses)}]];
FCClearScalarProducts[];
Return[ampMsq[3]]
]


(* ::Section:: *)
(*Tests*)


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})
fixConvention[arg_]:=symmetriseTU[arg/.{s->(-t-u)}]//Expand//Simplify//Expand


(* ::Subsection:: *)
(*VV -> VV*)


externalSignature={V[1],V[1]}->{V[1],V[1]};
(*externalSignature={V[2],V[2]}->{V[2],V[2]};*)


ampSq[1]=makeAmplitudeSquared[externalSignature,All,True];


ampSq[1]/.SUNF1->SUNF(*//InputForm*)//SUNSimplify


ampSq[1]/.{SUNN1->2}//SUNSimplify


resAMYSU3=16 dA CA^2 g^4 (3-(s u)/t^2-(s t)/u^2-(t u)/s^2)/.{dA->8,CA->3}
resAMYSU2=16 dA CA^2 g^4 (3-(s u)/t^2-(s t)/u^2-(t u)/s^2)/.{dA->3,CA->2}


fixConvention[ampSq[1]/.{SUNN1->3}]
fixConvention[ampSq[1]/.{SUNN1->2}]


(* should be zero - pure gauge amplitude unaffected by presence of scalars or fermions *)
fixConvention[ampSq[1]-resAMYSU3/.{SUNN1->3,g3->g}]
fixConvention[ampSq[1]-resAMYSU2/.{SUNN1->2,g3->g}]


(* ::Subsection:: *)
(*SS->SS*)


(*InsertFields[tops[[All]],externalSignature,InsertionLevel->{Classes},Model->modLoc,GenericModel->modLoc];
CreateFeynAmp[%,PreFactor->1,Truncated->False,
	GaugeRules->{GaugeXi[S[_]]->1,GaugeXi[V[_]]->1,FAGaugeXi[S[_]]->1,FAGaugeXi[V[_]]->1}]//InputForm*)


externalSignature={S[1],S[1]}->{S[1],S[1]};
ampSq[S]=makeAmplitudeSquared[externalSignature,All,True];


ampSq[S]/.{SUNN1->3}
ampSq[S]/.{SUNN1->5}


(* ::Subsection:: *)
(*FF->FF*)


externalSignature={F[3],-F[3]}->{F[3],-F[3]};
externalSignature={F[3],F[3]}->{F[3],F[3]};
externalSignature={F[2],-F[1]}->{-F[1],F[2]};
(*externalSignature={F[1],-F[1]}->{-F[1],F[1]};*)

(*InsertFields[tops[[All]],externalSignature,InsertionLevel->{Classes},Model->modLoc,GenericModel->modLoc];
CreateFeynAmp[%,PreFactor->1,Truncated->False,
	GaugeRules->{GaugeXi[S[_]]->1,GaugeXi[V[_]]->1,FAGaugeXi[S[_]]->1,FAGaugeXi[V[_]]->1}]*)


(*makeAmplitude[externalSignature,All]*)


ampSq["F"]=makeAmplitudeSquared[externalSignature,All,True]//Contract//SpinorChainEvaluate;
(*%//FullForm;*)
%//Contract//Expand;
Collect[%/.{SUNN1->3,SUNN2->5},{g3,g5,y}]


externalSignature={F[3],-F[4]}->{F[3],-F[4]};

InsertFields[tops[[All]],externalSignature,InsertionLevel->{Classes},Model->modLoc,GenericModel->modLoc];
CreateFeynAmp[%,PreFactor->1,Truncated->False,
	GaugeRules->{GaugeXi[S[_]]->1,GaugeXi[V[_]]->1,FAGaugeXi[S[_]]->1,FAGaugeXi[V[_]]->1}]


ampSq["F"]=makeAmplitudeSquared[externalSignature,All,True]//Contract//SUNSimplify//SpinorChainEvaluate
(*%//FullForm;*)
%//Contract//Expand
(*%/.{mPhi->0,mPsi->0,CA->2,CF->3/4}*)


externalSignature={F[2],-F[2]}->{F[2],-F[2]};
ampSq["F"]=makeAmplitudeSquared[externalSignature,All,True]//Contract//SUNSimplify//SpinorChainEvaluate
(*%//FullForm;*)
%//Contract//Expand
(*%/.{mPhi->0,mPsi->0,CA->2,CF->3/4}*)


(* ::Subsection:: *)
(*SS->FF*)


externalSignature={S[1],-S[1]}->{F[2],-F[2]};
externalSignature={S[1],-S[1]}->{-F[2],F[2]};
externalSignature={-S[1],S[1]}->{F[2],-F[2]};
externalSignature={-S[1],S[1]}->{-F[1],F[1]};


ampSq["SS->FF"]=makeAmplitudeSquared[externalSignature,All,True]//Contract//SUNSimplify//SpinorChainEvaluate;
(*%//FullForm;*)
%//Contract//Expand;
%/.{mPhi->0,mPsi->0,mChi->0,CA->2,CF->3/4}


signatures={
	{S[1],-S[1]}->{F[1],-F[1]},
	{S[1],-S[1]}->{-F[1],F[1]},
	{-S[1],S[1]}->{F[1],-F[1]},
	{-S[1],S[1]}->{-F[1],F[1]}	
};

ampSq["SS->FF"]=
Table[
	makeAmplitudeSquared[sig,All,True]//Contract//SUNSimplify//SpinorChainEvaluate,
	{sig,signatures}];
(*%//FullForm;*)
%//Contract//Expand;
%/.{mPhi->0,mPsi->0,mChi->0,CA->2,CF->3/4}
%//Total//Simplify
fixConvention[%]


(* ::Subsection:: *)
(*SF->VF*)


externalSignature={-S[1],F[1]}->{V[2],F[2]};
makeAmplitudeSquared[externalSignature,All,True]//Contract//SUNSimplify
(*%/.{Pair[x__]->0}*)

(*DoPolarizationSums[%,momenta[[3]],0,ExtraFactor->1]//Contract//SUNSimplify
SetMandelstam[s,t,u,p1,p2,-p3,-p4,masses/.List->Sequence];
%%//Simplify*)



(* ::Section:: *)
(*Results*)


removeIncomplete={Pair[x__]->0};
removeIncomplete={ };


particleList = {S[1],-S[1],V[1],V[2],F[1],-F[1],F[3],-F[3]};

(* squared summed matrix elements with a vector on leg 1 *)
createProcesses[particle_]:=Module[{processes},
	processes=makeProcesses[particle];
	DeleteCases[
	Table[
		makeMName[process]->Collect[makeAmplitudeSquared[process,All]/.{\[Lambda]->lam}/.removeIncomplete,
		{lam,g,y},Expand],
	{process,processes}],_->0]
];


results = Module[{i = 0, n = Length[particleList], results},
  results = Monitor[
    Table[
      i++;
      createProcesses[particle],
      {particle, particleList}
    ],
    Row[{"Processing ", particleList[[i]], " which is particle ", i, " of ", n, " (", NumberForm[100. (i-1)/n, {3, 1}], "%) "}]
  ] // Flatten[#, 1] &;
  results
];


feynAssociation=Association[results];


(* ::Section:: *)
(*Exporting to ascii for C++*)


(* writing as JSON *)
makeJsonObject[fields_, parameters_, results_] := 
 Module[{toString, replaceSpecials, getRelevantParameters, fieldsJson, matrixElementsJson},

  (* helpers *)
  toString[arg_] := If[StringQ[arg], arg, ToString[arg, InputForm]];
  replaceSpecials[arg_String] := 
    StringReplace[arg, {
    "mPsi" -> "mPsi",
    "Pi" -> "_pi",
    "s" -> "_s", "t" -> "_t", "u" -> "_u"}];
  getRelevantParameters[arg_] := 
    Select[parameters,Not[FreeQ[arg,#]]&];

  (* particle entries *)
  fieldsJson = Table[
    {"index" -> i - 1, "name" -> toString[fields[[i]]]},
    {i, Length[fields]}
  ];

  (* matrix element entries *)
  matrixElementsJson = Map[
    {
      "externalParticles" -> (#[[1]] /. M[a__] :> {a}),
     (* "parameters" -> Map[toString,getRelevantParameters[#[[2]]]],*)
      "parameters" -> (toString /@ getRelevantParameters[#[[2]]] ),
      "expression" -> replaceSpecials[toString[#[[2]]]]
    } &,
    results
  ];

  {
  "particles" -> fieldsJson,
  "matrixElements" -> matrixElementsJson
  }
]


(*
	here for the Yukawa model,where the matrix elements are stored in
	the object called MatrixElements
*)
toExportAsJSON=makeJsonObject[particleNames,Join[parameters,{SUNN,SUNN1,SUNN2}],results//SMPToSymbol];
Export[FileNameJoin[{NotebookDirectory[],StringJoin[model,".test.json"]}],toExportAsJSON];
