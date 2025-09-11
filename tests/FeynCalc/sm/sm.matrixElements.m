(* ::Package:: *)

(* ::Section:: *)
(*Matrix elements  in the Standard Model*)


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
model="SMQCD";
modLoc={SM,UnitarySM};
genericLoc={Lorentz,UnitaryLorentz};
InitializeModel["SMQCD"];


particleTypes={
	F[3,{3}](* t-quark *),
	-F[3,{3}](* anti-t-quark *),
	F[4,{3}](* b-quark *),
	-F[4,{3}](* anti-b-quark *),
	F[2,{3}](* tauon *),
	-F[2,{3}](* anti-tauon *),
	F[1,{3}](* tau-neutrino *),
	-F[1,{3}](* anti-tau-neutrino *),
	(*F[3,{2}](* charm *),
	-F[3,{2}](* anti-charm *),
	F[4,{2}](* strange *),
	-F[4,{2}](* anti-strange *),
	F[2,{2}](* muon *),
	-F[2,{2}](* anti-muon *),
	F[1,{2}](* mu-neutrino *),
	-F[1,{2}](* anti-mu-neutrino *),
	F[3,{1}](* up *),
	-F[3,{1}](* anti-up *),
	F[4,{1}](* down *),
	-F[4,{1}](* anti-down *),
	F[2,{1}](* electron *),
	-F[2,{1}](* anti-electron *),
	F[1,{1}](* electron-neutrino *),
	-F[1,{1}](* anti-electron-neutrino *),*)
	V[5](* g *),
	V[3](* W- *),
	-V[3](* W+ *),
	V[2](* Z *),
	V[1](* photon *),
	S[1](* physical Higgs *),
	S[2](* uncharged Goldstone *),
	S[3](* charged Goldstone*),
	-S[3](* charged anti-Goldstone*)
};

particleNames=Map[ToString[TheLabel[#]]<>If[Head[#]===Times,"bar",""]&,particleTypes]
makeProcesses[firstParticle_,particles_]:=DeleteDuplicates[Flatten[Table[{firstParticle,a}->{b,c},{a,particles},{b,particles},{c,particles}],2]];


masses={
	SMP["m_e"],SMP["m_mu"],SMP["m_tau"],
	SMP["m_u"],SMP["m_d"],SMP["m_c"],
	SMP["m_s"],SMP["m_t"],SMP["m_b"],
	SMP["m_H"],SMP["m_W"],SMP["m_Z"],
	SMP["m_q"],SMP["m_Q"],SMP["m_qu"],
	SMP["m_qd"],SMP["m_l"],SMP["m_pi"]
};
zeroMasses=Thread[masses->0]


(* getting expressions for parameters *)
SMP[];


(* this is where the FeynArts files made by FeynRules can be found, ending in .gen and .mod *)
modLoc=FileNameJoin[{NotebookDirectory[],model}];
genericLoc=FileNameJoin[{NotebookDirectory[],model}];
InitializeModel[modLoc,GenericModel->genericLoc];


Get[modLoc<>".pars"]
M$ExtParams//MatrixForm
parameters=Map[#[[1]]&,M$ExtParams]


(* sub out Weinberg angle *)
removeExplicitWeinbergAngle[arg_]:=arg/.{SMP["e"]->SMP["g'_W"]SMP["cos_W"]}/.{SMP["sin_W"]->SMP["e"]/SMP["g_W"],SMP["cos_W"]->SMP["e"]/SMP["g'_W"]}


(* sub out vev *)
removeVev[arg_]:=arg/.v->0(*MH/Sqrt[2lam]*);


(* remove hypercharge *)
removeHyperCharge[arg_]:=Expand[
removeExplicitWeinbergAngle[
Collect[
arg,
{SMP["m_W"],SMP["m_Z"]}]
]
/.{SMP["e"]->SMP["g'_W"](* holds to leading order in Subscript[\[Theta], W]*)}
]/.{SMP["g'_W"]->0}


(*mapParticleToInteger:={S[1]->0,F[1]->1,-F[1]->2};*)
mapParticleToInteger:=Thread[particleTypes->(Range[Length[particleTypes]]-1)];
mapParticleToInteger
makeMName[process_]:=M[process[[1,1]],process[[1,2]],process[[2,1]],process[[2,2]]]/.mapParticleToInteger


(* ::Section:: *)
(*Automating steps*)


tops=CreateTopologies[0,2->2];
massFn[particle_]:=0;


(* converting between particle index types *)
momenta={Momentum[p1,4],Momentum[p2,4],Momentum[p3,4],Momentum[p4,4]};
particleIndices={{1,1},{1,2},{2,1},{2,2}};
integerToParticleIndex[integer_]:=particleIndices[[integer]];
particleIndexToInteger[particleIndex_]:=Position[particleIndices,particleIndex][[1]][[1]];


(* polarization sums *)
polsums[x_,vec_,aux_,spinfac_]:=
	x//Collect2[#,Pair[_,Momentum[Polarization[vec,__]]]]&//Isolate[#,{Polarization[vec,__]}]&//DoPolarizationSums[#,vec,aux,ExtraFactor->spinfac]&//FixedPoint[ReleaseHold,#]&


insertCouplingsForMasses={
	SMP["m_e"]->0,
	SMP["m_mu"]->0,
	SMP["m_tau"]->0,
	SMP["m_u"]->0,
	SMP["m_d"]->0,
	SMP["m_c"]->0,
	SMP["m_s"]->0,
	SMP["m_t"]->yt v,
	SMP["m_b"]->0,
	SMP["m_H"]->Sqrt[2 lam]v,
	SMP["m_W"]->1/2 SMP["g_W"]v,
	SMP["m_Z"]->1/2 Sqrt[SMP["g_W"]^2+SMP["g'_W"]^2]v
};


ClearAll[makeAmplitude]
makeAmplitude[process_,channels_:All]:=Module[
	{
		diags,ampFA,ampFC,ampSq, ampMsq,ampMsqSummed,masses,transverse,i
	},
(* assumes 2->2 topology *)
(* create diagrams for process *)
(*diags=InsertFields[tops[[channels]],process,InsertionLevel->{Classes},Model->modLoc,GenericModel->modLoc]*)
If[model=="SM",
diags=InsertFields[tops[[channels]],process,InsertionLevel->{Particles},Model->modLoc,GenericModel->modLoc];
,
diags=InsertFields[tops[[channels]],process,InsertionLevel->{Particles},Model->"SMQCD"];
];

(* create FeynArts amplitude from diagrams *)
ampFA[1]=CreateFeynAmp[diags,PreFactor->1,
	Truncated->False,
	GaugeRules->{
(*		GaugeXi[S[___]]->1,
		GaugeXi[V[___]]->1,
		FAGaugeXi[S[___]]->1,
		FAGaugeXi[V[___]]->1*)
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
	ChangeDimension->4,
	List->False,
	SMP->True,
	Contract->True,
	DropSumOver->True,
	UndoChiralSplittings->True,
	TransversePolarizationVectors->transverse
];
ampFC[2]=ampFC[1](*/.PropagatorDenominator[x_,y_]\[RuleDelayed]FeynAmpDenominator[PropagatorDenominator[x,y]]*);
ampFC[3]=ampFC[2]/.Conjugate[PolarizationVector][_,a_,b_]->Conjugate[PolarizationVector[a,b]]/.PolarizationVector[_,a_,b_]->PolarizationVector[a,b];
ampFC[4]=ampFC[3]/.Spinor[a__]Spinor[b__]->Spinor[a] . Spinor[b];
ampFC[5]=removeVev[ampFC[4]/.insertCouplingsForMasses/.v->0];
ampFC[6]=removeHyperCharge[ampFC[5]];
Return[ampFC[6]//Contract]
]


ClearAll[makeAmplitudeSquared]
makeAmplitudeSquared[process_,channels_:All,dropUnphysicalPolarizations_:True]:=
Module[
	{
		ampM,ampSq, ampMsq,masses,stu,i,j,indexPair
	},
(* assumes 2->2 topology *)
masses={0,0,0,0};
stu=Plus@@(masses^2);
ampM[1]=makeAmplitude[process,channels]
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
ampSq[10]=SUNSimplify[ampSq[9],SUNNToCACF->False]/.SUNN->3;
(* FINAL SIMPLIFICATIONS *)
ampMsq[1]=PropagatorDenominatorExplicit[ampSq[10]];
(* simplify *)
ampMsq[3]=Expand[TrickMandelstam[ampMsq[1],{s,t,u,Plus@@(masses)}]];
ampMsq[4]=Simplify[ampMsq[3]/.insertCouplingsForMasses,vev>0]/.vev->0;
FCClearScalarProducts[];
Return[ampMsq[4]]
]


(* ::Section:: *)
(*Tests*)


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})
fixConvention[arg_]:=symmetriseTU[arg/.{s->(-t-u)}]//Expand//Simplify//Expand


makeAmplitude[{S[1],S[1]}->{S[1],S[1]}];
makeAmplitudeSquared[{S[1],S[1]}->{S[1],S[1]}]


makeAmplitude[{F[3,{1}],-F[3,{1}]}->{V[5],V[5]}];
makeAmplitudeSquared[{F[3,{1}],-F[3,{1}]}->{V[5],V[5]},All,True]
(* only the former agrees with the literature *)
makeAmplitudeSquared[{F[3,{1}],-F[3,{1}]}->{V[5],V[5]},All,False]


(* ::Subsection:: *)
(*VV -> VV*)


externalSignature={V[5],V[5]}->{V[5],V[5]};
makeAmplitude[externalSignature];
ampSq[1]=makeAmplitudeSquared[externalSignature]


ampSq[1]/.{SUNN->3}//SUNSimplify


resAMYSUN=16 dA CA^2 g^4 (3-(s u)/t^2-(s t)/u^2-(t u)/s^2);
resAMYSU3=resAMYSUN/.{dA->(CA^2-1)}/.{CA->3}


fixConvention[ampSq[1]/.{SUNN->3}]


(* should be zero - pure gauge amplitude unaffected by presence of scalars or fermions *)
fixConvention[ampSq[1]-resAMYSU3]/.{SMP["g_s"]->g}


(*WW -> WW*)
externalSignature={V[3],V[3]}->{V[3],V[3]};
makeAmplitude[externalSignature];
makeAmplitudeSquared[externalSignature]


(* ::Subsection:: *)
(*SS->SS*)


externalSignature={S[1],S[1]}->{S[1],S[1]};
makeAmplitude[externalSignature];
makeAmplitudeSquared[externalSignature]


(* ::Subsection:: *)
(*FF->FF*)


externalSignature={F[3,{3}],-F[3,{3}]}->{F[4,{3}],-F[4,{3}]};
makeAmplitude[externalSignature];
makeAmplitudeSquared[externalSignature]


(* ::Subsection:: *)
(*FF->VV*)


externalSignature={F[3,{1}],-F[3,{1}]}->{V[5],V[5]};
makeAmplitudeSquared[externalSignature,All,True]
(* only the former agrees with the literature *)
makeAmplitudeSquared[externalSignature,All,False]


(* ::Subsection:: *)
(*FV->FV*)


externalSignature={F[3,{1}],V[5]}->{F[3,{1}],V[5]};
makeAmplitude[externalSignature];
makeAmplitudeSquared[externalSignature]


(* ::Subsection:: *)
(*SS->FF*)


externalSignature={S[1],S[1]}->{F[3,{3}],-F[3,{3}]};
makeAmplitude[externalSignature];
makeAmplitudeSquared[externalSignature]


(* ::Subsection:: *)
(*SF->VF*)


externalSignature={F[3,{3}],-F[3,{3}]}->{S[1],V[2]};
makeAmplitude[externalSignature];
makeAmplitudeSquared[externalSignature]


particleNames
particleTypes


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
    "gs" -> "gs", "yt" -> "yt",
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


(* ::Section:: *)
(*Results*)


removeIncomplete={Pair[x__]->0};
removeIncomplete={ };


mapParticleToName[particleList_]:=Map[ToString[TheLabel[#]]<>If[Head[#]===Times,"bar",""]&,particleList]
mapParticleToInteger:=Thread[particleTypes->(Range[Length[particleTypes]]-1)];


(* squared summed matrix elements with a vector on leg 1 *)
createProcesses[particle_,ingoingParticleList_]:=Module[{processes},
	processes=makeProcesses[particle,ingoingParticleList];
	DeleteCases[
	Table[
		makeMName[process]->makeAmplitudeSquared[process,All,True]/.removeIncomplete,
	{process,processes}],_->0]
];

runProccesses[ingoingParticleList_]:=Module[{i = 0, n = Length[ingoingParticleList], results},
  results = Monitor[
    Table[
      i++;
      createProcesses[particle, ingoingParticleList],
      {particle, ingoingParticleList}
    ],
    Row[{"Processing ", ingoingParticleList[[i]], " which is particle ", i, " of ", n, " (", NumberForm[100. (i-1)/n, {3, 1}], "%) "}]
  ] // Flatten[#, 1] &;
  results
];


(* ::Subsection:: *)
(*SU(3) Yang-Mills (on external legs)*)


particleTypes = {
	V[5](* g *)
};
particleNames=mapParticleToName[particleTypes];
results = runProccesses[particleTypes];


results
(M[0,0,0,0]/.results/.{SMP["g_s"]->g})//Simplify//fixConvention
resAMYSU3-%//fixConvention


resultsExport=results/.{SMP["g_s"]->gs};
feynAssociation=Association[resultsExport];


(*
	here for the Yukawa model,where the matrix elements are stored in
	the object called MatrixElements
*)
exportParticles = "g";
exportParameters={gs};
toExportAsJSON=makeJsonObject[particleNames,Join[exportParameters,{SUNN}],results//SMPToSymbol];
Export[FileNameJoin[{NotebookDirectory[],StringJoin[model,".",exportParticles,".test.json"]}],toExportAsJSON];


(* ::Subsection:: *)
(*QCD with one generation of quarks (on external legs)*)


particleTypes = {
	F[3,{3}](* t-quark *),
	-F[3,{3}](* anti-t-quark *),
	F[4,{3}](* b-quark *),
	-F[4,{3}](* anti-b-quark *),
	V[5](* g *)
};
particleNames=mapParticleToName[particleTypes];
results = runProccesses[particleTypes];


resultsExport=results/.{SMP["g_W"]->gw,SMP["g_s"]->gs,SMP["y_t"]->yt};
feynAssociation=Association[resultsExport];


(*
	here for the Yukawa model,where the matrix elements are stored in
	the object called MatrixElements
*)
exportParticles = "tbg";
exportParameters = {gw,gs,yt};
toExportAsJSON=makeJsonObject[particleNames,Join[exportParameters,{SUNN}],resultsExport//SMPToSymbol];
Export[FileNameJoin[{NotebookDirectory[],StringJoin[model,".",exportParticles,".test.json"]}],toExportAsJSON];


(* ::Subsection:: *)
(*Higgs-tau sector (Subscript[y, tau]=0)*)


particleTypes = {
	(*F[3,{3}](* t-quark *),
	-F[3,{3}](* anti-t-quark *),
	F[4,{3}](* b-quark *),
	-F[4,{3}](* anti-b-quark *),*)
	F[2,{3}](* tauon *),
	-F[2,{3}](* anti-tauon *),
	(*F[1,{3}](* tau-neutrino *),
	-F[1,{3}](* anti-tau-neutrino *),*)
	S[1](* physical Higgs *),
	S[2](* uncharged Goldstone *),
	S[3](* charged Goldstone*),
	-S[3](* charged anti-Goldstone*)
};
results = runProccesses[particleTypes];


resultsExport=results/.{SMP["g_W"]->gw,SMP["g_s"]->gs,SMP["y_t"]->yt};
feynAssociation=Association[resultsExport];


(*
	here for the Yukawa model,where the matrix elements are stored in
	the object called MatrixElements
*)
particleNames=mapParticleToName[particleTypes];
particleNames=(* by hand *){"tau","taubar","H","G0","Gp","Gpbar"}
exportParticles = "HiggsTau";
exportParameters = {gw,gs,yt,lam};
toExportAsJSON=makeJsonObject[particleNames,Join[exportParameters,{SUNN}],resultsExport//SMPToSymbol];
Export[FileNameJoin[{NotebookDirectory[],StringJoin[model,".",exportParticles,".test.json"]}],toExportAsJSON];


(* ::Subsection:: *)
(*Higgs-top sector (Subscript[y, top]!=0)*)


particleTypes = {
	F[3,{3}](* t-quark *),
	-F[3,{3}](* anti-t-quark *),
	(*F[4,{3}](* b-quark *),
	-F[4,{3}](* anti-b-quark *),*)
	(*F[2,{3}](* tauon *),
	-F[2,{3}](* anti-tauon *),*)
	(*F[1,{3}](* tau-neutrino *),
	-F[1,{3}](* anti-tau-neutrino *),*)
	S[1](* physical Higgs *),
	S[2](* uncharged Goldstone *),
	S[3](* charged Goldstone*),
	-S[3](* charged anti-Goldstone*)
};

results = runProccesses[particleTypes];


resultsExport=results/.{SMP["g_W"]->gw,SMP["g_s"]->gs,SMP["y_t"]->yt};
feynAssociation=Association[resultsExport];


(*
	here for the Yukawa model,where the matrix elements are stored in
	the object called MatrixElements
*)
particleNames=mapParticleToName[particleTypes];
particleNames=(* by hand *){"t","tbar","H","G0","Gp","Gpbar"}
exportParticles = "HiggsTop";
exportParameters = {gw,gs,yt,lam};
toExportAsJSON=makeJsonObject[particleNames,Join[exportParameters,{SUNN}],resultsExport//SMPToSymbol];
Export[FileNameJoin[{NotebookDirectory[],StringJoin[model,".",exportParticles,".test.json"]}],toExportAsJSON];


(* ::Subsection:: *)
(*Electroweak bosons sector*)


particleTypes = {
	V[3](* W- *),
	-V[3](* W+ *),
	V[2](* Z *),
	V[1](* photon *),
	S[1](* physical Higgs *),
	S[2](* uncharged Goldstone *),
	S[3](* charged Goldstone*),
	-S[3](* charged anti-Goldstone*)
};
results = runProccesses[particleTypes];


resultsExport=results/.{SMP["g_W"]->gw,SMP["g_s"]->gs,SMP["y_t"]->yt};
feynAssociation=Association[resultsExport];


(*
	here for the Yukawa model,where the matrix elements are stored in
	the object called MatrixElements
*)
exportParticleNames=mapParticleToName[particleTypes]
exportParticleNames=(* by hand *){"W","Wbar","Z","gamma","H","G0","Gp","Gpbar"}
exportParticles = "EWbosons";
exportParameters = {gw,gs,yt,lam};
toExportAsJSON=makeJsonObject[exportParticleNames,Join[exportParameters,{SUNN}],resultsExport//SMPToSymbol];
Export[FileNameJoin[{NotebookDirectory[],StringJoin[model,".",exportParticles,".test.json"]}],toExportAsJSON];

