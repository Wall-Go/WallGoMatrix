(* ::Package:: *)

(* ::Section:: *)
(*Matrix elements  in a n-scalar model*)


Quit[];


(* ::Subsection:: *)
(*Loading packages*)


$LoadAddOns={"FeynArts"};
Get["FeynCalc`"];
$FAVerbose=0;
SetOptions[FourVector,FeynCalcInternal->False];


FCDisableTraditionalFormOutput[]


(* ::Subsection:: *)
(*Choosing model*)


(* choose model *)
model="n-scalars";
modLoc=FileNameJoin[{NotebookDirectory[],model}];
genericLoc=modLoc;
InitializeModel[modLoc,GenericModel->modLoc];


(* relevant processes for Boltzmann collision integrals with given particle in position p1 *)
particleTypes=DeleteCases[F$Particles,U[1]|-U[1]];
particleNames=Map[ToString[TheLabel[#]]<>If[Head[#]===Times,"bar",""]&,particleTypes]
makeProcesses[firstParticle_]:=DeleteDuplicates[Flatten[Table[{firstParticle,a}->{b,c},{a,particleTypes},{b,particleTypes},{c,particleTypes}],2]];


(* getting expressions for parameters *)
Get[modLoc<>".pars"]
parameters=Map[#[[1]]&,M$ExtParams]
parameters//Length


(* indices for particles *)
mapParticleToInteger:=Thread[particleTypes->(Range[Length[particleTypes]]-1)];
mapParticleToInteger
makeMName[process_]:=M[process[[1,1]],process[[1,2]],process[[2,1]],process[[2,2]]]/.mapParticleToInteger


(* ::Section:: *)
(*Automating steps*)


dimension=4;
tops=CreateTopologies[0,2->2];
momenta={Momentum[p1,dimension],Momentum[p2,dimension],Momentum[p3,dimension],Momentum[p4,dimension]};
indices={{1,1},{1,2},{2,1},{2,2}};


ClearAll[makeAmplitude]
makeAmplitude[process_,channels_:All]:=Module[{diags,ampFA,ampFC,ampSq, ampMsq,ampMsqSummed,masses,i,transverse},
(* assumes 2->2 topology *)
(* create diagrams for process *)
diags=InsertFields[tops[[channels]],process,InsertionLevel->{Particles},Model->modLoc,GenericModel->modLoc];
(* create FeynArts amplitude from diagrams *)
ampFA[1]=CreateFeynAmp[diags,PreFactor->1,Truncated->False,GaugeRules->{GaugeXi[S[_]]->1,GaugeXi[V[_]]->1,FAGaugeXi[S[_]]->1,FAGaugeXi[V[_]]->1}];
ampFA[2]=ampFA[1]/.PropagatorDenominator[x_,y_]:>FeynAmpDenominator[PropagatorDenominator[x,y]];
(* transverse momenta *)
transverse={};
For[i=1,i<=4,i++,
If[process[[indices[[i]]/.List->Sequence]]==V[1],
AppendTo[transverse,momenta[[i]]]
]]; 
(* converting to FeynCalc form *)
ampFC[1]=FCFAConvert[ampFA[2],
IncomingMomenta->{Momentum[p1,dimension],Momentum[p2,dimension]},
OutgoingMomenta->{Momentum[p3,dimension],Momentum[p4,dimension]},
LorentzIndexNames->{\[Mu],\[Nu],\[Rho],\[Sigma]},
ChangeDimension->dimension,
List->False,
SMP->False,
Contract->False,
DropSumOver->False,
UndoChiralSplittings->True,
TransversePolarizationVectors->transverse
];
ampFC[2]=ampFC[1](*/.PropagatorDenominator[x_,y_]\[RuleDelayed]FeynAmpDenominator[PropagatorDenominator[x,y]]*);
ampFC[3]=ampFC[2]/.Conjugate[PolarizationVector][_,a_,b_]->Conjugate[PolarizationVector[a,b]]/.PolarizationVector[_,a_,b_]->PolarizationVector[a,b];
ampFC[4]=ampFC[3]/.Spinor[a__]Spinor[b__]->Spinor[a] . Spinor[b];
Return[ampFC[4]//Contract]
]


(* ::Subsubsection:: *)
(*Amplitude squared*)


ClearAll[makeAmplitudeSquared]
makeAmplitudeSquared[process_,channels_:All,simplify_:True]:=Module[{ampM,ampSq, ampMsq,masses,stu,i},
(* assumes 2->2 topology *)
masses=Join[Map[TheMass,process[[1]]],Map[TheMass,process[[2]]]];
stu=Plus@@(masses^2);
ampM[1]=makeAmplitude[process,channels];
(* making amplitude squared *)
ampSq[1]=ComplexConjugate[ampM[1]]ampM[1];
ampSq[2]=FeynAmpDenominatorExplicit[ampSq[1]];
FCClearScalarProducts[];
(* insert Mandelstam variables *)
SetMandelstam[s,t,u,p1,p2,-p3,-p4,(masses/.List->Sequence)];
ampSq[2]=FermionSpinSum[ampSq[2]]/.DiracTrace[aaa__]->DiracTrace[aaa,Mandelstam->{s,t,u,stu}];
ampSq[3]=DiracSimplify[ampSq[2]];
For[i=1,i<=4,i++,
(* setting on-shell momenta *)
If[masses[[i]]==0,ScalarProduct[momenta[[i]],momenta[[i]]]:=0];
(* summing over polarizations *)
ampSq[i+3]=ampSq[i+2];
If[process[[indices[[i]]/.List->Sequence]]==V[1],
ampSq[i+3]=ampSq[i+2]//DoPolarizationSums[#,momenta[[i]],0]&];
];
ampMsq[1]=PropagatorDenominatorExplicit[ampSq[7]];
(* simplify *)
ampMsq[2]=Expand[TrickMandelstam[ampSq[1],{s,t,u,Plus@@(masses)}]];
FCClearScalarProducts[];
Return[ampMsq[1]]
]


(* squared summed matrix elements with a scalar on leg 1 *)
results=Flatten[Map[
DeleteCases[
Table[makeMName[process]->makeAmplitudeSquared[process,All,False],{process,Evaluate[makeProcesses[#]]}]
,_->0]&,
particleTypes]]/.Map[TheMass[#]->0&,particleTypes]/. {msq[indices__]:>msq[Sort[List[indices]]/. List->Sequence],g[indices__]:>g[Sort[List[indices]]/. List->Sequence],lam[indices__]:>lam[Sort[List[indices]]/. List->Sequence]};


results=results//Simplify


(* ::Section:: *)
(*Exporting to JSON*)


(* writing as JSON *)
makeJsonObject[fields_,parameters_,results_]:=Module[{fieldsJson,matrixElementsJson,toString,getRelevantParameters,replaceSpecials},(*construct an object which can then be exported to JSON*)toString[arg_]:=If[StringQ[arg],arg,ToString[arg,InputForm]];
replaceSpecials[arg_]:=StringReplace[arg,{"Pi"->"_pi","s"->"_s","t"->"_t","u"->"_u"}];
getRelevantParameters[arg_]:=Select[parameters,Not[FreeQ[arg,#]]&];
fieldsJson=Table[{"index"->i-1,"name"->toString[fields[[i]]]},{i,1,Length[fields]}];
matrixElementsJson=Map[{"externalParticles"->#[[1]]/.M[a__]->List[a],"parameters"->Map[toString,getRelevantParameters[#[[2]]]],"expression"->replaceSpecials[toString[#[[2]]]]}&,results];
Return[{"particles"->fieldsJson,"matrixElements"->matrixElementsJson}]]
(*here for the Yukawa model,where the matrix elements are stored in the object called MatrixElements*)
toExportAsJSON=makeJsonObject[particleNames,parameters,results//SMPToSymbol];
Export[FileNameJoin[{NotebookDirectory[],StringJoin[model,".json"]}],toExportAsJSON];
