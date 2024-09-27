(* ::Package:: *)

(* ::Title:: *)
(*Matrix elements  in a Yukawa Model*)


(* ::Chapter:: *)
(*Loading packages*)


(* ::Input:: *)
(*$LoadAddOns={"FeynArts"};*)
(*Get["FeynCalc`"];*)
(*$FAVerbose=0;*)
(*SetOptions[FourVector,FeynCalcInternal->False];*)


(* ::Chapter:: *)
(*Choosing model*)


(* ::Input:: *)
(*(* this is where the FeynArts files made by FeynRules can be found, ending in .gen and .mod *)*)
(*model="yukawa";*)
(*modLoc=FileNameJoin[{NotebookDirectory[],model}];*)
(*genericLoc=modLoc;*)
(*InitializeModel[modLoc,GenericModel->modLoc];*)


(* ::Input:: *)
(*(* relevant processes for Boltzmann collision integrals with given particle in position p1 *)*)
(*particleTypes=DeleteCases[F$Particles,U[1]|-U[1]];*)
(*particleNames=Map[ToString[TheLabel[#]]<>If[Head[#]===Times,"bar",""]&,particleTypes]*)
(*makeProcesses[firstParticle_]:=DeleteDuplicates[Flatten[Table[{firstParticle,a}->{b,c},{a,particleTypes},{b,particleTypes},{c,particleTypes}],2]];*)


(* ::Input:: *)
(*(* getting expressions for parameters *)*)
(*Get[modLoc<>".pars"]*)
(*parameters=Map[#[[1]]&,M$ExtParams]*)


(* ::Input:: *)
(*(* indices for particles *)*)
(*mapParticleToInteger:=Thread[particleTypes->(Range[Length[particleTypes]]-1)];*)
(*mapParticleToInteger*)
(*makeMName[process_]:=M[process[[1,1]],process[[1,2]],process[[2,1]],process[[2,2]]]/.mapParticleToInteger*)


(* ::Chapter:: *)
(*Computing matrix elements*)


(* ::Input:: *)
(*dimension=4;*)
(*tops=CreateTopologies[0,2->2];*)
(*momenta={Momentum[p1,dimension],Momentum[p2,dimension],Momentum[p3,dimension],Momentum[p4,dimension]};*)
(*indices={{1,1},{1,2},{2,1},{2,2}};*)


(* ::Input:: *)
(*ClearAll[makeAmplitude]*)
(*makeAmplitude[process_,channels_:All]:=Module[{diags,ampFA,ampFC,ampSq, ampMsq,ampMsqSummed,masses,i},*)
(*(* assumes 2->2 topology *)*)
(*(* create diagrams for process *)*)
(*diags=InsertFields[tops[[channels]],process,InsertionLevel->{Particles},Model->modLoc];*)
(*(* create FeynArts amplitude from diagrams *)*)
(*ampFA[1]=CreateFeynAmp[diags,PreFactor->1,Truncated->False]; *)
(*(* converting to FeynCalc form *)*)
(*ampFC[1]=FCFAConvert[ampFA[1],*)
(*IncomingMomenta->{Momentum[p1,dimension],Momentum[p2,dimension]},*)
(*OutgoingMomenta->{Momentum[p3,dimension],Momentum[p4,dimension]},*)
(*LorentzIndexNames->{\[Mu],\[Nu],\[Rho],\[Sigma]},*)
(*ChangeDimension->dimension,*)
(*List->False,*)
(*SMP->True,*)
(*Contract->True,*)
(*DropSumOver->True,*)
(*UndoChiralSplittings->True*)
(*];*)
(*Return[ampFC[1]//Contract]*)
(*]*)


(* ::Input:: *)
(*ClearAll[makeAmplitudeSquared]*)
(*makeAmplitudeSquared[process_,channels_:All,simplify_:True]:=Module[{ampM,ampSq, ampMsq,masses,stu,i},*)
(*(* assumes 2->2 topology *)*)
(*masses=Join[Map[TheMass,process[[1]]],Map[TheMass,process[[2]]]];*)
(*stu=Plus@@(masses^2);*)
(*ampM[1]=makeAmplitude[process,channels];*)
(*(* making amplitude squared *)*)
(*ampSq[1]=ComplexConjugate[ampM[1]]ampM[1];*)
(*ampSq[2]=FeynAmpDenominatorExplicit[ampSq[1]];*)
(*FCClearScalarProducts[];*)
(*(* insert Mandelstam variables *)*)
(*SetMandelstam[s,t,u,p1,p2,-p3,-p4,(masses/.List->Sequence)];*)
(*ampSq[2]=FermionSpinSum[ampSq[2]]/.DiracTrace[aaa__]->DiracTrace[aaa,Mandelstam->{s,t,u,stu}];*)
(*ampSq[3]=DiracSimplify[ampSq[2]];*)
(*(* simplify *)*)
(*ampMsq[1]=Expand[TrickMandelstam[ampSq[3],{s,t,u,Plus@@(masses)}]];*)
(*FCClearScalarProducts[];*)
(*Return[ampMsq[1]]*)
(*]*)


(* ::Input:: *)
(*(* squared summed matrix elements with a scalar on leg 1 *)*)
(*scalar=DeleteCases[*)
(*Table[makeMName[process]->makeAmplitudeSquared[process,All,False],{process,Evaluate[makeProcesses[S[1]]]}]*)
(*,_->0];*)
(*(* setting masses to zero *)*)
(*scalar0=scalar/.Thread[Map[TheMass,particleTypes]->0]*)


(* ::Input:: *)
(*(* squared summed matrix elements with a top or anti-top on leg 1 *)*)
(*fermion=DeleteCases[*)
(*Join[*)
(*Table[makeMName[process]->makeAmplitudeSquared[process,All,False],{process,Evaluate[makeProcesses[F[1]]]}],*)
(*Table[makeMName[process]->makeAmplitudeSquared[process,All,False],{process,Evaluate[makeProcesses[-F[1]]]}]*)
(*]*)
(*,_->0];*)
(*(* setting masses to zero *)*)
(*fermion0=fermion/.Thread[Map[TheMass,particleTypes]->0]*)


(* ::Input:: *)
(*results=Flatten[{scalar0,fermion0},1];*)


(* ::Input:: *)
(*results//MatrixForm*)


(* ::Input:: *)
(*feynAssociation=Association[results];*)


(* ::Chapter:: *)
(*Exporting to JSON*)


(* ::Input:: *)
(*(* writing as JSON *)*)
(*makeJsonObject[fields_,parameters_,results_]:=Module[{fieldsJson,matrixElementsJson,toString,getRelevantParameters,replaceSpecials},(*construct an object which can then be exported to JSON*)toString[arg_]:=If[StringQ[arg],arg,ToString[arg,InputForm]];*)
(*replaceSpecials[arg_]:=StringReplace[arg,{"Pi"->"_pi","s"->"_s","t"->"_t","u"->"_u"}];*)
(*getRelevantParameters[arg_]:=Select[parameters,Not[FreeQ[arg,#]]&];*)
(*fieldsJson=Table[{"index"->i-1,"name"->toString[fields[[i]]]},{i,1,Length[fields]}];*)
(*matrixElementsJson=Map[{"externalParticles"->#[[1]]/.M[a__]->List[a],"parameters"->Map[toString,getRelevantParameters[#[[2]]]],"expression"->replaceSpecials[toString[#[[2]]]]}&,results];*)
(*Return[{"particles"->fieldsJson,"matrixElements"->matrixElementsJson}]]*)
(*(*here for the Yukawa model,where the matrix elements are stored in the object called MatrixElements*)*)
(*toExportAsJSON=makeJsonObject[particleNames,parameters,results//SMPToSymbol];*)
(*Export[FileNameJoin[{NotebookDirectory[],StringJoin[model,".json"]}],toExportAsJSON];*)
