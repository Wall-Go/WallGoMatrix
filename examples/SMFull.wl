(* ::Package:: *)

Quit[];


SetDirectory[NotebookDirectory[]];
(*Put this if you want to create multiple model-files with the same kernel*)
$GroupMathMultipleModels=True;
$LoadGroupMath=True;
<<../DRalgo/DRalgo.m
<<../src/matrixElements.m


(* ::Chapter:: *)
(*QCD+W boson*)


(* ::Section:: *)
(*Model*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
HiggsDoublet={{{0,0},{1},1/2},"C"};
RepScalar={HiggsDoublet};
CouplingName={gs,gw,gY};


Rep1={{{1,0},{1},1/6},"L"};
Rep2={{{1,0},{0},2/3},"R"};
Rep3={{{1,0},{0},-1/3},"R"};
Rep4={{{0,0},{1},-1/2},"L"};
Rep5={{{0,0},{0},-1},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


(* ::Text:: *)
(*The first element is the vector self - interaction matrix :*)


InputInv={{1,1},{True,False}};
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;


VMass=m2*MassTerm1;


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


QuarticTerm1=MassTerm1^2;


VQuartic=lam1H*QuarticTerm1;


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


Ysff=-GradYukawa[yt*YukawaDoublet[[1]]];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];


(* ::Section:: *)
(*SM quarks + gauge bosons*)


(* ::Subsection:: *)
(*SymmetryBreaking*)


vev={0,v,0,0};
SymmetryBreaking[vev]


(* ::Subsection:: *)
(*UserInput*)


(*Third generation of fermions*)
ReptL=CreateParticle[{{1,1}},"F"];
RepbL=CreateParticle[{{1,2}},"F"];
ReptR=CreateParticle[{{2,1}},"F"];
RepbR=CreateParticle[{3},"F"];
RepTauL=CreateParticle[{4},"F"];
RepTauR=CreateParticle[{5},"F"];



(*Second generation of fermions*)
RepCharmStrangeL=CreateParticle[{6},"F"];
RepCharmR=CreateParticle[{7},"F"];
RepStrangeR=CreateParticle[{8},"F"];
RepMuonL=CreateParticle[{9},"F"];
RepMuonR=CreateParticle[{10},"F"];


(*First generation of fermions*)
RepUpDownL=CreateParticle[{11},"F"];
ReUpR=CreateParticle[{12},"F"];
RepDownR=CreateParticle[{13},"F"];
RepElectronL=CreateParticle[{14},"F"];
RepElectronR=CreateParticle[{15},"F"];


(*Vector bosons*)
RepGluon=CreateParticle[{1},"V"]; (*Gluons*)
RepW=CreateParticle[{{2,1}},"V"]; (*SU2 gauge bosons*)
RepB=CreateParticle[{3},"V"]; (*U1 gauge boson*)


(*Scalars bosons*)
RepHiggs=CreateParticle[{{1,2}},"S"];
RepGoldstone=CreateParticle[{{1,1}},"S"];


ParticleList={ReptL,RepbL,ReptR,RepbR,RepTauL,RepTauR,RepCharmStrangeL,RepCharmR,RepStrangeR,RepMuonL,RepMuonR
			,RepUpDownL,ReUpR,RepDownR,RepElectronL,RepElectronR,RepGluon,RepW,RepB,RepHiggs,RepGoldstone};


(*Defining various masses and couplings*)


VectorMass=Join[
	Table[mg2,{i,1,RepGluon[[1]]//Length}],
	Table[mw2,{i,1,RepW[[1]]//Length}],{mb2}]; (*mb2 is the mass of the U(1) gauge field*)


FermionMass=Table[mq2,{i,1,Length[gvff[[1]]]}];


LeptonIndices=Union[ParticleList[[5]][[1]],ParticleList[[6]][[1]],ParticleList[[10]][[1]],ParticleList[[11]][[1]],ParticleList[[15]][[1]],ParticleList[[16]][[1]]];


FermionMass[[LeptonIndices]]=ConstantArray[ml2,Length[LeptonIndices]];


ScalarMass={mG2,mH2,mG2,mG2};
ParticleMasses={VectorMass,FermionMass,ScalarMass};


(*
up to the user to make sure that the same order is given in the python code
*)
UserMasses={mq2,ml2,mg2,mw2,mb2,mG2,mH2}; 
UserCouplings=Variables@Normal@{Ysff,gvss,gvff,gvvv,\[Lambda]4,\[Lambda]3}//DeleteDuplicates;


ParticleList={ReptL,RepbL,ReptR,RepbR,RepTauL,RepTauR,RepCharmStrangeL,RepCharmR,RepStrangeR,RepMuonL,RepMuonR
			,RepUpDownL,ReUpR,RepDownR,RepElectronL,RepElectronR,RepGluon,RepW,RepB,RepHiggs,RepGoldstone};
ParticleName={"TopL","BotL","TopR","BotR","tauL","tauR","CharmStrangeL","CharmR","StrangeR","MuonL","MuonR",
				"UpDownL","UpR","DownR","ElectronL","ElectronR","Gluon","W","B","Higgs","Goldstone"};


(*
	output of matrix elements
*)
OutputFile="matrixElementsFull";
SetDirectory[NotebookDirectory[]];

RepOptional={};
MatrixElements=ExportMatrixElements[OutputFile,ParticleList,UserMasses,UserCouplings,ParticleName,ParticleMasses,RepOptional,Format->"json"];


MatrixElements


(*Formatting the matrix elements*)
	toExportJson=makeJsonMatrixElements[ParticleName,UserCouplings,MatrixElements]



(*Exporting the result*)
	exportJsonMatrixElements[StringJoin[OutputFile,".json"],toExportJson]	


NotebookDirectory[]


(* ::Subsubsection::Closed:: *)
(*json matrix elements functions*)


(* ::Input:: *)
(*makeJsonMatrixElements::usage="makeJsonMatrixElements[particles,parameters,results] converts a list of particle names {'Phi',...}, a list of particle parameters {g,...}, and a list of matrix elements results in the form {M[0,0,0,0]->g^4 s/t,...} to a JSON object in a standard format.";*)
(*makeJsonMatrixElements[particles_,parameters_,resultsI_]:=Module[{particlesJson,matrixElementsJson,toString,getRelevantParameters,replaceSpecials,results=resultsI},*)
(*toString[arg_]:=If[StringQ[arg],arg,ToString[arg,InputForm]];*)
(*replaceSpecials[arg_]:=StringReplace[arg,{"Pi"->"_pi","sReplace"->"_s","tReplace"->"_t","uReplace"->"_u"}];*)
(*getRelevantParameters[arg_]:=Select[parameters,Not[FreeQ[arg,#]]&];*)
(*particlesJson=Table[<|"index"->i-1,"name"->toString[particles[[i]]]|>,{i,1,Length[particles]}];*)
(*results=results/.{s->sReplace, t->tReplace,u->uReplace};*)
(*matrixElementsJson=Map[<|"externalParticles"->#[[1]]/.M[a__]->List[a],"parameters"->Map[toString,getRelevantParameters[#[[2]]]],"expression"->replaceSpecials[toString[#[[2]]]]|>&,results];*)
(*Return[<|"particles"->particlesJson,"matrixElements"->matrixElementsJson|>]];*)


(* ::Input:: *)
(*testJsonMatrixElements::usage="testJsonMatrixElements[json] tests if a JSON object is of the expected form for exporting matrix elements.";*)
(*testJsonMatrixElements[json_]:=Module[{testBool,returnString,nParticles,expectedForm},*)
(*testBool=True;*)
(*returnString="Json object matches expected schema";*)
(*(* checking head *)*)
(*If[Head[json]!=Association,*)
(*returnString="Not Association";testBool=False];*)
(*(* checking dimensions *)*)
(*If[Dimensions[json]!={2},*)
(*returnString="Dimensions not {2}";testBool=False];*)
(*(* checking top level keys *)*)
(*If[Keys[json]!={"particles","matrixElements"},*)
(*returnString="Top level keys not {'particles','matrixElements'}";testBool=False];*)
(*(* checking lower level keys *)*)
(*If[Keys[json["particles"][[1]]]!={"index","name"},*)
(*returnString="'particles' keys not {'index','name'}";testBool=False];*)
(*If[Keys[json["matrixElements"][[1]]]!={"externalParticles","parameters","expression"},*)
(*returnString="'matrixElements' keys not {'externalParticles','parameters','expressions'}";testBool=False];*)
(*(* returning results *)*)
(*{testBool, returnString}*)
(*]*)


(* ::Input:: *)
(*splitJsonMatrixElements::usage="splitJsonMatrixElements[json] splits a JSON object containing matrix elements into a list {particleNames,parameters,results}.";*)
(*splitJsonMatrixElements[json_]:=Module[{particles,matrixElements,particleIndices,particleNames,matrixElementIndices,matrixElementParameters,matrixElementExpressions,parameters,expressions,results},*)
(*particles=json["particles"];*)
(*particleIndices=Map[#["index"]&,json["particles"]];*)
(*particleNames=Map[#["name"]&,json["particles"]];*)
(*matrixElements=json["matrixElements"];*)
(*matrixElementIndices=Map[#["externalParticles"]&,json["matrixElements"]];*)
(*matrixElementParameters=Map[#["parameters"]&,json["matrixElements"]];*)
(*matrixElementExpressions=Map[#["expression"]&,json["matrixElements"]];*)
(*parameters=Map[ToExpression,DeleteDuplicates[Flatten[matrixElementParameters]]];*)
(*expressions=Map[ToExpression,StringReplace[matrixElementExpressions,{RegularExpression["(\\W)_s"]->"$1s",RegularExpression["(\\W)_t"]->"$1t",RegularExpression["(\\W)_u"]->"$1u"}]];*)
(*results=Thread[matrixElementIndices->expressions];*)
(*results=Map[M[#[[1]]/.List->Sequence]->#[[2]]&,results];*)
(*{particleNames,parameters,results}*)
(*];*)


ExportToJSON[MatrixElement_,ParticleName_,UserCouplings_,file_]:=Block[{toExportJson},
	
(*Formatting the matrix elements*)
	toExportJson=makeJsonMatrixElements[ParticleName,UserCouplings,MatrixElement];
(*Exporting the result*)
	exportJsonMatrixElements[StringJoin[file,".json"],toExportJson];		
]


(* ::Input:: *)
(*(* reading JSON matrix elements *)*)
(*importJSONMatrixElements::usage="importJSONMatrixElements[file] imports a JSON file of matrix elements into a JSON object.";*)
(*importJSONMatrixElements[file_]:=Import[fileImport,"RawJSON"];*)


(* ::Input:: *)
(*(* export JSONMatrixElements *)*)
(*exportJsonMatrixElements::usage="exportJsonMatrixElements[file,jsonMatrixElements] exports a JSON object of matrix elements into a JSON file.";*)
(*exportJsonMatrixElements[file_,jsonMatrixElements_]:=Module[{test},*)
(*If[Not[StringQ[file]],Print["File must be a string"];Return[]];*)
(*If[StringTake[fileExport,-5]!=".json",Print["File must end in .json"];Return[]];*)
(*Export[file,jsonMatrixElements]*)
(*];*)
