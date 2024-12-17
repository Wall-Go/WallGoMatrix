(* ::Package:: *)

(*Quit[];*)


If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
(*Put this if you want to create multiple model-files with the same kernel*)
WallGo`WallGoMatrix`$GroupMathMultipleModels=True;
WallGo`WallGoMatrix`$LoadGroupMath=True;
Check[
    Get["../Kernel/WallGoMatrix.m"],
    Message[Get::noopen, "WallGo`WallGoMatrix` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*Full Standard Model*)


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


ImportModel[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Section:: *)
(*SM quarks + gauge bosons*)


(* ::Subsection:: *)
(*SymmetryBreaking*)


vev={0,v,0,0};
SymmetryBreaking[vev]


(* ::Subsection:: *)
(*UserInput*)


(*Third generation of fermions*)
ReptL=CreateParticle[{{1,1}},"F",mq2,"TopL"];
RepbL=CreateParticle[{{1,2}},"F",mq2,"BotL"];
ReptR=CreateParticle[{{2,1}},"F",mq2,"TopR"];
RepbR=CreateParticle[{3},"F",mq2,"BotR"];
RepTauL=CreateParticle[{4},"F",ml2,"TauL"];
RepTauR=CreateParticle[{5},"F",ml2,"TauR"];


(*Second generation of fermions*)
RepCharmStrangeL=CreateParticle[{6},"F",mq2,"CharmStrangeL"];
RepCharmR=CreateParticle[{7},"F",mq2,"CharmR"];
RepStrangeR=CreateParticle[{8},"F",mq2,"StrangeR"];
RepMuonL=CreateParticle[{9},"F",ml2,"MuonL"];
RepMuonR=CreateParticle[{10},"F",ml2,"MuonR"];


(*First generation of fermions*)
RepUpDownL=CreateParticle[{11},"F",mq2,"UpDownL"];
ReUpR=CreateParticle[{12},"F",mq2,"UpR"];
RepDownR=CreateParticle[{13},"F",mq2,"DownR"];
RepElectronL=CreateParticle[{14},"F",ml2,"ElectronL"];
RepElectronR=CreateParticle[{15},"F",ml2,"ElectronR"];


(*Vector bosons*)
RepGluon=CreateParticle[{1},"V",mg2,"Gluon"]; (*Gluons*)
RepW=CreateParticle[{{2,1}},"V",mW2,"W"]; (*SU2 gauge bosons*)
RepB=CreateParticle[{3},"V",mB2,"B"]; (*U1 gauge boson*)


(*Scalars bosons*)
RepHiggs=CreateParticle[{{1,2}},"S",mH2,"Higgs"];
RepGoldstone=CreateParticle[{{1,1}},"S",mG2,"Goldstone"];


ParticleList={
	ReptL,RepbL,ReptR,RepbR,RepTauL,RepTauR,
	RepCharmStrangeL,RepCharmR,RepStrangeR,RepMuonL,RepMuonR,
	RepUpDownL,ReUpR,RepDownR,RepElectronL,RepElectronR,
	RepGluon,RepW,RepB,RepHiggs,RepGoldstone
	};
(*
Light particles are never incoming particles 
*)
LightParticleList={};


(*
	output of matrix elements
*)
OutputFile="output/matrixElementsFull";

MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{
		TruncateAtLeadingLog->False,
		Format->{"json","txt"}}];


(*
	output of matrix elements 
*)
OutputFile="sm.test";

MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{
		TruncateAtLeadingLog->False,
		Format->{"json","txt"}
	}];


MatrixElements


(* ::Chapter:: *)
(*Tests*)


(* ::Section:: *)
(*Importing results from WallGo*)


{particles,parameters,MatrixElements}=ImportMatrixElements["sm.test.json"];


(* ::Section:: *)
(*Importing results from FeynCalc*)


{particlesFeyn,parametersFeyn,MatrixElementsFeyn}=ImportMatrixElements["sm.feyncalc.test.json"];


(* ::Section:: *)
(*Comparison tests*)


insertCouplings={lam1H->lam,gY->0 (* TEMPORARILY TURNING OFF U(1)*),gw->gW,gs->gs};


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})


UserMasses={mq2,ml2,mg2,mW2,mB2,mH2,mG2};
fixConvention[arg_]:=symmetriseTU[arg/.Thread[UserMasses->0]/.{s->(-t-u)}/.insertCouplings/.v->0]//Expand//Simplify//Expand


removeMissing[arg_]:=arg/.M[__]->0/.Missing["KeyAbsent", _]->0


particlesIndicesSubset={0,1,2,3,16,17,18,19};
particlesIndicesFeyn=Map[Last[Last[#]]&,particlesFeyn];


(* ::Subsection:: *)
(*Test hard*)


testList={};


(* {Higgs} -> {Higgs} *)
AppendTo[testList,
TestCreate[
	M[19,19,19,19]/.MatrixElements//fixConvention//removeMissing,
	M[0,0,0,0]/.MatrixElementsFeyn//fixConvention//removeMissing
]];


(* {Higgs} -> {W,Z,B,g} *)
AppendTo[testList,
TestCreate[
	Sum[M[19,19,c,d],{c,{16,17,18}},{d,{16,17,18}}]/.MatrixElements//fixConvention//removeMissing,
	Sum[M[0,0,c,d],{c,{5,6,7,8}},{d,{5,6,7,8}}]/.MatrixElementsFeyn//fixConvention//removeMissing
]];


(* {Higgs} -> {t, b} *)
AppendTo[testList,
TestCreate[
	Sum[M[19,19,c,d],{c,{0,1,2,3}},{d,{0,1,2,3}}]/.MatrixElements//fixConvention//removeMissing,
	Sum[M[0,0,c,d],{c,{1,2,3,4}},{d,{1,2,3,4}}]/.MatrixElementsFeyn//fixConvention//removeMissing
]];


(* I don't expect there should be any 1/Subscript[N, dof] here because this is the physical Higgs on leg 1 *)
Sum[M[19,19,c,d],{c,{0,1,2,3}},{d,{0,1,2,3}}]/.MatrixElements//fixConvention//removeMissing
Sum[M[0,0,c,d],{c,{1,2,3,4}},{d,{1,2,3,4}}]/.MatrixElementsFeyn//fixConvention//removeMissing


(* {t, b} -> {t, b} *)
AppendTo[testList,
TestCreate[
	Sum[M[a,b,c,d],{a,{0,1,2,3}},{b,{0,1,2,3}},{c,{0,1,2,3}},{d,{0,1,2,3}}]/.MatrixElements//fixConvention//removeMissing,
	Sum[1/(2 * 3) M[a,b,c,d],{a,{1,2,3,4}},{b,{1,2,3,4}},{c,{1,2,3,4}},{d,{1,2,3,4}}]/.MatrixElementsFeyn//fixConvention//removeMissing
]];


A["WallGo"]=Sum[M[a,b,c,d],{a,{0,1,2,3}},{b,{0,1,2,3}},{c,{0,1,2,3}},{d,{0,1,2,3}}]/.MatrixElements//fixConvention//removeMissing
A["FeynCalc"]=Sum[1/(2 * 3) M[a,b,c,d],{a,{1,2,3,4}},{b,{1,2,3,4}},{c,{1,2,3,4}},{d,{1,2,3,4}}]/.MatrixElementsFeyn//fixConvention//removeMissing
(* doesn't cancel exactly, and the difference involves the Yukawa coupling *)
A["WallGo"]-A["FeynCalc"]


report=TestReport[testList]
report["ResultsDataset"]



