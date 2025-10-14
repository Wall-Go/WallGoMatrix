(* ::Package:: *)

(*Quit[];*)


If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
(*Put this if you want to create multiple model-files with the same kernel*)
WallGo`WallGoMatrix`$GroupMathMultipleModels=True;
Check[
    Get["WallGo`WallGoMatrix`"],
    (*Get["../Kernel/WallGoMatrix.m"],*)
    Message[Get::noopen, "WallGo`WallGoMatrix` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*Full Standard Model*)


(* ::Section::Closed:: *)
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
		Format->{"json","txt"}
	}
];


(*
	output of matrix elements 
*)
OutputFile="output/matrixElementsFull.LL";

MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{
		TruncateAtLeadingLog->True,
		Format->{"json","txt"}
	}
];


MatrixElements
