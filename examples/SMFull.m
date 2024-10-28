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
ReptL=CreateParticle[{{1,1}},"F","TopL"];
RepbL=CreateParticle[{{1,2}},"F","BotL"];
ReptR=CreateParticle[{{2,1}},"F","TopR"];
RepbR=CreateParticle[{3},"F","BotR"];
RepTauL=CreateParticle[{4},"F","TauL"];
RepTauR=CreateParticle[{5},"F","TauR"];


(*Second generation of fermions*)
RepCharmStrangeL=CreateParticle[{6},"F","CharmStrangeL"];
RepCharmR=CreateParticle[{7},"F","CharmR"];
RepStrangeR=CreateParticle[{8},"F","StrangeR"];
RepMuonL=CreateParticle[{9},"F","MuonL"];
RepMuonR=CreateParticle[{10},"F","MuonR"];


(*First generation of fermions*)
RepUpDownL=CreateParticle[{11},"F","UpDownL"];
ReUpR=CreateParticle[{12},"F","UpR"];
RepDownR=CreateParticle[{13},"F","DownR"];
RepElectronL=CreateParticle[{14},"F","ElectronL"];
RepElectronR=CreateParticle[{15},"F","ElectronR"];


(*Vector bosons*)
RepGluon=CreateParticle[{1},"V","Gluon"]; (*Gluons*)
RepW=CreateParticle[{{2,1}},"V","W"]; (*SU2 gauge bosons*)
RepB=CreateParticle[{3},"V","B"]; (*U1 gauge boson*)


(*Scalars bosons*)
RepHiggs=CreateParticle[{{1,2}},"S","Higgs"];
RepGoldstone=CreateParticle[{{1,1}},"S","Goldstone"];


ParticleList={
	ReptL,RepbL,ReptR,RepbR,RepTauL,RepTauR,
	RepCharmStrangeL,RepCharmR,RepStrangeR,RepMuonL,RepMuonR,
	RepUpDownL,ReUpR,RepDownR,RepElectronL,RepElectronR,
	RepGluon,RepW,RepB,RepHiggs,RepGoldstone
	};


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
	output of matrix elements
*)
OutputFile="output/matrixElementsFull";

MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	ParticleMasses,
	{
		TruncateAtLeadingLog->False,
		Format->{"json","txt"}}];


(*
	output of matrix elements 
*)
OutputFile="output/matrixElementsFull.LL";

MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	ParticleMasses,
	{
		TruncateAtLeadingLog->True,
		Format->{"json","txt"}}];


MatrixElements
