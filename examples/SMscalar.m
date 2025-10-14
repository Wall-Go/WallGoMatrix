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
(*SM quarks + gauge bosons*)


(* ::Section:: *)
(*Model*)


Group={"SU3","SU2"};
RepAdjoint={{1,1},{2}};
CouplingName={gs,gw};


Rep1={{{1,0},{1}},"L"};
Rep2={{{1,0},{0}},"R"};
Rep3={{{1,0},{0}},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3};


HiggsDoublet={{{0,0},{1}},"C"};
RepScalar={HiggsDoublet};


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


InputInv={{1,1},{True,False}};
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
VMass=m2*MassTerm1;
\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;
QuarticTerm1=MassTerm1^2;
VQuartic=lam1H*QuarticTerm1;
\[Lambda]4=GradQuartic[VQuartic];

InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;
Ysff=-GradYukawa[yt1*YukawaDoublet[[1]]];
YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt1>0}]];


ImportModel[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Section:: *)
(*MatrixElements*)


(*
In DRalgo fermions are Weyl.
So to create one Dirac we need
one left-handed and
one right-handed fermion
*)


(* ::Subsection:: *)
(*SymmetryBreaking*)


vev={0,v,0,0};
SymmetryBreaking[vev]


(* ::Subsection:: *)
(*UserInput*)


(*
	Reps 1-4 are quarks,
	reps 5,6 are vector bosons
*)
(*left-handed top-quark*)
ReptL=CreateParticle[{{1,1}},"F",mq2,"TopL"];

(*right-handed top-quark*)
ReptR=CreateParticle[{2},"F",mq2,"TopR"];

(*right-handed bottom-quark*)
RepbR=CreateParticle[{3},"F",mq2,"BotR"];

(*Vector bosons*)
RepGluon=CreateParticle[{1},"V",mg2,"Gluon"];
RepW=CreateParticle[{{2,1}},"V",mW2,"W"];

(*Scalar particles*)
RepHiggs=CreateParticle[{1},"S",ms2,"Higgs"];

(*Light particles*)
LightFermions=CreateParticle[{{1,2},4,5,6,7,8,9},"F",mq2,"LightFermions"];


(*
These particles do not have out-of-eq contributions
*)
ParticleList={ReptL,ReptR,RepbR,RepGluon,RepW,RepHiggs};
(*
Light particles are never incoming particles 
*)
LightParticleList={LightFermions};


(*
	output of matrix elements
*)
OutputFile="output/matrixElements.scalar";
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	Format->{"json","txt","hdf5"}
];


MatrixElements//Expand


Import[OutputFile<>".hdf5"]


Import[OutputFile<>".hdf5","CouplingInfo"]


Import[OutputFile<>".hdf5","ParticleInfo"]


Import[OutputFile<>".hdf5","CouplingInfo"]


Import[OutputFile<>".hdf5","ParticleInfo"]
