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
(*2HDM*)


(*See 2211.13142 for implementation details -- note our different normalization in the
quartic couplings*)


(* ::Section:: *)
(*Model*)


Group={"SU3","SU2"};
RepAdjoint={{1,1},{2},0};
HiggsDoublet1={{{0,0},{1}},"C"};
HiggsDoublet2={{{0,0},{1}},"C"};
RepScalar={HiggsDoublet1,HiggsDoublet2};
CouplingName={g3,gw};


Rep1={{{1,0},{1}},"L"};
Rep2={{{1,0},{0}},"R"};
Rep3={{{1,0},{0}},"R"};
Rep4={{{0,0},{1}},"L"};
Rep5={{{0,0},{0}},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};


(* ::Text:: *)
(*The input for the gauge interactions to DRalgo are then given by*)


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The first element is the vector self-interaction matrix:*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


InputInv={{1,1},{True,False}};
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;
InputInv={{2,2},{True,False}};
MassTerm2=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;
InputInv={{1,2},{True,False}};
MassTerm3=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;
InputInv={{2,1},{True,False}};
MassTerm4=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;


VMass=(
	+m1*MassTerm1
	+m2*MassTerm2
	);


\[Mu]ij=GradMass[VMass[[1]]]//Simplify//SparseArray;


QuarticTerm1=MassTerm1[[1]]^2;
QuarticTerm2=MassTerm2[[1]]^2;
QuarticTerm3=MassTerm1[[1]]*MassTerm2[[1]];
QuarticTerm4=MassTerm3[[1]]*MassTerm4[[1]];
QuarticTerm5=(MassTerm3[[1]]^2+MassTerm4[[1]]^2)//Simplify;


VQuartic=(
	+lam1H/2*QuarticTerm1
	+lam2H/2*QuarticTerm2
	+lam3H*QuarticTerm3
	+lam4H*QuarticTerm4
	+lam5H/2*QuarticTerm5
	);


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


Ysff=-GradYukawa[yt1*YukawaDoublet1[[1]]];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt1>0}]];


ImportModel[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Section:: *)
(*MatrixElements*)


(* ::Subsection:: *)
(*SymmetryBreaking*)


vev={0,v,0,0,0,0,0,0};
SymmetryBreaking[vev,VevDependentCouplings->True]
(*SymmetryBreaking[vev]*)


(* ::Subsection:: *)
(*Grouping representations*)


(* ::Text:: *)
(*We are interested in the matrix elements for top quarks, weak gauge bosons, and the heavy scalars A and H^+-.*)
(**)
(*We thus give the top quark its own representation, and group all other quarks together.*)
(*The Higgs and Goldstone bosons are not allowed out of equilibrium, and since they have the same *)
(*thermal mass, they are grouped together in one representation.*)
(**)
(*Since we neglect U(1), we treat the W and Z as identical, and group them into one single representation.*)
(**)
(*The second Higgs field is split into a representation for the heavy scalars, A, H^+-, which we treat as identical,*)
(*since we assume lambda4 = lambda5. *)
(*The CP-even scalar H will not be allowed out-of-equilibrium and  is grouped in its own representation.*)


(*
In DRalgo fermions are Weyl.
So to create one Dirac we need
one left-handed and
one right-handed fermoon
*)


(*Quarks*)
ReptL=CreateParticle[{{1,1}},"F", mq2, "TopL"]; (*left-handed top-quark*)
ReptR=CreateParticle[{{2,1}},"F", mq2, "TopR"]; (*right-handed top-quark*)
RepLightQ = CreateParticle[{{1,2},3,6,7,8,11,12,13},"F", mq2, "LightQuark"]; (*light quarks*)


(*Leptons*)
RepLepL = CreateParticle[{4,9,14},"F", ml2, "LepL"]; (*left-handed leptons*)
RepLepR = CreateParticle[{5,10,15},"F", ml2, "LepR"]; (*right-handed leptons -- these don't contribute*)


(*Vector bosons*)
RepGluon=CreateParticle[{1},"V",mg2,"Gluon"]; (*Gluons*)
RepW=CreateParticle[{{2,1}},"V",mW2,"W"]; (*SU2 gauge bosons, approximating the W and the Z as the same particle*)


(*Scalars bosons*)
RepHiggs = CreateParticle[{1},"S", mh2,"Higgs"];(*Higgs and Goldstones*)
RepH=CreateParticle[{{2,2}},"S",mH2,"H"]; (*CP-even inert scalar*)
RepA=CreateParticle[{{2,3},{2,1}},"S",mA2,"A"]; (*CP-odd inert and charged scalars. Note that when lambda4 = lambda5, they have the same mass *)


ParticleList={ReptL,ReptR,RepLightQ,RepGluon,RepW, RepHiggs, RepA};
(*
Light particles are never incoming particles 
*)
LightParticleList={RepLepL, RepLepR, RepH};


(*
	output of matrix elements
*)
OutputFile="output/matrixElements.idm";

MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{
		TruncateAtLeadingLog->True,
		Replacements->{lam1H->0,lam2H->0,lam4H->0,lam5H->0},
		Format->{"json","txt"},
		NormalizeWithDOF->False
	}
];

