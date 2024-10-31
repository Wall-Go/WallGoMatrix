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
(*QCD*)


(* ::Section:: *)
(*Model*)


(* ::Text:: *)
(*Using SU (3) in the representation with Dynkin index {1, 1} and representation R = 8 .*)


Group={"SU3"};
RepAdjoint={{1,1}};
CouplingName={gs};


(* ::Text:: *)
(*No scalars are present in the theory:*)


RepScalar={};


(* ::Text:: *)
(*Fermions are implemented as Weyl spinors.*)
(*Therefore, to create one Dirac fermion, one left-handed and one right-handed fermion is needed:*)


Rep1={{{1,0}},"L"};
Rep2={{{1,0}},"R"};
RepFermion1Gen={Rep1,Rep2};


(* ::Text:: *)
(*For QCD, Nf=6 fermions are needed:*)


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen,RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The input for the gauge interactions to DRalgo are then given by:*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


(* ::Text:: *)
(*To invoke the model, it needs to be imported :*)


ImportModel[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Section:: *)
(*A model with 6 quarks and 1 gluon*)


(* ::Subsection:: *)
(*UserInput*)


(* ::Text:: *)
(*Representation for top quark:*)


RepTop=CreateParticle[{1,2},"F",mq2,"Top"];


(* ::Text:: *)
(*Representation for gluon :*)


RepGluon=CreateParticle[{1},"V",mg2,"Gluon"];


(* ::Text:: *)
(*Representation for 5 light quarks :*)


LightQuarks=CreateParticle[{3,4,5,6,7,8,9,10,11,12},"F",mq2,"LightParticle"];


(* ::Text:: *)
(*Collecting out - of - equilibrium particles :*)


(*
These particles do not necessarily have to be out of equilibrium
*)

ParticleList={RepTop,RepGluon};


(* ::Text:: *)
(*Collecting light particles (these are never incoming) :*)


LightParticleList={LightQuarks};


(* ::Subsubsection:: *)
(*Generating the matrix elements*)


OutputFile="output/matrixElements.qcd";
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{
		TruncateAtLeadingLog->True,
		Format->{"json","txt"}}];
