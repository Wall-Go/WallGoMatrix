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


(* ::Title:: *)
(*SU(3) Higgs*)


(* ::Section:: *)
(*Model*)


Group={"SU3"};
CouplingName={g};
RepAdjoint={{1,1}};
Higgs={{{1,0}},"C"}; (* fundamental *)
RepScalar={Higgs};


RepFermion={};


(* ::Text:: *)
(*The first element is the vector self-interaction matrix:*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


InputInv={{1,1},{True,False}}; (*This specifies that we want a \[Phi]^+\[Phi] term*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;


VMass=msq*MassTerm1[[1]];(*This is the \[Phi]^+\[Phi] term written in component form*)


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


QuarticTerm1=MassTerm1[[1]]^2; (*Because MassTerm1=\[Phi]^+\[Phi], we can write (\[Phi]^+\[Phi])^2=MassTerm1^2*)


VQuartic=\[Lambda]*QuarticTerm1;


\[Lambda]4=GradQuartic[VQuartic];


ImportModel[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Section:: *)
(*User Input*)


(* ::Subsection:: *)
(*SymmetryBreaking*)


(*vev={v,0};
SymmetryBreaking[vev]*)


(* ::Subsection:: *)
(*UserInput*)


(* vector *)
RepVector=CreateParticle[{1},"V",mv,"Vector"];

(* fermion *)
RepFermion={};

(* scalar *)
RepScalar=CreateParticle[{1},"S",ms,"Phi"];



(*
These particles do not necessarily have to be out of equilibrium
the remainin particle content is set as light
*)
ParticleList={RepScalar,RepVector};
LightParticleList={}


(*Defining various masses and couplings*)
UserMasses={ms,mv};


(*
	output of matrix elements
*)
SetDirectory[NotebookDirectory[]];
OutputFile="output/matrixElements.su3higgs";
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{TruncateAtLeadingLog->False,Format->{"json","txt"},NormalizeWithDOF->False}
]


(* ::Chapter:: *)
(*Tests*)


(* ::Section:: *)
(*Importing results from WallGo*)


{particles,parameters,MatrixElements}=ImportMatrixElements["output/matrixElements.su3higgs.json"];


(* ::Section:: *)
(*Importing results from FeynCalc*)


{particlesFeyn,parametersFeyn,MatrixElementsFeyn}=ImportMatrixElements["testFiles/sunhiggs.test.json"];


(* ::Section:: *)
(*Comparison tests*)


insertCouplings={Global`g->g,\[Lambda]->lam,SUNN->3};


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})


fixConvention[arg_]:=symmetriseTU[arg/.Thread[UserMasses->0]/.{s->(-t-u)}/.insertCouplings]//Expand//Simplify//Expand


removeMissing[arg_]:=arg/.M[__]->0/.Missing["KeyAbsent", _]->0


testsRulesWallGo[arg_]:=arg/.Flatten[particles]/.MatrixElements//fixConvention//removeMissing;
testsRulesFeynCalc[arg_]:=arg/.Flatten[particlesFeyn]/.MatrixElementsFeyn//fixConvention//removeMissing
testWallGo[particlesA_,particlesB_,particlesC_,particlesD_]:=Sum[
	M[a,b,c,d],{a,particlesA},{b,particlesB},{c,particlesC},{d,particlesD}
]//testsRulesWallGo
testFeynCalc[particlesA_,particlesB_,particlesC_,particlesD_]:=Sum[
	M[a,b,c,d],{a,particlesA},{b,particlesB},{c,particlesC},{d,particlesD}
]//testsRulesFeynCalc


(* ::Subsection:: *)
(*Test hard*)


particles
particlesFeyn


testList={};


(* scalar-scalar scattering*)
test["WallGo"][1]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"Phi"},
	{"Phi"}
]
test["FeynCalc"][1]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Phi","Phibar"}
]
AppendTo[testList,TestCreate[test["WallGo"][1],test["FeynCalc"][1],TestID->"SS->SS"]];


(* scalar to vector *)
test["WallGo"][2]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"Vector"},
	{"Vector"}
]
test["FeynCalc"][2]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"A"},
	{"A"}
]
AppendTo[testList,TestCreate[test["WallGo"][2],test["FeynCalc"][2],TestID->"SS->VV"]];


(* vector to scalar *)
test["WallGo"][3]=testWallGo[
	{"Vector"},
	{"Vector"},
	{"Phi"},
	{"Phi"}
]
test["FeynCalc"][3]=testFeynCalc[
	{"A"},
	{"A"},
	{"Phi","Phibar"},
	{"Phi","Phibar"}
]
AppendTo[testList,TestCreate[test["WallGo"][3],test["FeynCalc"][3],TestID->"VV->SS"]];


(* scalar-vector scattering *)
test["WallGo"][4]=testWallGo[
	{"Phi"},
	{"Vector"},
	{"Phi","Vector"},
	{"Phi","Vector"}
]
test["FeynCalc"][4]=testFeynCalc[
	{"Phi","Phibar"},
	{"A"},
	{"Phi","Phibar","A"},
	{"Phi","Phibar","A"}
]
AppendTo[testList,TestCreate[test["WallGo"][4],test["FeynCalc"][4],TestID->"SV->SV"]];


(* vector-scalar scattering *)
test["WallGo"][5]=testWallGo[
	{"Vector"},
	{"Phi"},
	{"Phi","Vector"},
	{"Phi","Vector"}
]
test["FeynCalc"][5]=testFeynCalc[
	{"A"},
	{"Phi","Phibar"},
	{"Phi","Phibar","A"},
	{"Phi","Phibar","A"}
]
AppendTo[testList,TestCreate[test["WallGo"][5],test["FeynCalc"][5],TestID->"VS->VS"]];


(* vector-vector scattering*)
test["WallGo"][6]=testWallGo[
	{"Vector"},
	{"Vector"},
	{"Vector"},
	{"Vector"}
]
test["FeynCalc"][6]=testFeynCalc[
	{"A"},
	{"A"},
	{"A"},
	{"A"}
]
AppendTo[testList,TestCreate[test["WallGo"][6],test["FeynCalc"][6],TestID->"VV->VV"]];


report=TestReport[testList]
report["ResultsDataset"]



