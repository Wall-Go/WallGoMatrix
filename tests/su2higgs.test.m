(* ::Package:: *)

Quit[];


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
(*SU(2) Higgs*)


(* ::Section:: *)
(*Model*)


Group={"SU2","U1"}; (* I've added the U1 field because without it the Higgs is only given 2 degrees of freedom *)
CouplingName={g,gu1};
RepAdjoint={{2},0};
Higgs={{{1},0},"C"}; (* fundamental *)
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


(* ::Subsection::Closed:: *)
(*SymmetryBreaking*)


(*vev={v,0};
SymmetryBreaking[vev]*)


(* ::Subsection:: *)
(*UserInput*)


(* vector *)
RepVector=CreateParticle[{1},"V",mv,"Vector"];
RepB=CreateParticle[{2},"V",ma,"ExtraVector"];

(* fermion *)
RepFermion={};

(* scalar *)
RepScalar=CreateParticle[{1},"S",ms,"Phi"];



(*
These particles do not necessarily have to be out of equilibrium
the remainin particle content is set as light
*)
ParticleList={RepScalar,RepVector,RepB};
LightParticleList={}


(*Defining various masses and couplings*)
UserMasses={ms,mv,ma};


(*
	output of matrix elements
*)
SetDirectory[NotebookDirectory[]];
OutputFile="output/matrixElements.su2higgs";
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{TruncateAtLeadingLog->False,Format->{"json","txt"},NormalizeWithDOF->False}
]


(* ::Section:: *)
(*Tests*)


(* ::Section:: *)
(*Importing results from WallGo*)


{particles,parameters,MatrixElements}=ImportMatrixElements["output/matrixElements.su2higgs.json"];


(* ::Section:: *)
(*Importing results from FeynCalc*)


{particlesFeyn,parametersFeyn,MatrixElementsFeyn}=ImportMatrixElements["sunhiggs.test.json"];


(* ::Section:: *)
(*Comparison tests*)


insertCouplings={Global`g->g,\[Lambda]->lam,SUNN->2,gu1->0};


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
AppendTo[testList,TestCreate[test["WallGo"][1],test["FeynCalc"][1],TestID->"WallGo vs FeynCalc: SS->SS"]];


(* understanding the pure-scalar part of the result *)
Mabcd=(2\[Lambda])(KroneckerDelta[a,b]KroneckerDelta[c,d]+KroneckerDelta[a,c]KroneckerDelta[b,d]+KroneckerDelta[a,d]KroneckerDelta[b,c])
(* should give (6\[Lambda])^2 for a single d.o.f. *)
Sum[Mabcd Mabcd,{a,1},{b,1},{c,1},{d,1}]
Sum[Mabcd Mabcd,{a,2},{b,2},{c,2},{d,2}] (* WallGo looks like we have 2 real scalars *)
Sum[Mabcd Mabcd,{a,4},{b,4},{c,4},{d,4}] (* FeynCalc looks like we have 4 real scalars *)


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
AppendTo[testList,TestCreate[test["WallGo"][2],test["FeynCalc"][2],TestID->"WallGo vs FeynCalc: SS->VV"]];


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
AppendTo[testList,TestCreate[test["WallGo"][3],test["FeynCalc"][3],TestID->"WallGo vs FeynCalc: VV->SS"]];


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
AppendTo[testList,TestCreate[test["WallGo"][4],test["FeynCalc"][4],TestID->"WallGo vs FeynCalc: SV->SV"]];


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
AppendTo[testList,TestCreate[test["WallGo"][5],test["FeynCalc"][5],TestID->"WallGo vs FeynCalc: VS->VS"]];


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
AppendTo[testList,TestCreate[test["WallGo"][6],test["FeynCalc"][6],TestID->"WallGo vs FeynCalc: VV->VV"]];


(* vector-vector scattering versus AMY *)
test["WallGo"][6]=testWallGo[
	{"Vector"},
	{"Vector"},
	{"Vector"},
	{"Vector"}
]
test["AMY"][6]=fixConvention[16 dA CA^2 g^4 (3-(s u)/t^2-(s t)/u^2-(t u)/s^2)/.{dA->3,CA->2}]
AppendTo[testList,TestCreate[test["WallGo"][6],test["AMY"][6],TestID->"WallGo vs AMY: VV->VV"]];


report=TestReport[testList]
report["ResultsDataset"]





