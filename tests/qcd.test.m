(* ::Package:: *)

Quit[];


(* Check Mathematica version *)
If[$VersionNumber < 13.3,
  Print["The Mathematica testing framework requires Mathematica version ", requiredVersion," or higher. You are using version ", currentVersion, "."];
  Abort[]
];

SetDirectory[NotebookDirectory[]];
(*Put this if you want to create multiple model-files with the same kernel*)
$GroupMathMultipleModels=True;
$LoadGroupMath=True;
<<../DRalgo/DRalgo.m
<<../src/matrixElements.m


(* ::Chapter:: *)
(*QCD*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3"};
RepAdjoint={{1,1}};
RepScalar={};
CouplingName={gs};


Rep1={{{1,0}},"L"};
Rep2={{{1,0}},"R"};
RepFermion1Gen={Rep1,Rep2};


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen,RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


(* ::Section::Closed:: *)
(*User Input*)


(*
In DRalgo fermions are Weyl.
So to create one Dirac we need
one left-handed and
one right-handed fermion
*)


(*Below
rep 1-6 are quarks,
rep 7 is a gluon
*)
Rep1=CreateOutOfEq[{1,2},"F"];
RepGluon=CreateOutOfEq[{1},"V"];


ParticleList={Rep1,RepGluon};


(*Defining various masses and couplings*)


VectorMass=Table[mg2,{i,1,Length[gvff]}];
FermionMass=Table[mq2,{i,1,Length[gvff[[1]]]}];
ScalarMass={};
ParticleMasses={VectorMass,FermionMass,ScalarMass};
(*
up to the user to make sure that the same order is given in the python code
*)
UserMasses={mq2,mg2};
UserCouplings={gs};


(*
	output of matrix elements
*)
OutputFile="matrixElements.qcd";
SetDirectory[NotebookDirectory[]];
ParticleName={"Top","Gluon"};
RepOptional={};
MatrixElements=ExportMatrixElements[OutputFile,ParticleList,UserMasses,UserCouplings,ParticleName,ParticleMasses,RepOptional];


MatrixElements


(* ::Section:: *)
(*Tests*)


(* ::Subsection:: *)
(*Test hard*)


testList={};


(*5 light quarks*)
AppendTo[testList,
TestCreate[1/2*(M[0,2,0,2]+M[0,2,2,0])/.MatrixElements/.{c[0]->1}//Simplify,
	1/6 ((80 (s^2+u^2))/(t-msq[1])^2+(80 (s^2+t^2))/(u-msq[1])^2)//Simplify
]];


(*q1q1->gg*)
AppendTo[testList,
TestCreate[1/2*(M[0,0,1,1])/.MatrixElements/.{c[0]->1}//Simplify,
	1/12 (+(128/3) ((t u)/(-t+msq[0])^2+(t u)/(-u+msq[0])^2))//Simplify
]];


(*q1 g->q1 g*)
AppendTo[testList,
TestCreate[1/2*(M[0,1,0,1]+M[0,1,1,0])/.MatrixElements/.{c[0]->1}//Simplify,
	(-(64/9) ((s u)/(u-msq[0])^2)+(16 (s^2+u^2))/(t-msq[1])^2)//Simplify
]];


(*tt->tt*)
AppendTo[testList,
TestCreate[1/2*(M[0,0,0,0])/.MatrixElements/.{c[0]->1}//Simplify,
	8/3*((s^2+u^2)/(t-msq[1])^2+(s^2+t^2)/(u-msq[1])^2)//Simplify
]];


(*g g->g g*)
AppendTo[testList,
TestCreate[1/2*(M[1,1,1,1])/.MatrixElements/.{c[0]->1}//Simplify,
	9*(+((s-u)^2/(t-msq[1])^2)+(s-t)^2/(u-msq[1])^2)//Simplify
]];


report=TestReport[testList]
report["ResultsDataset"]



