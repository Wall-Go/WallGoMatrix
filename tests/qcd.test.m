(* ::Package:: *)

(*Quit[];*)


(* Check Mathematica version *)
If[$VersionNumber < 13.3,
  Print["The Mathematica testing framework requires Mathematica version ", requiredVersion," or higher. You are using version ", currentVersion, "."];
  Abort[]
];

SetDirectory[NotebookDirectory[]];
(*Put this if you want to create multiple model-files with the same kernel*)
$GroupMathMultipleModels=True;
$LoadGroupMath=True;
(*<<../DRalgo/DRalgo.m
<<../src/matrixElements.old.m*)
<<../src/WallgoMatrix.m


(* ::Chapter:: *)
(*QCD*)


(* ::Section:: *)
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


ImportModel[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Section:: *)
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
Rep1=CreateParticle[{1,2},"F"];
RepGluon=CreateParticle[{1},"V"];


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
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	UserMasses,
	UserCouplings,
	ParticleName,
	ParticleMasses,
	{TruncateAtLeadingLog->False,Format->{"json","txt"}}];


MatrixElements


(* ::Section:: *)
(*Tests*)


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})
fixConvention[arg_]:=symmetriseTU[arg/.{s->(-t-u)}]//Expand//Simplify//Expand
removeMissing[arg_]:=arg/.M[__]->0/.Missing["KeyAbsent", _]->0


(* ::Subsection:: *)
(*Test with [hep-ph/0209353]*)


testList={};
replaceSU3={Nf->3,CA->3,CF->4/3,dF->3,dA->8};


(*5 light quarks*)
AppendTo[testList,
TestCreate[
	1/2*(M[0,2,0,2]+M[0,2,2,0])/.MatrixElements/.{gs->1}//fixConvention//removeMissing,
	5/(2*CA)*(8*dF^2*CF^2/dA((s^2+u^2)/(t-mg2)^2))+5/(2*CA)*(8*dF^2*CF^2/dA((s^2+t^2)/(u-mg2)^2))/.replaceSU3//fixConvention//removeMissing
]];


(*tt->gg*)
AppendTo[testList,
TestCreate[
	1/2*(M[0,0,1,1])/.MatrixElements/.{gs->1}//fixConvention//removeMissing,
	1/(2*2*CA)*(8*dF*CF^2((t*u)/(-t+mq2)^2+(t*u)/(-u+mq2)^2)-8*dF*CF*CA((t^2+u^2)/(-s+mg2)^2))/.replaceSU3//fixConvention//removeMissing
]];


(*q1 g->q1 g*)
AppendTo[testList,
TestCreate[
	1/2*(M[0,1,0,1]+M[0,1,1,0])/.MatrixElements/.{gs->1}//fixConvention//removeMissing,
	1/(2*CA)*(-8*dF*CF^2((s*u)/(-s+mq2)^2 +s*u/(-u+mq2)^2)+8*dF*CF*CA((s^2+u^2)/(t-mg2)^2))/.replaceSU3//fixConvention//removeMissing
]];


(*tt->tt*)
AppendTo[testList,
TestCreate[
	1/2*(M[0,0,0,0])/.MatrixElements/.{gs->1}/.Thread[UserMasses->0]//fixConvention//removeMissing,
	1/(2*CA)/2(
		+8*dF^2*CF^2/dA((s^2+u^2)/(-t+mg2)^2+(s^2+t^2)/(-u+mg2)^2)+16*dF*CF(CF-CA/2)(s^2/(t*u))
		+2*(8*dF^2*CF^2/dA((s^2+u^2)/(-t+mg2)^2+(u^2+t^2)/s^2)+16*dF*CF(CF-CA/2)(u^2/(t*s)))
		)/.replaceSU3/.Thread[UserMasses->0]//fixConvention//removeMissing
]];


(*g g->g g*)
AppendTo[testList,
TestCreate[
	1/2*(M[1,1,1,1])/.MatrixElements/.{gs->1}/.Thread[UserMasses->0]//fixConvention//removeMissing,
	1/(2*2(CA^2-1))(16*dA*CA^2(3-s*u/(t-mg2)^2-s*t/(u-mg2)^2-t*u/s^2))/.replaceSU3/.Thread[UserMasses->0]//fixConvention//removeMissing
]];


report=TestReport[testList]
report["ResultsDataset"]


(* ::Subsection:: *)
(*Test hard*)


file=FileNameJoin[{NotebookDirectory[],"qcd.test.json"}];
{particleNames,parameters,FeynMatrixElements}=ImportMatrixElements[file];


testList={};


(*5 light quarks*)
AppendTo[testList,
TestCreate[
	1/2*(M[0,2,0,2]+M[0,2,2,0])/.MatrixElements//fixConvention//removeMissing,
	1/2*(M[0,2,0,2]+M[0,2,2,0])/.FeynMatrixElements//fixConvention//removeMissing
]];


(*tt->gg*)
AppendTo[testList,
TestCreate[
	1/2*(M[0,0,1,1])/.MatrixElements//fixConvention//removeMissing,
	1/2*(M[0,0,1,1])/.FeynMatrixElements//fixConvention//removeMissing
]];


(*q1 g->q1 g*)
AppendTo[testList,
TestCreate[
	1/2*(M[0,1,0,1]+M[0,1,1,0])/.MatrixElements//fixConvention//removeMissing,
	1/2*(M[0,1,0,1]+M[0,1,1,0])/.FeynMatrixElements//fixConvention//removeMissing
]];


(*tt->tt*)
AppendTo[testList,
TestCreate[
	1/2*(M[0,0,0,0])/.MatrixElements//fixConvention//removeMissing,
	1/2*(M[0,0,0,0])/.FeynMatrixElements//fixConvention//removeMissing
]];


(*g g->g g*)
AppendTo[testList,
TestCreate[
	1/2*(M[1,1,1,1])/.MatrixElements//fixConvention//removeMissing,
	1/2*(M[1,1,1,1])/.FeynMatrixElements//fixConvention//removeMissing
]];


report=TestReport[testList]
report["ResultsDataset"]



