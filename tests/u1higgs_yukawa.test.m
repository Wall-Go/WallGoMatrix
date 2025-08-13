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
(*Abelian-Higgs-Yukawa Model*)


(* ::Section:: *)
(*Model*)


Group={"U1"};
RepAdjoint={0};
RepScalar={{{1},"C"}, {{0},"R"}};
CouplingName={g1};


RepFermion={{{1},"L"},{{1},"R"},{{2},"L"},{{2},"R"}};


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


(* all cubic scalar terms *)
(* 1/3 b3 \[Chi]^3 *)
InputInv={{2,2,2},{True,True,True}};
Chi3=CreateInvariant[Group,RepScalar,InputInv][[1]];
VCubic=b3/3 CubicTerm;
(* a1/2 \[Chi] \[Phi]^2 *)
InputInv={{2,1,1},{True,True,False}};
ChiPhi2=CreateInvariant[Group,RepScalar,InputInv][[1]];
(* cubic terms together *)
VCubic=b3/3 Chi3 + a1/2 ChiPhi2;
\[Lambda]3=GradCubic[VCubic];


(* all quartic scalar terms *)
(* lam \[Phi]^4 *)
InputInv={{1,1,1,1},{True,False,True,False}};
Phi4=CreateInvariant[Group,RepScalar,InputInv][[1]];
(* 1/4 b4 \[Chi]^4 *)
InputInv={{2,2,2,2},{True,True,True,True}};
Chi4=CreateInvariant[Group,RepScalar,InputInv][[1]];
(* a2/2 \[Chi]^2 \[Phi]^2 *)
InputInv={{1,1,2,2},{True,False,True,True}};
Phi2Chi2=CreateInvariant[Group,RepScalar,InputInv][[1]];
(* quartic terms together *)
VQuartic=b4/4 Chi4 + lam Phi4 + a2/2 Phi2Chi2;
\[Lambda]4=GradQuartic[VQuartic];


(* y \[Phi]^*(Subscript[\[Psi], L]Subscript[\[Xi], R]^++Subscript[\[Psi], R]Subscript[\[Xi], L]^+)*)
InputInv={{1,1,4},{True,True,False}};
Yukawa1=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
InputInv={{1,2,3},{False,False,True}};
Yukawa2=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
InputInv={{1,1,4},{False,True,False}};
Yukawa1HC=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
InputInv={{1,2,3},{True,False,True}};
Yukawa2HC=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;


Ysff=y*GradYukawa[Yukawa1+Yukawa2];
YsffC=y*GradYukawa[Yukawa1HC+Yukawa2HC];


Ysff//MatrixForm
YsffC//MatrixForm


ImportModel[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Section:: *)
(*User Input*)


(* ::Subsection:: *)
(*UserInput*)


(*
In DRalgo fermions are Weyl.
So to create one Dirac we need
one left-handed and
one right-handed fermion
*)


vev={v,0,0};
SymmetryBreaking[vev]


(*Below
rep 1-2 are fermions,
(*rep 3 is a scalar*)
*)
(* scalar *)
RepPhi=CreateParticle[{1},"S",ms,"Phi"];
RepChi=CreateParticle[{2},"S",ms,"Chi"];

(* left-handed fermion *)
RepPsi=CreateParticle[{1,2},"F",mf,"Psi"];
RepXi=CreateParticle[{3,4},"F",mf,"Xi"];

(*Vector bosons*)
RepA=CreateParticle[{1},"V",mv,"VectorU1"];


(*
These particles do not necessarily have to be out of equilibrium
*)
ParticleList={RepPhi,RepChi,RepPsi,RepXi,RepA};
LightParticleList={};


(*Defining various masses and couplings*)
UserMasses={mv,mf,ms};


(*
	output of matrix elements
*)
OutputFile="output/matrixElements.u1_higgs_yukawa";
SetDirectory[NotebookDirectory[]];
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{
		TruncateAtLeadingLog->False,
		NormalizeWithDOF->False,
		Format->{"json","txt"}}];


MatrixElements


(* ::Section:: *)
(*Importing results from WallGo*)


{particles,parameters,MatrixElements}=ImportMatrixElements["output/matrixElements.u1_higgs_yukawa.json"];


(* ::Section:: *)
(*Importing results from FeynCalc*)


{particlesFeyn,parametersFeyn,MatrixElementsFeyn}=ImportMatrixElements["u1-higgs-yukawa.test.json"];


(* ::Section:: *)
(*Comparison tests*)


(* ::Subsection:: *)
(*Translate input*)


insertCouplings={Global`g1->g,\[Lambda]->lam,SUNN->2,gu1->0,mChi->0,mPhi->0};
customCouplings={ms2->mPhi^2};


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})


fixConvention[arg_]:=symmetriseTU[
	arg/.customCouplings/.Thread[UserMasses->0]/.insertCouplings/.{s->(-t-u)}
	]//Expand//Simplify//Expand


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


(* ::Subsubsection::Closed:: *)
(*VVtoVV*)


process="VV->VV"
test["WallGo"][process]=testWallGo[
	{"VectorU1"},
	{"VectorU1"},
	{"VectorU1"},
	{"VectorU1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"A"},
	{"A"},
	{"A"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*FFtoFF*)


process="F1F1->F1F1"
test["WallGo"][process]=testWallGo[
	{"Psi"},
	{"Psi"},
	{"Psi"},
	{"Psi"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Psi","Psibar"},
	{"Psi","Psibar"},
	{"Psi","Psibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* F2F2->F2F2 *)
process="F2F2->F2F2"
test["WallGo"][process]=testWallGo[
	{"Xi"},
	{"Xi"},
	{"Xi"},
	{"Xi"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Xi","Xibar"},
	{"Xi","Xibar"},
	{"Xi","Xibar"},
	{"Xi","Xibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F1F2->F1F2"
test["WallGo"][process]=testWallGo[
	{"Psi"},
	{"Xi"},
	{"Psi"},
	{"Xi"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Xi","Xibar"},
	{"Psi","Psibar"},
	{"Xi","Xibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F1F1->F2F2"
test["WallGo"][process]=testWallGo[
	{"Psi"},
	{"Psi"},
	{"Xi"},
	{"Xi"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Psi","Psibar"},
	{"Xi","Xibar"},
	{"Xi","Xibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*FFtoVV*)


process="F1F1->VV"
test["WallGo"][process]=testWallGo[
	{"Psi"},
	{"Psi"},
	{"VectorU1"},
	{"VectorU1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Psi","Psibar"},
	{"A"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*FVtoFV*)


process="F1V->F1V"
test["WallGo"][process]=testWallGo[
	{"Psi"},
	{"VectorU1"},
	{"Psi"},
	{"VectorU1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"A"},
	{"Psi","Psibar"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*SStoSS*)


process="SS->SS"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"Phi"},
	{"Phi"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Phi","Phibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="S2S2->S2S2"
test["WallGo"][process]=testWallGo[
	{"Chi"},
	{"Chi"},
	{"Chi"},
	{"Chi"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Chi"},
	{"Chi"},
	{"Chi"},
	{"Chi"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="S1S1->S2S2"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"Chi"},
	{"Chi"}
]//Simplify
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Chi"},
	{"Chi"}
]//Simplify
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*SStoFF*)


process="SS->F1F1"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"Psi"},
	{"Psi"}
]/.{flag1->1,flag2->1,flag3->1,flag4->1}
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Psi","Psibar"},
	{"Psi","Psibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="S1S1->F2F2"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"Xi"},
	{"Xi"}
]/.{flag1->1,flag2->1,flag3->1,flag4->1}
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Xi","Xibar"},
	{"Xi","Xibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*FStoFS*)


process="F1S->F1S"
test["WallGo"][process]=testWallGo[
	{"Psi"},
	{"Phi"},
	{"Psi"},
	{"Phi"}
]/.{flag1->1,flag2->1,flag3->1,flag4->1}
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"Psi","Psibar"},
	{"Phi","Phibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="SS->F1F1_flippingRule"
test["WallGo"][process]=-testWallGo[
	{"Phi"},
	{"Phi"},
	{"Psi"},
	{"Psi"}
]/.{s->s1,t->t1,u->uu1}/.{s1->t,t1->s,uu1->u}//fixConvention[#]&
process="F1S->F1S_flippingRule"
test["WallGo"][process]=testWallGo[
	{"Psi"},
	{"Phi"},
	{"Psi"},
	{"Phi"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"]["SS->F1F1_flippingRule"]//Evaluate,
		test["WallGo"]["F1S->F1S_flippingRule"],
		TestID->"WallGo vs WallGo: "<>process]];


process="SS->F1F1_flippingRule"
test["FeynCalc"][process]=-testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Psi","Psibar"},
	{"Psi","Psibar"}
]/.{s->s1,t->t1,u->uu1}/.{s1->t,t1->s,uu1->u}//fixConvention[#]&
process="F1S->F1S_flippingRule"
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"Psi","Psibar"},
	{"Phi","Phibar"}
]
AppendTo[testList,
	TestCreate[
		test["FeynCalc"]["SS->F1F1_flippingRule"]//Evaluate,
		test["FeynCalc"]["F1S->F1S_flippingRule"],
		TestID->"FeynCalc vs FeynCalc: "<>process]];


(* ::Subsubsection:: *)
(*FStoFV*)


process="F1S->F1V"
test["WallGo"][process]=testWallGo[
	{"Psi"},
	{"Phi"},
	{"Psi"},
	{"VectorU1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"Psi","Psibar"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F1S->F2V"
test["WallGo"][process]=testWallGo[
	{"Psi"},
	{"Phi"},
	{"Xi"},
	{"VectorU1"}
]/.{(*flag[x_]->1*)}
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"Xi","Xibar"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection:: *)
(*FFtoSV*)


process="F1F1->SV"
test["WallGo"][process]=testWallGo[
	{"Psi"},
	{"Psi"},
	{"Phi"},
	{"VectorU1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F1F2->SV"
test["WallGo"][process]=testWallGo[
	{"Psi"},
	{"Xi"},
	{"Phi"},
	{"VectorU1"}
]/.{(*flag[x_]->1*)}
test["FeynCalc"][process]=testFeynCalc[
	{"A"},
	{"Phi","Phibar"},
	{"Xi","Xibar"},
	{"Psi","Psibar"}	
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*SStoVV*)


process="SS->VV"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"VectorU1"},
	{"VectorU1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"A"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*SVtoSV*)


process="SV->SV"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"VectorU1"},
	{"Phi"},
	{"VectorU1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"A"},
	{"Phi","Phibar"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*SStoSV*)


process="SS->SV"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"Phi"},
	{"VectorU1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="SS->SV"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"Chi"},
	{"Phi"},
	{"VectorU1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"Chi"},
	{"Phi","Phibar"},
	{"A"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection:: *)
(*Test report*)


report=TestReport[testList]
report["ResultsDataset"]


(* ::Section:: *)
(*Test g^2 y^2*)


particles
particlesFeyn


expandNumber[expr_, num_, repl_List] := 
  expr /. M[args__] :> Module[{lst = ({args} /. num -> repl)},
    lst = Map[If[ListQ[#], #, {#}] &, lst]; (* ensure each slot is a list *)
    Total[M @@@ Tuples[lst]]
  ];

elements=Select[MatrixElements//Expand, ! FreeQ[#, g1^2 y^2] &][[All,1]]
test1=elements/.MatrixElements//fixConvention;
test2=elements/.{1->2,2->3,3->5,4->7};
test2=test2//expandNumber[#, 0, {0, 1}]&//expandNumber[#, 3, {3, 4}]&//expandNumber[#, 5, {5, 6}]&;
test2=test2/.MatrixElementsFeyn//removeMissing//fixConvention;

((s1*test1-s2*test2)//Simplify)/.(s1-s2)->0


(* ::Section:: *)
(*Full tests*)


(* ::Subsection::Closed:: *)
(*FeynMatrixElements hard copy*)


FeynMatrixElements={
 M[0, 0, 0, 0] -> 16*lam^2 + (a1^4*(t + u)^2)/(16*t^2*u^2) + (8*g^2*lam*(t^2 + u^2 - s*(t + u)))/(t*u) + (g^4*(t^2 + u^2 - s*(t + u))^2)/(t^2*u^2) + 
   a1^2*(2*lam*(t^(-1) + u^(-1)) + (g^2*(t + u)*(t^2 + u^2 - s*(t + u)))/(2*t^2*u^2)),
 M[0, 1, 0, 1] -> 16*lam^2 + (a1^4*(s + t)^2)/(16*s^2*t^2) + (8*g^2*lam*(s^2 + t*(t - u) - s*u))/(s*t) + 
   (g^4*(s^2 + t*(t - u) - s*u)^2)/(s^2*t^2) + a1^2*(2*lam*(s^(-1) + t^(-1)) + (g^2*(s + t)*(s^2 + t*(t - u) - s*u))/(2*s^2*t^2)), 
 M[0, 1, 1, 0] -> 16*lam^2 + (a1^4*(s + u)^2)/(16*s^2*u^2) + (8*g^2*lam*(s^2 - s*t + u*(-t + u)))/(s*u) + (g^4*(s^2 - s*t + u*(-t + u))^2)/(s^2*u^2) + 
   a1^2*(2*lam*(s^(-1) + u^(-1)) + (g^2*(s + u)*(s^2 - s*t + u*(-t + u)))/(2*s^2*u^2)),
 M[0, 1, 2, 2] -> a2^2 + (2*a1*a2*b3)/s + (a1^3*b3*(t + u))/(2*s*t*u) + (a1^4*(t + u)^2)/(16*t^2*u^2) + 
   a1^2*(b3^2/s^2 + (a2*(t + u))/(2*t*u)),
 M[0, 1, 2, 7] -> (a1^2*g^2*(s + t - 3*u)*(t + u))/(4*t*u^2),
 M[0, 1, 3, 4] -> (8*g^4*t*u)/s^2 - (8*g^2*t*y^2)/s + (2*t*y^4)/u, 
 M[0, 1, 4, 3] -> (8*g^4*t*u)/s^2 - (8*g^2*u*y^2)/s + (2*u*y^4)/t,
 M[0, 1, 5, 6] -> (32*g^4*t*u)/s^2 + (16*g^2*u*y^2)/s + (2*u*y^4)/t, M[0, 1, 6, 5] -> (32*g^4*t*u)/s^2 + (16*g^2*t*y^2)/s + (2*t*y^4)/u, 
 M[0, 1, 7, 2] -> (a1^2*g^2*(s - 3*t + u)*(t + u))/(4*t^2*u),
 M[0, 1, 7, 7] -> (g^4*(s^2 + t^2 + 18*t*u + u^2 + 2*s*(t + u)))/(2*t*u), 
 M[0, 2, 0, 2] -> a2^2 + (2*a1*a2*b3)/t + (a1^3*b3*(s + u))/(2*s*t*u) + (a1^4*(s + u)^2)/(16*s^2*u^2) + a1^2*(b3^2/t^2 + (a2*(s + u))/(2*s*u)),
 M[0, 2, 0, 7] -> (a1^2*g^2*(s + t - 3*u)*(s + u))/(4*s*u^2), 
 M[0, 2, 2, 0] -> a2^2 + (a1^4*(s + t)^2)/(16*s^2*t^2) + a1^2*((a2*(s + t))/(2*s*t) + b3^2/u^2) + (2*a1*a2*b3)/u + (a1^3*b3*(s + t))/(2*s*t*u),
 M[0, 2, 4, 5] -> (a1^2*y^2)/(2*s), M[0, 2, 5, 4] -> (a1^2*y^2)/(2*s), 
 M[0, 2, 7, 0] -> (a1^2*g^2*(s + t)*(s - 3*t + u))/(4*s*t^2),
 M[0, 3, 0, 3] -> (-8*g^4*s*u)/t^2 + (8*g^2*u*y^2)/t - (2*u*y^4)/s,
 M[0, 3, 2, 5] -> -1/2*(a1^2*y^2)/t, 
 M[0, 3, 3, 0] -> (-8*g^4*s*t)/u^2 + (8*g^2*t*y^2)/u - (2*t*y^4)/s,
 M[0, 3, 5, 2] -> -1/2*(a1^2*y^2)/u,
 M[0, 3, 5, 7] -> (-4*g^2*(t^3 + 2*s^2*u + 8*t^2*u + 5*t*u^2 + 2*u^3 + s*(t^2 + 8*t*u + 3*u^2))*y^2)/(s*t*u), 
 M[0, 3, 7, 5] -> (-4*g^2*(2*s^2*t + 2*t^3 + 5*t^2*u + 8*t*u^2 + u^3 + s*(3*t^2 + 8*t*u + u^2))*y^2)/(s*t*u),
 M[0, 4, 0, 4] -> (-8*g^4*s*u)/t^2 + (8*g^2*s*y^2)/t - (2*s*y^4)/u, 
 M[0, 4, 4, 0] -> (-8*g^4*s*t)/u^2 + (8*g^2*s*y^2)/u - (2*s*y^4)/t,
 M[0, 5, 0, 5] -> (-32*g^4*s*u)/t^2 - (16*g^2*s*y^2)/t - (2*s*y^4)/u,
 M[0, 5, 5, 0] -> (-32*g^4*s*t)/u^2 - (16*g^2*s*y^2)/u - (2*s*y^4)/t, 
 M[0, 6, 0, 6] -> (-32*g^4*s*u)/t^2 - (16*g^2*u*y^2)/t - (2*u*y^4)/s,
 M[0, 6, 2, 4] -> -1/2*(a1^2*y^2)/t, M[0, 6, 4, 2] -> -1/2*(a1^2*y^2)/u, 
 M[0, 6, 4, 7] -> (2*g^2*(t^3 - 16*s^2*u - 4*t^2*u - 7*t*u^2 - 4*u^3 + s*(t^2 - 19*t*u - 12*u^2))*y^2)/(s*t*u),
 M[0, 6, 6, 0] -> (-32*g^4*s*t)/u^2 - (16*g^2*t*y^2)/u - (2*t*y^4)/s, 
 M[0, 6, 7, 4] -> (-2*g^2*(16*s^2*t + 4*t^3 + 7*t^2*u + 4*t*u^2 - u^3 + s*(12*t^2 + 19*t*u - u^2))*y^2)/(s*t*u),
 M[0, 7, 0, 2] -> -1/4*(a1^2*g^2*(3*s - t - u)*(s + u))/(s^2*u), 
 M[0, 7, 0, 7] -> (g^4*(s^2 + (t + u)^2 + 2*s*(t + 9*u)))/(2*s*u), M[0, 7, 2, 0] -> -1/4*(a1^2*g^2*(s + t)*(3*s - t - u))/(s^2*t), 
 M[0, 7, 4, 5] -> (2*g^2*(4*s^3 + t^3 + t^2*u - 2*t*u^2 - 2*u^3 + s^2*(5*t + 14*u) + s*(4*t^2 + 19*t*u + 16*u^2))*y^2)/(s*t*u), 
 M[0, 7, 5, 4] -> (2*g^2*(4*s^3 - 2*t^3 - 2*t^2*u + t*u^2 + u^3 + s^2*(14*t + 5*u) + s*(16*t^2 + 19*t*u + 4*u^2))*y^2)/(s*t*u),
 M[0, 7, 7, 0] -> (g^4*(s^2 + (t + u)^2 + 2*s*(9*t + u)))/(2*s*t), 
 M[1, 0, 0, 1] -> 16*lam^2 + (a1^4*(s + u)^2)/(16*s^2*u^2) + (8*g^2*lam*(s^2 - s*t + u*(-t + u)))/(s*u) + (g^4*(s^2 - s*t + u*(-t + u))^2)/(s^2*u^2) + 
   a1^2*(2*lam*(s^(-1) + u^(-1)) + (g^2*(s + u)*(s^2 - s*t + u*(-t + u)))/(2*s^2*u^2)),
 M[1, 0, 1, 0] -> 16*lam^2 + (a1^4*(s + t)^2)/(16*s^2*t^2) + (8*g^2*lam*(s^2 + t*(t - u) - s*u))/(s*t) + 
   (g^4*(s^2 + t*(t - u) - s*u)^2)/(s^2*t^2) + a1^2*(2*lam*(s^(-1) + t^(-1)) + (g^2*(s + t)*(s^2 + t*(t - u) - s*u))/(2*s^2*t^2)), 
 M[1, 0, 2, 2] -> a2^2 + (2*a1*a2*b3)/s + (a1^3*b3*(t + u))/(2*s*t*u) + (a1^4*(t + u)^2)/(16*t^2*u^2) + a1^2*(b3^2/s^2 + (a2*(t + u))/(2*t*u)),
 M[1, 0, 2, 7] -> (a1^2*g^2*(s + t - 3*u)*(t + u))/(4*t*u^2), 
 M[1, 0, 3, 4] -> (8*g^4*t*u)/s^2 - (8*g^2*u*y^2)/s + (2*u*y^4)/t,
 M[1, 0, 4, 3] -> (8*g^4*t*u)/s^2 - (8*g^2*t*y^2)/s + (2*t*y^4)/u,
 M[1, 0, 5, 6] -> (32*g^4*t*u)/s^2 + (16*g^2*t*y^2)/s + (2*t*y^4)/u, 
 M[1, 0, 6, 5] -> (32*g^4*t*u)/s^2 + (16*g^2*u*y^2)/s + (2*u*y^4)/t,
 M[1, 0, 7, 2] -> (a1^2*g^2*(s - 3*t + u)*(t + u))/(4*t^2*u),
 M[1, 0, 7, 7] -> (g^4*(s^2 + t^2 + 18*t*u + u^2 + 2*s*(t + u)))/(2*t*u), 
 M[1, 1, 1, 1] -> 16*lam^2 + (a1^4*(t + u)^2)/(16*t^2*u^2) + (8*g^2*lam*(t^2 + u^2 - s*(t + u)))/(t*u) + (g^4*(t^2 + u^2 - s*(t + u))^2)/(t^2*u^2) + 
   a1^2*(2*lam*(t^(-1) + u^(-1)) + (g^2*(t + u)*(t^2 + u^2 - s*(t + u)))/(2*t^2*u^2)),
 M[1, 2, 1, 2] -> a2^2 + (2*a1*a2*b3)/t + (a1^3*b3*(s + u))/(2*s*t*u) + (a1^4*(s + u)^2)/(16*s^2*u^2) + 
   a1^2*(b3^2/t^2 + (a2*(s + u))/(2*s*u)),
 M[1, 2, 1, 7] -> (a1^2*g^2*(s + t - 3*u)*(s + u))/(4*s*u^2), 
 M[1, 2, 2, 1] -> a2^2 + (a1^4*(s + t)^2)/(16*s^2*t^2) + a1^2*((a2*(s + t))/(2*s*t) + b3^2/u^2) + (2*a1*a2*b3)/u + (a1^3*b3*(s + t))/(2*s*t*u),
 M[1, 2, 3, 6] -> (a1^2*y^2)/(2*s),
 M[1, 2, 6, 3] -> (a1^2*y^2)/(2*s), 
 M[1, 2, 7, 1] -> (a1^2*g^2*(s + t)*(s - 3*t + u))/(4*s*t^2),
 M[1, 3, 1, 3] -> (-8*g^4*s*u)/t^2 + (8*g^2*s*y^2)/t - (2*s*y^4)/u,
 M[1, 3, 3, 1] -> (-8*g^4*s*t)/u^2 + (8*g^2*s*y^2)/u - (2*s*y^4)/t, 
 M[1, 4, 1, 4] -> (-8*g^4*s*u)/t^2 + (8*g^2*u*y^2)/t - (2*u*y^4)/s,
 M[1, 4, 2, 6] -> -1/2*(a1^2*y^2)/t, M[1, 4, 4, 1] -> (-8*g^4*s*t)/u^2 + (8*g^2*t*y^2)/u - (2*t*y^4)/s,
 M[1, 4, 6, 2] -> -1/2*(a1^2*y^2)/u, 
 M[1, 4, 6, 7] -> (-4*g^2*(t^3 + 2*s^2*u + 8*t^2*u + 5*t*u^2 + 2*u^3 + s*(t^2 + 8*t*u + 3*u^2))*y^2)/(s*t*u),
 M[1, 4, 7, 6] -> (-4*g^2*(2*s^2*t + 2*t^3 + 5*t^2*u + 8*t*u^2 + u^3 + s*(3*t^2 + 8*t*u + u^2))*y^2)/(s*t*u), 
 M[1, 5, 1, 5] -> (-32*g^4*s*u)/t^2 - (16*g^2*u*y^2)/t - (2*u*y^4)/s,
 M[1, 5, 2, 3] -> -1/2*(a1^2*y^2)/t,
 M[1, 5, 3, 2] -> -1/2*(a1^2*y^2)/u, 
 M[1, 5, 3, 7] -> (2*g^2*(t^3 - 16*s^2*u - 4*t^2*u - 7*t*u^2 - 4*u^3 + s*(t^2 - 19*t*u - 12*u^2))*y^2)/(s*t*u),
 M[1, 5, 5, 1] -> (-32*g^4*s*t)/u^2 - (16*g^2*t*y^2)/u - (2*t*y^4)/s, 
 M[1, 5, 7, 3] -> (-2*g^2*(16*s^2*t + 4*t^3 + 7*t^2*u + 4*t*u^2 - u^3 + s*(12*t^2 + 19*t*u - u^2))*y^2)/(s*t*u),
 M[1, 6, 1, 6] -> (-32*g^4*s*u)/t^2 - (16*g^2*s*y^2)/t - (2*s*y^4)/u, 
 M[1, 6, 6, 1] -> (-32*g^4*s*t)/u^2 - (16*g^2*s*y^2)/u - (2*s*y^4)/t,
 M[1, 7, 1, 2] -> -1/4*(a1^2*g^2*(3*s - t - u)*(s + u))/(s^2*u),
 M[1, 7, 1, 7] -> (g^4*(s^2 + (t + u)^2 + 2*s*(t + 9*u)))/(2*s*u), 
 M[1, 7, 2, 1] -> -1/4*(a1^2*g^2*(s + t)*(3*s - t - u))/(s^2*t),
 M[1, 7, 3, 6] -> (2*g^2*(4*s^3 + t^3 + t^2*u - 2*t*u^2 - 2*u^3 + s^2*(5*t + 14*u) + s*(4*t^2 + 19*t*u + 16*u^2))*y^2)/(s*t*u), 
 M[1, 7, 6, 3] -> (2*g^2*(4*s^3 - 2*t^3 - 2*t^2*u + t*u^2 + u^3 + s^2*(14*t + 5*u) + s*(16*t^2 + 19*t*u + 4*u^2))*y^2)/(s*t*u),
 M[1, 7, 7, 1] -> (g^4*(s^2 + (t + u)^2 + 2*s*(9*t + u)))/(2*s*t), 
 M[2, 0, 0, 2] -> a2^2 + (a1^4*(s + t)^2)/(16*s^2*t^2) + a1^2*((a2*(s + t))/(2*s*t) + b3^2/u^2) + (2*a1*a2*b3)/u + (a1^3*b3*(s + t))/(2*s*t*u),
 M[2, 0, 0, 7] -> (a1^2*g^2*u)/(s*t), 
 M[2, 0, 2, 0] -> a2^2 + (2*a1*a2*b3)/t + (a1^3*b3*(s + u))/(2*s*t*u) + (a1^4*(s + u)^2)/(16*s^2*u^2) + a1^2*(b3^2/t^2 + (a2*(s + u))/(2*s*u)),
 M[2, 0, 4, 5] -> (a1^2*y^2)/(2*s),
 M[2, 0, 5, 4] -> (a1^2*y^2)/(2*s), 
 M[2, 0, 7, 0] -> (a1^2*g^2*t)/(s*u),
 M[2, 1, 1, 2] -> a2^2 + (a1^4*(s + t)^2)/(16*s^2*t^2) + a1^2*((a2*(s + t))/(2*s*t) + b3^2/u^2) + (2*a1*a2*b3)/u + (a1^3*b3*(s + t))/(2*s*t*u),
 M[2, 1, 1, 7] -> (a1^2*g^2*u)/(s*t), 
 M[2, 1, 2, 1] -> a2^2 + (2*a1*a2*b3)/t + (a1^3*b3*(s + u))/(2*s*t*u) + (a1^4*(s + u)^2)/(16*s^2*u^2) + a1^2*(b3^2/t^2 + (a2*(s + u))/(2*s*u)),
 M[2, 1, 3, 6] -> (a1^2*y^2)/(2*s),
 M[2, 1, 6, 3] -> (a1^2*y^2)/(2*s), 
 M[2, 1, 7, 1] -> (a1^2*g^2*t)/(s*u),
 M[2, 2, 0, 1] -> a2^2 + (2*a1*a2*b3)/s + (a1^3*b3*(t + u))/(2*s*t*u) + (a1^4*(t + u)^2)/(16*t^2*u^2) + a1^2*(b3^2/s^2 + (a2*(t + u))/(2*t*u)), 
 M[2, 2, 1, 0] -> a2^2 + (2*a1*a2*b3)/s + (a1^3*b3*(t + u))/(2*s*t*u) + (a1^4*(t + u)^2)/(16*t^2*u^2) + a1^2*(b3^2/s^2 + (a2*(t + u))/(2*t*u)), 
 M[2, 2, 2, 2] -> 36*b4^2 + 48*b3^2*b4*(s^(-1) + t^(-1) + u^(-1)) + (16*b3^4*(t*u + s*(t + u))^2)/(s^2*t^2*u^2),
 M[2, 3, 1, 5] -> -1/2*(a1^2*y^2)/t,
 M[2, 3, 5, 1] -> -1/2*(a1^2*y^2)/u,
 M[2, 4, 0, 6] -> -1/2*(a1^2*y^2)/t, 
 M[2, 4, 6, 0] -> -1/2*(a1^2*y^2)/u,
 M[2, 5, 0, 3] -> -1/2*(a1^2*y^2)/t,
 M[2, 5, 3, 0] -> -1/2*(a1^2*y^2)/u,
 M[2, 6, 1, 4] -> -1/2*(a1^2*y^2)/t,
 M[2, 6, 4, 1] -> -1/2*(a1^2*y^2)/u,
 M[2, 7, 0, 1] -> (a1^2*g^2*s)/(t*u), 
 M[2, 7, 1, 0] -> (a1^2*g^2*s)/(t*u),
 M[3, 0, 0, 3] -> (-8*g^4*s*t)/u^2 + (8*g^2*t*y^2)/u - (2*t*y^4)/s,
 M[3, 0, 2, 5] -> -1/2*(a1^2*y^2)/u,
 M[3, 0, 3, 0] -> (-8*g^4*s*u)/t^2 + (8*g^2*u*y^2)/t - (2*u*y^4)/s, 
 M[3, 0, 5, 2] -> -1/2*(a1^2*y^2)/t,
 M[3, 0, 5, 7] -> (2*g^2*(s^3 + 8*t^3 - 4*t^2*u - 24*t*u^2 - 4*u^3 + s^2*(-10*t + 4*u) - s*(t^2 + 26*t*u + u^2))*y^2)/(s*t*u), 
 M[3, 0, 7, 5] -> (2*g^2*(s^3 + 2*s^2*(2*t - 5*u) - s*(t^2 + 26*t*u + u^2) - 4*(t^3 + 6*t^2*u + t*u^2 - 2*u^3))*y^2)/(s*t*u),
 M[3, 1, 1, 3] -> (-8*g^4*s*t)/u^2 + (8*g^2*s*y^2)/u - (2*s*y^4)/t, 
 M[3, 1, 3, 1] -> (-8*g^4*s*u)/t^2 + (8*g^2*s*y^2)/t - (2*s*y^4)/u,
 M[3, 2, 1, 5] -> -1/2*(a1^2*y^2)/u, M[3, 2, 5, 1] -> -1/2*(a1^2*y^2)/t,
 M[3, 3, 3, 3] -> (8*g^4*(t^4 + u^4 + s^2*(t + u)^2))/(t^2*u^2), 
 M[3, 4, 0, 1] -> (8*g^4*t*u)/s^2 - (8*g^2*t*y^2)/s + (2*t*y^4)/u,
 M[3, 4, 1, 0] -> (8*g^4*t*u)/s^2 - (8*g^2*u*y^2)/s + (2*u*y^4)/t,
 M[3, 4, 3, 4] -> (8*g^4*(s^4 + t^4 + s^2*u^2 + 2*s*t*u^2 + t^2*u^2))/(s^2*t^2), 
 M[3, 4, 4, 3] -> (8*g^4*(s^4 + s^2*t^2 + 2*s*t^2*u + t^2*u^2 + u^4))/(s^2*u^2),
 M[3, 4, 5, 6] -> (32*g^4*(t^2 + u^2))/s^2 + (16*g^2*t*y^2)/s + 4*y^4,
 M[3, 4, 6, 5] -> (32*g^4*(t^2 + u^2))/s^2 + (16*g^2*u*y^2)/s + 4*y^4, 
 M[3, 4, 7, 7] -> (8*g^4*(-2*s^2 + 7*t^2 + 12*t*u + 7*u^2 + 4*s*(t + u)))/(t*u),
 M[3, 5, 3, 5] -> (32*g^4*(s^2 + u^2))/t^2 + (16*g^2*u*y^2)/t + 4*y^4,
 M[3, 5, 5, 3] -> (32*g^4*(s^2 + t^2))/u^2 + (16*g^2*t*y^2)/u + 4*y^4, 
 M[3, 6, 1, 2] -> (a1^2*y^2)/(2*s),
 M[3, 6, 1, 7] -> (-2*g^2*(6*s^3 + t^3 + 8*t^2*u + 11*t*u^2 + 4*u^3 - s^2*(3*t + 10*u) - 2*s*(5*t^2 + 13*t*u + 10*u^2))*y^2)/(s*t*u),
 M[3, 6, 2, 1] -> (a1^2*y^2)/(2*s), 
 M[3, 6, 3, 6] -> (32*g^4*(s^2 + u^2))/t^2 + (16*g^2*s*y^2)/t + 4*y^4,
 M[3, 6, 6, 3] -> (32*g^4*(s^2 + t^2))/u^2 + (16*g^2*s*y^2)/u + 4*y^4, 
 M[3, 6, 7, 1] -> (-2*g^2*(6*s^3 + 4*t^3 + 11*t^2*u + 8*t*u^2 + u^3 - s^2*(10*t + 3*u) - 2*s*(10*t^2 + 13*t*u + 5*u^2))*y^2)/(s*t*u), 
 M[3, 7, 1, 5] -> (2*g^2*(t^3 + s^2*(5*t - 22*u) - 10*t^2*u - 5*t*u^2 + 4*u^3 + 2*s*(3*t^2 - 13*t*u - 5*u^2))*y^2)/(s*t*u),
 M[3, 7, 3, 7] -> (-8*g^4*(6*s^2 - 2*t^2 + 2*t*u + 5*u^2 + 3*s*(t + 3*u)))/(s*u), 
 M[3, 7, 5, 1] -> (2*g^2*(4*t^3 - 5*t^2*u - 10*t*u^2 + u^3 + s^2*(-22*t + 5*u) - 2*s*(5*t^2 + 13*t*u - 3*u^2))*y^2)/(s*t*u),
 M[3, 7, 7, 3] -> (-8*g^4*(6*s^2 + 5*t^2 + 2*t*u - 2*u^2 + 3*s*(3*t + u)))/(s*t), 
 M[4, 0, 0, 4] -> (-8*g^4*s*t)/u^2 + (8*g^2*s*y^2)/u - (2*s*y^4)/t,
 M[4, 0, 4, 0] -> (-8*g^4*s*u)/t^2 + (8*g^2*s*y^2)/t - (2*s*y^4)/u,
 M[4, 1, 1, 4] -> (-8*g^4*s*t)/u^2 + (8*g^2*t*y^2)/u - (2*t*y^4)/s, 
 M[4, 1, 2, 6] -> -1/2*(a1^2*y^2)/u, M[4, 1, 4, 1] -> (-8*g^4*s*u)/t^2 + (8*g^2*u*y^2)/t - (2*u*y^4)/s,
 M[4, 1, 6, 2] -> -1/2*(a1^2*y^2)/t, 
 M[4, 1, 6, 7] -> (2*g^2*(s^3 + 8*t^3 - 4*t^2*u - 24*t*u^2 - 4*u^3 + s^2*(-10*t + 4*u) - s*(t^2 + 26*t*u + u^2))*y^2)/(s*t*u), 
 M[4, 1, 7, 6] -> (2*g^2*(s^3 + 2*s^2*(2*t - 5*u) - s*(t^2 + 26*t*u + u^2) - 4*(t^3 + 6*t^2*u + t*u^2 - 2*u^3))*y^2)/(s*t*u),
 M[4, 2, 0, 6] -> -1/2*(a1^2*y^2)/u,
 M[4, 2, 6, 0] -> -1/2*(a1^2*y^2)/t, 
 M[4, 3, 0, 1] -> (8*g^4*t*u)/s^2 - (8*g^2*u*y^2)/s + (2*u*y^4)/t,
 M[4, 3, 1, 0] -> (8*g^4*t*u)/s^2 - (8*g^2*t*y^2)/s + (2*t*y^4)/u,
 M[4, 3, 3, 4] -> (8*g^4*(s^4 + s^2*t^2 + 2*s*t^2*u + t^2*u^2 + u^4))/(s^2*u^2), 
 M[4, 3, 4, 3] -> (8*g^4*(s^4 + t^4 + s^2*u^2 + 2*s*t*u^2 + t^2*u^2))/(s^2*t^2),
 M[4, 3, 5, 6] -> (32*g^4*(t^2 + u^2))/s^2 + (16*g^2*u*y^2)/s + 4*y^4,
 M[4, 3, 6, 5] -> (32*g^4*(t^2 + u^2))/s^2 + (16*g^2*t*y^2)/s + 4*y^4, 
 M[4, 3, 7, 7] -> (8*g^4*(-2*s^2 + 7*t^2 + 12*t*u + 7*u^2 + 4*s*(t + u)))/(t*u),
 M[4, 4, 4, 4] -> (8*g^4*(t^4 + u^4 + s^2*(t + u)^2))/(t^2*u^2),
 M[4, 5, 0, 2] -> (a1^2*y^2)/(2*s), 
 M[4, 5, 0, 7] -> (-2*g^2*(6*s^3 + t^3 + 8*t^2*u + 11*t*u^2 + 4*u^3 - s^2*(3*t + 10*u) - 2*s*(5*t^2 + 13*t*u + 10*u^2))*y^2)/(s*t*u),
 M[4, 5, 2, 0] -> (a1^2*y^2)/(2*s), 
 M[4, 5, 4, 5] -> (32*g^4*(s^2 + u^2))/t^2 + (16*g^2*s*y^2)/t + 4*y^4,
 M[4, 5, 5, 4] -> (32*g^4*(s^2 + t^2))/u^2 + (16*g^2*s*y^2)/u + 4*y^4, 
 M[4, 5, 7, 0] -> (-2*g^2*(6*s^3 + 4*t^3 + 11*t^2*u + 8*t*u^2 + u^3 - s^2*(10*t + 3*u) - 2*s*(10*t^2 + 13*t*u + 5*u^2))*y^2)/(s*t*u),
 M[4, 6, 4, 6] -> (32*g^4*(s^2 + u^2))/t^2 + (16*g^2*u*y^2)/t + 4*y^4, 
 M[4, 6, 6, 4] -> (32*g^4*(s^2 + t^2))/u^2 + (16*g^2*t*y^2)/u + 4*y^4,
 M[4, 7, 0, 6] -> (2*g^2*(t^3 + s^2*(5*t - 22*u) - 10*t^2*u - 5*t*u^2 + 4*u^3 + 2*s*(3*t^2 - 13*t*u - 5*u^2))*y^2)/(s*t*u), 
 M[4, 7, 4, 7] -> (-8*g^4*(6*s^2 - 2*t^2 + 2*t*u + 5*u^2 + 3*s*(t + 3*u)))/(s*u),
 M[4, 7, 6, 0] -> (2*g^2*(4*t^3 - 5*t^2*u - 10*t*u^2 + u^3 + s^2*(-22*t + 5*u) - 2*s*(5*t^2 + 13*t*u - 3*u^2))*y^2)/(s*t*u), 
 M[4, 7, 7, 4] -> (-8*g^4*(6*s^2 + 5*t^2 + 2*t*u - 2*u^2 + 3*s*(3*t + u)))/(s*t),
 M[5, 0, 0, 5] -> (-32*g^4*s*t)/u^2 - (16*g^2*s*y^2)/u - (2*s*y^4)/t,
 M[5, 0, 5, 0] -> (-32*g^4*s*u)/t^2 - (16*g^2*s*y^2)/t - (2*s*y^4)/u, 
 M[5, 1, 1, 5] -> (-32*g^4*s*t)/u^2 - (16*g^2*t*y^2)/u - (2*t*y^4)/s,
 M[5, 1, 2, 3] -> -1/2*(a1^2*y^2)/u, M[5, 1, 3, 2] -> -1/2*(a1^2*y^2)/t, 
 M[5, 1, 3, 7] -> (-4*g^2*(s^3 - 4*t^3 - 4*t^2*u - u^3 + s^2*(8*t + u) - s*(t^2 - 7*t*u + u^2))*y^2)/(s*t*u),
 M[5, 1, 5, 1] -> (-32*g^4*s*u)/t^2 - (16*g^2*u*y^2)/t - (2*u*y^4)/s, 
 M[5, 1, 7, 3] -> (-4*g^2*(s^3 - t^3 - 4*t*u^2 - 4*u^3 + s^2*(t + 8*u) - s*(t^2 - 7*t*u + u^2))*y^2)/(s*t*u),
 M[5, 2, 0, 3] -> -1/2*(a1^2*y^2)/u,
 M[5, 2, 3, 0] -> -1/2*(a1^2*y^2)/t, 
 M[5, 3, 3, 5] -> (32*g^4*(s^2 + t^2))/u^2 + (16*g^2*t*y^2)/u + 4*y^4,
 M[5, 3, 5, 3] -> (32*g^4*(s^2 + u^2))/t^2 + (16*g^2*u*y^2)/t + 4*y^4,
 M[5, 4, 0, 2] -> (a1^2*y^2)/(2*s), 
 M[5, 4, 0, 7] -> (4*g^2*(-3*s^3 + t^3 - s^2*u + 2*t^2*u + 2*t*u^2 + u^3 + 2*s*(4*t^2 + 5*t*u + 2*u^2))*y^2)/(s*t*u),
 M[5, 4, 2, 0] -> (a1^2*y^2)/(2*s),
 M[5, 4, 4, 5] -> (32*g^4*(s^2 + t^2))/u^2 + (16*g^2*s*y^2)/u + 4*y^4, 
 M[5, 4, 5, 4] -> (32*g^4*(s^2 + u^2))/t^2 + (16*g^2*s*y^2)/t + 4*y^4,
 M[5, 4, 7, 0] -> (4*g^2*(-3*s^3 - s^2*t + t^3 + 2*t^2*u + 2*t*u^2 + u^3 + 2*s*(2*t^2 + 5*t*u + 4*u^2))*y^2)/(s*t*u), 
 M[5, 5, 5, 5] -> (128*g^4*(t^4 + u^4 + s^2*(t + u)^2))/(t^2*u^2),
 M[5, 6, 0, 1] -> (32*g^4*t*u)/s^2 + (16*g^2*u*y^2)/s + (2*u*y^4)/t,
 M[5, 6, 1, 0] -> (32*g^4*t*u)/s^2 + (16*g^2*t*y^2)/s + (2*t*y^4)/u, 
 M[5, 6, 3, 4] -> (32*g^4*(t^2 + u^2))/s^2 + (16*g^2*t*y^2)/s + 4*y^4,
 M[5, 6, 4, 3] -> (32*g^4*(t^2 + u^2))/s^2 + (16*g^2*u*y^2)/s + 4*y^4,
 M[5, 6, 5, 6] -> (128*g^4*(s^4 + t^4 + s^2*u^2 + 2*s*t*u^2 + t^2*u^2))/(s^2*t^2), 
 M[5, 6, 6, 5] -> (128*g^4*(s^4 + s^2*t^2 + 2*s*t^2*u + t^2*u^2 + u^4))/(s^2*u^2),
 M[5, 6, 7, 7] -> (128*g^4*(-2*s^2 + 7*t^2 + 12*t*u + 7*u^2 + 4*s*(t + u)))/(t*u), 
 M[5, 7, 0, 3] -> (-2*g^2*(s^2*(t + 4*u) + s*(3*t^2 + 17*t*u - 2*u^2) + 2*(t^3 + 8*t^2*u + t*u^2 - 2*u^3))*y^2)/(s*t*u), 
 M[5, 7, 3, 0] -> (-2*g^2*(s^2*(4*t + u) + s*(-2*t^2 + 17*t*u + 3*u^2) + 2*(-2*t^3 + t^2*u + 8*t*u^2 + u^3))*y^2)/(s*t*u),
 M[5, 7, 5, 7] -> (-128*g^4*(6*s^2 - 2*t^2 + 2*t*u + 5*u^2 + 3*s*(t + 3*u)))/(s*u), 
 M[5, 7, 7, 5] -> (-128*g^4*(6*s^2 + 5*t^2 + 2*t*u - 2*u^2 + 3*s*(3*t + u)))/(s*t),
 M[6, 0, 0, 6] -> (-32*g^4*s*t)/u^2 - (16*g^2*t*y^2)/u - (2*t*y^4)/s,
 M[6, 0, 2, 4] -> -1/2*(a1^2*y^2)/u,
 M[6, 0, 4, 2] -> -1/2*(a1^2*y^2)/t, 
 M[6, 0, 4, 7] -> (-4*g^2*(s^3 - 4*t^3 - 4*t^2*u - u^3 + s^2*(8*t + u) - s*(t^2 - 7*t*u + u^2))*y^2)/(s*t*u),
 M[6, 0, 6, 0] -> (-32*g^4*s*u)/t^2 - (16*g^2*u*y^2)/t - (2*u*y^4)/s, 
 M[6, 0, 7, 4] -> (-4*g^2*(s^3 - t^3 - 4*t*u^2 - 4*u^3 + s^2*(t + 8*u) - s*(t^2 - 7*t*u + u^2))*y^2)/(s*t*u),
 M[6, 1, 1, 6] -> (-32*g^4*s*t)/u^2 - (16*g^2*s*y^2)/u - (2*s*y^4)/t, 
 M[6, 1, 6, 1] -> (-32*g^4*s*u)/t^2 - (16*g^2*s*y^2)/t - (2*s*y^4)/u,
 M[6, 2, 1, 4] -> -1/2*(a1^2*y^2)/u,
 M[6, 2, 4, 1] -> -1/2*(a1^2*y^2)/t,
 M[6, 3, 1, 2] -> (a1^2*y^2)/(2*s), 
 M[6, 3, 1, 7] -> (4*g^2*(-3*s^3 + t^3 - s^2*u + 2*t^2*u + 2*t*u^2 + u^3 + 2*s*(4*t^2 + 5*t*u + 2*u^2))*y^2)/(s*t*u),
 M[6, 3, 2, 1] -> (a1^2*y^2)/(2*s), M[6, 3, 3, 6] -> (32*g^4*(s^2 + t^2))/u^2 + (16*g^2*s*y^2)/u + 4*y^4, 
 M[6, 3, 6, 3] -> (32*g^4*(s^2 + u^2))/t^2 + (16*g^2*s*y^2)/t + 4*y^4,
 M[6, 3, 7, 1] -> (4*g^2*(-3*s^3 - s^2*t + t^3 + 2*t^2*u + 2*t*u^2 + u^3 + 2*s*(2*t^2 + 5*t*u + 4*u^2))*y^2)/(s*t*u), 
 M[6, 4, 4, 6] -> (32*g^4*(s^2 + t^2))/u^2 + (16*g^2*t*y^2)/u + 4*y^4,
 M[6, 4, 6, 4] -> (32*g^4*(s^2 + u^2))/t^2 + (16*g^2*u*y^2)/t + 4*y^4,
 M[6, 5, 0, 1] -> (32*g^4*t*u)/s^2 + (16*g^2*t*y^2)/s + (2*t*y^4)/u, 
 M[6, 5, 1, 0] -> (32*g^4*t*u)/s^2 + (16*g^2*u*y^2)/s + (2*u*y^4)/t,
 M[6, 5, 3, 4] -> (32*g^4*(t^2 + u^2))/s^2 + (16*g^2*u*y^2)/s + 4*y^4,
 M[6, 5, 4, 3] -> (32*g^4*(t^2 + u^2))/s^2 + (16*g^2*t*y^2)/s + 4*y^4, 
 M[6, 5, 5, 6] -> (128*g^4*(s^4 + s^2*t^2 + 2*s*t^2*u + t^2*u^2 + u^4))/(s^2*u^2),
 M[6, 5, 6, 5] -> (128*g^4*(s^4 + t^4 + s^2*u^2 + 2*s*t*u^2 + t^2*u^2))/(s^2*t^2), 
 M[6, 5, 7, 7] -> (128*g^4*(-2*s^2 + 7*t^2 + 12*t*u + 7*u^2 + 4*s*(t + u)))/(t*u),
 M[6, 6, 6, 6] -> (128*g^4*(t^4 + u^4 + s^2*(t + u)^2))/(t^2*u^2), 
 M[6, 7, 1, 4] -> (-2*g^2*(s^2*(t + 4*u) + s*(3*t^2 + 17*t*u - 2*u^2) + 2*(t^3 + 8*t^2*u + t*u^2 - 2*u^3))*y^2)/(s*t*u), 
 M[6, 7, 4, 1] -> (-2*g^2*(s^2*(4*t + u) + s*(-2*t^2 + 17*t*u + 3*u^2) + 2*(-2*t^3 + t^2*u + 8*t*u^2 + u^3))*y^2)/(s*t*u),
 M[6, 7, 6, 7] -> (-128*g^4*(6*s^2 - 2*t^2 + 2*t*u + 5*u^2 + 3*s*(t + 3*u)))/(s*u), 
 M[6, 7, 7, 6] -> (-128*g^4*(6*s^2 + 5*t^2 + 2*t*u - 2*u^2 + 3*s*(3*t + u)))/(s*t),
 M[7, 0, 0, 2] -> (a1^2*g^2*(s + t)*(s^2 + s*(-2*t + u) + t*(t + u)))/(4*s^2*t^2),
 M[7, 0, 0, 7] -> (4*g^4*(u*(t + u) + s*(2*t + u)))/(s*t), 
 M[7, 0, 2, 0] -> (a1^2*g^2*(s + u)*(s^2 + s*(t - 2*u) + u*(t + u)))/(4*s^2*u^2),
 M[7, 0, 4, 5] -> 2*g^2*(9 + (8*t)/u + (2*u)/t + (t + u)/s)*y^2,
 M[7, 0, 5, 4] -> 2*g^2*(9 + (2*t)/u + (8*u)/t + (t + u)/s)*y^2, 
 M[7, 0, 7, 0] -> (4*g^4*(t*(t + u) + s*(t + 2*u)))/(s*u),
 M[7, 1, 1, 2] -> (a1^2*g^2*(s + t)*(s^2 + s*(-2*t + u) + t*(t + u)))/(4*s^2*t^2),
 M[7, 1, 1, 7] -> (4*g^4*(u*(t + u) + s*(2*t + u)))/(s*t), 
 M[7, 1, 2, 1] -> (a1^2*g^2*(s + u)*(s^2 + s*(t - 2*u) + u*(t + u)))/(4*s^2*u^2),
 M[7, 1, 3, 6] -> 2*g^2*(9 + (8*t)/u + (2*u)/t + (t + u)/s)*y^2,
 M[7, 1, 6, 3] -> 2*g^2*(9 + (2*t)/u + (8*u)/t + (t + u)/s)*y^2, 
 M[7, 1, 7, 1] -> (4*g^4*(t*(t + u) + s*(t + 2*u)))/(s*u),
 M[7, 2, 0, 1] -> (a1^2*g^2*(t + u)*((t - u)^2 + s*(t + u)))/(4*t^2*u^2),
 M[7, 2, 1, 0] -> (a1^2*g^2*(t + u)*((t - u)^2 + s*(t + u)))/(4*t^2*u^2), 
 M[7, 3, 1, 5] -> g^2*(-18 - (2*s)/t - (16*s)/u - (4*u)/s - (2*u)/t)*y^2,
 M[7, 3, 3, 7] -> (-4*g^4*(5*s^2 + 3*s*(2*t + u) + t*(5*t + 3*u)))/(s*t),
 M[7, 3, 5, 1] -> g^2*(-18 - (16*s)/t - (4*t)/s - (2*s)/u - (2*t)/u)*y^2, 
 M[7, 3, 7, 3] -> (-4*g^4*(5*s^2 + 3*s*(t + 2*u) + u*(3*t + 5*u)))/(s*u),
 M[7, 4, 0, 6] -> g^2*(-18 - (2*s)/t - (16*s)/u - (4*u)/s - (2*u)/t)*y^2,
 M[7, 4, 4, 7] -> (-4*g^4*(5*s^2 + 3*s*(2*t + u) + t*(5*t + 3*u)))/(s*t), 
 M[7, 4, 6, 0] -> g^2*(-18 - (16*s)/t - (4*t)/s - (2*s)/u - (2*t)/u)*y^2,
 M[7, 4, 7, 4] -> (-4*g^4*(5*s^2 + 3*s*(t + 2*u) + u*(3*t + 5*u)))/(s*u),
 M[7, 5, 0, 3] -> g^2*(-18 - (2*s)/t - (4*s)/u - (16*u)/s - (2*u)/t)*y^2, 
 M[7, 5, 3, 0] -> g^2*(-18 - (4*s)/t - (16*t)/s - (2*s)/u - (2*t)/u)*y^2,
 M[7, 5, 5, 7] -> (-64*g^4*(5*s^2 + 3*s*(2*t + u) + t*(5*t + 3*u)))/(s*t),
 M[7, 5, 7, 5] -> (-64*g^4*(5*s^2 + 3*s*(t + 2*u) + u*(3*t + 5*u)))/(s*u), 
 M[7, 6, 1, 4] -> g^2*(-18 - (2*s)/t - (4*s)/u - (16*u)/s - (2*u)/t)*y^2,
 M[7, 6, 4, 1] -> g^2*(-18 - (4*s)/t - (16*t)/s - (2*s)/u - (2*t)/u)*y^2,
 M[7, 6, 6, 7] -> (-64*g^4*(5*s^2 + 3*s*(2*t + u) + t*(5*t + 3*u)))/(s*t), 
 M[7, 6, 7, 6] -> (-64*g^4*(5*s^2 + 3*s*(t + 2*u) + u*(3*t + 5*u)))/(s*u),
 M[7, 7, 0, 1] -> (4*g^4*(s^2 + 2*t*u + s*(t + u)))/(t*u),
 M[7, 7, 1, 0] -> (4*g^4*(s^2 + 2*t*u + s*(t + u)))/(t*u), 
 M[7, 7, 3, 4] -> (4*g^4*(5*t^2 + 6*t*u + 5*u^2 + 3*s*(t + u)))/(t*u),
 M[7, 7, 4, 3] -> (4*g^4*(5*t^2 + 6*t*u + 5*u^2 + 3*s*(t + u)))/(t*u),
 M[7, 7, 5, 6] -> (64*g^4*(5*t^2 + 6*t*u + 5*u^2 + 3*s*(t + u)))/(t*u), 
 M[7, 7, 6, 5] -> (64*g^4*(5*t^2 + 6*t*u + 5*u^2 + 3*s*(t + u)))/(t*u)};


(* ::Subsection:: *)
(*Helper functions*)


feynAssociation=Association[feynResults];


insertCouplings:={g1->g};


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})


fixConvention[arg_]:=symmetriseTU[arg/.insertCouplings/.{s->(-t-u)}]//Expand//Simplify//Expand


removeMissing[arg_]:=arg/.M[__]->0/.Missing["KeyAbsent", _]->0


(* ::Subsection:: *)
(*Test hard*)


totalDRalgo=Sum[M[a,b,c,d],{a,0,4},{b,0,4},{c,0,4},{d,0,4}]/.MatrixElements/.Thread[UserMasses->0]//removeMissing//fixConvention;
totalFeyn=Sum[M[a,b,c,d],{a,0,7},{b,0,7},{c,0,7},{d,0,7}]/.FeynMatrixElements//removeMissing//fixConvention;


totalFeyn=Sum[M[a,b,c,d],{a,0,7},{b,0,7},{c,0,7},{d,0,7}]/.MatrixElementsFeyn//removeMissing//fixConvention;


Collect[totalFeyn,{g,y,a1,a2,b3,b4,lam}];


Collect[s1*totalDRalgo-s2*totalFeyn,{g,y,a1,a2,b3,b4,lam},Simplify[fixConvention[#]]&]/.{s1-s2->0}(*/.{s1->1,s2->1}*)(*//Simplify*)


testList={};


(* everything *)
AppendTo[testList,
TestCreate[
	totalDRalgo,
	totalFeyn
]];


report=TestReport[testList]
report["ResultsDataset"]



