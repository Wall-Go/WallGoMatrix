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


(* ::Title:: *)
(*SU(3) + SU(2) + Higgs Yukawa Model*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3","SU2","U1"}; (* I've added the U1 field because without it the Higgs is only given 2 degrees of freedom *)
CouplingName={g5,g3,gU1};
RepAdjoint={{1,1},{2},0};
Higgs1={{{0,0},{1},0},"C"}; (* fundamental *)
RepScalar={Higgs1};


su2Reps = RepsUpToDimN[SU2,3];
Grid[Prepend[{#,RepName[SU2,#]}&/@ su2Reps,{"Dynkin coefficients","Name"}],
Frame->All,FrameStyle->LightGray]
su3Reps = RepsUpToDimN[SU3,8];
Grid[Prepend[{#,RepName[SU3,#]}&/@ su3Reps,{"Dynkin coefficients","Name"}],
Frame->All,FrameStyle->LightGray]


Rep1={{{1,0},{1},0},"L"};
Rep2={{{1,0},{0},0},"R"};
RepFermion={Rep1,Rep2};


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


InputInv={{1,1},{True,False}}; (*This specifies that we want a \[Phi]^+\[Phi] term*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;


VMass=msq*MassTerm1[[1]];(*This is the \[Phi]^+\[Phi] term written in component form*)


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


QuarticTerm1=MassTerm1[[1]]^2; (*Because MassTerm1=\[Phi]^+\[Phi], we can write (\[Phi]^+\[Phi])^2=MassTerm1^2*)


VQuartic=\[Lambda]*QuarticTerm1;


\[Lambda]4=GradQuartic[VQuartic];


(* y \[Phi]^*(Subscript[\[Psi], L]Subscript[\[Xi], R]^++Subscript[\[Psi], R]Subscript[\[Xi], L]^+)*)
InputInv={{1,1,2},{False,False,True}};
Yukawa1=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;


Ysff=-y*GradYukawa[Yukawa1+Yukawa2];
YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{y>0}]];


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


vev={0,0,0,v};
SymmetryBreaking[vev]


(*Below
rep 1-2 are fermions,
(*rep 3 is a scalar*)
*)
(* scalar *)
RepPhi=CreateParticle[{{1,1},{1,2}},"S",ms2,"Phi"]

(* left-handed fermion *)
RepPsiL=CreateParticle[{{1,1}},"F",mf2,"PsiL"]
RepXiL=CreateParticle[{{1,2}},"F",mf2,"XiL"]
RepPsiR=CreateParticle[{2},"F",mf2,"PsiR"]

(*Vector bosons*)
RepVec2=CreateParticle[{1},"V",mv2,"Vector2"]
RepVec1=CreateParticle[{2},"V",mv2,"Vector1"]
RepB=CreateParticle[{3},"V",mv2,"VectorU1"]


(*
These particles do not necessarily have to be out of equilibrium
*)
ParticleList={RepPhi,RepPsiL,RepPsiR,RepXiL,RepVec1,RepVec2,RepB};
LightParticleList={};


(*Defining various masses and couplings*)
UserMasses={mv2,mf2,ms2};


Clear[ConjExpand]
ConjExpand[e_, ass_:True] :=
  Assuming[ass,
    Refine @ Distribute[Conjugate[e], Plus, Conjugate, Plus]
  ];


Element[{T$1888411,U$1888411,s,t,u,mv2},Reals]


SimplifyConjugates[expr_, assumptions_List] :=
  expr //. Conjugate[x_] :> FullSimplify[Conjugate[x],
     Assumptions -> assumptions];


Prop /: Conjugate[Prop[x__]] := Prop[x]
(*UpValues[Prop] = {};*)
expr=a+ b + c + Conjugate[T$1888411 (s-u) Prop[t,mv2]+(s-t) U$1888411 Prop[u,mv2]]
SimplifyConjugates[%,Thread[{T$1888411,U$1888411,s,t,u,mv2}>0]]
(*%/.Conjugate[x_]:>FullSimplify[Conjugate[x],Thread[{T$1888411,U$1888411,s,t,u,mv2}>0](*,Element[{T$1888411,U$1888411,s,t,u,mv2},Reals]*)]
*)(*Refine[%,Assumptions->{s \[Element] Reals, t \[Element] Reals, u \[Element] Reals, mv2 \[Element] Reals,
   T$1888411 \[Element] Reals, U$1888411 \[Element] Reals}]*)


(*expr+b*expr/.Conjugate[x_]:>ConjExpand[x,{s \[Element] Reals, t \[Element] Reals, u \[Element] Reals, mv2 \[Element] Reals,
   T$1888411 \[Element] Reals, U$1888411 \[Element] Reals}]*)


(*ConjExpand[x + y, (*{x \[Element] Ima, y \[Element] Reals}*)]*)


(*
	output of matrix elements
*)
OutputFile="output/matrixElements.su3su2_higgs_yukawa";
SetDirectory[NotebookDirectory[]];
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{
		TruncateAtLeadingLog->False,
		NormalizeWithDOF->False,
		Verbose->True,
		Format->{"json","txt"}}];


MatrixElements;


(* ::Chapter:: *)
(*Tests*)


(* ::Section:: *)
(*Importing results from WallGo*)


{particles,parameters,MatrixElements}=ImportMatrixElements["output/matrixElements.su3su2_higgs_yukawa.json"];


(* ::Section:: *)
(*Importing results from FeynCalc*)


{particlesFeyn,parametersFeyn,MatrixElementsFeyn}=ImportMatrixElements["testFiles/sum-sun-higgs-yukawa.test.json"];


(* ::Section:: *)
(*Comparison tests*)


(* ::Subsection:: *)
(*Translate input*)


groupFactors={SUNN1->2,SUNN2->3};
customParameters={ms2->mPhi^2,mf2->mPsi^2,g->g,\[Lambda]->lam};
setMassesToZero={Thread[UserMasses->0],Thread[{mChi,mPhi,mPsi}->0]}//Flatten;
comparisonReplacements={
	groupFactors,
	customParameters,
	setMassesToZero
	}//Flatten;


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})

fixConvention[arg_]:=symmetriseTU[
	arg/.comparisonReplacements/.{s->(-t-u)}
	]//Expand//Simplify//Expand

removeMissing[arg_]:=arg/.M[__]->0/.Missing["KeyAbsent", _]->0


generateWallGo[particlesA_,particlesB_,particlesC_,particlesD_]:=Sum[
	M[a,b,c,d],{a,particlesA},{b,particlesB},{c,particlesC},{d,particlesD}
]/.Flatten[particles]/.MatrixElements/.customParameters//removeMissing;

generateFeynCalc[particlesA_,particlesB_,particlesC_,particlesD_]:=Sum[
	M[a,b,c,d],{a,particlesA},{b,particlesB},{c,particlesC},{d,particlesD}
]/.Flatten[particlesFeyn]/.MatrixElementsFeyn/.groupFactors//removeMissing;

testWallGo[particlesA_,particlesB_,particlesC_,particlesD_]:=
	generateWallGo[particlesA,particlesB,particlesC,particlesD]//fixConvention
testFeynCalc[particlesA_,particlesB_,particlesC_,particlesD_]:=
	generateFeynCalc[particlesA,particlesB,particlesC,particlesD]//fixConvention


(* ::Subsection:: *)
(*Test hard*)


particles
particlesFeyn


testList={};


(* ::Subsubsection::Closed:: *)
(*VVtoVV*)


process="SU3: VV->VV"
test["WallGo"][process]=testWallGo[
	{"Vector1"},
	{"Vector1"},
	{"Vector1"},
	{"Vector1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"A1"},
	{"A1"},
	{"A1"},
	{"A1"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* vector-vector scattering versus AMY *)
process="SU3: VV->VV"
test["WallGo"][process]=testWallGo[
	{"Vector1"},
	{"Vector1"},
	{"Vector1"},
	{"Vector1"}
]
test["AMY"][process]=fixConvention[16 dA CA^2 g^4 (3-(s u)/t^2-(s t)/u^2-(t u)/s^2)/.{dA->CA^2-1}/.{CA->2}/.{g->g3}]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["AMY"][process],
		TestID->"WallGo vs AMY:"<>process]];


process="SU5: VV->VV"
test["WallGo"][process]=testWallGo[
	{"Vector2"},
	{"Vector2"},
	{"Vector2"},
	{"Vector2"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"A2"},
	{"A2"},
	{"A2"},
	{"A2"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* vector-vector scattering versus AMY *)
process="SU3: VV->VV"
test["WallGo"][process]=testWallGo[
	{"Vector2"},
	{"Vector2"},
	{"Vector2"},
	{"Vector2"}
]
test["AMY"][process]=fixConvention[16 dA CA^2 g^4 (3-(s u)/t^2-(s t)/u^2-(t u)/s^2)/.{dA->CA^2-1}/.{CA->3}/.{g->g5}]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["AMY"][process],
		TestID->"WallGo vs AMY:"<>process]];


(* ::Subsubsection::Closed:: *)
(*FFtoFF*)


process="F1F1->F1F1"
test["WallGo"][process]=testWallGo[
	{"PsiL","PsiR","XiL"},
	{"PsiL","PsiR","XiL"},
	{"PsiL","PsiR","XiL"},
	{"PsiL","PsiR","XiL"}
](*/.{flag1[1]->2,flag1[2]->0,flag1[3]->0,flag1[4]->0,flag1[5]->2,flag1[6]->2}*)
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar","Chi","Chibar"},
	{"Psi","Psibar","Chi","Chibar"},
	{"Psi","Psibar","Chi","Chibar"},
	{"Psi","Psibar","Chi","Chibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F1F1->F1F1massive"
test["WallGo"][process]=generateWallGo[
	{"PsiL","PsiR","XiL"},
	{"PsiL","PsiR","XiL"},
	{"PsiL","PsiR","XiL"},
	{"PsiL","PsiR","XiL"}
]/.{mv2->0}//FullSimplify//Collect[#,{g5,g3,y},Expand]&//Expand//Simplify//Expand
test["FeynCalc"][process]=generateFeynCalc[
	{"Psi","Psibar","Chi","Chibar"},
	{"Psi","Psibar","Chi","Chibar"},
	{"Psi","Psibar","Chi","Chibar"},
	{"Psi","Psibar","Chi","Chibar"}
]//FullSimplify//Collect[#,{g5,g3,y},Expand]&//Expand//Simplify//Expand
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F1F1->F1F1"
test["WallGo"][process]=testWallGo[
	{"PsiL","XiL"},
	{"PsiL","XiL"},
	{"PsiL","XiL"},
	{"PsiL","XiL"}
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


process="F1F1->F1F1"
test["WallGo"][process]=testWallGo[
	{"PsiR"},
	{"PsiR"},
	{"PsiR"},
	{"PsiR"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Chi","Chibar"},
	{"Chi","Chibar"},
	{"Chi","Chibar"},
	{"Chi","Chibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F1F2->F1F2"
test["WallGo"][process]=testWallGo[
	{"PsiL","XiL"},
	{"PsiR"},
	{"PsiL","XiL"},
	{"PsiR"}
](*/.{flag1[1]->2,flag1[2]->0,flag1[3]->0,flag1[4]->0,flag1[5]->2,flag1[6]->2}*)
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Chi","Chibar"},
	{"Psi","Psibar"},
	{"Chi","Chibar"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F1F1->F2F2"
test["WallGo"][process]=testWallGo[
	{"PsiL","XiL"},
	{"PsiL","XiL"},
	{"PsiR"},
	{"PsiR"}
](*/.{flag1[1]->2,flag1[2]->0,flag1[3]->0,flag1[4]->0,flag1[5]->2,flag1[6]->2}*)
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Psi","Psibar"},
	{"Chi","Chibar"},
	{"Chi","Chibar"}
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
	{"PsiL","XiL"},
	{"PsiL","XiL"},
	{"Vector2"},
	{"Vector2"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Psi","Psibar"},
	{"A2"},
	{"A2"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F2F2->VV"
test["WallGo"][process]=testWallGo[
	{"PsiR"},
	{"PsiR"},
	{"Vector2"},
	{"Vector2"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Chi","Chibar"},
	{"Chi","Chibar"},
	{"A2"},
	{"A2"}
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
	{"PsiL","XiL"},
	{"Vector2"},
	{"PsiL","XiL"},
	{"Vector2"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"A2"},
	{"Psi","Psibar"},
	{"A2"}
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


(* ::Subsubsection::Closed:: *)
(*SStoFF*)


process="SS->F1F1"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"PsiL","XiL"},
	{"PsiL","XiL"}
](*/.{flag1->1,flag2->1,flag3->1,flag4->1}*)
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


process="SS->F1F1massive"
test["WallGo"][process]=generateWallGo[
	{"Phi"},
	{"Phi"},
	{"PsiL","XiL"},
	{"PsiL","XiL"}
]/.{mv2->0}//Simplify//ApartSquareFree
test["FeynCalc"][process]=generateFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Psi","Psibar"},
	{"Psi","Psibar"}
]//Simplify//ApartSquareFree
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="SS->F1F1"
test["WallGo"][process]=testWallGo[
	{"Phi"},
	{"Phi"},
	{"PsiL","XiL"},
	{"PsiR"}
](*/.{flag1->1,flag2->1,flag3->1,flag4->1}*)
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Psi","Psibar"},
	{"Chi","Chibar"}
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
	{"PsiL","XiL"},
	{"Phi"},
	{"PsiL","XiL"},
	{"Phi"}
](*/.{flag1->1,flag2->1,flag3->1,flag4->1}*)
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


(* ::Subsubsection::Closed:: *)
(*FStoFV*)


process="F1S->F1V"
test["WallGo"][process]=testWallGo[
	{"PsiL","XiL"},
	{"Phi"},
	{"PsiL","XiL"},
	{"Vector2"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"Psi","Psibar"},
	{"A2"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F1S->F2V1"
test["WallGo"][process]=testWallGo[
	{"PsiL","XiL"},
	{"Phi"},
	{"PsiR"},
	{"Vector1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"Chi","Chibar"},
	{"A1"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(*process="F1S->F2V1massive"
test["WallGo"][process]=generateWallGo[
	{"PsiL","XiL"},
	{"Phi"},
	{"PsiR"},
	{"Vector1"}
]//Simplify
test["FeynCalc"][process]=generateFeynCalc[
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"Chi","Chibar"},
	{"A1"}
](*//FullSimplify*)
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];*)


process="F1S->F2V2"
test["WallGo"][process]=testWallGo[
	{"PsiL","XiL"},
	{"Phi"},
	{"PsiR"},
	{"Vector2"}
](*/.{flag[1]->1,flag[2]->1,flag[3]->1}*)
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"Chi","Chibar"},
	{"A2"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


(* ::Subsubsection::Closed:: *)
(*FFtoSV*)


process="F1F1->SV1"
test["WallGo"][process]=testWallGo[
	{"PsiL","XiL"},
	{"PsiL","XiL"},
	{"Phi"},
	{"Vector1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"A1"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F1F1->SV2"
test["WallGo"][process]=testWallGo[
	{"PsiL","XiL"},
	{"PsiL","XiL"},
	{"Phi"},
	{"Vector2"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Psi","Psibar"},
	{"Phi","Phibar"},
	{"A2"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F1F2->SV1"
test["WallGo"][process]=testWallGo[
	{"PsiL","XiL"},
	{"PsiR"},
	{"Phi"},
	{"Vector1"}
](*/.{flag[1]->1,flag[2]->1,flag[3]->1}*)
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Chi","Chibar"},
	{"Phi","Phibar"},
	{"A1"}
]
AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["FeynCalc"][process],
		TestID->"WallGo vs FeynCalc: "<>process]];


process="F1F2->SV2"
test["WallGo"][process]=testWallGo[
	{"PsiL","XiL"},
	{"PsiR"},
	{"Phi"},
	{"Vector2"}
](*/.{flag[1]->1,flag[2]->1,flag[3]->1}*)
test["FeynCalc"][process]=testFeynCalc[
	{"Psi","Psibar"},
	{"Chi","Chibar"},
	{"Phi","Phibar"},
	{"A2"}
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
	{"Vector1"},
	{"Vector1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"A1"},
	{"A1"}
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
	{"Vector1"},
	{"Phi"},
	{"Vector1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"A1"},
	{"Phi","Phibar"},
	{"A1"}
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
	{"Vector1"}
]
test["FeynCalc"][process]=testFeynCalc[
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"Phi","Phibar"},
	{"A1"}
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
(*Full test*)


totalWallGo=Sum[M[a,b,c,d],{a,0,6},{b,0,6},{c,0,6},{d,0,6}]/.MatrixElements/.Thread[UserMasses->0]//removeMissing//fixConvention;
totalFeyn=Sum[M[a,b,c,d],{a,0,13},{b,0,13},{c,0,13},{d,0,13}]/.MatrixElementsFeyn//removeMissing//fixConvention;


Collect[totalFeyn,{g,y,lam}];


Collect[s1*totalWallGo-s2*totalFeyn,{g,y,lam},Simplify[fixConvention[#]]&]/.{s1-s2->0}


testList={};


(* everything *)
AppendTo[testList,
TestCreate[
	totalWallGo,
	totalFeyn
]];


report=TestReport[testList]
report["ResultsDataset"]



