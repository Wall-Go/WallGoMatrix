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
(*SM+sr1*)


(*see 1506.04741 [hep-ph]*)


(* ::Section::Closed:: *)
(*Model*)


HypY={Yl,Ye,Yq,Yu,Yd,Y\[Phi],Y\[Eta]};
repY=Thread[HypY->{-1,-2,1/3,4/3,-(2/3),1,2}];


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
scalar1={{{0,0},{1},Y\[Phi]/2},"C"};
scalar2={{{0,0},{0},0},"R"};
RepScalar={scalar1,scalar2}/.repY;
CouplingName={g3,gw,g1};


Rep1={{{1,0},{1},Yq/2},"L"};
Rep2={{{1,0},{0},Yu/2},"R"};
Rep3={{{1,0},{0},Yd/2},"R"};
Rep4={{{0,0},{1},Yl/2},"L"};
Rep5={{{0,0},{0},Ye/2},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5}/.repY;


(* ::Text:: *)
(*The input for the gauge interactions to DRalgo are then given by*)


RepFermion1Gen={Rep1,Rep2,Rep3,Rep1,Rep2,Rep3,Rep1,Rep2,Rep3,Rep4,Rep5}/.repY;
RepFermion3Gen={RepFermion1Gen}//Flatten[#,1]&;
(*RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;*)


(* ::Text:: *)
(*The first element is the vector self-interaction matrix:*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar]/.repY;


InputInv={{1,1},{True,False}};
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;
InputInv={{2,2},{True,True}};
MassTerm2=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;


VMass=(
	+m1*MassTerm1
	+b2/2*MassTerm2
	);


\[Mu]ij=GradMass[VMass[[1]]]//Simplify//SparseArray;


QuarticTerm1=MassTerm1[[1]]^2;
QuarticTerm2=MassTerm2[[1]]^2;
QuarticTerm3=MassTerm1[[1]]*MassTerm2[[1]];


VQuartic=(
	+lam*QuarticTerm1
	+b4/4*QuarticTerm2
	+a2/2*QuarticTerm3
	);


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{True,False,True}};
CubicTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2,2},{True,True,True}};
CubicTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VCubic=(
	+a1/2*CubicTerm1
	+b3/3*CubicTerm2
	);


\[Lambda]3=GradCubic[VCubic];


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


Ysff=-yt*GradYukawa[YukawaDoublet1[[1]]];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];


ImportModel[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Section:: *)
(*MatrixElements*)


(*
In DRalgo fermions are Weyl.
So to create one Dirac we need
one left-handed and
one right-handed fermion
*)


PrintFieldRepPositions["Scalar"][[1]]


(* ::Subsection:: *)
(*TopL, TopR*)


vev={0,v,0,0,0};
SymmetryBreaking[vev];


(*left-handed top-quark*)
ReptL=CreateParticle[{{1,1}},"F",mq2,"TopL"];

(*right-handed top-quark*)
ReptR=CreateParticle[{2},"F",mq2,"TopR"];

(*join topL and topR into one rep*)
Rept=CreateParticle[{{1,1},2},"F",mq2,"Top"];

(*left-handed bottom-quark*)
RepbL=CreateParticle[{{1,2}},"F",mb2,"BottomL"];

(*right-handed bottom-quark*)
RepbR=CreateParticle[{3},"F",mb2,"BottomR"];

(*join bottomL and bottomR into one rep*)
Repb=CreateParticle[{{1,2},3},"F",mb2,"Bottom"];

(*scalar reps*)
Reph=CreateParticle[{{1,2}},"S",mh2,"Higgs"];
Rep\[Phi]0=CreateParticle[{{1,1}},"S",mG2,"Goldstone"];
Rep\[Phi]0={{4},"S",mG2,"Goldstone0"};
Rep\[Phi]p={{1},"S",mGp2,"GoldstonePlus"};
Rep\[Phi]m={{3},"S",mGm2,"GoldstoneMinus"};

Reps=CreateParticle[{2},"S",ms2,"Singlet"];

(*Vector bosons*)
RepGluon=CreateParticle[{1},"V",mg2,"Gluon"];
RepW=CreateParticle[{{2,1}},"V",mW2,"W"];
RepB=CreateParticle[{{3,1}},"V",mB2,"B"];

(*Light fermions*)
LightFermions=CreateParticle[{4,5,6,7,8,9,10,11},"F",mq2,"LightFermions"];


(*List of user masses*)
UserMasses={mq2,mt2,mb2,mg2,mW2,mB2,mG2,mh2,mGm2,mGp2,ms2};


ParticleList={
	Rept,Repb,
	RepGluon,RepW,RepB,
	Reph,Rep\[Phi]0,Rep\[Phi]p,Rep\[Phi]m,
	Reps
	};
LightParticleList={LightFermions};


ParticleList={{{1,3,5,7,8,9},"F",mq2,"Top"},{{2,4,6,10,11,12},"F",mb2,"Bottom"},{{1,2,3,4,5,6,7,8},"V",mg2,"Gluon"},{{9,10,11},"V",mW2,"W"},{{12},"V",mB2,"B"},{{2},"S",mh2,"Higgs"},{{4},"S",mG2,"Goldstone0"},{{1},"S",mGp2,"GoldstonePlus"},{{3},"S",mGm2,"GoldstoneMinus"},{{5},"S",ms2,"Singlet"}}
LightParticleList={{{13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39},"F",mq2,"LightFermions"}}


(*
	output of matrix elements
*)
OutputFile="output/matrixElements.xsm";
SetDirectory[NotebookDirectory[]];
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{
		Verbose->True,
		TruncateAtLeadingLog->True,
		Replacements->{gw->0,g1->0},
		NormalizeWithDOF->False,
		Format->{"json","txt"}
	}
];


TruncateAtLeadingLogarithm[MatrixElements_]:=Module[{MatrixElementsF,U,S,T},

	MatrixElementsF=MatrixElements/.Flatten[Map[{
			Power[+#[[1]] + msq_, n_?Negative]->#[[2]]^(-n)*Power[+#[[1]] + msq, n],
			Power[-#[[1]] + msq_, n_?Negative]->#[[2]]^(-n)*Power[-#[[1]] + msq, n]
		}&,{{s,S},{t,T},{u,U}}]];
	
	Print[MatrixElementsF];
	MatrixElementsF=Map[{
		Plus@@Table[
		+SeriesCoefficient[#[[1]]/.{T->xLarge*T,U->xLarge*U},{xLarge,Infinity,-i}]
		,
		{i,If[#[[2]][[3]]==#[[2]][[4]],2,1],2}],
		#[[2]]}&,
		MatrixElementsF];
	Print[MatrixElementsF];

	MatrixElementsF=Expand[MatrixElementsF]/.{S*T->0,S*U->0,T*U->0};
	MatrixElementsF=Collect[MatrixElementsF,{S,T,U},Simplify]/.Thread[{T,U,S}->1];
	MatrixElementsF=DeleteCases[MatrixElementsF, {0,{a__}}];
	
	Return[MatrixElementsF];
]


elems={{M[0,1,0,1]/.MatrixElements,{0,1,0,1}}}


TruncateAtLeadingLogarithm[elems]


dropMixedMonomials[expr_] := 
  Collect[expr, {S, T, U},
    Function[term,
      If[Length[Intersection[{S, T, U}, Variables[term]]] > 1,
        0,
        term
      ]
    ]
  ];


T^2*S/. { S^a_. T^b_. :> 0 /; a > 0 && b > 0,
     S^a_. U^b_. :> 0 /; a > 0 && b > 0,
     T^a_. U^b_. :> 0 /; a > 0 && b > 0 }


(* ::Section:: *)
(*Tests try*)


(* ::Subsubsection:: *)
(*O(g3^4)*)


1/2*M[0,0,2,2]/.MatrixElements/.{mq2->0}//Expand
1/2*M[0,2,0,2]/.MatrixElements/.{mq2->0}//Expand
1/2*(M[0,1,0,1]+M[0,10,0,10])/.MatrixElements/.{mq2->0(*,mGm2->0,mGp2->0*)}/.{(*,yt->0*)}//Expand


(* ::Subsubsection::Closed:: *)
(*O(g3^2*yt^2)*)


1/2*M[0,0,5,2]/.MatrixElements/.{mt2->0}//Expand
1/2*M[0,0,6,2]/.MatrixElements/.{mt2->0}//Expand


1/2*(M[0,1,7,2]+M[0,1,8,2])/.MatrixElements/.{mt2->0,mb2->0}//Expand


1/2*M[0,2,0,5]/.MatrixElements/.{mt2->0}//Expand
1/2*M[0,2,0,6]/.MatrixElements/.{mt2->0}//Expand


1/2*(M[0,2,1,7]+M[0,2,1,8])/.MatrixElements/.{mb2->0}//Expand


1/2*(M[0,7,1,2]+M[0,8,1,2])/.MatrixElements/.{mt2->0}//Expand


(* ::Subsubsection::Closed:: *)
(*O(yt^4)*)


1/2*M[0,0,5,5]/.MatrixElements/.{mt2->0}//Expand
1/2*M[0,0,6,6]/.MatrixElements/.{mt2->0}//Expand


1/2*(M[0,0,7,7]+M[0,0,8,8]+M[0,0,7,8]+M[0,0,8,7])/.MatrixElements/.{mb2->0}//Expand


1/2*M[0,0,5,6]/.MatrixElements/.{mt2->0}//Expand


1/2*(M[0,1,5,7]+M[0,1,5,8])/.MatrixElements/.{mt2->0}//Expand
1/2*(M[0,1,6,7]+M[0,1,6,8])/.MatrixElements/.{mt2->0}//Expand


1/2*M[0,5,5,0]/.MatrixElements/.{mt2->0}//Expand
1/2*M[0,5,6,0]/.MatrixElements/.{mt2->0}//Expand
1/2*M[0,6,5,0]/.MatrixElements/.{mt2->0}//Expand
1/2*M[0,6,6,0]/.MatrixElements/.{mt2->0}//Expand


1/2*(M[0,7,5,1]+M[0,8,5,1])/.MatrixElements/.{mt2->0}//Expand
1/2*(M[0,7,6,1]+M[0,8,6,1])/.MatrixElements/.{mt2->0}//Expand


1/2*(M[0,7,7,0]+M[0,8,8,0]+M[0,7,8,0]+M[0,8,7,0])/.MatrixElements/.{mb2->0}//Expand


(* ::Chapter:: *)
(*Tests*)


{particles,parameters,MatrixElements}=ImportMatrixElements["output/matrixElements.xsm.json"];


(* ::Subsection:: *)
(*Translate input*)


insertCouplings={g->g,\[Lambda]->lam,gu1->0,mChi->0,mPhi->0,mPsi->0};
customCouplings={ms2->mPhi^2};


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})


fixConvention[arg_]:=symmetriseTU[
	arg/.customCouplings/.Thread[UserMasses->0]/.{s->(-t-u)}/.insertCouplings
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
(*Test with table 1 of  https://arxiv.org/pdf/1506.04741*)
(*Note that there are no internal Higgs fields propagating*)


particles


testList={};


(* ::Subsubsection:: *)
(*O(g3^4)*)


process="tt->gg"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Top"},
	{"Gluon"},
	{"Gluon"}
]//Expand
test["reference"][process]=128/3 g3^4(u/t+t/u)//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


process="tg->tg"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Gluon"},
	{"Top"},
	{"Gluon"}
]//Expand
test["reference"][process]=-(128/3)g3^4 s/u+96*g3^4 (s^2+u^2)/t^2//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


(*
	has additional yt^2 g3^2 term than reference due to interference BLL
	since there the Higgs is not propagating
*)
process="tq->tq"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Bottom","LightFermions"},
	{"Top"},
	{"Bottom","LightFermions"}
]/.{yt->0}//Expand
test["reference"][process]=160*g3^4 (u^2+s^2)/t^2//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


(* ::Subsubsection::Closed:: *)
(*O(g3^2*yt^2)*)


process="tt->hg"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Top"},
	{"Higgs"},
	{"Gluon"}
]//Expand
test["reference"][process]=8*yt^2*g3^2(u/t+t/u)//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


process="tt->\[Phi]0 g"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Top"},
	{"Goldstone0"},
	{"Gluon"}
]//Expand
test["reference"][process]=8*yt^2*g3^2(u/t+t/u)//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


process="tt->\[Phi]0 g"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Bottom"},
	{"GoldstonePlus","GoldstoneMinus"},
	{"Gluon"}
]//Expand
test["reference"][process]=8*yt^2*g3^2(u/t+t/u)//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


process="tg->th"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Gluon"},
	{"Top"},
	{"Higgs"}
]//Expand
test["reference"][process]=-8*yt^2*g3^2*s/t//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


process="tg->t\[Phi]0"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Gluon"},
	{"Top"},
	{"Goldstone0"}
]//Expand
test["reference"][process]=-8*yt^2*g3^2*s/t//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


process="tg->b\[Phi]+"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Gluon"},
	{"Bottom"},
	{"GoldstonePlus","GoldstoneMinus"}
]//Expand
test["reference"][process]=-8*yt^2*g3^2*s/t//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


process="t\[Phi]-->bg"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"GoldstonePlus","GoldstoneMinus"},
	{"Bottom"},
	{"Gluon"}
]//Expand
test["reference"][process]=-8*yt^2*g3^2*s/t//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


(* ::Subsubsection:: *)
(*O(yt^4)*)


process="tt->hh"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Top"},
	{"Higgs"},
	{"Higgs"}
]//Expand
test["reference"][process]=(3 t yt^4)/(2 u)+(3 u yt^4)/(2 t)//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


process="tt->\[Phi]0\[Phi]0"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Top"},
	{"Goldstone0"},
	{"Goldstone0"}
]//Expand
test["reference"][process]=(3 t yt^4)/(2 u)+(3 u yt^4)/(2 t)//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


process="tt->\[Phi]^+\[Phi]^-"
test["WallGo"][process]=testWallGo[
	{"Top"},
	{"Top"},
	{"GoldstonePlus","GoldStoneMinus"},
	{"GoldstonePlus","GoldStoneMinus"}
]
test["reference"][process]=3*yt^4*u/t//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


process="tt->h\[Phi]0"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Top"},
	{"Higgs"},
	{"Goldstone0"}
]//Expand
test["reference"][process]=(3 t yt^4)/(2 u)+(3 u yt^4)/(2 t)//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


process="tt->h\[Phi]+"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Bottom"},
	{"Higgs","GoldstoneMinus"},
	{"Higgs","GoldstoneMinus"}
]//Expand
test["reference"][process]=3/2*yt^4*u/t//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


process="tt->\[Phi]0\[Phi]+"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Bottom"},
	{"Goldstone0","GoldstoneMinus"},
	{"Goldstone0","GoldstoneMinus"}
]//Expand
test["reference"][process]=3/2*yt^4*u/t//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


process="th->ht"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Higgs"},
	{"Higgs"},
	{"Top"}
]//Expand
test["reference"][process]=-3/2*yt^4*s/t//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


process="th->\[Phi]0t"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Higgs"},
	{"Goldstone0"},
	{"Top"}
]//Expand
test["reference"][process]=-3/2*yt^4*s/t//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


process="t\[Phi]0->ht"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Goldstone0"},
	{"Higgs"},
	{"Top"}
]//Expand
test["reference"][process]=-3/2*yt^4*s/t//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


process="t\[Phi]0->\[Phi]0t"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"Goldstone0"},
	{"Goldstone0"},
	{"Top"}
]//Expand
test["reference"][process]=-3/2*yt^4*s/t//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


(*contains yt^4 contribution in comparision to reference*)
process="t\[Phi]-->hb"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"GoldstoneMinus","GoldstonePlus"},
	{"Higgs"},
	{"Bottom"}
]/.{yt^4/u->yt^4/u,yt^4/t->yt^4/t,yt^4->0}//Expand
test["reference"][process]=3/2*yt^4*u/t//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


(*contains yt^4 contribution in comparision to reference*)
process="t\[Phi]-->\[Phi]0b"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"GoldstoneMinus","GoldstonePlus"},
	{"Goldstone0"},
	{"Bottom"}
]/.{yt^4/u->yt^4/u,yt^4/t->yt^4/t,yt^4->0}//Expand
test["reference"][process]=3/2*yt^4*u/t//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


(*contains yt^4 contribution in comparision to reference*)
process="t\[Phi]+->\[Phi]+t"
test["WallGo"][process]=1/2*testWallGo[
	{"Top"},
	{"GoldstoneMinus","GoldstonePlus"},
	{"GoldstoneMinus","GoldstonePlus"},
	{"Top"}
]/.{yt^4/u->yt^4/u,yt^4/t->yt^4/t,yt^4->0}//Expand
test["reference"][process]=3*yt^4*u/t//fixConvention//removeMissing

AppendTo[testList,
	TestCreate[
		test["WallGo"][process]//Evaluate,
		test["reference"][process],
		TestID->"WallGo vs reference: "<>process]];


(* ::Subsection:: *)
(*Report*)


report=TestReport[testList]
report["ResultsDataset"]

