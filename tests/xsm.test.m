(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
(*Put this if you want to create multiple model-files with the same kernel*)
$GroupMathMultipleModels=True;
$LoadGroupMath=True;
<<../src/WallGoMatrix.m


(* ::Chapter:: *)
(*SM+sr1*)


(*see 1506.04741 [hep-ph]*)


(* ::Section:: *)
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


(* ::Subsection:: *)
(*TopL, TopR*)


vev={0,v,0,0,0};
SymmetryBreaking[vev];


(*left-handed top-quark*)
ReptL=CreateParticle[{{1,1}},"F"];

(*right-handed top-quark*)
ReptR=CreateParticle[{2},"F"];

(*join topL and topR into one rep*)
Rept={Join[ReptL[[1]],ReptR[[1]]],"F"};

(*left-handed bottom-quark*)
RepbL=CreateParticle[{{1,2}},"F"];

(*right-handed bottom-quark*)
RepbR=CreateParticle[{3},"F"];

(*join bottomL and bottomR into one rep*)
Repb={Join[RepbL[[1]],RepbR[[1]]],"F"};

(*scalar reps*)
Reph=CreateParticle[{{1,2}},"S"];
Rep\[Phi]0=CreateParticle[{{1,1}},"S"];
Rep\[Phi]pm=CreateParticle[{{1,3}},"S"];
Rep\[Phi]0={{1},"S"};
Rep\[Phi]p={{4},"S"};
Rep\[Phi]m={{3},"S"};

Reps=CreateParticle[{2},"S"];

(*Vector bosons*)
RepGluon=CreateParticle[{1},"V"];
RepW=CreateParticle[{{2,1}},"V"];
RepB=CreateParticle[{{3,1}},"V"];


VectorMass=Join[
	Table[mg2,{i,1,RepGluon[[1]]//Length}],
	Table[mw2,{i,1,RepW[[1]]//Length}],
	Table[mb2,{i,1,RepB[[1]]//Length}]];
FermionMass=Table[mq2,{i,1,Length[gvff[[1]]]}];
ScalarMass={mG2,mh2,mGm2,mGp2,ms2};
ParticleMasses={VectorMass,FermionMass,ScalarMass};
(*
up to the user to make sure that the same order is given in the python code
*)
UserMasses={mq2,mg2,mw2,mb2,mG2,mh2,mGm2,mGp2,ms2}; 
UserCouplings=Variables@Normal@{Ysff,gvss,gvff,gvvv,\[Lambda]4,\[Lambda]3}//DeleteDuplicates;


ParticleList={
	Rept,Repb,
	RepGluon,RepW,RepB,
	Reph,Rep\[Phi]0,Rep\[Phi]p,Rep\[Phi]m,
	Reps
	};
ParticleName={
	"Top","Bottom",
	"Gluon","W","B",
	"Higgs","Goldstone0","GoldStoneMinus","GoldstonePlus",
	"Singlet"};


(*
	output of matrix elements
*)
OutputFile="matrixElements.xsm";
SetDirectory[NotebookDirectory[]];
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	UserMasses,
	UserCouplings,
	ParticleName,
	ParticleMasses,
	{
		Verbose->True,
		TruncateAtLeadingLog->True,
		Replacements->{gw->0,g1->0},
		NormalizeWithDOF->False,
		Format->{"json","txt"}}];


(* ::Section:: *)
(*Tests try*)


(* ::Subsubsection::Closed:: *)
(*O(g3^4)*)


1/2*M[0,0,2,2]/.MatrixElements/.{mq2->0}//Expand
1/2*M[0,2,0,2]/.MatrixElements/.{mq2->0}//Expand
1/2*(M[0,1,0,1]+M[0,10,0,10])/.MatrixElements/.{mq2->0,yt->0}//Expand


(* ::Subsubsection::Closed:: *)
(*O(g3^2*yt^2)*)


1/2*M[0,0,5,2]/.MatrixElements/.{mq2->0}//Expand
1/2*M[0,0,6,2]/.MatrixElements/.{mq2->0}//Expand


1/2*M[0,1,8,2]/.MatrixElements/.{mq2->0}//Expand


1/2*M[0,2,0,5]/.MatrixElements/.{mq2->0}//Expand
1/2*M[0,2,0,6]/.MatrixElements/.{mq2->0}//Expand


1/2*M[0,2,1,8]/.MatrixElements/.{mq2->0}//Expand


1/2*M[0,8,1,2]/.MatrixElements/.{mq2->0}//Expand


(* ::Subsubsection:: *)
(*O(yt^4)*)


1/2*M[0,0,5,5]/.MatrixElements/.{mq2->0}//Expand
1/2*M[0,0,6,6]/.MatrixElements/.{mq2->0}//Expand


1/2*(M[0,0,7,7]+M[0,0,8,8])/.MatrixElements/.{mq2->0}//Expand


1/2*M[0,0,5,6]/.MatrixElements/.{mq2->0}//Expand


1/2*M[0,1,5,8]/.MatrixElements/.{mq2->0}//Expand
1/2*M[0,1,6,7]/.MatrixElements/.{mq2->0}//Expand


1/2*M[0,5,5,0]/.MatrixElements/.{mq2->0}//Expand
1/2*M[0,5,6,0]/.MatrixElements/.{mq2->0}//Expand
1/2*M[0,6,5,0]/.MatrixElements/.{mq2->0}//Expand
1/2*M[0,6,6,0]/.MatrixElements/.{mq2->0}//Expand


1/2*M[0,8,5,1]/.MatrixElements/.{mq2->0}//Expand
1/2*M[0,7,6,1]/.MatrixElements/.{mq2->0}//Expand


1/2*M[0,7,7,0]/.MatrixElements/.{mq2->0}//Expand


(* ::Section:: *)
(*Tests*)


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})
fixConvention[arg_]:=symmetriseTU[arg/.{s->(-t-u)}]//Expand//Simplify//Expand
removeMissing[arg_]:=arg/.M[__]->0/.Missing["KeyAbsent", _]->0


(* ::Subsection:: *)
(*Test with [Hyperlink[1506.04741,{URL["https://arxiv.org/pdf/1506.04741"], None},Apply[Sequence, {ActiveStyle -> {"HyperlinkActive"}, BaseStyle -> {"Hyperlink"}, HyperlinkAction -> "Recycled"}]]]*)


testList={};


(* ::Subsubsection::Closed:: *)
(*O(g3^4)*)


(*tt->gg*)
AppendTo[testList,
TestCreate[
	1/2*M[0,0,2,2]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	128/3 g3^4(u/t+t/u)//fixConvention//removeMissing,
	TestID->"tt->gg"
]];


(*tg->tg*)
AppendTo[testList,
TestCreate[
	1/2*M[0,2,0,2]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	-(128/3)g3^4 s/u+96*g3^4 (s^2+u^2)/t^2//fixConvention//removeMissing,
	TestID->"tg->tg"
]];


(*tq->tq*)
(*has additional yt^2 g3^2 term than reference*)
AppendTo[testList,
TestCreate[
	1/2*(M[0,1,0,1]+M[0,10,0,10])/.MatrixElements/.{yt->0}/.Thread[UserMasses->0]//fixConvention//removeMissing,
	160*g3^4 (u^2+s^2)/t^2//fixConvention//removeMissing,
	TestID->"tq->tq"
]];


(* ::Subsubsection::Closed:: *)
(*O(g3^2*yt^2)*)


	1/2*M[0,0,5,2]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing
	8*yt^2*g3^2(u/t+t/u)//fixConvention//removeMissing


(*tt->hg*)
AppendTo[testList,
TestCreate[
	1/2*M[0,0,5,2]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	8*yt^2*g3^2(u/t+t/u)//fixConvention//removeMissing,
	TestID->"tt->hg"
]];
(*tt->\[Phi]0 g*)
AppendTo[testList,
TestCreate[
	1/2*M[0,0,6,2]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	8*yt^2*g3^2(u/t+t/u)//fixConvention//removeMissing,
	TestID->"tt->\[Phi]0g"
]];


(*tt->\[Phi]+ g*)
AppendTo[testList,
TestCreate[
	1/2*(M[0,1,2,8]+M[0,1,8,2])/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	8*yt^2*g3^2(u/t+t/u)//fixConvention//removeMissing,
	TestID->"tt->\!\(\*SuperscriptBox[\(\[Phi]\), \(+\)]\)g"
]];


(*tg->th*)
AppendTo[testList,
TestCreate[
	1/2*M[0,2,0,5]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	-8*yt^2*g3^2*s/t//fixConvention//removeMissing,
	TestID->"tg->th"
]];
(*tg->t\[Phi]0*)
AppendTo[testList,
TestCreate[
	1/2*M[0,2,0,6]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	-8*yt^2*g3^2*s/t//fixConvention//removeMissing,
	TestID->"tg->t\[Phi]0"
]];


(*tg->b\[Phi]+*)
AppendTo[testList,
TestCreate[
	1/2*(M[0,2,1,8]+M[0,2,8,1])/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	-8*yt^2*g3^2*s/t//fixConvention//removeMissing,
	TestID->"tg->\!\(\*SuperscriptBox[\(b\[Phi]\), \(+\)]\)"
]];


(*t\[Phi]-->bg*)
AppendTo[testList,
TestCreate[
	1/2*(M[0,8,1,2]+M[0,8,2,1])/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	-8*yt^2*g3^2*s/t//fixConvention//removeMissing,
	TestID->"\!\(\*SuperscriptBox[\(t\[Phi]\), \(-\)]\)->bg"
]];


(* ::Subsubsection:: *)
(*O(yt^4)*)


(*tt->hh*)
AppendTo[testList,
TestCreate[
	1/2*M[0,0,5,5]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	(3 t yt^4)/(2 u)+(3 u yt^4)/(2 t)//fixConvention//removeMissing,
	TestID->"tt->hh"
]];
(*tt->\[Phi]0\[Phi]0*)
AppendTo[testList,
TestCreate[
	1/2*M[0,0,6,6]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	(3 t yt^4)/(2 u)+(3 u yt^4)/(2 t)//fixConvention//removeMissing,
	TestID->"tt->\[Phi]0\[Phi]0"
]];


(*tt->\[Phi]^+\[Phi]^-*)
AppendTo[testList,
TestCreate[
	1/2*(M[0,0,7,7]+M[0,0,8,8])/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	3*yt^4*u/t//fixConvention//removeMissing,
	TestID->"tt->\!\(\*SuperscriptBox[\(\[Phi]\), \(+\)]\)\!\(\*SuperscriptBox[\(\[Phi]\), \(-\)]\)"
]];


(*tt->h\[Phi]0*)
AppendTo[testList,
TestCreate[
	1/2*M[0,0,5,6]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	(3 t yt^4)/(2 u)+(3 u yt^4)/(2 t)//fixConvention//removeMissing,
	TestID->"tt->h\[Phi]0"
]];


(*tt->h\[Phi]+*)
AppendTo[testList,
TestCreate[
	1/2*(M[0,1,5,8]+M[0,1,8,5])/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	3/2*yt^4*u/t//fixConvention//removeMissing,
	TestID->"tt->\!\(\*SuperscriptBox[\(h\[Phi]\), \(+\)]\)"
]];
(*tt->\[Phi]0\[Phi]+*)
AppendTo[testList,
TestCreate[
	1/2*(M[0,1,6,8]+M[0,1,8,6])/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	3/2*yt^4*u/t//fixConvention//removeMissing,
	TestID->"tt->\!\(\*SuperscriptBox[\(\[Phi]0\[Phi]\), \(+\)]\)"
]];


(*th->ht*)
AppendTo[testList,
TestCreate[
	1/2*M[0,5,5,0]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	-3/2*yt^4*s/t//fixConvention//removeMissing,
	TestID->"th->th"
]];
(*th->\[Phi]0t*)
AppendTo[testList,
TestCreate[
	1/2*M[0,5,6,0]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	-3/2*yt^4*s/t//fixConvention//removeMissing,
	TestID->"th->\[Phi]0t"
]];
(*t\[Phi]0->ht*)
AppendTo[testList,
TestCreate[
	1/2*M[0,6,5,0]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	-3/2*yt^4*s/t//fixConvention//removeMissing,
	TestID->"t\[Phi]0->ht"
]];
(*t\[Phi]0->\[Phi]0t*)
AppendTo[testList,
TestCreate[
	1/2*M[0,6,6,0]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	-3/2*yt^4*s/t//fixConvention//removeMissing,
	TestID->"t\[Phi]0->\[Phi]0t"
]];


(*t\[Phi]-->hb*)
AppendTo[testList,
TestCreate[
	1/2*M[0,8,5,1]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	3/2*yt^4*u/t//fixConvention//removeMissing,
	TestID->"\!\(\*SuperscriptBox[\(t\[Phi]\), \(-\)]\)->hb"
]];
(*t\[Phi]-->\[Phi]0b*)
AppendTo[testList,
TestCreate[
	1/2*M[0,8,6,1]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	3/2*yt^4*u/t//fixConvention//removeMissing,
	TestID->"\!\(\*SuperscriptBox[\(t\[Phi]\), \(-\)]\)->\[Phi]0b"
]];


(*t\[Phi]+->\[Phi]+t*)
AppendTo[testList,
TestCreate[
	1/2*M[0,7,7,0]/.MatrixElements/.Thread[UserMasses->0]//fixConvention//removeMissing,
	3*yt^4*u/t//fixConvention//removeMissing,
	TestID->"\!\(\*SuperscriptBox[\(t\[Phi]\), \(+\)]\)->\!\(\*SuperscriptBox[\(\[Phi]\), \(+\)]\)t"
]];


(* ::Subsection:: *)
(*Report*)


report=TestReport[testList]
report["ResultsDataset"]



