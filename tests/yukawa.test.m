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


(* ::Chapter:: *)
(*Yukawa Model*)


(* ::Section::Closed:: *)
(*Model*)


Group={"U1"};
RepAdjoint={0};
RepScalar={{{0},"R"}};
CouplingName={g1};


RepFermion={{{0},"L"},{{0},"R"}};


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


(* \[Sigma] \[Phi] *)
InputInv={{1},{True}};
LinearTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VLinear=0;(* turned off *)\[Sigma] LinearTerm
\[Lambda]1=GradTadpole[VLinear];


(* 1/2m^2\[Phi]^2 *)
InputInv={{1,1},{True,True}}; 
MassTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VMass=msq/2 MassTerm
\[Mu]ij=GradMass[VMass]//SparseArray;


(* 1/6\[Gamma]\[Phi]^3 *)
InputInv={{1,1,1},{True,True,True}};
CubicTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VCubic=\[Gamma]/6 CubicTerm
\[Lambda]3=GradCubic[VCubic];


(* 1/24\[Lambda]\[Phi]^4 *)
InputInv={{1,1,1,1},{True,True,True,True}};
QuarticTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VQuartic=\[Lambda]/24 MassTerm^2
\[Lambda]4=GradQuartic[VQuartic];


(* m(Subscript[\[Psi], R]Subscript[\[Psi], L]+Subscript[\[Psi]^+, L]Subscript[\[Psi], R]^+)*)
InputInv={{2,1},{True,True}}; (*Subscript[\[Psi], R]^+Subscript[\[Psi], L]*)
MassTerm1=CreateInvariantFermion[Group,RepFermion,InputInv][[1]]
InputInv={{1,2},{False,False}};  (*Subscript[\[Psi]^+, L]Subscript[\[Psi], R]*)
MassTerm2=CreateInvariantFermion[Group,RepFermion,InputInv][[1]]


\[Mu]IJ=m\[Psi]*GradMassFermion[MassTerm1];
\[Mu]IJC=m\[Psi]*GradMassFermion[MassTerm2];


(* y \[Phi](Subscript[\[Psi], R]Subscript[\[Psi], L]+Subscript[\[Psi]^+, L]Subscript[\[Psi], R]^+)*)
InputInv={{1,2,1},{True,True,True}};
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
InputInv={{1,1,2},{True,False,False}}; 
YukawaDoublet2=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;


Ysff= y*GradYukawa[YukawaDoublet1];
YsffC=y*GradYukawa[YukawaDoublet2];


ImportModel[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Section:: *)
(*User Input*)


(* ::Subsection:: *)
(*SymmetryBreaking*)


vev={v};
SymmetryBreaking[vev]


(* ::Subsection:: *)
(*UserInput*)


(*
In DRalgo fermions are Weyl.
So to create one Dirac we need
one left-handed and
one right-handed fermion
*)


(*Below
rep 1 is a scalar
rep 2-3 are fermions,
*)
(* scalar *)
RepScalar=CreateParticle[{1},"S",ms,"Phi"];

(* left-handed fermion *)
RepFermionL=CreateParticle[{1},"F",mf,"PsiL"];

(* right-handed fermion *)
RepFermionR=CreateParticle[{2},"F",mf,"PsiR"];

(*Vector bosons*)
RepZ=CreateParticle[{1},"V",mv,"LightParticle"];


(*
These particles do not necessarily have to be out of equilibrium
RepZ can later be assumed to be a light particle
*)
ParticleList={RepScalar,RepFermionL,RepFermionR};
LightParticleList={RepZ};


(*Defining various masses and couplings*)
UserMasses={ms,mf};


(*
	output of matrix elements
*)
OutputFile="output/matrixElements.yukawa";
SetDirectory[NotebookDirectory[]];
MatrixElements=ExportMatrixElements[
	OutputFile,
	ParticleList,
	LightParticleList,
	{
		TruncateAtLeadingLog->False,
		Format->{"json","txt"}}];


MatrixElements


(* ::Section:: *)
(*Tests*)


file=FileNameJoin[{NotebookDirectory[],"yukawa.test.json"}];
{particleNames,parameters,FeynMatrixElements}=ImportMatrixElements[file];


insertCouplings={\[Lambda]->lam,\[Gamma]->(g+lam v),y->y};


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})


fixConvention[arg_]:=symmetriseTU[arg/.Thread[UserMasses->0]/.{s->(-t-u)}/.insertCouplings/.v->0]//Expand//Simplify//Expand


removeMissing[arg_]:=arg/.M[__]->0/.Missing["KeyAbsent", _]->0


(* ::Subsection:: *)
(*Test hard*)


testList={};


(* scalar-scalar scattering*)
AppendTo[testList,
TestCreate[
	M[0,0,0,0]/.MatrixElements//fixConvention//removeMissing,
	M[0,0,0,0]/.FeynMatrixElements//fixConvention//removeMissing
]];


(* scalar to fermions *)
AppendTo[testList,
TestCreate[
	Sum[M[0,0,c,d],{c,1,2},{d,1,2}]/.MatrixElements//fixConvention//removeMissing,
	Sum[M[0,0,c,d],{c,1,2},{d,1,2}]/.FeynMatrixElements//fixConvention//removeMissing
]];


(* fermions to scalar *)
AppendTo[testList,
TestCreate[
	Sum[M[c,d,0,0],{c,1,2},{d,1,2}]/.MatrixElements//fixConvention//removeMissing,
	Sum[1/2 M[c,d,0,0],{c,1,2},{d,1,2}]/.FeynMatrixElements//fixConvention//removeMissing (* explicit 1/2 is due to average over leg 1 *)
]];


(* scalar-fermion scattering *)
AppendTo[testList,
TestCreate[
	Sum[M[0,c,d,0]+M[0,c,0,d],{c,1,2},{d,1,2}]/.MatrixElements//fixConvention//removeMissing,
	Sum[M[0,c,d,0]+M[0,c,0,d],{c,1,2},{d,1,2}]/.FeynMatrixElements//fixConvention//removeMissing
]];


(* fermion-scalar scattering *)
AppendTo[testList,
TestCreate[
	Sum[M[c,0,d,0]+M[c,0,0,d],{c,1,2},{d,1,2}]/.MatrixElements//fixConvention//removeMissing,
	Sum[1/2 (M[c,0,d,0]+M[c,0,0,d]),{c,1,2},{d,1,2}]/.FeynMatrixElements//fixConvention//removeMissing (* explicit 1/2 is due to average over leg 1 *)
]];


(* fermion-fermion scattering*)
AppendTo[testList,
TestCreate[
	Sum[M[a,b,c,d],{a,1,2},{b,1,2},{c,1,2},{d,1,2}]/.MatrixElements//fixConvention//removeMissing,
	Sum[1/2 M[a,b,c,d],{a,1,2},{b,1,2},{c,1,2},{d,1,2}]/.FeynMatrixElements//fixConvention//removeMissing (* explicit 1/2 is due to average over leg 1 *)
]];


report=TestReport[testList]
report["ResultsDataset"]
