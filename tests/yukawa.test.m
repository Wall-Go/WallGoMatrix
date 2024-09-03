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
(*Yukawa Model*)


(* ::Section:: *)
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
rep 1-2 are fermions,
(*rep 3 is a scalar*)
*)
(* scalar *)
RepScalar=CreateOutOfEq[{1},"S"];

(* left-handed fermion *)
RepFermionL=CreateOutOfEq[{1},"F"];

(* right-handed fermion *)
RepFermionR=CreateOutOfEq[{2},"F"];

(*Vector bosons*)
RepZ=CreateOutOfEq[{1},"V"];


(*
These particles do not necessarily have to be out of equilibrium
the remainin particle content is set as light
*)
ParticleList={RepScalar,RepFermionL,RepFermionR};


(*Defining various masses and couplings*)


VectorMass=Table[mv,{i,1,RepZ[[1]]//Length}];
FermionMass=Table[mf,{i,1,Length[gvff[[1]]]}];
ScalarMass=Table[ms,{i,1,Length[gvss[[1]]]}];
ParticleMasses={VectorMass,FermionMass,ScalarMass};
(*
up to the user to make sure that the same order is given in the python code
*)
UserMasses={ms,mf,mf};
UserCouplings={CouplingName,\[Lambda],\[Gamma],y}//Flatten;


(*
	output of matrix elements
*)
OutputFile="matrixElements.yukawa";
SetDirectory[NotebookDirectory[]];
ParticleName={"Phi","PsiL","PsiR"};
RepOptional={};
MatrixElements=ExportMatrixElements[OutputFile,ParticleList,UserMasses,UserCouplings,ParticleName,ParticleMasses,RepOptional];


MatrixElements


(* ::Section:: *)
(*Tests*)


FeynMatrixElements={M[0, 0, 0, 0] -> g^4/s^2 + g^4/t^2 + (2*g^4)/(s*t) + g^4/u^2 + 
   (2*g^4)/(s*u) + (2*g^4)/(t*u) + 
   lam*((2*g^2)/s + (2*g^2)/t + (2*g^2)/u + (4*g^3*v)/s^2 + 
     (4*g^3*v)/t^2 + (8*g^3*v)/(s*t) + (4*g^3*v)/u^2 + 
     (8*g^3*v)/(s*u) + (8*g^3*v)/(t*u)) + 
   lam^2*(1 + (4*g*v)/s + (4*g*v)/t + (4*g*v)/u + 
     (6*g^2*v^2)/s^2 + (6*g^2*v^2)/t^2 + (12*g^2*v^2)/(s*t) + 
     (6*g^2*v^2)/u^2 + (12*g^2*v^2)/(s*u) + 
     (12*g^2*v^2)/(t*u)) + lam^3*((2*v^2)/s + (2*v^2)/t + 
     (2*v^2)/u + (4*g*v^3)/s^2 + (4*g*v^3)/t^2 + 
     (8*g*v^3)/(s*t) + (4*g*v^3)/u^2 + (8*g*v^3)/(s*u) + 
     (8*g*v^3)/(t*u)) + lam^4*(v^4/s^2 + v^4/t^2 + 
     (2*v^4)/(s*t) + v^4/u^2 + (2*v^4)/(s*u) + (2*v^4)/(t*u)), 
 M[0, 0, 1, 2] -> ((2*g^2)/s + (4*g*lam*v)/s + (2*lam^2*v^2)/s)*
    y^2 + (-4 + (2*t)/u + (2*u)/t)*y^4, 
 M[0, 0, 2, 1] -> ((2*g^2)/s + (4*g*lam*v)/s + (2*lam^2*v^2)/s)*
    y^2 + (-4 + (2*t)/u + (2*u)/t)*y^4, 
 M[0, 1, 0, 1] -> ((2*g^2)/t + (4*g*lam*v)/t + (2*lam^2*v^2)/t)*
    y^2 + (-4 + (2*s)/u + (2*u)/s)*y^4, 
 M[0, 1, 1, 0] -> ((2*g^2)/u + (4*g*lam*v)/u + (2*lam^2*v^2)/u)*
    y^2 + (-4 + (2*s)/t + (2*t)/s)*y^4, 
 M[0, 2, 0, 2] -> ((2*g^2)/t + (4*g*lam*v)/t + (2*lam^2*v^2)/t)*
    y^2 + (-4 + (2*s)/u + (2*u)/s)*y^4, 
 M[0, 2, 2, 0] -> ((2*g^2)/u + (4*g*lam*v)/u + (2*lam^2*v^2)/u)*
    y^2 + (-4 + (2*s)/t + (2*t)/s)*y^4, 
 M[1, 0, 0, 1] -> ((2*g^2)/u + (4*g*lam*v)/u + (2*lam^2*v^2)/u)*
    y^2 + ((4*s)/t + (4*t)/s - (2*u^2)/(s*t))*y^4, 
 M[1, 0, 1, 0] -> ((2*g^2)/t + (4*g*lam*v)/t + (2*lam^2*v^2)/t)*
    y^2 + ((4*s)/u - (2*t^2)/(s*u) + (4*u)/s)*y^4, 
 M[1, 1, 1, 1] -> (8 + (2*s^2)/(t*u) - (2*t)/u - (2*u)/t)*y^4, 
 M[1, 2, 0, 0] -> ((2*g^2)/s + (4*g*lam*v)/s + (2*lam^2*v^2)/s)*
    y^2 + ((-2*s^2)/(t*u) + (4*t)/u + (4*u)/t)*y^4, 
 M[1, 2, 1, 2] -> (8 - (2*s)/t - (2*t)/s + (2*u^2)/(s*t))*y^4, 
 M[1, 2, 2, 1] -> (8 - (2*s)/u + (2*t^2)/(s*u) - (2*u)/s)*y^4, 
 M[2, 0, 0, 2] -> ((2*g^2)/u + (4*g*lam*v)/u + (2*lam^2*v^2)/u)*
    y^2 + ((4*s)/t + (4*t)/s - (2*u^2)/(s*t))*y^4, 
 M[2, 0, 2, 0] -> ((2*g^2)/t + (4*g*lam*v)/t + (2*lam^2*v^2)/t)*
    y^2 + ((4*s)/u - (2*t^2)/(s*u) + (4*u)/s)*y^4, 
 M[2, 1, 0, 0] -> ((2*g^2)/s + (4*g*lam*v)/s + (2*lam^2*v^2)/s)*
    y^2 + ((-2*s^2)/(t*u) + (4*t)/u + (4*u)/t)*y^4, 
 M[2, 1, 1, 2] -> (8 - (2*s)/u + (2*t^2)/(s*u) - (2*u)/s)*y^4, 
 M[2, 1, 2, 1] -> (8 - (2*s)/t - (2*t)/s + (2*u^2)/(s*t))*y^4, 
 M[2, 2, 2, 2] -> (8 + (2*s^2)/(t*u) - (2*t)/u - (2*u)/t)*y^4};


feynAssociation=Association[feynResults];


insertCouplings={c[1]->lam,c[2]->(g+lam v),c[3]->y};


symmetriseTU[arg_]:=1/2 (arg)+1/2 (arg/.{t->tt}/.{u->t, tt->u})


fixConvention[arg_]:=symmetriseTU[arg/.msq[i_]->0/.{s->(-t-u)}/.insertCouplings]//Expand//Simplify//Expand


removeMissing[arg_]:=arg/.M[__]->0/.Missing["KeyAbsent", _]->0


(* ::Subsection:: *)
(*Test hard*)


testList={};


(* scalar-scalar scattering*)
AppendTo[testList,
TestCreate[M[0,0,0,0]/.MatrixElements//fixConvention//removeMissing,
	M[0,0,0,0]/.FeynMatrixElements//fixConvention//removeMissing
]];


(* scalar to fermions *)
AppendTo[testList,
TestCreate[Sum[M[0,0,c,d],{c,1,2},{d,1,2}]/.MatrixElements//fixConvention//removeMissing,
	Sum[M[0,0,c,d],{c,1,2},{d,1,2}]/.FeynMatrixElements//fixConvention//removeMissing
]];


(* fermions to scalar *)
AppendTo[testList,
TestCreate[Sum[M[c,d,0,0],{c,1,2},{d,1,2}]/.MatrixElements//fixConvention//removeMissing,
	Sum[1/2 M[c,d,0,0],{c,1,2},{d,1,2}]/.FeynMatrixElements//fixConvention//removeMissing (* explicit 1/2 is due to average over leg 1 *)
]];


(* scalar-fermion scattering *)
AppendTo[testList,
TestCreate[Sum[M[0,c,d,0]+M[0,c,0,d],{c,1,2},{d,1,2}]/.MatrixElements//fixConvention//removeMissing,
	Sum[M[0,c,d,0]+M[0,c,0,d],{c,1,2},{d,1,2}]/.FeynMatrixElements//fixConvention//removeMissing
]];


(* fermion-scalar scattering *)
AppendTo[testList,
TestCreate[Sum[M[c,0,d,0]+M[c,0,0,d],{c,1,2},{d,1,2}]/.MatrixElements//fixConvention//removeMissing,
	Sum[1/2 (M[c,0,d,0]+M[c,0,0,d]),{c,1,2},{d,1,2}]/.FeynMatrixElements//fixConvention//removeMissing (* explicit 1/2 is due to average over leg 1 *)
]];


(* fermion-fermion scattering*)
AppendTo[testList,
TestCreate[Sum[M[a,b,c,d],{a,1,2},{b,1,2},{c,1,2},{d,1,2}]/.MatrixElements//fixConvention//removeMissing,
	Sum[1/2 M[a,b,c,d],{a,1,2},{b,1,2},{c,1,2},{d,1,2}]/.FeynMatrixElements//fixConvention//removeMissing (* explicit 1/2 is due to average over leg 1 *)
]];


report=TestReport[testList]
report["ResultsDataset"]







