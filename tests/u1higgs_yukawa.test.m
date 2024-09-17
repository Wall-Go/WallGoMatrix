(* ::Package:: *)

Quit[];


(* Check Mathematica version *)
If[$VersionNumber < 13.3,
  Print["The Mathematica testing framework requires Mathematica version ", requiredVersion," or higher. You are using version ", currentVersion, "."];
  Abort[]
];

SetDirectory[NotebookDirectory[]];
(*Put this if you want to create multiple model-files with the same kernel*)
Global`$GroupMathMultipleModels=True;
Global`$LoadGroupMath=True;
(*Needs["DRalgo`","../DRalgo/DRalgo.m"]*)
(*<<DRalgo`*)
Needs["matrixElements`","../src/matrixElements.m"]


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
(* a2/2 \[Phi]^2 \[Chi]^2 *)
InputInv={{1,1,2,2},{True,False,True,True}};
Phi2Chi2=CreateInvariant[Group,RepScalar,InputInv][[1]];
(* quartic terms together *)
VQuartic= lam Phi4 + b4/4 Chi4 + a2/2 Phi2Chi2;
\[Lambda]4=GradQuartic[VQuartic];


(* y (\[Phi]^*)(\[Psi]L (\[Xi]R^*)+ \[Psi]R (\[Xi]L^*)) *)
InputInv={{1,3,2},{False,True,False}};
Yukawa1=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
InputInv={{1,4,1},{False,True,False}}; 
Yukawa2=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
(* Hermitian conjugates *)
InputInv={{1,3,2},{True,False,True}};
Yukawa1HC=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
InputInv={{1,4,1},{True,False,True}}; 
Yukawa2HC=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;


Ysff=y*GradYukawa[Yukawa1+Yukawa2];
YsffC=y*GradYukawa[Yukawa1HC+Yukawa2HC];


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


(* scalars *)
RepPhi=CreateOutOfEq[{1},"S"];
RepChi=CreateOutOfEq[{2},"S"];

(* fermions *)
RepPsi=CreateOutOfEq[{1,2},"F"];
RepXi=CreateOutOfEq[{3,4},"F"];

(* vector bosons *)
RepA=CreateOutOfEq[{1},"V"];


(*
These particles do not necessarily have to be out of equilibrium
the remainin particle content is set as light
*)
ParticleList={RepPhi,RepChi,RepPsi,RepXi,RepA};


(*Defining various masses and couplings*)


VectorMass=Table[mv,{i,1,RepA[[1]]//Length}];
FermionMass=Table[mf,{i,1,Length[gvff[[1]]]}];
ScalarMass=Table[ms,{i,1,Length[gvss[[1]]]}];
ParticleMasses={VectorMass,FermionMass,ScalarMass}
(*
up to the user to make sure that the same order is given in the python code
*)
UserMasses={0};
UserCouplings={g1,b3,b4,a1,a2,lam,y}//Flatten;


(*
	output of matrix elements
*)
OutputFile="matrixElements.u1_higgs_yukawa";
SetDirectory[NotebookDirectory[]];
ParticleName={"Phi","Chi","Psi","Xi","A"};
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



