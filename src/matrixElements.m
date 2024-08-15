(* ::Package:: *)

BeginPackage["matrixElements`"]


(*List of public functions*)

CreateOutOfEq::usage=
"CreateOutOfEq[{{\!\(\*SubscriptBox[\(r\), \(1\)]\),\!\(\*SubscriptBox[\(m\), \(1\)]\)},...,{\!\(\*SubscriptBox[\(r\), \(n\)]\),\!\(\*SubscriptBox[\(m\), \(n\)]\)}},\"R\"] is a function that groups n particles into one representation for
Fermions (R=F),
Vector Bosons (R=V), and
Scalars (R=S)."
SymmetryBreaking::usage=
"Classify different scalar, fermion, and vector representations into respective particles using VEV-induced masses.
As an output particle i are given as {\!\(\*SubscriptBox[\(r\), \(i\)]\),\!\(\*SubscriptBox[\(m\), \(i\)]\)} where \!\(\*SubscriptBox[\(r\), \(i\)]\) is the label of the representation and \!\(\*SubscriptBox[\(m\), \(i\)]\) is the label of the particle mass in that representation"


(*
	Functions from groupmath are used to create the model.
*)
Get["DRalgo`"];
Print["DRalgo is an independent package"];
Print["Please Cite DRalgo: Comput.Phys.Commun. 288 (2023) 108725 \[Bullet] e-Print: 2205.08815 [hep-th]
"];


(* ::Title:: *)
(*Matrix elements*)


(* ::Section::Closed:: *)
(*Help functions*)


DiagonalTensor2[s_SparseArray,a_Integer,b_Integer] := With[
    {
    s1=Flatten[s,{{a},{b}}]
    },
     Table[i,{i,s1}]//Table[#[[i,i]],{i,1,Length[#]}]&//SparseArray
    ]


RangeToIndices[ListI_]:=Block[{},
(*Creates a list of numbers between i and j when ListI=i;;j*)
	Table[i,{i,ListI[[1]],ListI[[2]]}]
];


RepToIndices[ListI_]:=Block[{},
(*Converts a list of range of numbers to a flat array*)
	Table[RangeToIndices[x],{x,ListI}]//Flatten//Sort
];


CreateOutOfEq[Indices_,Type_]:=Block[{PosScalar,PosVector,PosFermion,temp,Particle},
(*Combines selected representations and obtain their indices*)
	If[Type=="F",
		Particle={};
		PosFermion=PrintFermionRepPositions[];
		Do[
			If[Length[i]>=2,
					AppendTo[Particle,{FermionMassiveReps[[i[[1]]]][[i[[2]]]][[1]]}]
				,
					AppendTo[Particle,{RepToIndices[{PosFermion[[i]]}]}]
			];
		,{i,Indices}];
		Return[{Flatten[Particle],Type}]
	];


	If[Type=="V",
		Particle={};
		PosVector=PrintGaugeRepPositions[];
		Do[
			If[Length[i]>=2,
					AppendTo[Particle,{GaugeMassiveReps[[i[[1]]]][[i[[2]]]][[1]]}]
				,
					AppendTo[Particle,{RepToIndices[{PosVector[[i]]}]}]
			];
		,{i,Indices}];
		Return[{Flatten[Particle],Type}]
	];
	
	If[Type=="S",
		Particle={};
		PosScalar=PrintScalarRepPositions[];
		Do[
			If[Length[i]>=2,
					AppendTo[Particle,{ScalarMassiveReps[[i[[1]]]][[i[[2]]]][[1]]}]
				,
					AppendTo[Particle,{RepToIndices[{PosScalar[[i]]}]}]
			];
		,{i,Indices}];
		Return[{Flatten[Particle],Type}]
	];

];


SymmetryBreaking[vev_] :=Block[{PosVector,PosFermion,PosScalar,count},
(*
	
*)
	PosVector=PrintGaugeRepPositions[];
	
	GaugeMassiveReps=Table[SymmetryBreakingGauge[i,vev],{i,1,Length[PosVector]}];
	
	Do[
		If[!NumericQ[Total[GaugeMassiveReps[[i]][[;;,2]],-1]],
			Print[Style[StringJoin["Gauge rep ",ToString[i]," splits into particles with mass squared:"],Bold]];
			count=0;
			Do[
				count++;
				Print[{i,count},":   ",particle[[2]] ];
				,{particle,GaugeMassiveReps[[i]]}]
		]
	,{i,1,Length[PosVector]}];


	PosFermion=PrintFermionRepPositions[];
	
	FermionMassiveReps=Table[SymmetryBreakingFermion[i,vev],{i,1,Length[PosFermion]}];
	
	Do[
		If[!NumericQ[Total[FermionMassiveReps[[i]][[;;,2]],-1]],
			Print[Style[StringJoin["Fermion rep ",ToString[i]," splits into particles with mass:"],Bold]];
			count=0;
			Do[
				count++;
				Print[{i,count},":   ",particle[[2]] ];
				,{particle,FermionMassiveReps[[i]]}]
		]
	,{i,1,Length[PosFermion]}];
	

	PosScalar=PrintScalarRepPositions[];
	
	ScalarMassiveReps=Table[SymmetryBreakingScalar[i,vev],{i,1,Length[PosScalar]}];
	
	Do[
		If[!NumericQ[Total[ScalarMassiveReps[[i]][[;;,2]],-1]],
			Print[Style[StringJoin["Scalar rep ",ToString[i]," splits into particles with mass squared:"],Bold]];
			count=0;
			Do[
				count++;
				Print[{i,count},":   ",particle[[2]] ];
				,{particle,ScalarMassiveReps[[i]]}]
		]
	,{i,1,Length[PosScalar]}];
	
]


SymmetryBreakingGauge[Indices_,vev_] :=Block[{PosVector,Habij,massV,gaugeInd,posHeavy,posLight,rep,val,val2,pos,pos2},
(*
	Finds the masses for the gauge-representation Indices.
*)
	PosVector=PrintGaugeRepPositions[];
	
(*Gauge bosons*)
	Habij=Contract[gvss,gvss,{{3,5}}]//Transpose[#,{1,3,2,4}]&//SparseArray;
	massV=-Activate@TensorContract[Inactive@TensorProduct[Habij,vev,vev],{{3,5},{4,6}}]//SparseArray;
	gaugeInd=Delete[DiagonalTensor2[massV,1,2]//Normal//ArrayRules,-1]/.(({i1_}->x_)->i1);
	
	posHeavy=Intersection[RangeToIndices[PosVector[[Indices]]],gaugeInd];
	posLight=Complement[RangeToIndices[PosVector[[Indices]]],gaugeInd];
		
	If[Length[DiagonalTensor2[massV,1,2][[posHeavy]]]==0,
			rep={{posLight,0}};
		,
			rep={};
			val=DiagonalTensor2[massV,1,2][[posHeavy]]["NonzeroValues"]//DeleteDuplicates;
			val2=DiagonalTensor2[massV,1,2][[posHeavy]]["NonzeroValues"];
			pos=Table[i,{i,Length[posHeavy]}];
	
			rep={};
	
			Do[
				pos2=Table[posHeavy[[pos[[a]][[1]]]],{a,Position[val2,a]}];
				AppendTo[rep,{pos2,a}];
			,{a,val}];

			AppendTo[rep,{posLight,0}];
			If[posLight=={},rep=Drop[rep,-1]];
	];
	rep
]


SymmetryBreakingScalar[Indices_,vev_] :=Block[{PosScalar,scalarInd,massS,posHeavy,posLight,rep,val,val2,pos,pos2},
(*
	Finds the masses for the fermion-representation Indices.
*)
	PosScalar=PrintScalarRepPositions[];
	
(*Scalars*)
	massS=\[Mu]ij + (vev . \[Lambda]4 . vev/2)//Normal//SparseArray;
	scalarInd=Delete[massS//ArrayRules,-1]/.(({i1_,i2_}->x_)->i1);

	posHeavy=Intersection[RangeToIndices[PosScalar[[Indices]]],scalarInd];
	posLight=Complement[RangeToIndices[PosScalar[[Indices]]],scalarInd];
	
	If[Length[massS[[posHeavy,;;]]]==0,
			rep={{posLight,0}};
		,	
			val=massS[[posHeavy,;;]]["NonzeroValues"]//DeleteDuplicates;
			pos=Table[i,{i,Length[posHeavy]}];
	
			rep={};
	
			Do[
				pos2=Table[posHeavy[[pos[[a]][[1]]]],{a,Position[massS[[posHeavy,;;]]["NonzeroValues"],a]}];
				AppendTo[rep,{pos2,a}];
			,{a,val}];

			AppendTo[rep,{posLight,0}];
			If[posLight=={},rep=Drop[rep,-1]];
	];
	rep
]


SymmetryBreakingFermion[Indices_,vev_] :=Block[{PosFermion,fermionInd,massF,posHeavy,posLight,rep,val,val2,pos,pos2},
(*
	Finds the masses for the fermion-representation Indices.
*)
	PosFermion=PrintFermionRepPositions[];
	
(*Fermions*)
	massF=\[Mu]IJ + (vev . Ysff);
	fermionInd=Delete[massF//ArrayRules,-1]/.(({i1_,i2_}->x_)->i1);
	
	posHeavy=Intersection[RangeToIndices[PosFermion[[Indices]]],fermionInd];
	posLight=Complement[RangeToIndices[PosFermion[[Indices]]],fermionInd];
	
	If[Length[massF[[posHeavy,;;]]]==0,
			rep={{posLight,0}};
		,	
			val=massF[[posHeavy,;;]]["NonzeroValues"]//DeleteDuplicates;
			pos=Table[i,{i,Length[posHeavy]}];
	
			rep={};
	
			Do[
				pos2=Table[posHeavy[[pos[[a]][[1]]]],{a,Position[massF[[posHeavy,;;]]["NonzeroValues"],a]}];
				AppendTo[rep,{pos2,a}];
			,{a,val}];

			AppendTo[rep,{posLight,0}];
			If[posLight=={},rep=Drop[rep,-1]];
	];
	rep
]


Contract[tensor1_,tensor2_,indices_]:=Activate @ TensorContract[
        Inactive[TensorProduct][tensor1,tensor2], indices]


Contract[tensor1_,tensor2_,tensor3_,indices_]:=Activate @ TensorContract[
        Inactive[TensorProduct][tensor1,tensor2,tensor3], indices]


Contract[tensor1_,tensor2_,tensor3_,tensor4_,indices_]:=Activate @ TensorContract[
        Inactive[TensorProduct][tensor1,tensor2,tensor3,tensor4], indices]


(* ::Section:: *)
(*Matrix elements*)


(* ::Subsubsection::Closed:: *)
(*V1V2toV3V4*)


CreateMatrixElementV1V2toV3V4[particle1_,particle2_,particle3_,particle4_,vectorMass_]:=
Block[{s,t,u,gTensor,leadingLog},
(*
In QCD this process depends on up to
two Lorentz structures if the particles are all the same; and
one Lorentz structure if they are all different.
*)
leadingLog=False;
If[
	particle1[[2]]!="V"||
	particle2[[2]]!="V"||
	particle3[[2]]!="V"||
	particle4[[2]]!="V",
	Return[0];
,
(*Coupling constants that we will need*)
	gTensor=Table[gvvv[[;;,Particle1[[1]],Particle2[[1]]]],
		{Particle1,{particle1,particle2,particle3,particle4}},
		{Particle2,{particle1,particle2,particle3,particle4}}];

(*Group invariants that multiply various Lorentz Structures*)
	C1=Table[
		Tr[gTensor[[1,3]][[a]] . Transpose[gTensor[[1,3]][[b]]]]
		Tr[gTensor[[2,4]][[a]] . Transpose[gTensor[[2,4]][[b]]]],
		{a,1,Length[gTensor[[1,3]]]},{b,1,Length[gTensor[[1,3]]]}];
	C2=Table[
		Tr[gTensor[[1,4]][[a]] . Transpose[gTensor[[1,4]][[b]]]]
		Tr[gTensor[[2,3]][[a]] . Transpose[gTensor[[2,3]][[b]]]],
		{a,1,Length[gTensor[[1,3]]]},{b,1,Length[gTensor[[1,3]]]}];
	C3=Table[
		Tr[gTensor[[1,2]][[a]] . Transpose[gTensor[[1,2]][[b]]]]
		Tr[gTensor[[3,4]][[a]] . Transpose[gTensor[[3,4]][[b]]]],
		{a,1,Length[gTensor[[1,3]]]},{b,1,Length[gTensor[[1,3]]]}];
		
(*Vector propagators*)
	vectorPropT=Table[1/(t-i),{i,vectorMass}];
	vectorPropU=Table[1/(u-i),{i,vectorMass}];

(*Since there are two diagrams there can be 3 Lorentz structures after squaring and summing over spins*)
(*They are just hardcoded for now*)
	A1=16(-1/4)(s-u)^2; (*squared t-channel diagram*)
	A2=16(-1/4)(s-t)^2; (*squared u-channel diagram*)
	A3=16*t*u/s^2; (*squared s-channel diagram*)

(*We now generate the matrix element for Q1+Q2->Q3+Q4*)
	Res1=A1*DiagonalMatrix[vectorPropT] . C1 . DiagonalMatrix[vectorPropT];
	Res2=A2*DiagonalMatrix[vectorPropU] . C2 . DiagonalMatrix[vectorPropU];
	Res3=A3*C3;
(*add terms independent on the mandelstam variables*)
	Res4=-16*((C1+C2)*(1-1/4)+C3);
	
(*Factor 4 from anti-particle contributions*)
	Return[-Total[Res1+Res2+If[leadingLog,Res3+Res4,0],-1]]
]
];


(* ::Subsubsection::Closed:: *)
(*Q1Q2toQ3Q4*)


CreateMatrixElementQ1Q2toQ3Q4Pre[particle1_,particle2_,particle3_,particle4_,vectorMass_,scalarMass_]:=Block[{},
	Return[1/2*(
		+CreateMatrixElementQ1Q2toQ3Q4[particle1,particle2,particle3,particle4,vectorMass,scalarMass]
		+CreateMatrixElementQ1Q2toQ3Q4[particle1,particle2,particle4,particle3,vectorMass,scalarMass]
		)];
];


CreateMatrixElementQ1Q2toQ3Q4[particle1_,particle2_,particle3_,particle4_,vectorMass_,scalarMass_]:=
Block[{s,t,u,gTensor, YTensor,YTensorC,scalarPropT,scalarPropU,vectorPropU,vectorPropT,
		C1,C2,C1Y,C2Y,A1,A2,Res1,Res2},
(*
In QCD this process depends on up to
two Lorentz structures if the particles are all the same; and
one Lorentz structure if they are all different.
*)
If[
	particle1[[2]]!="F"||
	particle2[[2]]!="F"||
	particle3[[2]]!="F"||
	particle4[[2]]!="F",
	Return[0];
,
(*Coupling constants that we will need*)
	gTensor=Table[gvff[[;;,Particle1[[1]],Particle2[[1]]]],
		{Particle1,{particle1,particle2,particle3,particle4}},
		{Particle2,{particle1,particle2,particle3,particle4}}];
	
	YTensor=Table[Ysff[[;;,Particle1[[1]],Particle2[[1]]]],
		{Particle1,{particle1,particle2,particle3,particle4}},
		{Particle2,{particle1,particle2,particle3,particle4}}];
		
	YTensorC=Table[YsffC[[;;,Particle1[[1]],Particle2[[1]]]],
		{Particle1,{particle1,particle2,particle3,particle4}},
		{Particle2,{particle1,particle2,particle3,particle4}}];		
				
(*Group invariants that multiply various Lorentz Structures*)
	C1=Table[
		Tr[gTensor[[1,3]][[a]] . gTensor[[3,1]][[b]]]
		Tr[gTensor[[2,4]][[a]] . gTensor[[4,2]][[b]]],
		{a,1,Length[gTensor[[1,3]]]},{b,1,Length[gTensor[[2,4]]]}];
	C2=Table[
		Tr[gTensor[[1,4]][[a]] . gTensor[[4,1]][[b]]]
		Tr[gTensor[[2,3]][[a]] . gTensor[[3,2]][[b]]],
		{a,1,Length[gTensor[[1,3]]]},{b,1,Length[gTensor[[1,3]]]}];

	C1Y=Table[
		Tr[YTensor[[1,3]][[a]] . YTensorC[[3,1]][[b]]]
		Tr[YTensor[[2,4]][[a]] . YTensorC[[4,2]][[b]]],
		{a,1,Length[YTensor[[1,3]]]},{b,1,Length[YTensor[[2,4]]]}];
	C2Y=Table[
		Tr[YTensor[[1,4]][[a]] . YTensorC[[4,1]][[b]]]
		Tr[YTensor[[2,3]][[a]] . YTensorC[[3,2]][[b]]],
		{a,1,Length[YTensor[[1,3]]]},{b,1,Length[YTensor[[1,3]]]}];
		

(*Vector propagators*)
	vectorPropT=Table[1/(t-i),{i,vectorMass}];
	vectorPropU=Table[1/(u-i),{i,vectorMass}];

(*Scalar propagators*)
	scalarPropT=Table[1/(t-i),{i,scalarMass}];
	scalarPropU=Table[1/(u-i),{i,scalarMass}];
	
(*
Since there are two diagrams there can be
3 Lorentz structures after squaring and summing over spins
*)
(*They are just hardcoded for now*)
	A1=2(s^2+u^2); (*squared t-channel diagram*)
	A2=2(s^2+t^2); (*squared u-channel diagram*)

(*We now generate the matrix element for Q1+Q2->Q3+Q4*)
	(*Vector contribution*)
	Res1=A1*Tr[DiagonalMatrix[vectorPropT] . C1 . DiagonalMatrix[vectorPropT]]; 
	Res2=A2*Tr[DiagonalMatrix[vectorPropU] . C2 . DiagonalMatrix[vectorPropU]];
	
	(*Yukawa contribution; the 2* takes into account anti-particles*)
	If[Length[scalarPropU]>0&&Length[scalarPropT]>0,
		Res1+=Tr[DiagonalMatrix[scalarPropT] . C1Y . DiagonalMatrix[scalarPropT]];
		Res2+=Tr[DiagonalMatrix[scalarPropU] . C2Y . DiagonalMatrix[scalarPropU]];
	];
	
	Return[4*(Res1+Res2)] (*factor of 4 from adding anti-quark contributions*)
]
];


(* ::Subsubsection::Closed:: *)
(*Q1V1toQ1V1*)


CreateMatrixElementQ1V1toQ1V1[particle1_,particle2_,particle3_,particle4_,vectorMass_,fermionMass_]:=
Block[{s,t,u,gTensor,leadingLog},
(*
In QCD this process depends on up to
two Lorentz structures if the particles are all the same; and
one Lorentz structure if they are all different.
*)
leadingLog=False;
If[ (
	(particle1[[2]]=="F"&&particle2[[2]]=="V")||
	(particle1[[2]]=="V"&&particle2[[2]]=="F")
	)&&(
	(particle3[[2]]=="F"&&particle4[[2]]=="V")||
	(particle3[[2]]=="V"&&particle4[[2]]=="F")),

	p1=particle1[[1]];
	p2=particle2[[1]];
	p3=particle3[[1]];
	p4=particle4[[1]];
(*Changing the order of the particles so that it is always QV->QV*)	
	If[particle1[[2]]=="V"&&particle2[[2]]=="F",
		temp=p1;
		p1=p2;
		p2=temp;
	];
	If[particle3[[2]]=="V"&&particle4[[2]]=="F",
		temp=p3;
		p3=p4;
		p4=temp;
	];
(*Coupling constants that we will need*)
	gTensor[1,2]=gvff[[p2,p1,;;]];
	gTensor[2,1]=gvff[[p2,;;,p1]];
	
	gTensor[3,4]=gvff[[p4,p3,;;]];
	gTensor[4,3]=gvff[[p4,;;,p3]];
	
	gTensor[1,4]=gvff[[p4,p1,;;]];
	gTensor[4,1]=gvff[[p4,;;,p1]];
	
	gTensor[2,3]=gvff[[p2,p3,;;]];
	gTensor[3,2]=gvff[[p2,;;,p3]];

(*Vector propagators*)
	vectorPropT=Table[1/(t-i),{i,vectorMass}];
	
(*Fermion propagators*)
	fermionPropU=Table[1/(u-i),{i,fermionMass}];
	fermionPropS=Table[1/(s),{i,fermionMass}];
	
(*Group invariants that multiply various Lorentz Structures*)
	C1=Tr[
		Sum[gTensor[2,1][[a]] . gTensor[1,2][[a]],{a,Length[gTensor[1,2]]}] .
		Sum[gTensor[4,3][[b]] . gTensor[3,4][[b]],{b,Length[gTensor[3,4]]}]];
	
	C2=Tr[
		DiagonalMatrix[fermionPropU] . Sum[gTensor[4,1][[a]] . gTensor[1,4][[a]],{a,Length[gTensor[1,4]]}] .
		DiagonalMatrix[fermionPropU] . Sum[gTensor[3,2][[b]] . gTensor[2,3][[b]],{b,Length[gTensor[2,3]]}]];
	
	C3=Sum[vectorPropT[[a]]vectorPropT[[b]]
		Tr[gTensor[4,1][[b]] . gTensor[1,2][[a]] . gTensor[4,3][[b]] . gTensor[2,3][[a]]],
		{a,Length[gTensor[1,2]]},{b,Length[gTensor[1,4]]}];
	C3+=Sum[vectorPropT[[a]]vectorPropT[[b]]
		Tr[gTensor[2,1][[b]] . gTensor[1,4][[a]] . gTensor[3,2][[b]] . gTensor[3,4][[a]]],
		{a,Length[gTensor[1,4]]},{b,Length[gTensor[3,2]]}];
	
	C4=Tr[
		Sum[vectorPropT[[a]]gTensor[2,1][[a]] . gTensor[1,2][[a]],{a,Length[gTensor[1,2]]}] .
		Sum[vectorPropT[[b]]gTensor[4,3][[b]] . gTensor[3,4][[b]],{b,Length[gTensor[3,4]]}]];
	C5=Tr[
		Sum[vectorPropT[[a]]gTensor[4,1][[a]] . gTensor[1,4][[a]],{a,Length[gTensor[1,4]]}] .
		Sum[vectorPropT[[b]]gTensor[3,2][[b]] . gTensor[2,3][[b]],{b,Length[gTensor[2,3]]}]];
	
(*Since there are two diagrams there can be 3 Lorentz structures after squaring and summing over spins*)
(*The interference diagram does however not contribute at leading-log, so I omit it*)
(*They are just hardcoded for now*)
	A1=-4*u/s; (*squared t-channel diagram*)
	A2=-4*s*u; (*squared u-channel diagram*)
	A3=4*(s^2+u^2);  (*squared s-channel diagram*)

(*We now generate the matrix element for Q1+Q2->Q3+Q4*)
	Res1=A1*C1;
	Res2=A2*C2;
	Res3=-(A3*C3-8*(C4*s^2+C5*u^2))//Simplify;

	Return[2(If[leadingLog,Res1,0]+Res2+Res3)] (*Factor of 2 from anti-quark contribution*)
,
	Return[0]
]	
];


(* ::Subsubsection::Closed:: *)
(*Q1Q2toV1V2*)


CreateMatrixElementQ1Q2toV1V2[particle1_,particle2_,particle3_,particle4_,fermionMass_]:=
Block[{s,t,u,gTensor,gTensorT,leadingLog},
(*
In QCD this process depends on up to
two Lorentz structures if the particles are all the same; and
one Lorentz structure if they are all different.
*)
leadingLog=False;
If[
	(particle1[[2]]!="F"||particle2[[2]]!="F"||particle3[[2]]!="V"||particle4[[2]]!="V")&&
	(particle1[[2]]!="V"||particle2[[2]]!="V"||particle3[[2]]!="F"||particle4[[2]]!="F"),
	Return[0];
,
	p1=particle1[[1]];
	p2=particle2[[1]];
	p3=particle3[[1]];
	p4=particle4[[1]];
If[
	particle1[[2]]=="V"&&
	particle2[[2]]=="V"&&
	particle3[[2]]=="F"&&
	particle4[[2]]=="F",
(*Just changing the order of the particles so that it is always QQ->VV*)	
		temp=p1;
		p1=p3;
		p3=temp;
		
		temp=p2;
		p2=p4;
		p4=temp;
	];
(*Coupling constants that we will need*)
	gTensor[1,3]=gvff[[p3,p1,;;]];
	gTensor[3,1]=gvff[[p3,;;,p1]];
	
	gTensor[2,4]=gvff[[p4,p2,;;]];
	gTensor[4,2]=gvff[[p4,;;,p2]];

	gTensor[1,4]=gvff[[p4,p1,;;]];
	gTensor[4,1]=gvff[[p4,;;,p1]];
	
	gTensor[2,3]=gvff[[p3,p2,;;]];
	gTensor[3,2]=gvff[[p3,;;,p2]];

(*Fermion propagators*)
	fermionPropU=Table[1/(u-i),{i,fermionMass}];
	fermionPropT=Table[1/(t-i),{i,fermionMass}];

(*Group invariants that multiply various Lorentz Structures*)
	C1=Tr[
		DiagonalMatrix[fermionPropT] . Sum[gTensor[3,1][[a]] . gTensor[1,3][[a]],{a,Length[gTensor[1,3]]}] .
		DiagonalMatrix[fermionPropT] . Sum[gTensor[4,2][[b]] . gTensor[2,4][[b]],{b,Length[gTensor[2,4]]}]];
	C2=Tr[
		DiagonalMatrix[fermionPropU] . Sum[gTensor[4,1][[a]] . gTensor[1,4][[a]],{a,Length[gTensor[1,4]]}] .
		DiagonalMatrix[fermionPropU] . Sum[gTensor[3,2][[b]] . gTensor[2,3][[b]],{b,Length[gTensor[2,3]]}]];
	C3=Sum[
		Tr[gTensor[4,1][[b]] . gTensor[1,3][[a]] . gTensor[4,2][[b]] . gTensor[2,3][[a]]],
		{a,Length[gTensor[1,3]]},{b,Length[gTensor[4,2]]}];
	C3+=Sum[
		Tr[gTensor[3,1][[b]] . gTensor[1,4][[a]] . gTensor[3,2][[b]] . gTensor[2,4][[a]]],
		{a,Length[gTensor[1,4]]},{b,Length[gTensor[3,2]]}];

(*Multiplying with propagators for each particle*)

(*Since there are two diagrams there can be 3 Lorentz structures after squaring and summing over spins*)
(*The interference diagram does however not contribute at leading-log, so I omit it*)
(*They are hardcoded for now*)
	A1=4*t*u; (*squared t-channel diagram*)
	A2=4*t*u; (*squared u-channel diagram*)
	A3=4*(t^2+u^2)/s^2; (*squared s-channel diagram*)

(*Generate the matrix element for Q1+Q2->Q3+Q4*)
	Res1=A1*C1;
	Res2=A2*C2;
	Res3=A3*C3
		-8*t^2/s^2(t^2 C1/.(#->0&/@fermionMass))
		-8*u^2/s^2(u^2 C2/.(#->0&/@fermionMass))//Simplify;

	Return[2*(Res1+Res2+If[leadingLog,Res3,0])](*Factor of 2 from anti-particles*)
]	
];


(* ::Subsubsection::Closed:: *)
(*S1S2toS3S4*)


CreateMatrixElementS1S2toS3S4[particle1_,particle2_,particle3_,particle4_,vectorMass_,scalarMass_]:=
Block[{s,t,u,gTensor,leadingLog,\[Lambda]3Tensor,vectorPropT,scalarPropT,vectorPropU,scalarPropU,
		CSS,CTT,CUU,CSS\[Lambda],CTT\[Lambda],CUU\[Lambda],ASS,ATT,AUU,ResLL,ResQuartic},
(*
	There is no trilinear scalar interaction taken into account as the symmetric phase is assumed.
	Thus only vector trillinear processes and quartic scalar processes are taken into account.
*)
leadingLog=True;
If[
	particle1[[2]]!="S"||
	particle2[[2]]!="S"||
	particle3[[2]]!="S"||
	particle4[[2]]!="S",
	Return[0];
,
(*Coupling constants that we will need*)
	gTensor=Table[gvss[[;;,Particle1[[1]],Particle2[[1]]]],
		{Particle1,{particle1,particle2,particle3,particle4}},
		{Particle2,{particle1,particle2,particle3,particle4}}];

	\[Lambda]3Tensor=Table[\[Lambda]3[[;;,Particle1[[1]],Particle2[[1]]]],
		{Particle1,{particle1,particle2,particle3,particle4}},
		{Particle2,{particle1,particle2,particle3,particle4}}];
		
(*Group invariants that multiply various Lorentz Structures*)
	CSS=Table[
		Tr[gTensor[[1,2]][[a]] . Transpose[gTensor[[1,2]][[b]]]]
		Tr[gTensor[[3,4]][[a]] . Transpose[gTensor[[3,4]][[b]]]],
		{a,1,Length[gvss]},
		{b,1,Length[gvss]}];
	CTT=Table[
		Tr[gTensor[[1,3]][[a]] . Transpose[gTensor[[1,3]][[b]]]]
		Tr[gTensor[[2,4]][[a]] . Transpose[gTensor[[2,4]][[b]]]],
		{a,1,Length[gvss]},
		{b,1,Length[gvss]}];
	CUU=Table[
		Tr[gTensor[[1,4]][[a]] . Transpose[gTensor[[1,4]][[b]]]]
		Tr[gTensor[[2,3]][[a]] . Transpose[gTensor[[2,3]][[b]]]],
		{a,1,Length[gvss]},
		{b,1,Length[gvss]}];

				
	CSS\[Lambda]=Table[
		Tr[\[Lambda]3Tensor[[1,2]][[a]] . Transpose[\[Lambda]3Tensor[[1,2]][[b]]]]
		Tr[\[Lambda]3Tensor[[3,4]][[a]] . Transpose[\[Lambda]3Tensor[[3,4]][[b]]]],
		{a,1,Length[\[Lambda]3]},
		{b,1,Length[\[Lambda]3]}];
	CTT\[Lambda]=Table[
		Tr[\[Lambda]3Tensor[[1,3]][[a]] . Transpose[\[Lambda]3Tensor[[1,3]][[b]]]]
		Tr[\[Lambda]3Tensor[[2,4]][[a]] . Transpose[\[Lambda]3Tensor[[2,4]][[b]]]],
		{a,1,Length[\[Lambda]3]},
		{b,1,Length[\[Lambda]3]}];
	CUU\[Lambda]=Table[
		Tr[\[Lambda]3Tensor[[1,4]][[a]] . Transpose[\[Lambda]3Tensor[[1,4]][[b]]]]
		Tr[\[Lambda]3Tensor[[2,3]][[a]] . Transpose[\[Lambda]3Tensor[[2,3]][[b]]]],
		{a,1,Length[\[Lambda]3]},
		{b,1,Length[\[Lambda]3]}];		
	
		
(*Vector propagators*)
	vectorPropT=Table[1/(t-i),{i,vectorMass}];
	vectorPropU=Table[1/(u-i),{i,vectorMass}];

(*Scalar propagators*)
	scalarPropT=Table[1/(t-i),{i,scalarMass}];
	scalarPropU=Table[1/(u-i),{i,scalarMass}];
		
(*Lorentz structures*)
	ASS=(t-u)^2/4; (*squared s-channel diagram*)
	ATT=(s-u)^2/4; (*squared t-channel diagram*)
	AUU=(t-s)^2/4; (*squared u-channel diagram*)


(*Leading-log terms*)
	ResLL=(
		+ATT*vectorPropT . CTT . vectorPropT
		+AUU*vectorPropU . CUU . vectorPropU);
		
	ResLL+=(
		+ scalarPropT . CTT\[Lambda] . scalarPropT
		+ scalarPropU . CUU\[Lambda] . scalarPropU);

	ResQuartic=Total[\[Lambda]4[[particle1[[1]],particle2[[1]],particle3[[1]],particle4[[1]]]]^2,-1];

	Return[ResLL]
]
];


(* ::Subsubsection::Closed:: *)
(*S1S2toF1F2*)


CreateMatrixElementS1S2toF1F2[particle1_,particle2_,particle3_,particle4_,fermionMass_]:=
Block[{s,t,u,gTensor,gTensorT,gTensorC,gTensorCT,gTensorVF,gTensorVFC,gTensorVS,
		Res,CTT,CUU,ATT,AUU},
(*
		There is no trilinear scalar interaction taken into account as the symmetric phase is assumed.
		The full Yukawa contribution is included, as well as the s-channel exchange of gauge bosons.
		No mixing between gauge and Yukawa terms was added.
*)
If[
	(particle1[[2]]!="S"||particle2[[2]]!="S"||particle3[[2]]!="F"||particle4[[2]]!="F")&&
	(particle1[[2]]!="F"||particle2[[2]]!="F"||particle3[[2]]!="S"||particle4[[2]]!="S"),
	Return[0];
,
	p1=particle1[[1]];
	p2=particle2[[1]];
	p3=particle3[[1]];
	p4=particle4[[1]];
If[
	particle1[[2]]=="F"&&
	particle2[[2]]=="F"&&
	particle3[[2]]=="S"&&
	particle4[[2]]=="S",
(*Just changing the order of the particles so that it is always QQ->VV*)	
		temp=p1;
		p1=p3;
		p3=temp;
		
		temp=p2;
		p2=p4;
		p4=temp;
	];
(*Coupling constants that we will need*)
	gTensor[1,3]=Ysff[[p1,p3,;;]];
	gTensor[2,4]=Ysff[[p2,p4,;;]];
	gTensor[1,4]=Ysff[[p1,p4,;;]];
	gTensor[2,3]=Ysff[[p2,p3,;;]];

	gTensorT[1,3]=Ysff[[p1,;;,p3]];
	gTensorT[2,4]=Ysff[[p2,;;,p4]];
	gTensorT[1,4]=Ysff[[p1,;;,p4]];
	gTensorT[2,3]=Ysff[[p2,;;,p3]];
		
	gTensorC[1,3]=YsffC[[p1,p3,;;]];
	gTensorC[2,4]=YsffC[[p2,p4,;;]];
	gTensorC[1,4]=YsffC[[p1,p4,;;]];
	gTensorC[2,3]=YsffC[[p2,p3,;;]];

	gTensorCT[1,3]=YsffC[[p1,;;,p3]];
	gTensorCT[2,4]=YsffC[[p2,;;,p4]];
	gTensorCT[1,4]=YsffC[[p1,;;,p4]];
	gTensorCT[2,3]=YsffC[[p2,;;,p3]];

	gTensorVS=gvss[[;;,p1,p2]];
	gTensorVF=gvff[[;;,p3,p4]];
	gTensorVFC=gvff[[;;,p4,p3]];

(*Fermion propagators*)
	fermionPropT=Table[1/(t-i),{i,fermionMass}];
	fermionPropU=Table[1/(u-i),{i,fermionMass}];

(*Coupling factors that appear*)
	CTT=Table[
		Tr[gTensor[1,3][[;;,;;,a]] . Transpose[gTensorC[1,3][[;;,;;,b]]]]
		Tr[gTensor[2,4][[;;,;;,b]] . Transpose[gTensorC[2,4][[;;,;;,a]]]],
		{a,1,Length[YsffC[[1,;;]]]},
		{b,1,Length[YsffC[[1,;;]]]}
	];
	
	CUU=Table[
		Tr[gTensor[1,4][[;;,;;,a]] . Transpose[gTensorC[1,4][[;;,;;,b]]]]
		Tr[gTensor[2,3][[;;,;;,b]] . Transpose[gTensorC[2,3][[;;,;;,a]]]],
		{a,1,Length[YsffC[[1,;;]]]},
		{b,1,Length[YsffC[[1,;;]]]}
	];
	
(*Lorentz structures*)
	ATT=t*u;
	AUU=t*u;
	
(*The result*)
	Res=(
		+ATT*fermionPropT . CTT . fermionPropT
		+AUU*fermionPropU . CUU . fermionPropU);
	
(*The full result*)
	Return[2*Res] (*factor of 2 from anti-particles*)
]	
];


(* ::Subsubsection::Closed:: *)
(*Q1S1toQ1S1*)


CreateMatrixElementQ1S1toQ1S1[particle1_,particle2_,particle3_,particle4_,vectorMass_,fermionMass_]:=
Block[{s,t,u,gTensor,leadingLog,gTensorC,gTensorVF,gTensorVFC,gTensorVS,symfac,Res,vectorPropT,fermionPropU,
		CUU,CTT,AUU,ATT,p1,p2,p3,p4,temp},
(*
	
*)
leadingLog=False;
If[ (
	(particle1[[2]]=="F"&&particle2[[2]]=="S")||
	(particle1[[2]]=="S"&&particle2[[2]]=="F")
	)&&(
	(particle3[[2]]=="F"&&particle4[[2]]=="S")||
	(particle3[[2]]=="S"&&particle4[[2]]=="F")),

	p1=particle1[[1]];
	p2=particle2[[1]];
	p3=particle3[[1]];
	p4=particle4[[1]];
	symfac=1;
(*Changing the order of the particles so that it is always QV->QV*)	
	If[particle1[[2]]=="S"&&particle2[[2]]=="F",
		temp=p1;
		p1=p2;
		p2=temp;
		symfac=2; (*accounting for anti-particles*)
	];
	If[particle3[[2]]=="S"&&particle4[[2]]=="F",
		temp=p3;
		p3=p4;
		p4=temp;
	];
	
(*Coupling constants that we will need*)
	gTensor[1,2]=Ysff[[p2,p1,;;]];
	gTensor[3,4]=Ysff[[p4,p3,;;]];

	gTensor[1,4]=Ysff[[p4,p1,;;]];
	gTensor[3,2]=Ysff[[p2,p3,;;]];
	
	gTensorC[1,2]=YsffC[[p2,p1,;;]];
	gTensorC[3,4]=YsffC[[p4,p3,;;]];

	gTensorC[1,4]=YsffC[[p4,p1,;;]];
	gTensorC[3,2]=YsffC[[p2,p3,;;]];
	
	gTensorVS=gvss[[;;,p2,p4]];
	gTensorVF=gvff[[;;,p1,p3]];
	gTensorVFC=gvff[[;;,p3,p1]];
	
(*Vector propagators*)
	vectorPropT=Table[1/(t-i),{i,vectorMass}];
	
(*Fermion propagators*)
	fermionPropU=Table[1/(u-i),{i,fermionMass}];
	
(*Group invariants that multiply various Lorentz Structures*)
	CUU=Table[
		Tr[gTensor[1,4][[;;,;;,a]] . Transpose[gTensorC[1,4][[;;,;;,b]]]]
		Tr[gTensorC[3,2][[;;,;;,a]] . Transpose[gTensor[3,2][[;;,;;,b]]]],
		{a,Length[Ysff[[1]]]},
		{b,Length[Ysff[[1]]]}];

	CTT=Table[
		Tr[gTensorVS[[a]] . Transpose[gTensorVS[[b]]]]
		Tr[gTensorVF[[a]] . gTensorVFC[[b]]],
		{a,Length[gvss]},{b,Length[gvss]}];
		
(*Lorentz structures that appear*)
(*The first three are the fermion channels, the third is the vector channel*)
	AUU=s*u;
	ATT=4*u*s; (*Note to self: double check the sign here*)
	
(*The result*)
	Res=(
		AUU*fermionPropU . CUU . fermionPropU
		+ATT*vectorPropT . CTT . vectorPropT);
	
	Return[2*Res] (*Factor of 2 from anti-quark contribution*)
,
	Return[0]
]	
];


(* ::Subsubsection::Closed:: *)
(*Q1S1toQ1V1*)


SortQ1S1toQ1V1[L_]:=Block[{helpList,ordering},
	helpList=L[[;;,2]];
	ordering={1,2,3,4};
	If[helpList[[1]]=="V"||helpList[[2]]=="V",ordering={3,4,1,2};];
	If[helpList[[ordering[[1]]]]=="S",ordering[[{1,2}]]=ordering[[{2,1}]]];
	If[helpList[[ordering[[3]]]]=="V",ordering[[{3,4}]]=ordering[[{4,3}]]];
	
	Return[L[[ordering]]]
]


CreateMatrixElementQ1S1toQ1V1[particle1_,particle2_,particle3_,particle4__,fermionMass_]:=
Block[{s,t,u,gTensor,leadingLog,gTensorC,gTensorVFC,gTensorVF,gTensorVS,SortHelp},
(*
	Only terms from the fermion propagators are taken into account.
*)
leadingLog=False;
If[ (
	(particle1[[2]]=="F"&&particle2[[2]]=="S")||
	(particle1[[2]]=="S"&&particle2[[2]]=="F")
	)&&(
	(particle3[[2]]=="F"&&particle4[[2]]=="V")||
	(particle3[[2]]=="V"&&particle4[[2]]=="F"))||(
	(particle1[[2]]=="F"&&particle2[[2]]=="V")||
	(particle1[[2]]=="V"&&particle2[[2]]=="F")
	)&&(
	(particle3[[2]]=="F"&&particle4[[2]]=="S")||
	(particle3[[2]]=="S"&&particle4[[2]]=="F")),
	
(*Changing the order of the particles so that it is always QS->QV*)	
	SortHelp=Join[{particle1},{particle2},{particle3},{particle4}]//SortQ1S1toQ1V1[#]&;
	p1=SortHelp[[1]][[1]];
	p2=SortHelp[[2]][[1]];
	p3=SortHelp[[3]][[1]];
	p4=SortHelp[[4]][[1]];

(*Coupling constants that we will need*)
	gTensor[1,2]=Ysff[[p2,p1,;;]];
	gTensorC[1,2]=YsffC[[p2,p1,;;]];

	gTensor[3,2]=Ysff[[p2,p3,;;]];
	gTensorC[3,2]=YsffC[[p2,p3,;;]];
		
	gTensorVF[3,4]=gvff[[p4,p3,;;]];
	gTensorVFC[3,4]=gvff[[p4,;;,p3]];
	
	gTensorVF[1,4]=gvff[[p4,p1,;;]];
	gTensorVFC[1,4]=gvff[[p4,;;,p1]];
	
(*Fermion propagators*)
	fermionPropU=Table[1/(u-i),{i,fermionMass}];
	
(*Group structures that appear*)
	CSS=Sum[
		Tr[gTensor[1,2][[;;,;;,a]] . Transpose[gTensorC[1,2][[;;,;;,b]]]]
		Tr[gTensorVF[3,4][[;;,;;,a]] . Transpose[gTensorVFC[3,4][[;;,b,;;]]]],
		{a,Length[Ysff[[1]]]},{b,Length[Ysff[[1]]]}];
		
	CUU=Table[
		Tr[gTensor[3,2][[;;,;;,a]] . Transpose[gTensorC[3,2][[;;,;;,b]]]]
		Tr[gTensorVF[1,4][[;;,;;,a]] . Transpose[gTensorVFC[1,4][[;;,b,;;]]]],
		{a,Length[Ysff[[1]]]},{b,Length[Ysff[[1]]]}];
	
(*Lorentz structures that appear*)
	ASS=2;
	AUU=2;
	
(*Results*)
	ResSS=CSS*ASS;
	ResUU=AUU*fermionPropU . CUU . fermionPropU;
	Return[ResSS+ResUU]
,
	Return[0]
]	
];


(* ::Section:: *)
(*Getting the matrix elements for out-of-Equilibrium particles*)


degreeOfFreedom[particle_]:=Block[{dof},

	dof=Length[particle[[1]]];
	
	(*factor of 2 from anti-particles*)
	If[particle[[2]]=="F",dof*=2];
	
	(*Factor of 2 from spins*)
	If[particle[[2]]=="V",dof*=2]; 
	
	If[particle[[2]]=="S",dof];
	
	Return[dof];
]


(* ::Subsubsection::Closed:: *)
(*Q1Q2toQ3Q4*)


ExtractOutOfEqElementQ1Q2toQ3Q4[particleList_,LightParticles_,ParticleMasses_]:=
Block[{OutOfEqParticles,MatrixElements,Elements,CollElements,VectorMass,FermionMass,ScalarMass},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into
	Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

	VectorMass=ParticleMasses[[1]];
	FermionMass=ParticleMasses[[2]];
	ScalarMass=ParticleMasses[[3]];

(*Generate all matrix elements*)
	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];

(*Divide the incoming particle by its degree's of freedom*)
	MatrixElements=Table[
		1/degreeOfFreedom[particleList[[a]]]
		CreateMatrixElementQ1Q2toQ3Q4Pre[particleList[[a]],b,c,d,VectorMass,ScalarMass],
		{a,OutOfEqParticles},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
(*This is a list of all non-zero matrix elements*)
	Elements=MatrixElements["NonzeroPositions"];  
(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,

	MatrixElements=Table[Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}]; 
(*We now select all elements where deltaFparticle is amongst the scattered particles*)
(*deltaF is here a list of 1 and 0s*)
(*Gives all light-particle the label specified by the first lightparticle*)
	deltaF=Elements/. x_?NumericQ /;MemberQ[ LightParticles,x ] -> LightParticles[[1]]; 
	
(*Now we create the full list of distinct collision elements*)
	CollElements=Table[{MatrixElements[[i]],deltaF[[i]]},{i,1,Length[Elements]}];

(*We now add all elements with the same deltaF list*)
	symmetries=Gather[CollElements,#1[[2]]==#2[[2]]&];
	CollElements=Table[{symmetries[[i]][[;;,1]]//Total[#,-1]&,symmetries[[i]][[1,2]]},{i,1,Length[symmetries]}];
	
	Return[CollElements]
]
];


(* ::Subsubsection::Closed:: *)
(*Q1V1toQ1V1*)


ExtractOutOfEqElementQ1V1toQ1V1[particleList_,LightParticles_,ParticleMasses_]:=
Block[{OutOfEqParticles,MatrixElements,Elements,CollElements,VectorMass,FermionMass,ScalarMass},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into
	Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

(*particleList is the complete list of particles*)

	VectorMass=ParticleMasses[[1]];
	FermionMass=ParticleMasses[[2]];
	ScalarMass=ParticleMasses[[3]];
(*First we generate all matrix elements*)
	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];

(*Divide the incoming particle by its degree's of freedom*)
	MatrixElements=Table[
		1/degreeOfFreedom[particleList[[a]]]
		CreateMatrixElementQ1V1toQ1V1[particleList[[a]],b,c,d,VectorMass,FermionMass],
		{a,OutOfEqParticles},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
	(*This is a list of all non-zero matrix elements*)
	Elements=MatrixElements["NonzeroPositions"];


(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,
	MatrixElements=Table[Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}];

(*We now select all elements where deltaFparticle is amongst the scattered particles*)
(*deltaF is here a list of 1 and 0s*)
(*Gives all light-particle the label specified by the first lightparticle*)
	deltaF=Elements/. x_?NumericQ /;MemberQ[ LightParticles,x ] -> LightParticles[[1]];
	
(*Now we create the full list of distinct collision elements*)
	CollElements=Table[{MatrixElements[[i]],deltaF[[i]]},{i,1,Length[Elements]}];

(*We now add all elements with the same deltaF list*)
	symmetries=Gather[CollElements,#1[[2]]==#2[[2]]&];
	CollElements=Table[{symmetries[[i]][[;;,1]]//Total[#,-1]&,symmetries[[i]][[1,2]]},{i,1,Length[symmetries]}];
	
	Return[CollElements]
	
]	
];


(* ::Subsubsection::Closed:: *)
(*Q1Q2toV1V2*)


ExtractOutOfEqElementQ1Q2toV1V2[particleList_,LightParticles_,ParticleMasses_]:=
Block[{OutOfEqParticles,MatrixElements,Elements,CollElements,VectorMass,FermionMass,ScalarMass},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into
	Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

(*particleList is the complete list of particles*)
	VectorMass=ParticleMasses[[1]];
	FermionMass=ParticleMasses[[2]];
	ScalarMass=ParticleMasses[[3]];
(*First we generate all matrix elements*)
	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];

(*Divide the incoming particle by its degree's of freedom*)
	MatrixElements=Table[
		1/degreeOfFreedom[particleList[[a]]]
		CreateMatrixElementQ1Q2toV1V2[particleList[[a]],b,c,d,FermionMass],
		{a,OutOfEqParticles},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
	(*This is a list of all non-zero matrix elements*)
	Elements=MatrixElements["NonzeroPositions"];

(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,

	MatrixElements=Table[Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}];

(*We now select all elements where deltaFparticle is amongst the scattered particles*)
(*deltaF is here a list of 1 and 0s*)
(*Gives all light-particle the label specified by the first lightparticle*)
	deltaF=Elements/. x_?NumericQ /;MemberQ[ LightParticles,x ] -> LightParticles[[1]];
	
(*Now we create the full list of distinct collision elements*)
	CollElements=Table[{MatrixElements[[i]],deltaF[[i]]},{i,1,Length[Elements]}];

(*We now add all elements with the same deltaF list*)
	symmetries=Gather[CollElements,#1[[2]]==#2[[2]]&];
	CollElements=Table[{symmetries[[i]][[;;,1]]//Total[#,-1]&,symmetries[[i]][[1,2]]},{i,1,Length[symmetries]}];
	
	Return[CollElements]	
]	
];


(* ::Subsubsection::Closed:: *)
(*V1V2toV3V4*)


ExtractOutOfEqElementV1V2toV3V4[particleList_,LightParticles_,ParticleMasses_]:=
Block[{OutOfEqParticles,MatrixElements,Elements,CollElements,VectorMass,FermionMass,ScalarMass},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into Sum_deltaFparticle
	\[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

(*particleList is the complete list of particles*)
	VectorMass=ParticleMasses[[1]];
	FermionMass=ParticleMasses[[2]];
	ScalarMass=ParticleMasses[[3]];
(*Generate all matrix elements*)
	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];
	
(*Divide the incoming particle by its degree's of freedom*)
	MatrixElements=Table[
		1/degreeOfFreedom[particleList[[a]]]
		CreateMatrixElementV1V2toV3V4[particleList[[a]],b,c,d,VectorMass],
		{a,OutOfEqParticles},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
	(*This is a list of all non-zero matrix elements*)
	Elements=MatrixElements["NonzeroPositions"];

(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,

	MatrixElements=Table[Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}];
(*We now select all elements where deltaFparticle is amongst the scattered particles*)
(*deltaF is here a list of 1 and 0s*)
(*Gives all light-particle the label specified by the first lightparticle*)
	deltaF=Elements/. x_?NumericQ /;MemberQ[ LightParticles,x ] -> LightParticles[[1]];
	
(*Now we create the full list of distinct collision elements*)
	CollElements=Table[{MatrixElements[[i]],deltaF[[i]]},{i,1,Length[Elements]}];

(*We now add all elements with the same deltaF list*)
	symmetries=Gather[CollElements,#1[[2]]==#2[[2]]&];
	CollElements=Table[{symmetries[[i]][[;;,1]]//Total[#,-1]&,symmetries[[i]][[1,2]]},{i,1,Length[symmetries]}];
	
	Return[CollElements]	
]	
];


(* ::Subsubsection::Closed:: *)
(*S1S2toS3S4*)


ExtractOutOfEqElementS1S2toS3S4[particleList_,LightParticles_,ParticleMasses_]:=
Block[{OutOfEqParticles,MatrixElements,Elements,CollElements,VectorMass,FermionMass,ScalarMass},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into Sum_deltaFparticle
	\[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

(*particleList is the complete list of particles*)
	VectorMass=ParticleMasses[[1]];
	FermionMass=ParticleMasses[[2]];
	ScalarMass=ParticleMasses[[3]];
(*Generate all matrix elements*)
	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];
	
(*Divide the incoming particle by its degree's of freedom*)
	MatrixElements=Table[
		1/degreeOfFreedom[particleList[[a]]]
		CreateMatrixElementS1S2toS3S4[particleList[[a]],b,c,d,VectorMass,ScalarMass],
		{a,OutOfEqParticles},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
	(*This is a list of all non-zero matrix elements*)
	Elements=MatrixElements["NonzeroPositions"];

(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,

	MatrixElements=Table[Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}];

(*We now select all elements where deltaFparticle is amongst the scattered particles*)
(*deltaF is here a list of 1 and 0s*)
(*Gives all light-particle the label specified by the first lightparticle*)
	deltaF=Elements/. x_?NumericQ /;MemberQ[ LightParticles,x ] -> LightParticles[[1]];
	
(*Now we create the full list of distinct collision elements*)
	CollElements=Table[{MatrixElements[[i]],deltaF[[i]]},{i,1,Length[Elements]}];

(*We now add all elements with the same deltaF list*)
	symmetries=Gather[CollElements,#1[[2]]==#2[[2]]&];
	CollElements=Table[{symmetries[[i]][[;;,1]]//Total[#,-1]&,symmetries[[i]][[1,2]]},{i,1,Length[symmetries]}];
	
	Return[CollElements]	
]	
];


(* ::Subsubsection::Closed:: *)
(*S1S2toF1F2*)


ExtractOutOfEqElementS1S2toF1F2[particleList_,LightParticles_,ParticleMasses_]:=
Block[{OutOfEqParticles,MatrixElements,Elements,CollElements,VectorMass,FermionMass,ScalarMass},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into
	Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

(*particleList is the complete list of particles*)
	VectorMass=ParticleMasses[[1]];
	FermionMass=ParticleMasses[[2]];
	ScalarMass=ParticleMasses[[3]];
(*First we generate all matrix elements*)
	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];

(*Divide the incoming particle by its degree's of freedom*)
	MatrixElements=Table[
		1/degreeOfFreedom[particleList[[a]]]
		CreateMatrixElementS1S2toF1F2[particleList[[a]],b,c,d,FermionMass],
		{a,OutOfEqParticles},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
	(*This is a list of all non-zero matrix elements*)
	Elements=MatrixElements["NonzeroPositions"];

(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,

	MatrixElements=Table[Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}];
(*We now select all elements where deltaFparticle is amongst the scattered particles*)
(*deltaF is here a list of 1 and 0s*)
(*Gives all light-particle the label specified by the first lightparticle*)
	deltaF=Elements/. x_?NumericQ /;MemberQ[ LightParticles,x ] -> LightParticles[[1]];
	
(*Now we create the full list of distinct collision elements*)
	CollElements=Table[{MatrixElements[[i]],deltaF[[i]]},{i,1,Length[Elements]}];

(*We now add all elements with the same deltaF list*)
	symmetries=Gather[CollElements,#1[[2]]==#2[[2]]&];
	CollElements=Table[{symmetries[[i]][[;;,1]]//Total[#,-1]&,symmetries[[i]][[1,2]]},{i,1,Length[symmetries]}];
	
	Return[CollElements]	
]	
];


(* ::Subsubsection::Closed:: *)
(*Q1S1toQ1S1*)


ExtractOutOfEqElementQ1S1toQ1S1[particleList_,LightParticles_,ParticleMasses_]:=
Block[{OutOfEqParticles,MatrixElements,Elements,CollElements,VectorMass,FermionMass,ScalarMass},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into
	Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

(*particleList is the complete list of particles*)
	VectorMass=ParticleMasses[[1]];
	FermionMass=ParticleMasses[[2]];
	ScalarMass=ParticleMasses[[3]];
	
(*First we generate all matrix elements*)
	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];

(*Divide the incoming particle by its degree's of freedom*)
	MatrixElements=Table[
		1/degreeOfFreedom[particleList[[a]]]
		CreateMatrixElementQ1S1toQ1S1[particleList[[a]],b,c,d,VectorMass,FermionMass],
		{a,OutOfEqParticles},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
	(*This is a list of all non-zero matrix elements*)
	Elements=MatrixElements["NonzeroPositions"];


(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,

	MatrixElements=Table[Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}];
(*We now select all elements where deltaFparticle is amongst the scattered particles*)
(*deltaF is here a list of 1 and 0s*)
(*Gives all light-particle the label specified by the first lightparticle*)
	deltaF=Elements/. x_?NumericQ /;MemberQ[ LightParticles,x ] -> LightParticles[[1]];
	
(*Now we create the full list of distinct collision elements*)
	CollElements=Table[{MatrixElements[[i]],deltaF[[i]]},{i,1,Length[Elements]}];

(*We now add all elements with the same deltaF list*)
	symmetries=Gather[CollElements,#1[[2]]==#2[[2]]&];
	CollElements=Table[{symmetries[[i]][[;;,1]]//Total[#,-1]&,symmetries[[i]][[1,2]]},{i,1,Length[symmetries]}];
	
	Return[CollElements]
	
]	
];


(* ::Subsubsection::Closed:: *)
(*Q1S1toQ1V1*)


ExtractOutOfEqElementQ1S1toQ1V1[particleList_,LightParticles_,ParticleMasses_]:=
Block[{OutOfEqParticles,MatrixElements,Elements,CollElements,VectorMass,FermionMass,ScalarMass},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into
	Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

(*particleList is the complete list of particles*)
	VectorMass=ParticleMasses[[1]];
	FermionMass=ParticleMasses[[2]];
	ScalarMass=ParticleMasses[[3]];
	
(*First we generate all matrix elements*)
	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];

(*Divide the incoming particle by its degree's of freedom*)
	MatrixElements=Table[
		1/degreeOfFreedom[particleList[[a]]]
		CreateMatrixElementQ1S1toQ1V1[particleList[[a]],b,c,d,FermionMass],
		{a,OutOfEqParticles},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
	(*This is a list of all non-zero matrix elements*)
	Elements=MatrixElements["NonzeroPositions"];


(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,

	MatrixElements=Table[ Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}];

(*We now select all elements where deltaFparticle is amongst the scattered particles*)
(*deltaF is here a list of 1 and 0s*)
(*Gives all light-particle the label specified by the first lightparticle*)
	deltaF=Elements/. x_?NumericQ /;MemberQ[ LightParticles,x ] -> LightParticles[[1]];
	
(*Now we create the full list of distinct collision elements*)
	CollElements=Table[{MatrixElements[[i]],deltaF[[i]]},{i,1,Length[Elements]}];

(*We now add all elements with the same deltaF list*)
	symmetries=Gather[CollElements,#1[[2]]==#2[[2]]&];
	CollElements=Table[{symmetries[[i]][[;;,1]]//Total[#,-1]&,symmetries[[i]][[1,2]]},{i,1,Length[symmetries]}];
	
	Return[CollElements]
	
]	
];


(* ::Subsubsection::Closed:: *)
(*Extract elements*)


ExtractOutOfEqElement[particleList_,LightParticles_,ParticleMasses_]:=
Block[{CollEllQ1Q2toQ3Q4,CollEllQ1V1toQ1V1,CollEllQ1Q2toV1V2,CollEllV1V2toV3V4,CollEllTotal,CollEllS1S2toS3S4,CollEllS1S2toF1F2,
CollEllQ1S1toQ1S1,CollEllQ1S1toQ1V1},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*Essentially this is generates the elements going into Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]*)

(*particleList is the complete list of particles*)


(*First we extract the result for all subprocesses*)
CollEllQ1Q2toQ3Q4=ExtractOutOfEqElementQ1Q2toQ3Q4[particleList,LightParticles,ParticleMasses];
CollEllQ1V1toQ1V1=ExtractOutOfEqElementQ1V1toQ1V1[particleList,LightParticles,ParticleMasses];
CollEllQ1Q2toV1V2=ExtractOutOfEqElementQ1Q2toV1V2[particleList,LightParticles,ParticleMasses];
CollEllV1V2toV3V4=ExtractOutOfEqElementV1V2toV3V4[particleList,LightParticles,ParticleMasses];
CollEllS1S2toS3S4=ExtractOutOfEqElementS1S2toS3S4[particleList,LightParticles,ParticleMasses];
CollEllS1S2toF1F2=ExtractOutOfEqElementS1S2toF1F2[particleList,LightParticles,ParticleMasses];
CollEllQ1S1toQ1S1=ExtractOutOfEqElementQ1S1toQ1S1[particleList,LightParticles,ParticleMasses];


CollEllTotal=Join[
	CollEllQ1Q2toQ3Q4,
	CollEllQ1V1toQ1V1,
	CollEllQ1Q2toV1V2,
	CollEllV1V2toV3V4,
	CollEllS1S2toS3S4,
	CollEllS1S2toF1F2,
	CollEllQ1S1toQ1S1
	];


Return[CollEllTotal]

];


(* ::Section::Closed:: *)
(*Exporting to C++*)


ExportMatrixElements[file_,particleList_,UserMasses_,UserCouplings_,ParticleName_,ParticleMasses_,RepOptional_:{}]:=
Block[{ExportTXT,ExportH5,
	Cij,ParticleInfo,LightParticles,particleListFull,CouplingInfo,MatrixElements,
	OutOfEqParticles,RepMasses,RepCouplings,C
	},

(*
Extracting the out-of-eq particles
*)
		OutOfEqParticles=Table[i,{i,1,Length[particleList]}];
		particleListFull=particleList; (*This list includes the light particles which are added below*)
		
(*Adding all remaining particles as light*)
	posFermions=Position[particleList, _?(# =="F" &)][[;;,1]];
	NonEqFermions=Table[particleList[[i]][[1]],{i,posFermions}]//Flatten[#]&;
	LightFermions={Complement[RepToIndices[PrintFermionRepPositions[]],NonEqFermions],"F"};
	
	posVectors=Position[particleList, _?(# =="V" &)][[;;,1]];
	NonEqVectors=Table[particleList[[i]][[1]],{i,posVectors}]//Flatten[#]&;
	LightVectors={Complement[RepToIndices[PrintGaugeRepPositions[]],NonEqVectors],"V"};
	
	posScalars=Position[particleList, _?(# =="S" &)][[;;,1]];
	NonEqScalars=Table[particleList[[i]][[1]],{i,posScalars}]//Flatten[#]&;
	LightScalars={Complement[RepToIndices[PrintScalarRepPositions[]],NonEqScalars],"S"};
	
	If[Length[LightVectors[[1]]]>0,AppendTo[particleListFull,LightVectors]];
	If[Length[LightFermions[[1]]]>0,AppendTo[particleListFull,LightFermions]];
	If[Length[LightScalars[[1]]]>0,AppendTo[particleListFull,LightScalars]];

	LightParticles=Table[i,{i,Length[particleList]+1,Length[particleListFull]}];
	
	particleListFull=particleListFull;
	LightParticles=LightParticles;
(*
Replacement rules for converting to the format used by the c++ code
*)	
	RepMasses=Table[UserMasses[[i]]->msq[i-1],{i,1,Length[UserMasses]}];
	RepCouplings=Table[UserCouplings[[i]]->Symbol["c"][i-1],{i,1,Length[UserCouplings]}];
	
	
(*
Loop over all out-of-eq particles and extracting the matrix elements
*)	
	
	MatrixElements=ExtractOutOfEqElement[particleListFull,LightParticles,ParticleMasses];
	
(*
Extract various C^{ij} components
*)	
	Cij=ConstantArray[0,{Length[OutOfEqParticles],Length[OutOfEqParticles]}];

	Do[
		Elem=Extract[MatrixElements,Position[MatrixElements[[;;,2]],{i,___}]];
			Do[
				Cij[[i,j]]=Extract[Elem,#]&/@Union[Position[Elem[[;;,2]],{_,j,__}],Position[Elem[[;;,2]],{_,_,j,_}],Position[Elem[[;;,2]],{_,__,j}]];
				,{j,OutOfEqParticles}],
		{i,OutOfEqParticles}
	];

(*Metadata*)
	ParticleInfo=Table[{ToString[OutOfEqParticles[[i]]-1],ParticleName[[i]]},{i,Length[OutOfEqParticles]}];
	AppendTo[ParticleInfo,{ToString[Length[OutOfEqParticles]],"LightParticle"}];
	
	CouplingInfo=Table[{ToString[UserCouplings[[i]]],ToString[Symbol["c"][i-1]]},{i,1,Length[UserCouplings]}];
	
(*
Rewrite the matrix element to an export format
*)
	MatrixElements=Table[MatrixElemToC@i/.RepCouplings/.RepMasses/.RepOptional,{i,MatrixElements}];
	
	Table[
		Cij[[i,j]]=Table[MatrixElemToC@k/.RepCouplings/.RepMasses,{k,Cij[[i,j]]}];,
		{i,OutOfEqParticles},{j,OutOfEqParticles}];
	
	ExportTXT=MatrixElements;
	
(*Adding metadata*)
	PrependTo[ExportTXT,ParticleInfo];
	PrependTo[ExportTXT,CouplingInfo];	
	
(*
In the text-file all matrix elements are directly listed
*)
	Export[StringJoin[file,".txt"],ExportTXT];
	
(*In the hdf5 file we separate them into Cij components*)
	ExportH5=Reap[Do[
		CijName=StringJoin["MatrixElements",ParticleName[[i]],ParticleName[[j]]];
		Sow[
			writeData=Table[{ToString[FortranForm[a[[1]]]],ToString[FortranForm[a[[2]]]]},{a,Cij[[i,j]]}];
			If[Length[Cij[[i,j]]]==0,writeData=""];
			CijName -> {"Data" -> writeData}
			];
		,{i,OutOfEqParticles},{j,OutOfEqParticles}]];
	
	ExportH5=Flatten[ExportH5[[2]][[1]]];

(*Adding metadata*)
	AppendTo[ExportH5,"ParticleInfo"->{"Data"->ParticleInfo}];
	AppendTo[ExportH5,"CouplingInfo"->{"Data"->CouplingInfo}];
	
(*Writing the final result to file*)
	Export[StringJoin[file,".hdf5"],ExportH5];

	Return[MatrixElements]
];


MatrixElemToC[MatrixElem_]:=Block[{Ind},
	Ind=MatrixElem[[2]];
	
	Return[M[Ind[[1]]-1,Ind[[2]]-1,Ind[[3]]-1,Ind[[4]]-1]->MatrixElem[[1]]]
]


Begin["`Private`"]
End[]


EndPackage[]







