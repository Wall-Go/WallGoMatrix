(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
(*Put this if you want to create multiple model-files with the same kernel*)
$GroupMathMultipleModels=True;
$LoadGroupMath=True;
<<../DRalgo/DRalgo.m


(* ::Chapter:: *)
(*QCD+W boson*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3","SU2"};
RepAdjoint={{1,1},{2}};
CouplingName={gs,gw};


Rep1={{{1,0},{1}},"L"};
Rep2={{{1,0},{0}},"R"};
Rep3={{{1,0},{0}},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3};


HiggsDoublet={{{0,0},{1}},"C"};
RepScalar={HiggsDoublet};


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


InputInv={{1,1},{True,False}};
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
VMass=m2*MassTerm1;
\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;
QuarticTerm1=MassTerm1^2;
VQuartic=\[Lambda]1H*QuarticTerm1;
\[Lambda]4=GradQuartic[VQuartic];

InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;
Ysff=-GradYukawa[yt*YukawaDoublet[[1]]];
YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];


(* ::Title:: *)
(*Matrix elements*)


(*move this to .m file*)


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


SymmetryBreaking[vev_] :=Block[{PosVector,PosFermion,PosScalar},
(*
	
*)
	PosVector=PrintGaugeRepPositions[];
	
	GaugeMassiveReps=Table[SymmetryBreakingGauge[i,vev],{i,1,Length[PosVector]}];
	
	Do[
		If[!NumericQ[Total[GaugeMassiveReps[[i]][[;;,2]],-1]],
			Print[Style[StringJoin["Gauge rep ",ToString[i]," splits into:"],Bold]];
			Do[
				Print["One particle with mass squared ",particle[[2]] ];
				,{particle,GaugeMassiveReps[[i]]}]
		]
	,{i,1,Length[PosVector]}];


	PosFermion=PrintFermionRepPositions[];
	
	FermionMassiveReps=Table[SymmetryBreakingFermion[i,vev],{i,1,Length[PosFermion]}];
	
	Do[
		If[!NumericQ[Total[FermionMassiveReps[[i]][[;;,2]],-1]],
			Print[Style[StringJoin["Fermion rep ",ToString[i]," splits into:"],Bold]];
			Do[
				Print["One particle with mass  ",particle[[2]] ];
				,{particle,FermionMassiveReps[[i]]}]
		]
	,{i,1,Length[PosFermion]}];
	

	PosScalar=PrintScalarRepPositions[];
	
	ScalarMassiveReps=Table[SymmetryBreakingScalar[i,vev],{i,1,Length[PosScalar]}];
	
	Do[
		If[!NumericQ[Total[ScalarMassiveReps[[i]][[;;,2]],-1]],
			Print[Style[StringJoin["Scalar rep ",ToString[i]," splits into:"],Bold]];
			Do[
				Print["One particle with mass squared ",particle[[2]] ];
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
	gaugeInd=Delete[DiagonalTensor2[massV,1,2]//Normal//ArrayRules,-1]/.(({a_}->x_)->a);
	
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
	
(*Fermions*)
	massS=\[Mu]ij+vev . \[Lambda]4 . vev/2;
	scalarInd=Delete[massS//ArrayRules,-1]/.(({a_,b_}->x_)->a);

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
	massF=\[Mu]IJ+vev . Ysff;
	fermionInd=Delete[massF//ArrayRules,-1]/.(({a_,b_}->x_)->a);
	
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


(* ::Section::Closed:: *)
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


CreateMatrixElementQ1Q2toQ3Q4[particle1_,particle2_,particle3_,particle4_,vectorMass_]:=
Block[{s,t,u,gTensor,leadingLog},
(*
In QCD this process depends on up to
two Lorentz structures if the particles are all the same; and
one Lorentz structure if they are all different.
*)
leadingLog=False;
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
		Tr[gTensor[[3,1]][[a]] . gTensor[[1,4]][[b]] . gTensor[[4,2]][[a]] . gTensor[[2,3]][[b]]],
		{a,1,Length[gTensor[[1,3]]]},{b,1,Length[gTensor[[1,3]]]}];
	C3+=Table[
		Tr[gTensor[[1,3]][[a]] . gTensor[[3,2]][[b]] . gTensor[[2,4]][[a]] . gTensor[[4,1]][[b]]],
		{a,1,Length[gTensor[[1,3]]]},{b,1,Length[gTensor[[1,3]]]}];
	
	C4=Table[
		Tr[gTensor[[1,2]][[a]] . Transpose[gTensor[[1,2]][[b]]]]
		Tr[gTensor[[3,4]][[a]] . Transpose[gTensor[[3,4]][[b]]]],
		{a,1,Length[gTensor[[1,2]]]},{b,1,Length[gTensor[[1,2]]]}];
	
	C5=Table[
		Tr[gTensor[[1,3]][[a]] . gTensor[[3,4]][[b]] . gTensor[[4,2]][[a]] . gTensor[[2,1]][[b]]],
		{a,1,Length[gTensor[[1,3]]]},{b,1,Length[gTensor[[1,3]]]}];
	C5+=Table[
		Tr[gTensor[[3,1]][[a]] . gTensor[[1,2]][[b]] . gTensor[[2,4]][[a]] . gTensor[[4,3]][[b]]],
		{a,1,Length[gTensor[[1,3]]]},{b,1,Length[gTensor[[1,3]]]}];	
	
(*Vector propagators*)
	vectorPropT=Table[1/(t-i),{i,vectorMass}];
	vectorPropU=Table[1/(u-i),{i,vectorMass}];
	vectorPropS=Table[1/(s-i),{i,vectorMass}];

(*
Since there are two diagrams there can be
3 Lorentz structures after squaring and summing over spins
*)
(*They are just hardcoded for now*)
	A1=2(s^2+u^2); (*squared t-channel diagram*)
	A2=2(s^2+t^2); (*squared u-channel diagram*)
	A3=4 s^2; (*Interference between u and t channel diagrams*)
	A4=2(t^2+u^2)/s^2; (*Squared s-channel*)
	A5=4 u^2; (*Interference between s and t channel diagrams*)

(*We now generate the matrix element for Q1+Q2->Q3+Q4*)
	Res1=A1*DiagonalMatrix[vectorPropT] . C1 . DiagonalMatrix[vectorPropT]; 
	Res2=1/2*(A2*DiagonalMatrix[vectorPropU] . C2 . DiagonalMatrix[vectorPropU]);
	Res3=1/2*(A3*DiagonalMatrix[vectorPropU] . C3 . DiagonalMatrix[vectorPropT])/.(#->0&/@VectorMass);
	Res4=A4*C4;
	Res5=(A5*DiagonalMatrix[vectorPropS] . C5 . DiagonalMatrix[vectorPropT])/.(#->0&/@VectorMass);

	Return[ Total[Res1+Res2+Res3+If[leadingLog,Res4+Res5,0],-1]]
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
	gTensor[1,3]=gvff[[;;,p1,p3]];
	gTensor[1,2]=gvff[[p2,p1,;;]]//Flatten[#,{{3},{2},{1}}]&;
	gTensor[1,4]=gvff[[p4,p1,;;]]//Flatten[#,{{3},{2},{1}}]&;
	
	gTensor[3,4]=gvff[[p4,;;,p3]]//Flatten[#,{{2},{3},{1}}]&;
	gTensor[2,4]=gvvv[[p2,p4,;;]]//Flatten[#,{{3},{1},{2}}]&;
	gTensor[2,3]=gvff[[p2,;;,p3]]//Flatten[#,{{2},{3},{1}}]&;

	LV=Length[gTensor[1,3]];
	LF=Length[gTensor[1,2]];

(*Vector propagators*)
	vectorPropT=Table[1/(t-i),{i,vectorMass}];
	
(*Fermion propagators*)
	fermionPropU=Table[1/(u-i),{i,fermionMass}];
	fermionPropS=Table[1/(s),{i,fermionMass}];
	
(*Group invariants that multiply various Lorentz Structures*)
	C1=Sum[vectorPropT[[c]]vectorPropT[[d]]
		Tr[gTensor[1,3][[c]] . Transpose[gTensor[1,3][[d]]]]
		Tr[gTensor[2,4][[c]] . Transpose[gTensor[2,4][[d]]]],{c,1,LV},{d,1,LV}];
	C2=Sum[fermionPropU[[L]] fermionPropU[[K]]
		Tr[gTensor[1,4][[L]] . Transpose[gTensor[1,4][[K]]]]
		Tr[gTensor[2,3][[L]] . Transpose[gTensor[2,3][[K]]]],{L,1,LF},{K,1,LF}];
	C3=Sum[fermionPropS[[L]] fermionPropS[[K]]
		Tr[gTensor[1,2][[L]] . Transpose[gTensor[1,2][[K]]]]
		Tr[gTensor[3,4][[L]] . Transpose[gTensor[3,4][[K]]]],{L,1,LF},{K,1,LF}];
	
(*Since there are two diagrams there can be 3 Lorentz structures after squaring and summing over spins*)
(*The interference diagram does however not contribute at leading-log, so I omit it*)
(*They are just hardcoded for now*)
	A1=16(s^2+u^2); (*squared t-channel diagram*)
	A2=-4*16*s*u; (*squared u-channel diagram*)
	A3=-4*16*u*s;  (*squared s-channel diagram*)

(*We now generate the matrix element for Q1+Q2->Q3+Q4*)
	Res1=A1*C1;
	Res2=A2*C2;
	Res3=A3*C3;

(*Factor of 2 from anti-particle contribution*)
	Return[(Res1+Res2+If[leadingLog,Res3,0])]
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
	symmFac=1;
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
		symmFac=2; (*for quark final states we have to add vv-> q qbar+vv->qbar q; this overcounting is always compensated for later*)
	];
(*Coupling constants that we will need*)
	gTensor[1,3]=gvff[[p3,p1,;;]];
	gTensor[2,4]=gvff[[p4,;;,p2]];
	gTensorT[1,3]=gvff[[p3,;;,p1]];
	gTensorT[2,4]=gvff[[p4,p2,;;]];

	gTensor[1,4]=gvff[[p4,p1,;;]];
	gTensor[2,3]=gvff[[p3,;;,p2]];
	gTensorT[1,4]=gvff[[p4,;;,p1]];
	gTensorT[2,3]=gvff[[p3,p2,;;]];
	
	gTensor[1,2]=gvff[[;;,p1,p2]];
	gTensor[3,4]=gvvv[[;;,p3,p4]];
	
	LV=Length[gTensor[1,3]];
	F=Length[gTensor[1,2][[1]]];

(*Fermion propagators*)
	fermionPropU=Table[1/(u-i),{i,fermionMass}];
	fermionPropT=Table[1/(t-i),{i,fermionMass}];

(*Group invariants that multiply various Lorentz Structures*)
(*Multiplying with propagators for each particle*)
	gTensorT[1,3]=DiagonalTensor2[TensorProduct[gTensorT[1,3],fermionPropT],2,4]//Flatten[#,{{2},{1},{3}}]&;
	gTensor[1,3]=DiagonalTensor2[TensorProduct[gTensor[1,3],fermionPropT],3,4]//Flatten[#,{{2},{3},{1}}]&;
	
	C1=Sum[
		Tr[gTensor[1,3][[a]] . gTensor[2,4][[b]] . gTensorT[2,4][[b]] . gTensorT[1,3][[a]]],
		{a,Length[gTensor[1,3]]},{b,Length[gTensor[2,4]]}];

	gTensorT[1,4]=DiagonalTensor2[TensorProduct[gTensorT[1,4],fermionPropU],2,4]//Flatten[#,{{2},{1},{3}}]&;
	gTensor[1,4]=DiagonalTensor2[TensorProduct[gTensor[1,4],fermionPropU],3,4]//Flatten[#,{{2},{3},{1}}]&;

	C2=Sum[
		Tr[gTensor[1,4][[a]] . gTensor[2,3][[b]] . gTensorT[2,3][[b]] . gTensorT[1,4][[a]]],
		{a,Length[gTensor[1,4]]},{b,Length[gTensor[2,3]]}];
	C3=-Sum[
		Tr[gTensor[1,2][[a]] . Transpose[gTensor[1,2][[b]]]]
		Tr[gTensor[3,4][[a]] . Transpose[gTensor[3,4][[b]]]],
		{a,Length[gTensor[3,4]]},{b,Length[gTensor[3,4]]}];
(*Since there are two diagrams there can be 3 Lorentz structures after squaring and summing over spins*)
(*The interference diagram does however not contribute at leading-log, so I omit it*)
(*They are hardcoded for now*)
	A1=4*t*u; (*squared t-channel diagram*)
	A2=4*t*u; (*squared u-channel diagram*)
	A3=16*(t^2+u^2)/s^2; (*squared s-channel diagram*)

(*Generate the matrix element for Q1+Q2->Q3+Q4*)
	Res1=A1*C1;
	Res2=A2*C2;
	Res3=A3*C3;

(*Factor of 2 from anti-particles*)
	Return[symmFac*(Res1+Res2+If[leadingLog,Res3,0])]
]	
];


(* ::Subsubsection::Closed:: *)
(*S1S2toS3S4*)


CreateMatrixElementS1S2toS3S4[particle1_,particle2_,particle3_,particle4_,vectorMass_]:=
Block[{s,t,u,gTensor,leadingLog},
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

(*Group invariants that multiply various Lorentz Structures*)
	CSS=Table[
		Tr[gTensor[[1,2]][[a]] . Transpose[gTensor[[1,2]][[b]]]]
		Tr[gTensor[[3,4]][[a]] . Transpose[gTensor[[3,4]][[b]]]],
		{a,1,Length[gvss]},{b,1,Length[gvss]}];
	CTT=Table[
		Tr[gTensor[[1,3]][[a]] . Transpose[gTensor[[1,3]][[b]]]]
		Tr[gTensor[[2,4]][[a]] . Transpose[gTensor[[2,4]][[b]]]],
		{a,1,Length[gvss]},{b,1,Length[gvss]}];
	CUU=Table[
		Tr[gTensor[[1,4]][[a]] . Transpose[gTensor[[1,4]][[b]]]]
		Tr[gTensor[[2,3]][[a]] . Transpose[gTensor[[2,3]][[b]]]],
		{a,1,Length[gvss]},{b,1,Length[gvss]}];
		
	CMixST=Table[
		Tr[gTensor[[1,2]][[a]] . gTensor[[2,4]][[b]] . Transpose[gTensor[[3,4]][[a]]] . Transpose[gTensor[[1,3]][[b]]]],
		{a,1,Length[gvss]},{b,1,Length[gvss]}];
		
	CMixSU=Table[
		Tr[gTensor[[1,2]][[a]] . gTensor[[2,3]][[b]] . gTensor[[3,4]][[a]] . Transpose[gTensor[[1,4]][[b]]]],
		{a,1,Length[gvss]},{b,1,Length[gvss]}];
	CMixTU=Table[
		Tr[gTensor[[1,3]][[a]] . Transpose[gTensor[[2,3]][[b]]] . gTensor[[2,4]][[a]] . Transpose[gTensor[[1,4]][[b]]]],
		{a,1,Length[gvss]},{b,1,Length[gvss]}];
		
(*Vector propagators*)
	vectorPropT=Table[1/(t-i),{i,vectorMass}];
	vectorPropU=Table[1/(u-i),{i,vectorMass}];
	
(*Lorentz structures*)
	ASS=(t-u)^2; (*squared u-channel diagram*)
	ATT=(s-u)^2; (*squared t-channel diagram*)
	AUU=(t-s)^2; (*squared u-channel diagram*)
	
	ASU=-2(t-u)*(u-s);
	AST=-2(t-u)*(t-s);
	ATU=2(u-s)(t-s);

(*Leading-log terms*)
	ResLL=ATT*vectorPropT . CTT . vectorPropT+AUU vectorPropU . CUU . vectorPropU;
	
(*non-divergent terms*)
	ResBLL=Total[ASU*CMixSU/(s*u)+AST*CMixST/(s*t)+ATU*CMixTU/(t*u),-1];
	
	ResQuartic=Total[\[Lambda]4[[particle1[[1]],particle2[[1]],particle3[[1]],particle4[[1]]]]^2,-1];
	
	
	Return[ResLL+If[leadingLog,ResBLL+ResQuartic,0]]
]
];


CreateMatrixElementS1S2toF1F2[particle1_,particle2_,particle3_,particle4_,fermionMass_]:=
Block[{s,t,u,gTensor,gTensorT,gTensorC,gTensorCT,leadingLog},
(*

*)
leadingLog=False;
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

(*Fermion propagators*)
	fermionPropT=Table[1/(t-i),{i,FermionMass}];
	fermionPropU=Table[1/(u-i),{i,FermionMass}];
	
(*Coupling factors that appear*)
	CTT=Table[
	Tr[gTensor[1,3][[;;,;;,a]] . Transpose[gTensorC[1,3][[;;,;;,a]]]]Tr[gTensor[2,4][[;;,;;,b]] . Transpose[gTensorC[2,4][[;;,;;,b]]]],{a,1,Length[gvff[[1,;;]]]},{b,1,Length[gvff[[1,;;]]]}
	];
	CUU=Table[
	Tr[gTensor[1,4][[;;,;;,a]] . Transpose[gTensorC[1,4][[;;,;;,a]]]]Tr[gTensor[2,3][[;;,;;,b]] . Transpose[gTensorC[2,3][[;;,;;,b]]]],{a,1,Length[gvff[[1,;;]]]},{b,1,Length[gvff[[1,;;]]]}
	];
	
	CUT=Sum[Tr[gTensor[1,3][[;;,;;,a]] . Transpose[gTensorCT[2,3][[;;,b,;;]]] . gTensorT[2,4][[;;,a,;;]] . Transpose[gTensorC[1,4][[;;,;;,b]]]],{a,1,Length[gvff[[1,;;]]]},{b,1,Length[gvff[[1,;;]]]}
	];
(*Lorentz structures*)

	ATT=t u;
	AUU=t u;
	
	AUT=-2;
	
(*The result*)
	ResLL=ATT fermionPropT . CTT . fermionPropT+AUU fermionPropU . CUU . fermionPropU;
	
	ResBLL=AUT*CUT;
(*The full result*)
	Return[ResLL+ResBLL]
]	
];


(* ::Subsubsection::Closed:: *)
(*S1S2toF1F2*)


CreateMatrixElementS1S2toF1F2[particle1_,particle2_,particle3_,particle4_,fermionMass_]:=
Block[{s,t,u,gTensor,gTensorT,gTensorC,gTensorCT,leadingLog,gTensorVF,gTensorVFC,gTensorVS},
(*
		There is no trilinear scalar interaction taken into account as the symmetric phase is assumed.
		The full Yukawa contribution is included, as well as the s-channel exchange of gauge bosons.
		No mixing between gauge and Yukawa terms was added.
*)
leadingLog=False;
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
	fermionPropT=Table[1/(t-i),{i,FermionMass}];
	fermionPropU=Table[1/(u-i),{i,FermionMass}];

(*Coupling factors that appear*)
	CTT=1/4Table[
	Tr[gTensor[1,3][[;;,;;,a]] . Transpose[gTensorC[1,3][[;;,;;,a]]]]Tr[gTensor[2,4][[;;,;;,b]] . Transpose[gTensorC[2,4][[;;,;;,b]]]],{a,1,Length[YsffC[[1,;;]]]},{b,1,Length[YsffC[[1,;;]]]}
	];
	CUU=1/4Table[
	Tr[gTensor[1,4][[;;,;;,a]] . Transpose[gTensorC[1,4][[;;,;;,a]]]]Tr[gTensor[2,3][[;;,;;,b]] . Transpose[gTensorC[2,3][[;;,;;,b]]]],{a,1,Length[YsffC[[1,;;]]]},{b,1,Length[YsffC[[1,;;]]]}
	];
	
	CUT=1/2 Sum[Tr[gTensor[1,3][[;;,;;,a]] . Transpose[gTensorCT[2,3][[;;,b,;;]]] . gTensorT[2,4][[;;,a,;;]] . Transpose[gTensorC[1,4][[;;,;;,b]]]],{a,1,Length[YsffC[[1,;;]]]},{b,1,Length[YsffC[[1,;;]]]}
	];
	
	CSS=1/2Sum[Tr[gTensorVS[[a]] . gTensorVS[[b]]]Tr[gTensorVF[[a]] . gTensorVFC[[b]]],{a,1,Length[gvff]},{b,1,Length[gvff]}];
(*Lorentz structures*)

	ATT=t u;
	AUU=t u;
	AUT=-2;
	ASS=u t /s^2;
(*The result*)
	ResLL=ATT fermionPropT . CTT . fermionPropT+AUU fermionPropU . CUU . fermionPropU;
	
	ResBLL=AUT*CUT+ASS CSS;
(*The full result*)
	Return[2(ResLL+ResBLL)]
]	
];


(* ::Subsubsection::Closed:: *)
(*Q1S1toQ1S1*)


CreateMatrixElementQ1S1toQ1S1[particle1_,particle2_,particle3_,particle4_,vectorMass_,fermionMass_]:=
Block[{s,t,u,gTensor,leadingLog,gTensorC,gTensorVF,gTensorVFC,gTensorVS},
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
(*Changing the order of the particles so that it is always QV->QV*)	
	If[particle1[[2]]=="S"&&particle2[[2]]=="F",
		temp=p1;
		p1=p2;
		p2=temp;
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
	CSS=1/4Sum[Tr[gTensor[1,2][[;;,;;,a]] . Transpose[gTensorC[1,2][[;;,;;,b]]]]Tr[gTensorC[3,4][[;;,;;,a]] . Transpose[gTensor[3,4][[;;,;;,b]]]],{a,Length[Ysff[[1]]]},{b,Length[Ysff[[1]]]}];

	CUU=1/4 Table[Tr[gTensor[1,4][[;;,;;,a]] . Transpose[gTensorC[1,4][[;;,;;,b]]]]Tr[gTensorC[3,2][[;;,;;,a]] . Transpose[gTensor[3,2][[;;,;;,b]]]],{a,Length[Ysff[[1]]]},{b,Length[Ysff[[1]]]}];

	CTT=1/2 Table[Tr[gTensorVS[[a]] . gTensorVS[[b]]]Tr[gTensorVF[[a]] . gTensorVFC[[b]]],{a,Length[gvss]},{b,Length[gvss]}];
(*Lorentz structures that appear*)
(*The first three are the fermion channels, the third is the vector channel*)
	ASS=u/s;
	AUU=s u;
	
	ATT=u s;
	
	ResSS=CSS ASS;
	ResUU=fermionPropU . CUU . fermionPropU;
	
	ResTT=vectorPropT . CTT . vectorPropT;
	

	Return[2*(ResSS+ResUU+ResTT)]
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
	CSS=1/4 Sum[Tr[gTensor[1,2][[;;,;;,a]] . Transpose[gTensorC[1,2][[;;,;;,b]]]]Tr[gTensorVF[3,4][[;;,;;,a]] . Transpose[gTensorVFC[3,4][[;;,b,;;]]]],{a,Length[Ysff[[1]]]},{b,Length[Ysff[[1]]]}];
	CUU=1/4 Table[Tr[gTensor[3,2][[;;,;;,a]] . Transpose[gTensorC[3,2][[;;,;;,b]]]]Tr[gTensorVF[1,4][[;;,;;,a]] . Transpose[gTensorVFC[1,4][[;;,b,;;]]]],{a,Length[Ysff[[1]]]},{b,Length[Ysff[[1]]]}];
	
(*Lorentz structures that appear*)
	ASS=4*u/s;
	AUU=4*s;
	
(*Results*)
	ResSS=CSS*ASS;
	ResUU=AUU fermionPropU . CUU . fermionPropU;
	Return[ResSS+ResUU]
,
	Return[0]
]	
];


(* ::Section::Closed:: *)
(*Getting the matrix elements for out-of-Equilibrium particles*)


degreeOfFreedom[particle_]:=Block[{dof},

	dof=Length[particle[[1]]];
	
	(*Factor of 2 from helicites, factor of 2 from anti-particles*)
	If[particle[[2]]=="F",dof];
	
	(*Factor of 2 from spins*)
	If[particle[[2]]=="V",dof*=2]; 
	
	If[particle[[2]]=="S",dof];
	
	Return[dof];
]


(* ::Subsubsection::Closed:: *)
(*Q1Q2toQ3Q4*)


ExtractOutOfEqElementQ1Q2toQ3Q4[particleList_,LightParticles_]:=
Block[{OutOfEqParticles,MatrixElements,Elements,CollElements},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into
	Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

(*particleList is the complete list of particles*)

(*Generate all matrix elements*)
	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];

(*Divide the incoming particle by its degree's of freedom*)
	MatrixElements=Table[
		1/degreeOfFreedom[particleList[[a]]]
		CreateMatrixElementQ1Q2toQ3Q4[particleList[[a]],b,c,d,VectorMass],
		{a,OutOfEqParticles},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
(*This is a list of all non-zero matrix elements*)
	Elements=MatrixElements["NonzeroPositions"];  

(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,

(*Now we add identical contributions*)
(*The Q1Q2->Q1Q2 process is the same as the Q1Q2->Q2Q1 process*)
	symmetries=Gather[Elements,(#1[[1;;2]]==#2[[1;;2]])&];
	symmetries=Gather[#,#1[[3;;4]]==Sort[#2[[3;;4]]]&]&/@symmetries//Flatten[#,1]&;
	MultiPlicity=Table[Length[a],{a,symmetries}];
	Elements=Table[i[[1]],{i,symmetries}];
	(*MatrixElements is now just a list*)
	MatrixElements=Table[MultiPlicity[[i]] Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}]; 
	
(*We also divide by the number of degrees of freedom of incomingParticle*)
(*The 2 is from the anti-particle contribution. I.e q1 q2->q1 q2 +q1 q2Bar->q1 q2Bar*)
(*Factor of 2 from over-counting ab->cd +ab->dc, or if c=d the 1/2 is a symmetry factor*)
	MatrixElements=2*MatrixElements/2;

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


ExtractOutOfEqElementQ1V1toQ1V1[particleList_,LightParticles_]:=
Block[{OutOfEqParticles,MatrixElements,Elements,CollElements},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into
	Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

(*particleList is the complete list of particles*)

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

(*Now we add identical contributions*)
(*The Q1Q2->Q1Q2 process is the same as the Q1Q2->Q2Q1 process*)
	symmetries=Gather[Elements,(#1[[1;;2]]==#2[[1;;2]])&];
	symmetries=Gather[#,#1[[3;;4]]==Sort[#2[[3;;4]]]&]&/@symmetries//Flatten[#,1]&;

	MultiPlicity=Table[Length[a],{a,symmetries}];
	Elements=Table[i[[1]],{i,symmetries}];
	(*MatrixElements is now just a list*)
	MatrixElements=Table[MultiPlicity[[i]] Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}];

	MatrixElements=MatrixElements/2;(*Factor of 2 from over-counting ab->cd +ab->dc, or if c=d the 1/2 is a symmetry factor*)

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


ExtractOutOfEqElementQ1Q2toV1V2[particleList_,LightParticles_]:=
Block[{OutOfEqParticles,MatrixElements,Elements,CollElements},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into
	Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

(*particleList is the complete list of particles*)

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

(*Now we add identical contributions*)
(*The Q1Q2->Q1Q2 process is the same as the Q1Q2->Q2Q1 process*)
	symmetries=Gather[Elements,(#1[[1;;2]]==#2[[1;;2]])&];
	symmetries=Gather[#,#1[[3;;4]]==Sort[#2[[3;;4]]]&]&/@symmetries//Flatten[#,1]&;
	MultiPlicity=Table[Length[a],{a,symmetries}];
	Elements=Table[i[[1]],{i,symmetries}];
	(*MatrixElements is now just a list*)
	MatrixElements=Table[MultiPlicity[[i]] Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}];

(*We also divide by the number of degrees of freedom of incomingParticle*)
(*Factor of 2 from over-counting ab->cd +ab->dc, or if c=d the 1/2 is a symmetry factor*)
	MatrixElements=MatrixElements/2;

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


ExtractOutOfEqElementV1V2toV3V4[particleList_,LightParticles_]:=
Block[{OutOfEqParticles,MatrixElements,Elements,CollElements},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into Sum_deltaFparticle
	\[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

(*particleList is the complete list of particles*)

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

(*Now we add identical contributions*)
(*The Q1Q2->Q1Q2 process is the same as the Q1Q2->Q2Q1 process*)
	symmetries=Gather[Elements,(#1[[1;;2]]==#2[[1;;2]])&];
	symmetries=Gather[#,#1[[3;;4]]==Sort[#2[[3;;4]]]&]&/@symmetries//Flatten[#,1]&;
	MultiPlicity=Table[Length[a],{a,symmetries}];
	Elements=Table[i[[1]],{i,symmetries}];
	(*MatrixElements is now just a list*)
	MatrixElements=Table[MultiPlicity[[i]] Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}];

(*We also divide by the number of degrees of freedom of incomingParticle*)
(*Factor of 2 from over-counting ab->cd +ab->dc, or if c=d the 1/2 is a symmetry factor*)
	MatrixElements=MatrixElements/2;

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


ExtractOutOfEqElementS1S2toS3S4[particleList_,LightParticles_]:=
Block[{OutOfEqParticles,MatrixElements,Elements,CollElements},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into Sum_deltaFparticle
	\[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

(*particleList is the complete list of particles*)

(*Generate all matrix elements*)
	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];
	
(*Divide the incoming particle by its degree's of freedom*)
	MatrixElements=Table[
		1/degreeOfFreedom[particleList[[a]]]
		CreateMatrixElementS1S2toS3S4[particleList[[a]],b,c,d,VectorMass],
		{a,OutOfEqParticles},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
	(*This is a list of all non-zero matrix elements*)
	Elements=MatrixElements["NonzeroPositions"];

(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,

(*Now we add identical contributions*)
	symmetries=Gather[Elements,(#1[[1;;2]]==#2[[1;;2]])&];
	symmetries=Gather[#,#1[[3;;4]]==Sort[#2[[3;;4]]]&]&/@symmetries//Flatten[#,1]&;
	MultiPlicity=Table[Length[a],{a,symmetries}];
	Elements=Table[i[[1]],{i,symmetries}];
	(*MatrixElements is now just a list*)
	MatrixElements=Table[MultiPlicity[[i]] Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}];

(*We also divide by the number of degrees of freedom of incomingParticle*)
(*Factor of 2 from over-counting ab->cd +ab->dc, or if c=d the 1/2 is a symmetry factor*)
	MatrixElements=MatrixElements/2;

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


ExtractOutOfEqElementS1S2toF1F2[particleList_,LightParticles_]:=
Block[{OutOfEqParticles,MatrixElements,Elements,CollElements},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into
	Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

(*particleList is the complete list of particles*)

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

(*Now we add identical contributions*)
(*The Q1Q2->Q1Q2 process is the same as the Q1Q2->Q2Q1 process*)
	symmetries=Gather[Elements,(#1[[1;;2]]==#2[[1;;2]])&];
	symmetries=Gather[#,#1[[3;;4]]==Sort[#2[[3;;4]]]&]&/@symmetries//Flatten[#,1]&;
	MultiPlicity=Table[Length[a],{a,symmetries}];
	Elements=Table[i[[1]],{i,symmetries}];
	(*MatrixElements is now just a list*)
	MatrixElements=Table[MultiPlicity[[i]] Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}];

(*We also divide by the number of degrees of freedom of incomingParticle*)
(*Factor of 2 from over-counting ab->cd +ab->dc, or if c=d the 1/2 is a symmetry factor*)
	MatrixElements=2*MatrixElements/2;

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


ExtractOutOfEqElementQ1S1toQ1S1[particleList_,LightParticles_]:=
Block[{OutOfEqParticles,MatrixElements,Elements,CollElements},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into
	Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

(*particleList is the complete list of particles*)

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

(*Now we add identical contributions*)
(*The Q1Q2->Q1Q2 process is the same as the Q1Q2->Q2Q1 process*)
	symmetries=Gather[Elements,(#1[[1;;2]]==#2[[1;;2]])&];
	symmetries=Gather[#,#1[[3;;4]]==Sort[#2[[3;;4]]]&]&/@symmetries//Flatten[#,1]&;

	MultiPlicity=Table[Length[a],{a,symmetries}];
	Elements=Table[i[[1]],{i,symmetries}];
	(*MatrixElements is now just a list*)
	MatrixElements=Table[MultiPlicity[[i]] Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}];

	MatrixElements=2 MatrixElements/2;(*Factor of 2 from over-counting ab->cd +ab->dc, or if c=d the 1/2 is a symmetry factor*)

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


ExtractOutOfEqElementQ1S1toQ1V1[particleList_,LightParticles_]:=
Block[{OutOfEqParticles,MatrixElements,Elements,CollElements},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into
	Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

(*particleList is the complete list of particles*)

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

(*Now we add identical contributions*)
(*The Q1Q2->Q1Q2 process is the same as the Q1Q2->Q2Q1 process*)
	symmetries=Gather[Elements,(#1[[1;;2]]==#2[[1;;2]])&];
	symmetries=Gather[#,#1[[3;;4]]==Sort[#2[[3;;4]]]&]&/@symmetries//Flatten[#,1]&;

	MultiPlicity=Table[Length[a],{a,symmetries}];
	Elements=Table[i[[1]],{i,symmetries}];
	(*MatrixElements is now just a list*)
	MatrixElements=Table[MultiPlicity[[i]] Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}];

	MatrixElements=MatrixElements/2;(*Factor of 2 from over-counting ab->cd +ab->dc, or if c=d the 1/2 is a symmetry factor*)

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


ExtractOutOfEqElement[particleList_,LightParticles_]:=
Block[{CollEllQ1Q2toQ3Q4,CollEllQ1V1toQ1V1,CollEllQ1Q2toV1V2,CollEllV1V2toV3V4,CollEllTotal,CollEllS1S2toS3S4,CollEllS1S2toF1F2,
CollEllQ1S1toQ1S1,CollEllQ1S1toQ1V1},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*Essentially this is generates the elements going into Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]*)

(*particleList is the complete list of particles*)


(*First we extract the result for all subprocesses*)
CollEllQ1Q2toQ3Q4=ExtractOutOfEqElementQ1Q2toQ3Q4[particleList,LightParticles];
CollEllQ1V1toQ1V1=ExtractOutOfEqElementQ1V1toQ1V1[particleList,LightParticles];
CollEllQ1Q2toV1V2=ExtractOutOfEqElementQ1Q2toV1V2[particleList,LightParticles];
CollEllV1V2toV3V4=ExtractOutOfEqElementV1V2toV3V4[particleList,LightParticles];
CollEllS1S2toS3S4=ExtractOutOfEqElementS1S2toS3S4[particleList,LightParticles];
CollEllS1S2toF1F2=ExtractOutOfEqElementS1S2toF1F2[particleList,LightParticles];
CollEllQ1S1toQ1S1=ExtractOutOfEqElementQ1S1toQ1S1[particleList,LightParticles];
CollEllQ1S1toQ1V1=ExtractOutOfEqElementQ1S1toQ1V1[particleList,LightParticles];

CollEllTotal=Join[CollEllQ1Q2toQ3Q4,CollEllQ1V1toQ1V1,CollEllQ1Q2toV1V2,CollEllV1V2toV3V4,CollEllS1S2toS3S4
				,CollEllS1S2toF1F2,CollEllQ1S1toQ1S1,CollEllQ1S1toQ1V1];


Return[CollEllTotal]

];


(* ::Section::Closed:: *)
(*Exporting to C++*)


ExportMatrixElements[file_,particleList__,UserMasses_,UserCouplings_,ParticleName_]:=
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
	
	particleListFull2=particleListFull;
	LightParticles2=LightParticles;
(*
Replacement rules for converting to the format used by the c++ code
*)	
	RepMasses=Table[UserMasses[[i]]->msq[i-1],{i,1,Length[UserMasses]}];
	RepCouplings=Table[UserCouplings[[i]]->Symbol["c"][i-1],{i,1,Length[UserCouplings]}];
	
(*
Loop over all out-of-eq particles and extracting the matrix elements
*)	
	MatrixElements=ExtractOutOfEqElement[particleListFull,LightParticles];
	
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
	MatrixElements=Table[MatrixElemToC@i/.RepCouplings/.RepMasses,{i,MatrixElements}];
	
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


(* ::Title:: *)
(*SM quarks + gauge bosons*)


(* ::Subtitle:: *)
(*SymmetryBreaking*)


vev={0,v,0,0};


SymmetryBreaking[vev]


(* ::Subtitle:: *)
(*UserInput*)


(*
In DRalgo fermions are Weyl.
So to create one Dirac we need
one left-handed and
one right-handed fermoon
*)


(*
	Reps 1-4 are quarks,
	reps 5,6 are vector bosons
*)
(*left-handed top-quark*)
ReptL=CreateOutOfEq[{{1,1}},"F"];

(*right-handed top-quark*)
ReptR=CreateOutOfEq[{2},"F"];

(*right-handed bottom-quark*)
RepbR=CreateOutOfEq[{3},"F"];

(*Vector bosons*)
RepGluon=CreateOutOfEq[{1},"V"];
RepW=CreateOutOfEq[{{2,1}},"V"];

(*Scalar particles*)
RepHiggs=CreateOutOfEq[{1},"S"];


ParticleList={ReptL,ReptR,RepbR,RepGluon,RepW,RepHiggs};


(*Defining various masses and couplings*)


VectorMass=Join[
	Table[mg2,{i,1,RepGluon[[1]]//Length}],
	Table[mw2,{i,1,RepW[[1]]//Length}]];
FermionMass=Table[mq2,{i,1,Length[gvff[[1]]]}];
(*
up to the user to make sure that the same order is given in the python code
*)
UserMasses={mq2,mq2,mg2,mw2}; 
UserCouplings=CouplingName;


OutputFile="matrixElements.scalar";
SetDirectory[NotebookDirectory[]];
ParticleName={"TopL","TopR","BotR","Gluon","W","Higgs"};
MatrixElements=ExportMatrixElements[OutputFile,ParticleList,UserMasses,UserCouplings,ParticleName];


MatrixElements//Expand


(*g g->g g*)
M[0,3,0,3]/.MatrixElements
(*g g->g g*)
M[3,3,3,3]/.MatrixElements
(*t g->t g*)
M[1,3,1,3]/.MatrixElements
(*t q->t q*)
5/4*M[1,5,1,5]/.MatrixElements


Import["MatrixElem.hdf5"]


Import["MatrixElem.hdf5","CouplingInfo"]


Import["MatrixElem.hdf5","ParticleInfo"]


Import["MatrixElem.hdf5","CouplingInfo"]


Import["MatrixElem.hdf5","ParticleInfo"]


(* ::Title:: *)
(*Scalar processes*)
