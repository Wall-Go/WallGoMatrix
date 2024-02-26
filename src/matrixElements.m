(* ::Package:: *)

(* ::Title:: *)
(*Matrix elements*)


(* ::Section:: *)
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


CreateOutOfEq[Indices_,Type_]:=Block[{PosScalar,PosVector,PosFermion},
(*Combines selected representations and obtain their indices*)
	PosScalar=PrintScalarRepPositions[];
	PosVector=PrintGaugeRepPositions[];
	PosFermion=PrintFermionRepPositions[];
	If[Type=="F",
		temp=PosFermion[[Indices]];
		Return[{RepToIndices[temp],Type}]
	];

	If[Type=="S",
		temp=PosScalar[[Indices]];
		Return[{RepToIndices[temp],Type}]
	];

	If[Type=="V",
		temp=PosVector[[Indices]];
		Return[{RepToIndices[temp],Type}]
	];

];


SymmetryBreaking[Indices_,vev_] :=Block[{PosScalar,PosVector,PosFermion},
(*
Routine that extracts representations given their vev.
This way one can identify which e.g. parts of multiplet fermions are light or heavy
*)
	PosScalar=PrintScalarRepPositions[];
	PosVector=PrintGaugeRepPositions[];
	PosFermion=PrintFermionRepPositions[];
	
(*Fermions*)
	massF=vev . Ysff;
	fermionInd=Delete[massF//ArrayRules,-1]/.(({a_,b_}->x_)->a);
	
	posHeavy=Intersection[RangeToIndices[PosFermion[[Indices]]],fermionInd];
	posLight=Complement[RangeToIndices[PosFermion[[Indices]]],fermionInd];
		
	val=massF[[posHeavy,;;]]["NonzeroValues"]//DeleteDuplicates;
	pos=Table[i,{i,Length[posHeavy]}];
	
	rep={};
	
	Do[
	pos2=Table[posHeavy[[pos[[a]][[1]]]],{a,Position[massF[[posHeavy,;;]]["NonzeroValues"],a]}];
	AppendTo[rep,{pos2,a}];
	,{a,val}];

	AppendTo[rep,{posLight,0}];
	
	rep
]


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


(* ::Subsubsection:: *)
(*Q1Q2toQ3Q4*)


CreateMatrixElementQ1Q2toQ3Q4Pre[particle1_,particle2_,particle3_,particle4_,vectorMass_]:=Block[{},
	Return[
		1/2*(CreateMatrixElementQ1Q2toQ3Q4[particle1,particle2,particle3,particle4,vectorMass]
		+CreateMatrixElementQ1Q2toQ3Q4[particle1,particle2,particle4,particle3,vectorMass])
	];
];


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
		Tr[gTensor[[1,3]][[a]] . gTensor[[3,1]][[b]]]
		Tr[gTensor[[2,4]][[a]] . gTensor[[4,2]][[b]]],
		{a,1,Length[gTensor[[1,3]]]},{b,1,Length[gTensor[[2,4]]]}];
	C2=Table[
		Tr[gTensor[[1,4]][[a]] . gTensor[[4,1]][[b]]]
		Tr[gTensor[[2,3]][[a]] . gTensor[[3,2]][[b]]],
		{a,1,Length[gTensor[[1,3]]]},{b,1,Length[gTensor[[1,3]]]}];
	
	C3=Table[
		Tr[gTensor[[1,3]][[a]] . gTensor[[3,2]][[b]] . gTensor[[2,4]][[a]] . gTensor[[4,1]][[b]]],
		{a,1,Length[gTensor[[1,3]]]},{b,1,Length[gTensor[[2,4]]]}];
	C3+=Table[
		Tr[gTensor[[3,1]][[a]] . gTensor[[1,4]][[b]] . gTensor[[4,2]][[a]] . gTensor[[2,3]][[b]]],
		{a,1,Length[gTensor[[1,4]]]},{b,1,Length[gTensor[[2,3]]]}];
	
	C4=Table[
		Tr[gTensor[[1,2]][[a]] . gTensor[[2,1]][[b]]]
		Tr[gTensor[[3,4]][[a]] . gTensor[[4,3]][[b]]],
		{a,1,Length[gTensor[[1,2]]]},{b,1,Length[gTensor[[3,4]]]}];
	
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
	A3=4*s^2; (*Interference between u and t channel diagrams*)
	A4=2(t^2+u^2)/s^2; (*Squared s-channel*)
	A5=4*u^2; (*Interference between s and t channel diagrams*)
	A6=4*s^2; (*Interference between s and u channel diagrams*)

(*We now generate the matrix element for Q1+Q2->Q3+Q4*)
	Res1=A1*DiagonalMatrix[vectorPropT] . C1 . DiagonalMatrix[vectorPropT]; 
	Res2=A2*DiagonalMatrix[vectorPropU] . C2 . DiagonalMatrix[vectorPropU];
	Res3=A3*DiagonalMatrix[vectorPropU] . C3 . DiagonalMatrix[vectorPropT]/.(#->0&/@VectorMass);
	Res4=A4*C4;
	Res5=1/2*(A5*DiagonalMatrix[vectorPropS] . C5 . DiagonalMatrix[vectorPropT])/.(#->0&/@VectorMass);
	Res5+=1/2*(A6*DiagonalMatrix[vectorPropS] . C5 . DiagonalMatrix[vectorPropU])/.(#->0&/@VectorMass);


	Return[ Total[Res1+Res2+Res3+If[leadingLog,Res4+Res5,0],-1]]
]
];


(* ::Subsubsection:: *)
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
		Tr[gTensor[4,1][[b]] . gTensor[1,4][[a]] . gTensor[4,3][[b]] . gTensor[2,3][[a]]],
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

	Return[(Res1+Res2+If[leadingLog,Res3,0])]
,
	Return[0]
]	
];


(* ::Subsubsection:: *)
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

(*Factor of 2 from anti-particles*)
	Return[symmFac*(Res1+Res2+If[leadingLog,Res3,0])]
]	
];


(* ::Section:: *)
(*Getting the matrix elements for out-of-Equilibrium particles*)


degreeOfFreedom[particle_]:=Block[{dof},

	dof=Length[particle[[1]]];
	
	(*Factor of 2 from helicites, factor of 2 from anti-particles*)
	If[particle[[2]]=="F",dof];
	
	(*Factor of 2 from spins*)
	If[particle[[2]]=="V",dof*=2]; 
	
	Return[dof];
]


(* ::Subsubsection:: *)
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
		CreateMatrixElementQ1Q2toQ3Q4Pre[particleList[[a]],b,c,d,VectorMass],
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


(* ::Subsubsection:: *)
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


(* ::Subsubsection:: *)
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


(* ::Subsubsection:: *)
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


(* ::Subsubsection:: *)
(*Extract elements*)


ExtractOutOfEqElement[particleList_,LightParticles_]:=
Block[{CollEllQ1Q2toQ3Q4,CollEllQ1V1toQ1V1,CollEllQ1Q2toV1V2,CollEllV1V2toV3V4,CollEllTotal},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*Essentially this is generates the elements going into Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]*)

(*particleList is the complete list of particles*)


(*First we extract the result for all subprocesses*)
CollEllQ1Q2toQ3Q4=ExtractOutOfEqElementQ1Q2toQ3Q4[particleList,LightParticles];
CollEllQ1V1toQ1V1=ExtractOutOfEqElementQ1V1toQ1V1[particleList,LightParticles];
CollEllQ1Q2toV1V2=ExtractOutOfEqElementQ1Q2toV1V2[particleList,LightParticles];
CollEllV1V2toV3V4=ExtractOutOfEqElementV1V2toV3V4[particleList,LightParticles];

CollEllTotal=Join[CollEllQ1Q2toQ3Q4,CollEllQ1V1toQ1V1,CollEllQ1Q2toV1V2,CollEllV1V2toV3V4];


(*I have so far not added all the terms with the same deltaF indices, which we should of course do*)
(*This is trivial to do, and I'll do it after we have agreed what format we want*)

Return[CollEllTotal]

];


(* ::Section:: *)
(*Exporting to C++*)


ExportMatrixElements[file_,particleList_,LightParticles_,UserMasses_,UserCouplings_,ParticleName_]:=
Block[{ExportTXT,ExportH5,
	Cij,ParticleInfo,CouplingInfo,MatrixElements,
	OutOfEqParticles,RepMasses,RepCouplings,C
	},

(*
Extracting the out-of-eq particles
*)
	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];

(*
Replacement rules for converting to the format used by the c++ code
*)	
	RepMasses=Table[UserMasses[[i]]->msq[i-1],{i,1,Length[UserMasses]}];
	RepCouplings=Table[UserCouplings[[i]]->Symbol["c"][i-1],{i,1,Length[UserCouplings]}];
	
(*
Loop over all out-of-eq particles and extracting the matrix elements
*)	
	MatrixElements=ExtractOutOfEqElement[particleList,LightParticles];
	
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
