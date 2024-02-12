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


(* ::Section:: *)
(*Matrix elements*)


(* ::Subsubsection:: *)
(*V1V2toV3V4*)


CreateMatrixElementV1V2toV3V4[particle1_,particle2_,particle3_,particle4_,vectorMass_]:=
Block[{s,t,u},
(*
In QCD this process depends on up to
two Lorentz structures if the particles are all the same; and
one Lorentz structure if they are all different.
*)
If[
	particle1[[2]]!="V"||
	particle2[[2]]!="V"||
	particle3[[2]]!="V"||
	particle4[[2]]!="V",
	Return[0];
,
(*Coupling constants that we will need*)
	g13=gvvv[[;;,particle1[[1]],particle3[[1]]]];
	g24=gvvv[[;;,particle2[[1]],particle4[[1]]]];
	g14=gvvv[[;;,particle1[[1]],particle4[[1]]]];
	g23=gvvv[[;;,particle2[[1]],particle3[[1]]]];

(*Group invariants that multiply various Lorentz Structures*)
	C1=Table[Tr[g13[[a]] . Transpose[g13[[b]]]]Tr[g24[[a]] . Transpose[g24[[b]]]],{a,1,Length[g13]},{b,1,Length[g13]}];
	C2=Table[Tr[g14[[a]] . Transpose[g14[[b]]]]Tr[g23[[a]] . Transpose[g23[[b]]]],{a,1,Length[g13]},{b,1,Length[g13]}];
	C3=Table[Tr[g13[[a]] . Transpose[g23[[b]]] . g24[[a]] . Transpose[g14[[b]]]],{a,1,Length[g13]},{b,1,Length[g13]}];
	C3+=Table[Tr[g14[[a]] . Transpose[g24[[b]]] . g23[[a]] . Transpose[g13[[b]]]],{a,1,Length[g13]},{b,1,Length[g13]}];
		
(*Vector propagators*)
	vectorPropT=Table[1/(t-i),{i,vectorMass}];
	vectorPropU=Table[1/(u-i),{i,vectorMass}];

(*Since there are two diagrams there can be 3 Lorentz structures after squaring and summing over spins*)
(*They are just hardcoded for now*)
	A1=16(-1/4)(s-u)^2; (*squared t-channel diagram*)
	A2=16(-1/4)(s-t)^2; (*squared u-channel diagram*)

(*We now generate the matrix element for Q1+Q2->Q3+Q4*)
	Res1=A1*DiagonalMatrix[vectorPropT] . C1 . DiagonalMatrix[vectorPropT];
	Res2=A2*DiagonalMatrix[vectorPropU] . C2 . DiagonalMatrix[vectorPropU];
	
(*Factor 4 from anti-particle contributions*)
	Return[-Total[Res1+Res2,-1]]
]
	
];


(* ::Subsubsection:: *)
(*Q1Q2toQ3Q4*)


CreateMatrixElementQ1Q2toQ3Q4[particle1_,particle2_,particle3_,particle4_,vectorMass_]:=
Block[{s,t,u},
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
	g13=gvff[[;;,particle1[[1]],particle3[[1]]]];
	g24=gvff[[;;,particle2[[1]],particle4[[1]]]];
	g14=gvff[[;;,particle1[[1]],particle4[[1]]]];
	g23=gvff[[;;,particle2[[1]],particle3[[1]]]];

(*Group invariants that multiply various Lorentz Structures*)
	C1=Table[Tr[g13[[a]] . Transpose[g13[[b]]]]Tr[g24[[a]] . Transpose[g24[[b]]]],{a,1,Length[g13]},{b,1,Length[g13]}];
	C2=Table[Tr[g14[[a]] . Transpose[g14[[b]]]]Tr[g23[[a]] . Transpose[g23[[b]]]],{a,1,Length[g13]},{b,1,Length[g13]}];
	C3=Table[Tr[g13[[a]] . Transpose[g23[[b]]] . g24[[a]] . Transpose[g14[[b]]]],{a,1,Length[g13]},{b,1,Length[g13]}];
	C3+=Table[Tr[g14[[a]] . Transpose[g24[[b]]] . g23[[a]] . Transpose[g13[[b]]]],{a,1,Length[g13]},{b,1,Length[g13]}];
		
(*Vector propagators*)
	vectorPropT=Table[1/(t-i),{i,vectorMass}];
	vectorPropU=Table[1/(u-i),{i,vectorMass}];

(*
Since there are two diagrams there can be
3 Lorentz structures after squaring and summing over spins
*)
(*They are just hardcoded for now*)
	A1=2(s^2+u^2); (*squared t-channel diagram*)
	A2=2(s^2+t^2); (*squared u-channel diagram*)
	A3=4 s^2; (*Interference between u and t channel diagrams*)

(*We now generate the matrix element for Q1+Q2->Q3+Q4*)
	Res1=A1*DiagonalMatrix[vectorPropT] . C1 . DiagonalMatrix[vectorPropT];
	Res2=A2*DiagonalMatrix[vectorPropU] . C2 . DiagonalMatrix[vectorPropU];
	Res3=A3*DiagonalMatrix[vectorPropU] . C3 . DiagonalMatrix[vectorPropT];
	
(*Factor 4 from anti-particle contributions*)
	Return[ Total[Res1+Res2+Res3,-1]]
]
];


(* ::Subsubsection:: *)
(*Q1V1toQ1V1*)


CreateMatrixElementQ1V1toQ1V1[particle1_,particle2_,particle3_,particle4_,vectorMass_,fermionMass_]:=
Block[{s,t,u},
(*
In QCD this process depends on up to
two Lorentz structures if the particles are all the same; and
one Lorentz structure if they are all different.
*)
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
	g13=gvff[[;;,p1,p3]];
	g12=gvff[[p2,p1,;;]]//Flatten[#,{{3},{2},{1}}]&;
	g14=gvff[[p4,p1,;;]]//Flatten[#,{{3},{2},{1}}]&;
	
	g34=gvff[[p4,;;,p3]]//Flatten[#,{{2},{3},{1}}]&;
	g24=gvvv[[p2,p4,;;]]//Flatten[#,{{3},{1},{2}}]&;
	g23=gvff[[p2,;;,p3]]//Flatten[#,{{2},{3},{1}}]&;

	LV=Length[g13];
	LF=Length[g12];

(*Vector propagators*)
	vectorPropT=Table[1/(t-i),{i,vectorMass}];
	
(*Fermion propagators*)
	fermionPropU=Table[1/(u-i),{i,fermionMass}];
	
(*Group invariants that multiply various Lorentz Structures*)
	C1=Sum[vectorPropT[[c]]vectorPropT[[d]]Tr[g13[[c]] . Transpose[g13[[d]]]]Tr[g24[[c]] . Transpose[g24[[d]]]],{c,1,LV},{d,1,LV}];
	C2=Sum[fermionPropU[[L]] fermionPropU[[K]]Tr[g14[[L]] . Transpose[g14[[K]]]]Tr[g23[[L]] . Transpose[g23[[K]]]],{L,1,LF},{K,1,LF}];
	
(*Since there are two diagrams there can be 3 Lorentz structures after squaring and summing over spins*)
(*The interference diagram does however not contribute at leading-log, so I omit it*)
(*They are just hardcoded for now*)
	A1=16(s^2+u^2); (*squared t-channel diagram*)
	A2=-4*16*s*u; (*squared u-channel diagram*)

(*We now generate the matrix element for Q1+Q2->Q3+Q4*)
	Res1=A1*C1;
	Res2=A2*C2;

(*Factor of 2 from anti-particle contribution*)
	Return[(Res1+Res2)]
,
	Return[0]
]	
];


(* ::Subsubsection:: *)
(*Q1Q2toV1V2*)


CreateMatrixElementQ1Q2toV1V2[particle1_,particle2_,particle3_,particle4_,fermionMass_]:=
Block[{s,t,u},
(*
In QCD this process depends on up to
two Lorentz structures if the particles are all the same; and
one Lorentz structure if they are all different.
*)
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
	g13=gvff[[p3,p1,;;]];
	g24=gvff[[p4,;;,p2]];
	g13T=gvff[[p3,;;,p1]];
	g24T=gvff[[p4,p2,;;]];

	g14=gvff[[p4,p1,;;]];
	g23=gvff[[p3,;;,p2]];
	g14T=gvff[[p4,;;,p1]];
	g23T=gvff[[p3,p2,;;]];
	
	LV=Length[g13];
	F=Length[g12[[1]]];

(*Fermion propagators*)
	fermionPropU=Table[1/(u-i),{i,fermionMass}];
	fermionPropT=Table[1/(t-i),{i,fermionMass}];

(*Group invariants that multiply various Lorentz Structures*)
(*Multiplying with propagators for each particle*)
	g13T=DiagonalTensor2[TensorProduct[g13T,fermionPropT],2,4]//Flatten[#,{{2},{1},{3}}]&;
	g13=DiagonalTensor2[TensorProduct[g13,fermionPropT],3,4]//Flatten[#,{{2},{3},{1}}]&;
	
	C1=Sum[Tr[g13[[a]] . g24[[b]] . g24T[[b]] . g13T[[a]]],{a,Length[g13]},{b,Length[g24]}];

	g14T=DiagonalTensor2[TensorProduct[g14T,fermionPropU],2,4]//Flatten[#,{{2},{1},{3}}]&;
	g14=DiagonalTensor2[TensorProduct[g14,fermionPropU],3,4]//Flatten[#,{{2},{3},{1}}]&;

	C2=Sum[Tr[g14[[a]] . g23[[b]] . g23T[[b]] . g14T[[a]]],{a,Length[g14]},{b,Length[g23]}];
	
(*Since there are two diagrams there can be 3 Lorentz structures after squaring and summing over spins*)
(*The interference diagram does however not contribute at leading-log, so I omit it*)
(*They are hardcoded for now*)
	A1=4 t u; (*squared t-channel diagram*)
	A2=4 t u; (*squared u-channel diagram*)

(*Generate the matrix element for Q1+Q2->Q3+Q4*)
	Res1=A1*C1;
	Res2=A2*C2;

(*Factor of 2 from anti-particles*)
	Return[(Res1+Res2)]
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
(*Essentially this is generates the elements going into Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]*)

(*particleList is the complete list of particles*)

(*Generate all matrix elements*)
	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];

(*Divide the incoming particle by its degree's of freedom*)
	MatrixElements=Table[
		1/degreeOfFreedom[particleList[[a]]]
		CreateMatrixElementQ1Q2toQ3Q4[particleList[[a]],b,c,d,GluonMass],
		{a,OutOfEqParticles},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
(*This is a list of all non-zero matrix elements*)
	Elements=MatrixElements["NonzeroPositions"];  

(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,

(*Now we add identical contributions*)
(*The Q1Q2->Q1Q2 process is the same as the Q1Q2->Q2Q1 process*)
	symmetries=Split[Elements,#1[[3;;4]]==Sort[#2[[3;;4]]]&];
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
	symmetries=Split[CollElements,#1[[2]]==#2[[2]]&];
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
(*Essentially this is generates the elements going into Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]*)

(*particleList is the complete list of particles*)

(*First we generate all matrix elements*)
	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];

(*Divide the incoming particle by its degree's of freedom*)
	MatrixElements=Table[
		1/degreeOfFreedom[particleList[[a]]]
		CreateMatrixElementQ1V1toQ1V1[particleList[[a]],b,c,d,GluonMass,QuarkMass],
		{a,OutOfEqParticles},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
	(*This is a list of all non-zero matrix elements*)
	Elements=MatrixElements["NonzeroPositions"];


(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,

(*Now we add identical contributions*)
(*The Q1Q2->Q1Q2 process is the same as the Q1Q2->Q2Q1 process*)
	symmetries=Split[Elements,#1[[3;;4]]==Sort[#2[[3;;4]]]&];
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
	symmetries=Split[CollElements,#1[[2]]==#2[[2]]&];
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
(*Essentially this is generates the elements going into Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]*)

(*particleList is the complete list of particles*)

(*First we generate all matrix elements*)
	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];

(*Divide the incoming particle by its degree's of freedom*)
	MatrixElements=Table[
		1/degreeOfFreedom[particleList[[a]]]
		CreateMatrixElementQ1Q2toV1V2[particleList[[a]],b,c,d,QuarkMass],
		{a,OutOfEqParticles},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
	(*This is a list of all non-zero matrix elements*)
	Elements=MatrixElements["NonzeroPositions"];

(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,

(*Now we add identical contributions*)
(*The Q1Q2->Q1Q2 process is the same as the Q1Q2->Q2Q1 process*)
	symmetries=Split[Elements,#1[[3;;4]]==Sort[#2[[3;;4]]]&];
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
	symmetries=Split[CollElements,#1[[2]]==#2[[2]]&];
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
(*Essentially this is generates the elements going into Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]*)

(*particleList is the complete list of particles*)

(*Generate all matrix elements*)
	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];
	
(*Divide the incoming particle by its degree's of freedom*)
	MatrixElements=Table[
		1/degreeOfFreedom[particleList[[a]]]
		CreateMatrixElementV1V2toV3V4[particleList[[a]],b,c,d,GluonMass],
		{a,OutOfEqParticles},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
	(*This is a list of all non-zero matrix elements*)
	Elements=MatrixElements["NonzeroPositions"];

(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,

(*Now we add identical contributions*)
(*The Q1Q2->Q1Q2 process is the same as the Q1Q2->Q2Q1 process*)
	symmetries=Split[Elements,#1[[3;;4]]==Sort[#2[[3;;4]]]&];
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
	symmetries=Split[CollElements,#1[[2]]==#2[[2]]&];
	CollElements=Table[{symmetries[[i]][[;;,1]]//Total[#,-1]&,symmetries[[i]][[1,2]]},{i,1,Length[symmetries]}];
	
	Return[CollElements]	
]	
];


(* ::Subsubsection::Closed:: *)
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
	OutOfEqParticles,RepMasses,RepCouplings
	},

(*
Extracting the out-of-eq particles
*)
	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];

(*
Replacement rules for converting to the format used by the c++ code
*)	
	RepMasses=Table[UserMasses[[i]]->msq[i-1],{i,1,Length[UserMasses]}];
	RepCouplings=Table[UserCouplings[[i]]->Coupling[i-1],{i,1,Length[UserCouplings]}];
	
(*
Loop over all out-of-eq particles and extracting the matrix elements
*)	
	MatrixElements=ExtractOutOfEqElement[particleList,LightParticles]//FullSimplify;
	
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
	
	CouplingInfo=Table[{ToString[UserCouplings[[i]]],ToString[Coupling[i-1]]},{i,1,Length[UserCouplings]}];
	
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
