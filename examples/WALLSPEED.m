(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
<<../DRalgo.m


(* ::Chapter:: *)
(*QCD*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3"};
RepAdjoint={{1,1}};
RepScalar={};
CouplingName={gs};


Rep1={{{1,0}},"L"};
Rep2={{{1,0}},"R"};
RepFermion1Gen={Rep1,Rep2};


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen,RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


(* ::Section::Closed:: *)
(*Help Functions*)


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


CreateOutOfEq[Indices_,Type_]:=Block[{},
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


(* ::Section::Closed:: *)
(*Matrix elements*)


CreateMatrixElementQ1Q2toQ3Q4[particle1_,particle2_,particle3_,particle4_,vectorMass_]:=Block[{},
(*In QCD this process depends on up to two Lorentz structures if the particles are all the same; and one if they are all different.*)
If[particle1[[2]]!="F"||particle2[[2]]!="F"||particle3[[2]]!="F"||particle4[[2]]!="F",
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

(*Since there are two diagrams there can be 3 Lorentz structures after squaring and summing over spins*)
(*They are just hardcoded for now*)
	A1=2(s^2+u^2); (*squared t-channel diagram*)
	A2=2(s^2+t^2); (*squared u-channel diagram*)
	A3=4 s^2; (*Interference between u and t channel diagrams*)

(*We now generate the matrix element for Q1+Q2->Q3+Q4*)
	Res1=A1*DiagonalMatrix[vectorPropT] . C1 . DiagonalMatrix[vectorPropT];
	Res2=A2*DiagonalMatrix[vectorPropU] . C2 . DiagonalMatrix[vectorPropU];
	Res3=A3*DiagonalMatrix[vectorPropU] . C3 . DiagonalMatrix[vectorPropT];
	
	(*the 4 comes from anti-particle contributions*)
	Return[4 Total[Res1+Res2+Res3,-1]]
]
	
];


CreateMatrixElementQ1V1toQ1V1[particle1_,particle2_,particle3_,particle4_,vectorMass_,fermionMass_]:=Block[{},
(*In QCD this process depends on up to two Lorentz structures if the particles are all the same; and one if they are all different.*)
If[((particle1[[2]]=="F"&&particle2[[2]]=="V")||(particle1[[2]]=="V"&&particle2[[2]]=="F"))
	&&((particle3[[2]]=="F"&&particle4[[2]]=="V")||(particle3[[2]]=="V"&&particle4[[2]]=="F")),

	p1=particle1[[1]];
	p2=particle2[[1]];
	p3=particle3[[1]];
	p4=particle4[[1]];
(*Just changing the order of the particles so that it is always QV->QV*)	
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
	C1=Sum[vectorPropT[[c]]vectorPropT[[d]]Tr[g13[[c]] . Transpose[g13[[d]]]]Tr[g24[[c]] . g24[[d]]],{c,1,LV},{d,1,LV}];
	C2=Sum[fermionPropU[[L]] fermionPropU[[K]]Tr[g14[[L]] . Transpose[g23[[L]]]g23[[K]] . Transpose[g14[[K]]]],{L,1,LF},{K,1,LF}];

(*Since there are two diagrams there can be 3 Lorentz structures after squaring and summing over spins*)
(*The interference diagram does however not contribute at leading-log, so I omit it*)
(*They are just hardcoded for now*)
	A1=-16(s^2+t^2); (*squared t-channel diagram*)
	A2=-4 s u; (*squared u-channel diagram*)

(*We now generate the matrix element for Q1+Q2->Q3+Q4*)
	Res1=A1*C1;
	Res2=A2*C2;

	
	Return[2(Res1+Res2)](*Factor of 2 from anti-particle contribution*)
,
	Return[0]
]	
	
];


CreateMatrixElementQ1Q2toV1V2[particle1_,particle2_,particle3_,particle4_,fermionMass_]:=Block[{},
(*In QCD this process depends on up to two Lorentz structures if the particles are all the same; and one if they are all different.*)

If[(particle1[[2]]!="F"||particle2[[2]]!="F"||particle3[[2]]!="V"||particle4[[2]]!="V")&&(particle1[[2]]!="V"||particle2[[2]]!="V"||particle3[[2]]!="F"||particle4[[2]]!="F"),
	Return[0];
,

	p1=particle1[[1]];
	p2=particle2[[1]];
	p3=particle3[[1]];
	p4=particle4[[1]];
If[particle1[[2]]=="V"&&particle2[[2]]=="V"&&particle3[[2]]=="F"&&particle4[[2]]=="F",
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
	
	C1=Sum[Tr[g13[[a]] . g24[[b]] . g24T[[b]] . g13T[[a]]],{a,LV},{b,LV}];
	
	g14T=DiagonalTensor2[TensorProduct[g14T,fermionPropU],2,4]//Flatten[#,{{2},{1},{3}}]&;
	g14=DiagonalTensor2[TensorProduct[g14,fermionPropU],3,4]//Flatten[#,{{2},{3},{1}}]&;
	
	C2=Sum[Tr[g14[[a]] . g23[[b]] . g23T[[b]] . g14T[[a]]],{a,LV},{b,LV}];
	
(*Since there are two diagrams there can be 3 Lorentz structures after squaring and summing over spins*)
(*The interference diagram does however not contribute at leading-log, so I omit it*)
(*They are just hardcoded for now*)
	A1=4 t u; (*squared t-channel diagram*)
	A2=4 t u; (*squared u-channel diagram*)

(*We now generate the matrix element for Q1+Q2->Q3+Q4*)
	Res1=A1*C1;
	Res2=A2*C2;
	
	Return[2*(Res1+Res2)](*Factor of 2 from anti-particles*)
]
	
];


(* ::Section:: *)
(*Getting the matrix elements for out - of - eq particles*)


degreeOfFreedom[particle_]:=Block[{},

	dof=Length[particle[[1]]];

	If[particle[[2]]=="F",dof*=2*2]; (*Factor of 2 from helicites, factor of 2 from anti-particles*)
	If[particle[[2]]=="V",dof*=2]; (*Factor of 2 from spins*)
	
	Return[dof];
]


ExtractOutOfEqElementQ1Q2toQ3Q4[incomingParticle_,deltaFparticle_,particleList_]:=Block[{},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*Essentially this is generates the elements going into Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]*)

(*particleList is the complete list of particles*)

(*First we generate all matrix elements*)
	MatrixElements=Table[CreateMatrixElementQ1Q2toQ3Q4[a,b,c,d,GluonMass],{a,particleList},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
	Elements=MatrixElements["NonzeroPositions"];  (*This is a list of all non-zero matrix elements*)

(*We now extract all matrix elements where the incoming particle is of "incomingParticle" type *)
	PosIncoming=Position[Elements,{incomingParticle,__}];
	Elements=Extract[Elements,PosIncoming];

(*We now demand that, at least, one of the scattered particles need to be a deltaFparticle*)	
	PosDeltaF=Union[Position[Elements,{_,deltaFparticle,__}],Position[Elements,{_,_,deltaFparticle,_}],Position[Elements,{_,__,deltaFparticle}]];
	Elements=Extract[Elements,PosDeltaF];
	
(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,

(*Now we add identical contributions*)
(*The Q1Q2->Q1Q2 process is the same as the Q1Q2->Q2Q1 process*)
	symmetries=Split[Elements,#1[[3;;4]]==Sort[#2[[3;;4]]]&];
	MultiPlicity=Table[Length[a],{a,symmetries}];
	Elements=Table[i[[1]],{i,symmetries}];
	MatrixElements=Table[MultiPlicity[[i]] Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}]; (*MatrixElements is now just a list*)

(*We also divide by the number of degrees of freedom of incomingParticle*)
	MatrixElements=MatrixElements/degreeOfFreedom[particleList[[incomingParticle]]];

(*We now select all elements where deltaFparticle is amongst the scattered particles*)
(*deltaF is here a list of 1 and 0s*)
	deltaF=(Elements/. x_?NumericQ /; Not@MatchQ[x, deltaFparticle ] -> 0 )/. x_?NumericQ /; MatchQ[x, deltaFparticle ] -> 1;
	
(*Now we create the full list of distinct collision elements*)
	CollElements=Table[{MatrixElements[[i]],deltaF[[i]]},{i,1,Length[Elements]}];

(*We now add all elements with the same deltaF list*)
	symmetries=Split[CollElements,#1[[2]]==#2[[2]]&];
	CollElements=Table[{symmetries[[i]][[;;,1]]//Total[#,-1]&,symmetries[[i]][[1,2]]},{i,1,Length[symmetries]}];
	
	Return[CollElements]
	
]
	
];


ExtractOutOfEqElementQ1V1toQ1V1[incomingParticle_,deltaFparticle_,particleList_]:=Block[{},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*Essentially this is generates the elements going into Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]*)

(*particleList is the complete list of particles*)

(*First we generate all matrix elements*)
	MatrixElements=Table[CreateMatrixElementQ1V1toQ1V1[a,b,c,d,GluonMass,QuarkMass],{a,particleList},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
	Elements=MatrixElements["NonzeroPositions"];  (*This is a list of all non-zero matrix elements*)

(*We now extract all matrix elements where the incoming particle is of "incomingParticle" type *)
	PosIncoming=Position[Elements,{incomingParticle,__}];
	Elements=Extract[Elements,PosIncoming];
	
(*We now demand that, at least, one of the scattered particles need to be a deltaFparticle*)	
	PosDeltaF=Union[Position[Elements,{_,deltaFparticle,__}],Position[Elements,{_,_,deltaFparticle,_}],Position[Elements,{_,__,deltaFparticle}]];
	Elements=Extract[Elements,PosDeltaF];
	
(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,

(*Now we add identical contributions*)
(*The Q1Q2->Q1Q2 process is the same as the Q1Q2->Q2Q1 process*)
	symmetries=Split[Elements,#1[[3;;4]]==Sort[#2[[3;;4]]]&];
	MultiPlicity=Table[Length[a],{a,symmetries}];
	Elements=Table[i[[1]],{i,symmetries}];
	MatrixElements=Table[MultiPlicity[[i]] Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}]; (*MatrixElements is now just a list*)

(*We also divide by the number of degrees of freedom of incomingParticle*)
	MatrixElements=MatrixElements/degreeOfFreedom[particleList[[incomingParticle]]];

(*We now select all elements where deltaFparticle is amongst the scattered particles*)
(*deltaF is here a list of 1 and 0s*)
	deltaF=(Elements/. x_?NumericQ /; Not@MatchQ[x, deltaFparticle ] -> 0 )/. x_?NumericQ /; MatchQ[x, deltaFparticle ] -> 1;
	
(*Now we create the full list of distinct collision elements*)
	CollElements=Table[{MatrixElements[[i]],deltaF[[i]]},{i,1,Length[Elements]}];

(*We now add all elements with the same deltaF list*)
	symmetries=Split[CollElements,#1[[2]]==#2[[2]]&];
	CollElements=Table[{symmetries[[i]][[;;,1]]//Total[#,-1]&,symmetries[[i]][[1,2]]},{i,1,Length[symmetries]}];
	
	Return[CollElements]
	
]
	
];


ExtractOutOfEqElementQ1Q2toV1V2[incomingParticle_,deltaFparticle_,particleList_]:=Block[{},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*Essentially this is generates the elements going into Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]*)

(*particleList is the complete list of particles*)

(*First we generate all matrix elements*)
	MatrixElements=Table[CreateMatrixElementQ1Q2toV1V2[a,b,c,d,QuarkMass],{a,particleList},{b,particleList},{c,particleList},{d,particleList}]//SparseArray;
	Elements=MatrixElements["NonzeroPositions"];  (*This is a list of all non-zero matrix elements*)

(*We now extract all matrix elements where the incoming particle is of "incomingParticle" type *)
	PosIncoming=Position[Elements,{incomingParticle,__}];
	Elements=Extract[Elements,PosIncoming];

(*We now demand that, at least, one of the scattered particles need to be a deltaFparticle*)	
	PosDeltaF=Union[Position[Elements,{_,deltaFparticle,__}],Position[Elements,{_,_,deltaFparticle,_}],Position[Elements,{_,__,deltaFparticle}]];
	Elements=Extract[Elements,PosDeltaF];
	
(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,

(*Now we add identical contributions*)
(*The Q1Q2->Q1Q2 process is the same as the Q1Q2->Q2Q1 process*)
	symmetries=Split[Elements,#1[[3;;4]]==Sort[#2[[3;;4]]]&];
	MultiPlicity=Table[Length[a],{a,symmetries}];
	Elements=Table[i[[1]],{i,symmetries}];
	MatrixElements=Table[MultiPlicity[[i]] Extract[MatrixElements,Elements[[i]]],{i,1,Length[Elements]}]; (*MatrixElements is now just a list*)

(*We also divide by the number of degrees of freedom of incomingParticle*)
	MatrixElements=MatrixElements/degreeOfFreedom[particleList[[incomingParticle]]];

(*We now select all elements where deltaFparticle is amongst the scattered particles*)
(*deltaF is here a list of 1 and 0s*)
	deltaF=(Elements/. x_?NumericQ /; Not@MatchQ[x, deltaFparticle ] -> 0 )/. x_?NumericQ /; MatchQ[x, deltaFparticle ] -> 1;
	
(*Now we create the full list of distinct collision elements*)
	CollElements=Table[{MatrixElements[[i]],deltaF[[i]]},{i,1,Length[Elements]}];

(*We now add all elements with the same deltaF list*)
	symmetries=Split[CollElements,#1[[2]]==#2[[2]]&];
	CollElements=Table[{symmetries[[i]][[;;,1]]//Total[#,-1]&,symmetries[[i]][[1,2]]},{i,1,Length[symmetries]}];
	
	Return[CollElements]
	
]
	
];


ExtractOutOfEqElement[incomingParticle_,deltaFparticle_,particleList_]:=Block[{},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*Essentially this is generates the elements going into Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]*)

(*particleList is the complete list of particles*)


(*First we extract the result for all subprocesses*)
CollEllQ1Q2toQ3Q4=ExtractOutOfEqElementQ1Q2toQ3Q4[incomingParticle,deltaFparticle,particleList];
CollEllQ1V1toQ1V1=ExtractOutOfEqElementQ1V1toQ1V1[incomingParticle,deltaFparticle,particleList];
CollEllQ1Q2toV1V2=ExtractOutOfEqElementQ1Q2toV1V2[incomingParticle,deltaFparticle,particleList];

CollEllTotal=Join[CollEllQ1Q2toQ3Q4,CollEllQ1V1toQ1V1,CollEllQ1Q2toV1V2];


(*I have so far not added all the terms with the same deltaF indices, which we should of course do*)
(*This is trivial to do, and I'll do it after we have agreed what format we want*)

Return[CollEllTotal]

];


(* ::Title:: *)
(*A model with 6 quarks and 1 gluon*)


(*Creating the reps*)


PosScalar=PrintScalarRepPositions[];
PosVector=PrintGaugeRepPositions[];
PosFermion=PrintFermionRepPositions[];


(*In DRalgo fermion are Weyl. So to create one DIrac we need one left-handed and one right-handed fermoon*)


(*Below Reps 1-6 are quarks, and rep 7 is a gluon*)


Rep1=CreateOutOfEq[{1,2},"F"];
Rep2=CreateOutOfEq[{3,4},"F"];
Rep3=CreateOutOfEq[{5,6},"F"];
Rep4=CreateOutOfEq[{7,8},"F"];
Rep5=CreateOutOfEq[{9,10},"F"];
Rep6=CreateOutOfEq[{11,12},"F"];
RepV=CreateOutOfEq[{1},"V"];


ParticleList={Rep1,Rep2,Rep3,Rep4,Rep5,Rep6,RepV};


(*Defining various masses*)


GluonMass=Table[mg2,{i,1,Length[gvff]}];


QuarkMass=Table[mq2,{i,1,Length[gvff[[1]]]}];


(*The deltaC[top,top] contributions*)


ExtractOutOfEqElement[1,1,ParticleList]


(*The deltaC[top,gluon] contributions*)


ExtractOutOfEqElement[1,7,ParticleList]


(*The deltaC[gluon,gluon] contributions*)


ExtractOutOfEqElement[7,7,ParticleList]


(*The deltaC[gluon,top] contributions*)


ExtractOutOfEqElement[7,1,ParticleList]


(*Examples of the scattering elements from various processes. With the top as an external particle, and the top deltaF terms*)


ExtractOutOfEqElementQ1Q2toQ3Q4[1,1,ParticleList]//Expand


ExtractOutOfEqElementQ1V1toQ1V1[1,1,ParticleList]//Expand


ExtractOutOfEqElementQ1Q2toV1V2[1,1,ParticleList]//Expand
