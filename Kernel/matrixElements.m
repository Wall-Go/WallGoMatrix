(* ::Package:: *)

(* ::Title:: *)
(*Matrix elements*)


(* ::Section::Closed:: *)
(*Help functions*)


(*
	Trick to simplify tensors because Mathematica 13 sucks.
*)
SimplifySparse[s_SparseArray] := With[
    {
    elem =Simplify[s["NonzeroValues"]],
    pos=s["NonzeroPositions"],
    dim = Dimensions[s]
    },
SparseArray[pos->elem,dim,0]
    ]



(*
	Creates tensors used in intermediate steps
*)
CreateHelpTensors[]:=Module[{},
	If[verbose,Print["Creating Help Tensors"]];

(*Tensors that are built from two scalar-vector trillinear couplings*)
		Habij=Contract[gvss,gvss,{{3,5}}]//Transpose[#,{1,3,2,4}]&//SimplifySparse;
		Hg=TensorContract[Habij,{{3,4}}]//SimplifySparse;
		\[CapitalLambda]g=TensorContract[Habij,{{1,2}}]//SimplifySparse;            
		HabijV=Habij+Transpose[Habij,{2,1,3,4}]//SparseArray//SimplifySparse;
		
(*Tensor that is built from two structure constants*)
		GabcdV=gvvv . gvvv//SparseArray//SimplifySparse;
	
(*Tensor that is built from two fermion-vector trillinear couplings*)
		HabIJF=Contract[gvff,gvff,{{3,5}}]//Transpose[#,{1,3,2,4}]&//SimplifySparse;
        
(*Tensors that is built from two Yukawa couplings*)
		Ysij=Contract[Ysff,YsffC,{{2, 5},{3,6}}]//SimplifySparse;
		YsijC=Contract[YsffC,Ysff,{{2, 5},{3,6}}]//SimplifySparse;
        
	If[mode>=1,
(*Tensor that is built from two scalar quartics*)
		\[CapitalLambda]\[Lambda] =Flatten[\[Lambda]4,{{1},{2},{3,4}}] . Flatten[\[Lambda]4,{1,2}]//SparseArray//SimplifySparse;

(*Invariant tensors built from two Yukawa couplings*)
		YTemp=Ysff . Transpose[YsffC]//Transpose[#,{1,3,2,4}]&//Transpose[#,{1,2,4,3}]&//SimplifySparse;
		YTempC=YsffC . Transpose[Ysff]//Transpose[#,{1,3,2,4}]&//Transpose[#,{1,2,4,3}]&//SimplifySparse;
        
		Yhelp=Flatten[YTemp,{{1},{2},{3,4}}] . Flatten[YTemp,{4,3}]//SparseArray//SimplifySparse;     
		YhelpC=Flatten[YTempC,{{1},{2},{3,4}}] . Flatten[YTempC,{4,3}]//SparseArray//SimplifySparse;
	];


];


(*Sums all the elements of a tensor and simplifies with the assumption that all variables are real*)
TotalConj[s_] := Simplify[Total[s,-1],Assumptions->VarAsum] 


DiagonalTensor2[s_SparseArray,a_Integer,b_Integer] := With[
    {
    s1=Flatten[s,{{a},{b}}]
    },
     Table[i,{i,s1}]//Table[#[[i,i]],{i,1,Length[#]}]&//SparseArray
    ]


(*Permutes a rank 4 tensor via: Subscript[T, i1 i2 i3 i4]->Subscript[T, a b c d]*)
OrderArray[s_SparseArray,a_Integer,b_Integer,c_Integer,d_Integer] := Flatten[s,{{a},{b},{c},{d}}]


Contract[tensor1_,tensor2_,indices_]:=Activate @ TensorContract[
        Inactive[TensorProduct][tensor1,tensor2], indices] (*Performs tensor contractions on two tensors*)


Contract[tensor1_,tensor2_,tensor3_,indices_]:=Activate @ TensorContract[
        Inactive[TensorProduct][tensor1,tensor2,tensor3], indices](*Performs tensor contractions on three tensors*)


Contract[tensor1_,tensor2_,tensor3_,tensor4_,indices_]:=Activate @ TensorContract[
        Inactive[TensorProduct][tensor1,tensor2,tensor3,tensor4], indices] (*Performs tensor contractions on four tensors*)


ListToMat[list1_] := DiagonalMatrix[list1]//SparseArray (*Takes a 1-d list and turns it into a rank 2 sparse array*)


RangeToIndices[ListI_]:=Block[{},
(*Creates a list of numbers between i and j when ListI=i;;j*)
	Table[i,{i,ListI[[1]],ListI[[2]]}]
];


RepToIndices[ListI_]:=Block[{},
(*Converts a list of range of numbers to a flat array*)
	Table[RangeToIndices[x],{x,ListI}]//Flatten//Sort
];


(* ::Subsubsection::Closed:: *)
(*CreateParticle*)


CreateParticle[Indices_,Type_,Mass_,Name_]:=Block[
{
	Field,FieldPosition,temp,Particle
},
(*Combines selected representations and obtain their indices*)
	Particle={};
	Field=Switch[Type,
		"V","Vector",
		"F","Fermion",
		"S","Scalar"];

	FieldPosition=PrintFieldRepPositions[Field];
	Do[
		If[Length[i]>=2,
			(*if there is more than 1 argument the SymmetryBreaking["X"] split is assumed*)
			AppendTo[Particle,{MassiveReps[Field][[i[[1]]]][[i[[2]]]][[1]]}],
			AppendTo[Particle,{RepToIndices[{FieldPosition[[i]]}]}]
		];
	,{i,Indices}];
	Return[{Flatten[Particle],Type,Mass,Name}]
];


(* ::Subsubsection::Closed:: *)
(*SymmetryBreaking*)


Options[SymmetryBreaking]={
	VevDependentCouplings->False
	};


SymmetryBreakingField["Vector"][Indices_,vev_] :=Module[
{
	PosVector,Habij,massV,gaugeInd,posHeavy,posLight,
	rep,val,val2,pos,pos2
},
(*
	Finds the masses for the gauge-representation Indices.
*)
	PosVector=PrintFieldRepPositions["Vector"];
	
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


SymmetryBreakingField["Scalar"][Indices_,vev_] :=Module[
{
	PosScalar,scalarInd,massS,posHeavy,posLight,rep,val,val2,pos,pos2
},
(*
	Finds the masses for the fermion-representation Indices.
*)
	PosScalar=PrintFieldRepPositions["Scalar"];
	
(*Scalars*)
	massS=(
		+ \[Mu]ij
		+ (vev . \[Lambda]4 . vev/2)
		+ (vev . \[Lambda]3)
		)//Normal//SparseArray;
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
				pos2=Table[posHeavy[[pos[[b]][[1]]]],{b,Position[massS[[posHeavy,posHeavy]]["NonzeroValues"],a]}];
				AppendTo[rep,{pos2,a}];
			,{a,val}];

			AppendTo[rep,{posLight,0}];
			If[posLight=={},rep=Drop[rep,-1]];
	];
	rep
]


SymmetryBreakingField["Fermion"][Indices_,vev_] :=Module[
{
	PosFermion,fermionInd,massF,posHeavy,posLight,rep,val,val2,pos,pos2
},
(*
	Finds the masses for the fermion-representation Indices.
*)
	PosFermion=PrintFieldRepPositions["Fermion"];
	
(*Fermions*)
	massF=(
		+ \[Mu]IJF
		+ (vev . Ysff));
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
	Return[rep]
]


SymmetryBreaking[vev_,OptionsPattern[]]:=Module[
{
	FieldPosition,PosVector,PosFermion,PosScalar,count,(*,
	GaugeMassiveReps,FermionMassiveReps,ScalarMassiveReps*)
	lambda
},	
	(* reset cubics *)
	\[Lambda]3=(\[Lambda]3//Normal)/.Thread[vev->0]//SparseArray;
	
	Do[
		FieldPosition[Field]=PrintFieldRepPositions[Field];
		(* TODO Gauge fix massive reps as internal variables*)
		MassiveReps[Field]=Table[SymmetryBreakingField[Field][i,vev],{i,1,Length[FieldPosition[Field]]}];
		Do[
			If[!NumericQ[Total[MassiveReps[Field][[i]][[;;,2]],-1]],
				Print[Style[StringJoin[Field<>" rep ",ToString[i]," splits into particles with mass squared:"],Bold]];
				count=0;
				Do[
					count++;
					Print[{i,count},":   ", PrintNonPrivate[particle[[2]]] ];
					,{particle,MassiveReps[Field][[i]]}]
			]
		,{i,1,Length[FieldPosition[Field]]}];
	,{Field,{"Vector","Fermion","Scalar"}}
	];
	
	(* 
		Additional terms are added to scalar trillinear coupligs if 
		the VevDependentCouplings option is used
	*)
	If[OptionValue[VevDependentCouplings]==True,
		\[Lambda]3=(
			+\[Lambda]3
			+\[Lambda]4 . vev
		);
	]
]


(* ::Section:: *)
(*Matrix elements*)


(* ::Subsubsection::Closed:: *)
(*V1V2toV3V4*)


CreateMatrixElement["V1V2toV3V4"][particle1_,particle2_,particle3_,particle4_,particleMass_]:=
Block[{
	s,t,u,gTensor,C1,C2,C3,
	vectorMass,
	vectorPropT,vectorPropS,vectorPropU,
	A1,A2,A3,
	totRes
},
(*
	This module returns the squared matrix element of VV->VV summed over all quantum numbers of the incoming particles.
	If the incoming particles are not of the form FV->FV the routine returns 0.
*)
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
	
(*Propagator masses*)
	vectorMass=particleMass[[1]];
				
(*Vector propagators*)
	vectorPropT=Table[Prop[t,i],{i,vectorMass}];
	vectorPropU=Table[Prop[u,i],{i,vectorMass}];
	vectorPropS=Table[Prop[s,i],{i,vectorMass}];

(*Lorentz structures that appear*)
(*
	Here we use the convention that a term like su/t^2->1/4-1/4 (s-u)^2/(t-msq)^2. See hep-ph/0302165.
*)
	A1=16(-1/4)(s-u)^2; (*squared t-channel diagram*)
	A2=16(-1/4)(s-t)^2; (*squared u-channel diagram*)
	A3=16(-1/4)(u-t)^2; (*squared s-channel diagram*)

(* === assembling the full result === *)
	totRes=A1*DiagonalMatrix[vectorPropT] . C1 . DiagonalMatrix[vectorPropT];
	totRes+=A2*DiagonalMatrix[vectorPropU] . C2 . DiagonalMatrix[vectorPropU];
	totRes+=A3*DiagonalMatrix[vectorPropS] . C3 . DiagonalMatrix[vectorPropS];

(*add terms independent on the mandelstam variables*)
	totRes+=-16*((C1+C2+C3)*(1-1/4));
	
	Return[-Total[totRes,-1]]
]
];


(* ::Subsubsection::Closed:: *)
(*F1F2toF3F4*)


CreateMatrixElement["F1F2toF3F4"][particle1_,particle2_,particle3_,particle4_,particleMass_]:=
Block[{
	s,t,u,gTensor,
	partInt,
	CSy,CTy,CUy,CSyC,CTyC,CUyC,CSyCC,CTyCC,CUyCC,
	CS,CT,CU,CSC,CTC,CUC,CSCC,CTCC,CUCC,
	CTrS,CTrT,CTrU,A3,A4,A5,A6,
	Tensor,YTensorC,
	vectorMass,scalarMass,
	scalarPropT,scalarPropU,vectorPropU,vectorPropT,
	C5,C1Y,C2Y,A1,A2,vectorPropS,totRes,scalarPropS,YTensor
},
(*
	This module returns the squared matrix element of FF->FF summed over all quantum numbers of the incoming particles.
	If the incoming particles are not of the form FF->FF the routine returns 0.
*)
If[
	particle1[[2]]!="F"||
	particle2[[2]]!="F"||
	particle3[[2]]!="F"||
	particle4[[2]]!="F",
	Return[0];
,
	Table[partInt[i]={particle1,particle2,particle3,particle4}[[i,1]], {i,1,4}];

(*Propagator masses*)
	vectorMass=particleMass[[1]];
	scalarMass=particleMass[[3]];

(*Vector propagators*)
	vectorPropT=Table[Prop[t,i],{i,vectorMass}]//ListToMat;
	vectorPropU=Table[Prop[u,i],{i,vectorMass}]//ListToMat;
	vectorPropS=Table[Prop[s,i],{i,vectorMass}]//ListToMat;

(*Scalar propagators*)
	scalarPropT=Table[Prop[t,i],{i,scalarMass}]//ListToMat;
	scalarPropU=Table[Prop[u,i],{i,scalarMass}]//ListToMat;
	scalarPropS=Table[Prop[s,i],{i,scalarMass}]//ListToMat;
	
(*Coupling constants that we will need*)
	gTensor=Table[gvff[[;;,Particle1,Particle2]],
		{Particle1,{partInt[1],partInt[2],partInt[3],partInt[4]}},
		{Particle2,{partInt[1],partInt[2],partInt[3],partInt[4]}}];
		
	YTensor=Table[Ysff[[;;,Particle1,Particle2]],
		{Particle1,{partInt[1],partInt[2],partInt[3],partInt[4]}},
		{Particle2,{partInt[1],partInt[2],partInt[3],partInt[4]}}];
		
	YTensorC=Table[YsffC[[;;,Particle1,Particle2]],
		{Particle1,{partInt[1],partInt[2],partInt[3],partInt[4]}},
		{Particle2,{partInt[1],partInt[2],partInt[3],partInt[4]}}];
	
	
(*Group invariants from scalar diagrams*)
	CSy=Contract[YTensor[[1,2]],scalarPropS . YTensor[[3,4]],{{1,4}}]//OrderArray[#,1,2,3,4]&;
	CTy=Contract[YTensor[[1,3]],scalarPropT . YTensor[[2,4]],{{1,4}}]//OrderArray[#,1,3,2,4]&;
	CUy=Contract[YTensor[[1,4]],scalarPropU . YTensor[[2,3]],{{1,4}}]//OrderArray[#,1,3,4,2]&;
	
	CSyC=Contract[YTensor[[1,2]],scalarPropS . YTensorC[[3,4]],{{1,4}}]//OrderArray[#,1,2,3,4]&;
	CTyC=Contract[YTensor[[1,3]],scalarPropT . YTensorC[[2,4]],{{1,4}}]//OrderArray[#,1,3,2,4]&;
	CUyC=Contract[YTensor[[1,4]],scalarPropU . YTensorC[[2,3]],{{1,4}}]//OrderArray[#,1,3,4,2]&;

	CSyCC=Contract[YTensorC[[1,2]],scalarPropS . YTensorC[[3,4]],{{1,4}}]//OrderArray[#,1,2,3,4]&;
	CTyCC=Contract[YTensorC[[1,3]],scalarPropT . YTensorC[[2,4]],{{1,4}}]//OrderArray[#,1,3,2,4]&;
	CUyCC=Contract[YTensorC[[1,4]],scalarPropU . YTensorC[[2,3]],{{1,4}}]//OrderArray[#,1,3,4,2]&;
	
(*Group invariants from vector diagrams*)
	CS=Contract[gTensor[[1,2]],vectorPropS . gTensor[[3,4]],{{1,4}}]//OrderArray[#,1,2,3,4]&;
	CT=Contract[gTensor[[1,3]],vectorPropT . gTensor[[2,4]],{{1,4}}]//OrderArray[#,1,3,2,4]&;
	CU=Contract[gTensor[[1,4]],vectorPropU . gTensor[[2,3]],{{1,4}}]//OrderArray[#,1,3,4,2]&;
	
	CTrS=Contract[gTensor[[2,1]],vectorPropS . gTensor[[4,3]],{{1,4}}]//OrderArray[#,2,1,4,3]&;
	CTrT=Contract[gTensor[[3,1]],vectorPropT . gTensor[[4,2]],{{1,4}}]//OrderArray[#,2,4,1,3]&;
	CTrU=Contract[gTensor[[4,1]],vectorPropU . gTensor[[3,2]],{{1,4}}]//OrderArray[#,2,4,3,1]&;

	CSC=Contract[gTensor[[1,2]],vectorPropS . gTensor[[4,3]],{{1,4}}]//OrderArray[#,1,2,4,3]&;
	CTC=Contract[gTensor[[1,3]],vectorPropT . gTensor[[4,2]],{{1,4}}]//OrderArray[#,1,4,2,3]&;
	CUC=Contract[gTensor[[1,4]],vectorPropU . gTensor[[3,2]],{{1,4}}]//OrderArray[#,1,4,3,2]&;
	
	CSCC=Contract[gTensor[[2,1]],vectorPropS . gTensor[[3,4]],{{1,4}}]//OrderArray[#,2,1,3,4]&;
	CTCC=Contract[gTensor[[3,1]],vectorPropT . gTensor[[2,4]],{{1,4}}]//OrderArray[#,2,3,1,4]&;
	CUCC=Contract[gTensor[[4,1]],vectorPropU . gTensor[[2,3]],{{1,4}}]//OrderArray[#,2,3,4,1]&;
	
	
	C5=Table[
		Tr[gTensor[[1,3]][[a]] . gTensor[[3,4]][[b]] . gTensor[[4,2]][[a]] . gTensor[[2,1]][[b]]],
		{a,1,Length[gTensor[[1,3]]]},{b,1,Length[gTensor[[1,3]]]}];
	C5+=Table[
		Tr[gTensor[[3,1]][[a]] . gTensor[[1,2]][[b]] . gTensor[[2,4]][[a]] . gTensor[[4,3]][[b]]],
		{a,1,Length[gTensor[[1,3]]]},{b,1,Length[gTensor[[1,3]]]}];		
			
				
(*Lorentz structures that appear in squared vector-exchange diagrams*)
	A1=2(s^2+u^2); (*squared t-channel diagram*)
	A2=2(s^2+t^2); (*squared u-channel diagram*)
	A3=4*s^2; (*Interference between u and t channel diagrams*)
	A4=2(t^2+u^2); (*Squared s-channel*)
	A5=4*u^2; (*Interference between s and t channel diagrams*)
	A6=4*t^2; (*Interference between s and u channel diagrams*)		
	
(* === assembling the full result === *)		
(*Result for squared vector-exchange diagrams*)
	totRes=A1*Total[CT CTrT,-1]; 
	totRes+=A2*Total[CU CTrU,-1]; 
	totRes+=1/2A3*Total[CT CTrU+CU CTrT,-1];
	totRes+=A4*Total[CS CTrS,-1]; 
	totRes+=1/2*A5 Total[(vectorPropS . C5 . vectorPropT),-1];
	totRes+=1/2*A6 Total[(vectorPropS . C5 . vectorPropU),-1];

(*Result for squared scalar-exchange diagrams*)	
	totRes+=  t^2/4*TotalConj[CTy Conjugate[CTy]];
	totRes+=  u^2/4*TotalConj[CUy Conjugate[CUy]];
	totRes+=  s^2/4*TotalConj[CSy Conjugate[CSy]];

	totRes+=  t^2/2*TotalConj[CTyC Conjugate[CTyC]];
	totRes+=  u^2/2*TotalConj[CUyC Conjugate[CUyC]];
	totRes+=  s^2/2*TotalConj[CSyC Conjugate[CSyC]];	
	
	totRes+=  t^2/4*TotalConj[CTyCC Conjugate[CTyCC]];
	totRes+=  u^2/4*TotalConj[CUyCC Conjugate[CUyCC]];
	totRes+=  s^2/4*TotalConj[CSyCC Conjugate[CSyCC]];
	
	
	totRes+=-  1/8(t^2-s^2+u^2)*TotalConj[CUy Conjugate[CTy]  +CTy Conjugate[CUy]];
	totRes+=-  1/8(s^2-t^2+u^2)*TotalConj[CUy Conjugate[CSy]  +CSy Conjugate[CUy]];
	totRes+=-  1/8(t^2-u^2+s^2)*TotalConj[CSy Conjugate[CTy]  +CTy Conjugate[CSy]];

	totRes+=-  1/8(t^2-s^2+u^2)*TotalConj[CUyCC Conjugate[CTyCC]  +CTyCC Conjugate[CUyCC]];
	totRes+=-  1/8(s^2-t^2+u^2)*TotalConj[CUyCC Conjugate[CSyCC]  +CSyCC Conjugate[CUyCC]];
	totRes+=-  1/8(t^2-u^2+s^2)*TotalConj[CSyCC Conjugate[CTyCC]  +CTyCC Conjugate[CSyCC]];	
	
(*
	Result for interfaces between vector and scalar diagrams---
	Only cross terms between diagrams can contribute
*)	
	totRes+=(-2*t*t)*1/2*TotalConj[CS*Conjugate[CTyC] + CTyC*Conjugate[CS]];
	totRes+=(-2*t*t)*1/2*TotalConj[CUC*Conjugate[CTyC] + CTyC*Conjugate[CUC]];

	totRes+=(-2*u*u)*1/2*TotalConj[CSC*Conjugate[CUyC] + CUyC*Conjugate[CSC]];
	totRes+=(-2*u*u)*1/2*TotalConj[CTC*Conjugate[CUyC] + CUyC*Conjugate[CTC]];

	totRes+=(-2*s*s)*1/2*TotalConj[CT*Conjugate[CSyC] + CSyC*Conjugate[CT]];
	totRes+=(-2*s*s)*1/2*TotalConj[CU*Conjugate[CSyC] + CSyC*Conjugate[CU]];
	
	Return[Refine[4*totRes,Assumptions->VarAsum]]
]
];


(* ::Subsubsection::Closed:: *)
(*F1F2toV1V2*)


CreateMatrixElement["F1F2toV1V2"][particle1_,particle2_,particle3_,particle4_,particleMass_]:=
Block[{
	s,t,u,gTensorV,gTensorF,gTensorF2,
	vectorMass,fermionMass,
	particleNull,
	vectorPropS,fermionPropU,fermionPropT,
	CS,CTrS,
	CSV,CTF,CUF,
	fermionIdentity,
	C1,C2,C3,C4,C5,C6,partInt,
	A1,A2,A3,
	Res1,Res2,Res3,Res4
},
(*
	This module returns the squared matrix element of FF->VV summed over all quantum numbers of the incoming particles.
	If the incoming particles are not of the form FF->VV the routine returns 0.
*)

If[
	(particle1[[2]]!="F"||particle2[[2]]!="F"||particle3[[2]]!="V"||particle4[[2]]!="V")&&
	(particle1[[2]]!="V"||particle2[[2]]!="V"||particle3[[2]]!="F"||particle4[[2]]!="F"),
	Return[0];
,
	Table[partInt[i]={particle1,particle2,particle3,particle4}[[i,1]], {i,1,4}];
If[
	particle1[[2]]=="V"&&
	particle2[[2]]=="V"&&
	particle3[[2]]=="F"&&
	particle4[[2]]=="F",
	
(*Just changing the order of the particles so that it is always FF->VV*)
	{partInt[1],partInt[3]} = {partInt[3],partInt[1]};
	{partInt[2],partInt[4]} = {partInt[4],partInt[2]};
];

	particleNull={};
	(*Simplify the ordering of particles*)
	(*Coupling constants that we will need*)
	gTensorV=Table[gvvv[[;;,Particle1,Particle2]],
		{Particle1,{particleNull,particleNull,partInt[3],partInt[4]}},
		{Particle2,{particleNull,particleNull,partInt[3],partInt[4]}}
		];
	
	gTensorF=Table[gvff[[;;,Particle1,Particle2]],
		{Particle1,{partInt[1],partInt[2],particleNull,particleNull}},
		{Particle2,{partInt[1],partInt[2],particleNull,particleNull}}
		];
		
	gTensorF2=Table[gvff[[Particle1,Particle2,;;]],
		{Particle1,{partInt[3],partInt[4],particleNull,particleNull}},
		{Particle2,{particleNull,particleNull,partInt[1],partInt[2]}}
		];
	
(*Propagator masses*)
	vectorMass=particleMass[[1]];
	fermionMass=particleMass[[2]];
		
(*Fermion propagators*)
	fermionPropU=Table[Prop[u,i],{i,fermionMass}]//ListToMat;
	fermionPropT=Table[Prop[t,i],{i,fermionMass}]//ListToMat;
	
(*Vector propagators*)
	vectorPropS=Table[Prop[s,i],{i,vectorMass}]//ListToMat;
	
	fermionIdentity=Table[1,{i,fermionMass}]; 
	
(*	Print[gTensorF[[1,2]]//MatrixForm];
	Print[gTensorV[[3,4]]//MatrixForm];*)
	
(*Group invariants from vector diagrams*)
	CSV=Contract[gTensorF[[1,2]], vectorPropS . gTensorV[[3,4]],{{1,4}}]//OrderArray[#,1,2,3,4]&;
	
(*Group invariants from Fermion diagrams*)	
	CTF=Contract[gTensorF2[[1,3]] . fermionPropT , gTensorF2[[2,4]],{{3,6}}]//OrderArray[#,1,3,2,4]&;
	CUF=Contract[gTensorF2[[1,4]] . fermionPropU , gTensorF2[[2,3]],{{3,6}}]//OrderArray[#,1,3,4,2]&;

(*Since there are two diagrams there can be 3 Lorentz structures after squaring and summing over spins,
	so naturally the final result contains 4 Lorentz structures for consistency.*)
	A1=4*t*u; (*squared t-channel diagram*)
	A2=4*t*u; (*squared u-channel diagram*)
	A3=4*(t^2+u^2); (*squared s-channel diagram*)

(* === assembling the full result === *)
	Res1=+A1*TotalConj[CTF*Conjugate[CTF]];
	Res2=+A2*TotalConj[CUF*Conjugate[CUF]];
	Res3=-A3*TotalConj[CSV*Conjugate[CSV]];
	
	Return[2*(Res1+Res2+Res3)](*Factor of 2 from anti-particles*)
]	
];


CreateMatrixElementOld["F1F2toV1V2"][particle1_,particle2_,particle3_,particle4_,particleMass_]:=
Block[{
	s,t,u,gTensor,gTensorT,
	vectorMass,fermionMass,
	vectorPropS,gVTensor,fermionPropU,fermionPropT,CS,CTrS,
	fermionIdentity,
	C1,C2,C3,C4,C5,C6,partInt,
	A1,A2,A3,
	Res1,Res2,Res3,Res4
},
(*
	This module returns the squared matrix element of FF->VV summed over all quantum numbers of the incoming particles.
	If the incoming particles are not of the form FF->VV the routine returns 0.
*)

If[
	(particle1[[2]]!="F"||particle2[[2]]!="F"||particle3[[2]]!="V"||particle4[[2]]!="V")&&
	(particle1[[2]]!="V"||particle2[[2]]!="V"||particle3[[2]]!="F"||particle4[[2]]!="F"),
	Return[0];
,
	Table[partInt[i]={particle1,particle2,particle3,particle4}[[i,1]], {i,1,4}];
If[
	particle1[[2]]=="V"&&
	particle2[[2]]=="V"&&
	particle3[[2]]=="F"&&
	particle4[[2]]=="F",
	
(*Just changing the order of the particles so that it is always FF->VV*)
	{partInt[1],partInt[3]} = {partInt[3],partInt[1]};
	{partInt[2],partInt[4]} = {partInt[4],partInt[2]};
];

(*Coupling constants that we will need*)
	gTensor[1,3]=gvff[[partInt[3],partInt[1],;;]];
	gTensor[3,1]=gvff[[partInt[3],;;,partInt[1]]];
	
	gTensor[2,4]=gvff[[partInt[4],partInt[2],;;]];
	gTensor[4,2]=gvff[[partInt[4],;;,partInt[2]]];

	gTensor[1,4]=gvff[[partInt[4],partInt[1],;;]];
	gTensor[4,1]=gvff[[partInt[4],;;,partInt[1]]];
	
	gTensor[2,3]=gvff[[partInt[3],partInt[2],;;]];
	gTensor[3,2]=gvff[[partInt[3],;;,partInt[2]]];
	
(*Propagator masses*)
	vectorMass=particleMass[[1]];
	fermionMass=particleMass[[2]];
		
(*Fermion propagators*)
	fermionPropU=Table[Prop[u,i],{i,fermionMass}];
	fermionPropT=Table[Prop[t,i],{i,fermionMass}];
	
(*Vector propagators*)
	vectorPropS=Table[Prop[s,i],{i,vectorMass}];
	
	fermionIdentity=Table[1,{i,fermionMass}]; 

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

	C4=Tr[
		DiagonalMatrix[fermionIdentity] . Sum[vectorPropS[[a]]gTensor[3,1][[a]] . gTensor[1,3][[a]],{a,Length[gTensor[1,3]]}] .
		DiagonalMatrix[fermionIdentity] . Sum[vectorPropS[[b]]gTensor[4,2][[b]] . gTensor[2,4][[b]],{b,Length[gTensor[2,4]]}]];
	C5=Tr[
		DiagonalMatrix[fermionIdentity] . Sum[vectorPropS[[a]] gTensor[4,1][[a]] . gTensor[1,4][[a]],{a,Length[gTensor[1,4]]}] .
		DiagonalMatrix[fermionIdentity] . Sum[vectorPropS[[b]]gTensor[3,2][[b]] . gTensor[2,3][[b]],{b,Length[gTensor[2,3]]}]];
		
	C6=Sum[
		Tr[vectorPropS[[a]]vectorPropS[[b]]gTensor[4,1][[b]] . gTensor[1,3][[a]] . gTensor[4,2][[b]] . gTensor[2,3][[a]]],
		{a,Length[gTensor[1,3]]},{b,Length[gTensor[4,2]]}];
	C6+=Sum[
		Tr[vectorPropS[[a]]vectorPropS[[b]]gTensor[3,1][[b]] . gTensor[1,4][[a]] . gTensor[3,2][[b]] . gTensor[2,4][[a]]],
		{a,Length[gTensor[1,4]]},{b,Length[gTensor[3,2]]}];

(*Since there are two diagrams there can be 3 Lorentz structures after squaring and summing over spins,
	so naturally the final result contains 4 Lorentz structures for consistency.*)
	A1=4*t*u; (*squared t-channel diagram*)
	A2=4*t*u; (*squared u-channel diagram*)
	A3=4*(t^2+u^2); (*squared s-channel diagram*)

(* === assembling the full result === *)
	Res1=A1*C1;
	Res2=A2*C2;
	Res3=A3*C6;
	Res3+=-8*t^2 C4;
	Res3+=-8*u^2 C5;
	
	If[
	(Length[gTensor[1,3]]== 3 && Length[gTensor[4,2]]==3) && (Length[gTensor[1,4]]==3 && Length[gTensor[3,2]]==3),
	Print[{particleMass,vectorMass,vectorPropS}];
	Print["Lengths gTensor ([1,3],[4,2],[1,4],[3,2]) ", Length[gTensor[1,3]],Length[gTensor[4,2]],Length[gTensor[1,4]],Length[gTensor[3,2]]];
	Print[Res3];
	];

	
	
	Return[2*(Res1+Res2+Res3)](*Factor of 2 from anti-particles*)
]	
];


(* ::Subsubsection::Closed:: *)
(*F1V1toF1V1*)


CreateMatrixElement["F1V1toF1V1"][particle1_,particle2_,particle3_,particle4_,particleMass_]:=
Block[{
	s,t,u,gTensor,resTot,
	t1,s1,u1,partInt,temp,kinFlip
},
(*
	This module returns the squared matrix element of FV->FV summed over all quantum numbers of the incoming particles.
	If the incoming particles are not of the form FV->FV the routine returns 0.
*)
	kinFlip=0; (*records if the external or incoming legs were swapped, and does a t<->u transformation for each swap*)
If[ (
	(particle1[[2]]=="F"&&particle2[[2]]=="V")||
	(particle1[[2]]=="V"&&particle2[[2]]=="F")
	)&&(
	(particle3[[2]]=="F"&&particle4[[2]]=="V")||
	(particle3[[2]]=="V"&&particle4[[2]]=="F")),

	Table[partInt[i]={particle1,particle2,particle3,particle4}[[i]], {i,1,4}];

(*Changing the order of the particles so that it is always FV->FV*)	
	If[particle1[[2]]=="V"&&particle2[[2]]=="F",
		{partInt[1],partInt[2]} = {partInt[2],partInt[1]};
		kinFlip+=1;
	];
	If[particle3[[2]]=="V"&&particle4[[2]]=="F",
		{partInt[3],partInt[4]} = {partInt[4],partInt[3]};
		kinFlip+=1;
	];
	
(*The result for FV->FV is the same as  FF->VV if we do some renaming of the mandelstam variables*)
(*Minus sign as we reversed the momentum of a fermion*)
	resTot=-CreateMatrixElement["F1F2toV1V2"][partInt[1],partInt[3],partInt[2],partInt[4],particleMass]/.{s->s1,t->t1,u->u1};
	resTot=resTot/.{s1->t,t1->s,u1->u};
	
	If[Mod[kinFlip,2]==1,resTot=resTot/.{t->t1,u->u1}/.{t1->u,u1->t};];
	Return[resTot]	
,
	Return[0]
]	
];


(* ::Subsubsection::Closed:: *)
(*S1S2toS3S4*)


CreateMatrixElement["S1S2toS3S4"][particle1_,particle2_,particle3_,particle4_,particleMass_]:=
Block[{
	s,t,u,gTensor,leadingLog,\[Lambda]3Tensor,\[Lambda]4Tensor,
	vectorMass,scalarMass,
	vectorPropT,scalarPropT,vectorPropU,scalarPropU,Particle3,Particle4,
	AS,AU,AT,scalarPropS,vectorPropS,
	TotRes,CS,CT,CU,CS\[Lambda],CT\[Lambda],CU\[Lambda],C\[Lambda]
},
(*
	This module returns the squared matrix element of SS->SS summed over all quantum numbers of the incoming particles.
	If the incoming particles are not of the form SS->SS the routine returns 0.
*)
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
	
	\[Lambda]4Tensor=Table[\[Lambda]4[[Particle1[[1]],Particle2[[1]],Particle3[[1]],Particle4[[1]]]],
		{Particle1,{particle1,particle2,particle3,particle4}},
		{Particle2,{particle1,particle2,particle3,particle4}},
		{Particle3,{particle1,particle2,particle3,particle4}},
		{Particle4,{particle1,particle2,particle3,particle4}}];
		
(*Propagator masses*)
	vectorMass=particleMass[[1]];
	scalarMass=particleMass[[3]];		

(*Vector propagators*)
	vectorPropT=Table[Prop[t,i],{i,vectorMass}]//ListToMat;
	vectorPropU=Table[Prop[u,i],{i,vectorMass}]//ListToMat;
	vectorPropS=Table[Prop[s,i],{i,vectorMass}]//ListToMat;

(*Scalar propagators*)
	scalarPropT=Table[Prop[t,i],{i,scalarMass}]//ListToMat;
	scalarPropU=Table[Prop[u,i],{i,scalarMass}]//ListToMat;
	scalarPropS=Table[Prop[s,i],{i,scalarMass}]//ListToMat;

(*Lorentz structures for vector exchanges*)
	AS=(t-u); (* s-channel diagram*)
	AT=(s-u); (* t-channel diagram*)
	AU=(s-t); (* u-channel diagram*)
	
(*Group invariants from vector diagrams*)
	CS=AS*Contract[gTensor[[1,2]],vectorPropS . gTensor[[3,4]],{{1,4}}]//OrderArray[#,1,2,3,4]&;
	CT=AT*Contract[gTensor[[1,3]],vectorPropT . gTensor[[2,4]],{{1,4}}]//OrderArray[#,1,3,2,4]&;
	CU=AU*Contract[gTensor[[1,4]],vectorPropU . gTensor[[2,3]],{{1,4}}]//OrderArray[#,1,3,4,2]&;
	
(*Group invariants from scalar diagrams. Kos, some say Kosm...*)	

	CS\[Lambda]= Contract[\[Lambda]3Tensor[[1,2]] , scalarPropS . \[Lambda]3Tensor[[3,4]],{{1,4}}]//OrderArray[#,1,2,3,4]&;
	CT\[Lambda]= Contract[\[Lambda]3Tensor[[1,3]] , scalarPropT . \[Lambda]3Tensor[[2,4]],{{1,4}}]//OrderArray[#,1,3,2,4]&;
	CU\[Lambda]= Contract[\[Lambda]3Tensor[[1,4]] , scalarPropU . \[Lambda]3Tensor[[2,3]],{{1,4}}]//OrderArray[#,1,3,4,2]&;
	C\[Lambda]=\[Lambda]4[[particle1[[1]],particle2[[1]],particle3[[1]],particle4[[1]]]];

(*Total contribution; as the external legs have no helicities we can just sum and square the amplitudes*)		
	TotRes=Total[(CS\[Lambda]+CT\[Lambda]+CU\[Lambda]+C\[Lambda]+CS+CT+CU) Conjugate[(CS\[Lambda]+CT\[Lambda]+CU\[Lambda]+C\[Lambda]+CS+CT+CU)],-1]; (*total amplitude-squared contribution from all diagrams*)
	
	Return[Refine[TotRes,Assumptions->VarAsum]]
]
];


(* ::Subsubsection::Closed:: *)
(*S1S2toF1F2*)


CreateMatrixElement["S1S2toF1F2"][particle1_,particle2_,particle3_,particle4_,particleMass_]:=
Block[{
	types,
	s,t,u,gTensor,gTensorT,gTensorC,gTensorCT,gTensorVF,gTensorVFC,gTensorVS,
	Res,CTT,CUU,ATT,AUU,p1,p2,p3,p4,gTensorS,gTensorF,YTensor,YTensorC,particleNull,
	\[Lambda]3Tensor,YTensor2,YTensor2C,
	vectorMass,fermionMass,scalarMass,
	vectorPropS,scalarPropS,fermionPropT,fermionPropU,
	CSV,CSVC,
	CSy,CTy,CUy,CTyC,CUyC,CTyCC,CUyCC,
	ASSY,A,TotRes,temp
},
(*
	This module returns the squared matrix element of SS->FF summed over all quantum numbers of the incoming particles.
	If the incoming particles are not of the form SS->FF the routine returns 0.
*)
	types = #[[2]] & /@ {particle1, particle2, particle3, particle4};
	{p1, p2, p3, p4} = #[[1]] & /@ {particle1, particle2, particle3, particle4};
  
	If[!(
		types === {"S", "S", "F", "F"} ||
	    types === {"F", "F", "S", "S"}),
		Return[0];
	];

	(* Change the order of the particles so that it is always SS->FF. *)	
	If[types === {"F", "F", "S", "S"},
		{p1, p3} = {p3, p1};
		{p2, p4} = {p4, p2};
	];
	
	particleNull={};
	(*Simplify the ordering of particles*)
	(*Coupling constants that we will need*)
	gTensorS=Table[gvss[[;;,Particle1,Particle2]],
		{Particle1,{p1,p2,particleNull,particleNull}},
		{Particle2,{p1,p2,particleNull,particleNull}}];
	
	gTensorF=Table[gvff[[;;,Particle1,Particle2]],
		{Particle1,{particleNull,particleNull,p3,p4}},
		{Particle2,{particleNull,particleNull,p3,p4}}];
		
	YTensor=Table[Ysff[[;;,Particle1,Particle2]],
		{Particle1,{particleNull,particleNull,p3,p4}},
		{Particle2,{particleNull,particleNull,p3,p4}}];
		
	YTensorC=Table[YsffC[[;;,Particle1,Particle2]],
		{Particle1,{particleNull,particleNull,p3,p4}},
		{Particle2,{particleNull,particleNull,p3,p4}}];
		
	\[Lambda]3Tensor=Table[\[Lambda]3[[;;,Particle1,Particle2]],
		{Particle1,{p1,p2,particleNull,particleNull}},
		{Particle2,{p1,p2,particleNull,particleNull}}];
		
	YTensor2=Table[Ysff[[Particle1,Particle2,;;]],
		{Particle1,{p1,p2,particleNull,particleNull}},
		{Particle2,{particleNull,particleNull,p3,p4}}];

	YTensor2C=Table[YsffC[[Particle1,Particle2,;;]],
		{Particle1,{p1,p2,particleNull,particleNull}},
		{Particle2,{particleNull,particleNull,p3,p4}}];

(*Propagator masses*)
	vectorMass=particleMass[[1]];
	fermionMass=particleMass[[2]];
	scalarMass=particleMass[[3]];

(*Vector propagators*)
	vectorPropS=Table[Prop[s,i],{i,vectorMass}]//ListToMat;

(*Scalar propagators*)
	scalarPropS=Table[Prop[s,i],{i,scalarMass}]//ListToMat;

(*Fermion propagators*)
	fermionPropT=Table[Prop[t,i],{i,fermionMass}]//ListToMat;
	fermionPropU=Table[Prop[u,i],{i,fermionMass}]//ListToMat;
	
(*Group invariants from vector diagrams*)
	CSV=Contract[gTensorS[[1,2]], vectorPropS . gTensorF[[3,4]],{{1,4}}]//OrderArray[#,1,2,3,4]&;
	CSVC=Contract[gTensorS[[1,2]], vectorPropS . gTensorF[[4,3]],{{1,4}}]//OrderArray[#,1,2,4,3]&;
	
(*Group invariants from Yukawa diagrams*)
	CSy=Contract[\[Lambda]3Tensor[[1,2]], scalarPropS . YTensor[[3,4]],{{1,4}}]//OrderArray[#,1,2,3,4]&; 
	
	CTy=Contract[YTensor2[[1,3]] . fermionPropT , YTensor2[[2,4]],{{3,6}}]//OrderArray[#,1,3,2,4]&;
	CUy=Contract[YTensor2[[1,4]] . fermionPropU , YTensor2[[2,3]],{{3,6}}]//OrderArray[#,1,3,4,2]&;
	
	CTyC=Contract[YTensor2[[1,3]] . fermionPropT , YTensor2C[[2,4]],{{3,6}}]//OrderArray[#,1,3,2,4]&;
	CUyC=Contract[YTensor2[[1,4]] . fermionPropU , YTensor2C[[2,3]],{{3,6}}]//OrderArray[#,1,3,4,2]&;
	
	CTyCC=Contract[YTensor2C[[1,3]] . fermionPropT , YTensor2C[[2,4]],{{3,6}}]//OrderArray[#,1,3,2,4]&;
	CUyCC=Contract[YTensor2C[[1,4]] . fermionPropU , YTensor2C[[2,3]],{{3,6}}]//OrderArray[#,1,3,4,2]&;

(*Lorentz structures for vector exchanges*)
(*
	Due to angular momentum conservation the SS->S->FF diagram does not mix with the other channels.
*)
	ASSY=s; (*Squared SS->S->FF diagram*)
	A=t*u; (*All other squared diagrams are proportional to u*t *)
	
(* === assembling the full result === *)
	(*SS->S->FF squared*)
	TotRes=ASSY * TotalConj[CSy*Conjugate[CSy]];
	
	(*Squared Yukawa diagrams with fermion exchanges*)
	TotRes+=+ A * TotalConj[CTy*Conjugate[CTy]]; (*squared t-channel diagram*)
	TotRes+=+ A * TotalConj[CUy*Conjugate[CUy]]; (*squared u-channel diagram*)
	TotRes+=- A * TotalConj[CUy*CTy+CTy*CUy]; (*mixed u & t-channel diagrams*)
	
	(*Squared s-channel diagram with a vector boson*)
	TotRes+= 4*A*TotalConj[CSV*Conjugate[CSV]]; 
	
	(*Mix between vector- and fermion-exchange diagrams*)
	TotRes+=+2*t*u*I*TotalConj[CSV*Conjugate[CTyC] - CTyC*Conjugate[CSV]];
	TotRes+=+2*t*u*I*TotalConj[CSVC*Conjugate[CUyC] - CUyC*Conjugate[CSVC]];

(*The full result*)
	Return[2*Refine[TotRes,Assumptions->VarAsum]] (*factor of 2 from anti-particles*)
];


(* ::Subsubsection::Closed:: *)
(*F1S1toF1S1*)


CreateMatrixElement["F1S1toF1S1"][particle1_,particle2_,particle3_,particle4_,particleMass_]:=
Block[{s,t,u,partInt,temp,resTot,u1,s1,t1,kinFlip},
(*
	This module returns the squared matrix element of FV->FV summed over all quantum numbers of the incoming particles.
	If the incoming particles are not of the form FV->FV the routine returns 0.
*)
	kinFlip=0; (*records if the external or incoming legs were swapped, and does a t<->u transformation for each swap*)
If[ (
	(particle1[[2]]=="F"&&particle2[[2]]=="S")||
	(particle1[[2]]=="S"&&particle2[[2]]=="F")
	)&&(
	(particle3[[2]]=="F"&&particle4[[2]]=="S")||
	(particle3[[2]]=="S"&&particle4[[2]]=="F")),

	Table[partInt[i]={particle1,particle2,particle3,particle4}[[i]], {i,1,4}];
(*
	Changing the order of the particles so that it is always SF->SF.
	Incidentally the opposite ordering as given in the name of the routine.
*)	
	If[particle1[[2]]=="F"&&particle2[[2]]=="S",
		{partInt[1], partInt[2]} = {partInt[2], partInt[1]};
		kinFlip+=1;
	];
	If[particle3[[2]]=="F"&&particle4[[2]]=="S",
		{partInt[3], partInt[4]} = {partInt[4], partInt[3]};
		kinFlip+=1;
	];
	
	(*The result for SF->SF is the same as  SS->FF if we do some renaming of the mandelstam variables*)
	(*Minus sign as we reveresed the momentum of a fermion*)
	resTot=-CreateMatrixElement["S1S2toF1F2"][partInt[1],partInt[3],partInt[2],partInt[4],particleMass]/.{s->s1,t->t1,u->u1};
	resTot=resTot/.{s1->t,t1->s,u1->u};
	
	If[Mod[kinFlip,2]==1,resTot=resTot/.{t->t1,u->u1}/.{t1->u,u1->t};];
	Return[Refine[resTot,Assumptions->VarAsum]]
,
	Return[0]
]	
];


(* ::Subsubsection::Closed:: *)
(*F1S1toF1V1*)


SortF1S1toF1V1[L_]:=
Block[{
	helpList,ordering,kinFlip
},
	helpList=L[[;;,2]];
	ordering={1,2,3,4};
	kinFlip=0;
	If[helpList[[1]]=="V"||helpList[[2]]=="V",ordering={3,4,1,2};];
	If[helpList[[ordering[[1]]]]=="S",
		ordering[[{1,2}]]=ordering[[{2,1}]];
		kinFlip+=1;
		];
	If[helpList[[ordering[[3]]]]=="V",
		ordering[[{3,4}]]=ordering[[{4,3}]];
		kinFlip+=1;
		];
	
	Return[{L[[ordering]],kinFlip}]
]


CreateMatrixElement["F1S1toF1V1"][particle1_,particle2_,particle3_,particle4_,particleMass_]:=
Block[{
	s,t,u,gTensor,gTensorC,gTensorVFC,gTensorVF,gTensorVS,SortHelp,
	fermionMass,
	fermionPropS,particleNull,gTensorF,YTensor,CS,resTot,
	YTensor2,gTensorF2,fermionPropU,CU,kinFlip,t2,u2,
	partInt
},
(*
	This module returns the squared matrix element of FS->FV summed over all quantum numbers of the incoming particles.
	If the incoming particles are not of the form FSFVFF the routine returns 0.
*)
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
	
(*Changing the order of the particles so that it is always FS->FV*)	
	{SortHelp,kinFlip}=Join[{particle1},{particle2},{particle3},{particle4}]//SortF1S1toF1V1[#]&;
	partInt[1]=SortHelp[[1]][[1]];
	partInt[2]=SortHelp[[2]][[1]];
	partInt[3]=SortHelp[[3]][[1]];
	partInt[4]=SortHelp[[4]][[1]];

(* Coupling constants needed for the calculation *)
	YTensor=Ysff[[partInt[2],partInt[1],;;]];
	YTensor2=Ysff[[partInt[2],partInt[3],;;]];
	
	gTensorF=gvff[[partInt[4],partInt[3],;;]];
	gTensorF2=gvff[[partInt[4],partInt[1],;;]];
	
(*Propagator masses*)
	fermionMass=particleMass[[2]];
	
(*Fermion propagators*)
	fermionPropS=Table[Prop[s,i],{i,fermionMass}]//ListToMat;
	fermionPropU=Table[Prop[u,i],{i,fermionMass}]//ListToMat;

(*Group structures*)
	CS=Contract[YTensor . fermionPropS, gTensorF,{{3,6}}]//OrderArray[#,2,1,4,3]&;
	CU=Contract[YTensor2 . fermionPropU, gTensorF2,{{3,6}}]//OrderArray[#,4,1,2,3]&;
	(* TODO include scalar exchange in t-channel. *)

(* === assembling the full result === *)
	
	(*Squared s-channel*)
	resTot=+2*s*u* TotalConj[CS*Conjugate[CS]];
	(*Squared u-channel*)
	resTot+=+2*s*u* TotalConj[CU*Conjugate[CU]];
	(*Mixed s & u channel*)
	resTot+=-2*s*u* TotalConj[CS*Conjugate[CU]+CU*Conjugate[CS]];
	
	If[
		Mod[kinFlip,2]==1,
		resTot=resTot/.{t->t2,u->u2}/.{t2->u,u2->t};
	];
	
	Return[-2*Refine[resTot,Assumptions->VarAsum]]
,
	Return[0]
]	
];


(* ::Subsubsection::Closed:: *)
(*F1F1toS1V1*)


CreateMatrixElement["F1F2toS1V1"][particle1_,particle2_,particle3_,particle4_,particleMass_]:=
Block[{s,t,u,s1,t1,u1,resTot},
(*
	This module returns the squared matrix element of FF->SV summed over all quantum numbers of the incoming particles.
	If the incoming particles are not of the form FF->SV the routine returns 0.
*)
If[ 
	(particle3[[2]]=="F"&&particle4[[2]]=="F")&&((particle1[[2]]=="S"&&particle2[[2]]=="V")||(particle1[[2]]=="V"&&particle2[[2]]=="S"))||
	(particle1[[2]]=="F"&&particle2[[2]]=="F")&&((particle3[[2]]=="S"&&particle4[[2]]=="V")||(particle3[[2]]=="V"&&particle4[[2]]=="S"))
	,
	
	(*The result for FF->SV is the same as  FS->VS if we do some renaming of the mandelstam variables*)
	(*Minus sign as we reversed the momentum of a fermion*)
	resTot=-CreateMatrixElement["F1S1toF1V1"][particle1,particle3,particle2,particle4,particleMass]/.{s->s1,t->t1,u->u1};

	resTot=resTot/.{t1->s,s1->t,u1->u}//Refine;	
	Return[resTot]
,
	Return[0]
]	
];


(* ::Subsubsection:: *)
(*S1S2toV1V2*)


CreateMatrixElement["S1S2toV1V2"][particle1_,particle2_,particle3_,particle4_,particleMass_]:=
Block[{s,t,u,
		partInt,temp,
		particleNull,gTensorS,gTensorV,
		vectorMass,scalarMass,
		vectorPropS,gTensorS2,
		ASS,A,ASU,AST,CT,CU,CTU,CV,totRes,
		ATU,scalarPropT,scalarPropU},
(*
	This module returns the squared matrix element of SS->VV summed over all quantum numbers of the incoming particles.
	If the incoming particles are not of the form SS->VV the routine returns 0.
*)
If[
	(particle1[[2]]!="S"||particle2[[2]]!="S"||particle3[[2]]!="V"||particle4[[2]]!="V")&&
	(particle1[[2]]!="V"||particle2[[2]]!="V"||particle3[[2]]!="S"||particle4[[2]]!="S"),
	Return[0];
,
	(* internal particles *)
	Table[partInt[i]={particle1,particle2,particle3,particle4}[[i,1]], {i,1,4}];
	
	(* reorder particles so that it is always SS -> VV *)
	If[
	  {particle1[[2]], particle2[[2]], particle3[[2]], particle4[[2]]} === {"V", "V", "S", "S"},
	  {partInt[1], partInt[3]} = {partInt[3], partInt[1]};
	  {partInt[2], partInt[4]} = {partInt[4], partInt[2]};
	];
		
	(* Just a trick to make the ordering of particles simpler *)
	particleNull={}; 
	
(*Coupling constants that we will need; and some some that we will not.*)
	gTensorS=Table[gvss[[;;,Particle1,Particle2]],
		{Particle1,{partInt[1],partInt[2],particleNull,particleNull}},
		{Particle2,{partInt[1],partInt[2],particleNull,particleNull}}];
		
	gTensorS2=Table[gvss[[Particle2,Particle1,;;]],
		{Particle1,{partInt[1],partInt[2],particleNull,particleNull}},
		{Particle2,{particleNull,particleNull,partInt[3],partInt[4]}}];	
		
	gTensorV=Table[gvvv[[Particle1,Particle2,;;]],
		{Particle1,{particleNull,particleNull,partInt[3],partInt[4]}},
		{Particle2,{particleNull,particleNull,partInt[3],partInt[4]}}];	
		
(*Propagator masses*)
	vectorMass=particleMass[[1]];
	scalarMass=particleMass[[3]];

(*Vector propagators*)
	vectorPropS=Table[Prop[s,i],{i,vectorMass}]//ListToMat;

(*Scalar propagators*)
	scalarPropT=Table[Prop[t,i],{i,scalarMass}]//ListToMat;
	scalarPropU=Table[Prop[u,i],{i,scalarMass}]//ListToMat;
	

(* === group invariants from vector diagrams === *)
(*
	CT and CU come from s- and u-channel scalar propagator contractions,
	CV comes from s-channel vector propagator contraction. 
	The indices in OrderArray specify the desired ordering for later use.
*)
	CT=-Contract[gTensorS2[[1,3]] . scalarPropT, gTensorS2[[2,4]],{{3,6}}]//OrderArray[#,2,4,1,3]&;
	CU=-Contract[gTensorS2[[1,4]] . scalarPropU, gTensorS2[[2,3]],{{3,6}}]//OrderArray[#,2,4,3,1]&;
	CV=gTensorV[[3,4]] . vectorPropS . gTensorS[[1,2]]//OrderArray[#,3,4,1,2]&;

(* === Lorentz structures that appear === *)
(*
	The coefficients here are the kinematic factors multiplying the 
	group invariants. These are simplified from more general forms 
*)
	ASS=-1/2(32*t*u);
	A=4*s;
	ASU=+4*s(t+2*u);
	AST=-4*s(u+2*t);
	ATU=+8*s^2;
	
(* === assembling the full result === *)
	totRes=ASS*Total[CV^2,-1];
	totRes+=ASU*Total[CV*CU,-1];
	totRes+=AST*Total[CV*CT,-1];
	totRes+=ATU*Total[CT*CU,-1];
	totRes+=A*Total[(t*CT+u*CU)*CU,-1];
	totRes+=A*Total[(t*CT+u*CU)*CT,-1];
	totRes+=4*Total[(t*CT+u*CU)^2,-1];

	Return[Refine[totRes/2,Assumptions->VarAsum]] 
]	
];


(* ::Subsubsection::Closed:: *)
(*S1V1toS1V1*)


CreateMatrixElement["S1V1toS1V1"][particle1_,particle2_,particle3_,particle4_,particleMass_]:=
Block[{
	s,t,u,partInt,temp,resTot,
	u1,s1,t1,kinFlip},
(*
	This module returns the squared matrix element of SV->SV summed over all quantum numbers of the incoming particles.
	If the incoming particles are not of the form SV->SV the routine returns 0.
*)
	kinFlip=0; (*records if the external or incoming legs were swapped, and does a t<->u transformation for each swap*)
If[ (
	(particle1[[2]]=="V"&&particle2[[2]]=="S")||
	(particle1[[2]]=="S"&&particle2[[2]]=="V")
	)&&(
	(particle3[[2]]=="V"&&particle4[[2]]=="S")||
	(particle3[[2]]=="S"&&particle4[[2]]=="V")),

	(* internal particles *)
	Table[partInt[i]={particle1,particle2,particle3,particle4}[[i]], {i,1,4}];

(*Changing the order of the particles so that it is always SV->SV*)	
	If[particle1[[2]]=="V"&&particle2[[2]]=="S",
		{partInt[1],partInt[2]} = {partInt[2],partInt[1]};
		kinFlip+=1;
	];
	If[particle3[[2]]=="V"&&particle4[[2]]=="S",
		{partInt[3],partInt[4]} = {partInt[4],partInt[3]};
		kinFlip+=1;
	];
	
	(*The result for SV->SV is identical to SS->VV (some would say twice as identical) if we do some renaming of the mandelstam variables*)
	resTot=CreateMatrixElement["S1S2toV1V2"][partInt[1],partInt[3],partInt[2],partInt[4],particleMass]/.{s->s1,t->t1,u->u1};
	resTot=resTot/.{s1->t,t1->s,u1->u};
	
	If[Mod[kinFlip,2]==1,resTot=resTot/.{t->t1,u->u1}/.{t1->u,u1->t};];
	Return[resTot]
,
	Return[0]
]	
];


(* ::Subsubsection::Closed:: *)
(*S1S2toS3V1*)


CreateMatrixElement["S1S2toS3V1"][particle1_,particle2_,particle3_,particle4_,particleMass_]:=
Block[{
	s,t,u,partInt,resTot,u1,s1,t1,temp1,temp2,
	particleNull,gTensor,\[Lambda]3Tensor,\[Lambda]4Tensor,
	scalarPropT,scalarPropU,A,CT,CU,CS,kinFlip,
	scalarMass,
	scalarPropS
},
(*
	This module returns the squared matrix element of SS->SV summed over all quantum numbers of the incoming particles.
	If the incoming particles are not of the form SS->SV the routine returns 0.
	
	I maintain that only a lunatic will actually have this matrix element appearing in their model, and if you are actually 
	reading this you have but proved my point.
*)
	kinFlip=0; (*records if the external or incoming legs were swapped, and does a t<->u transformation for each swap*)
If[ (
	(particle1[[2]]=="V"&&particle2[[2]]=="S")||
	(particle1[[2]]=="S"&&particle2[[2]]=="V")
	)&&(
	(particle3[[2]]=="S"&&particle4[[2]]=="S"))||
	(
	(particle3[[2]]=="V"&&particle4[[2]]=="S")||
	(particle3[[2]]=="S"&&particle4[[2]]=="V")
	)&&(
	(particle1[[2]]=="S"&&particle2[[2]]=="S")),

	Table[partInt[i]={particle1,particle2,particle3,particle4}[[i]], {i,1,4}];

(*Changing the order of the particles so that it is always SS->SV*)	
	If[particle1[[2]]=="V"||particle2[[2]]=="V",
		{partInt[1],partInt[2],partInt[3],partInt[4]} = {partInt[3],partInt[4],partInt[1],partInt[2]};
	];
	If[partInt[3][[2]]=="V",
		{partInt[3],partInt[4]} = {partInt[4],partInt[3]};
		kinFlip+=1;
	];
	
	particleNull={{}};
	(*Just a trick to not confuse the ordering of particles*)
	(*Coupling constants that we will need*)
	
	gTensor=Table[gvss[[partInt[4][[1]],Particle1[[1]],;;]],
		{Particle1,{partInt[1],partInt[2],partInt[3],particleNull}}];
	
	\[Lambda]3Tensor=Table[\[Lambda]3[[;;,Particle1[[1]],Particle2[[1]]]],
		{Particle1,{partInt[1],partInt[2],partInt[3]}},
		{Particle2,{partInt[1],partInt[2],partInt[3]}}];

(*Propagator masses*)
	scalarMass=particleMass[[3]];

(*Scalar propagators*)
	scalarPropT=Table[Prop[t,i],{i,scalarMass}]//ListToMat;
	scalarPropU=Table[Prop[u,i],{i,scalarMass}]//ListToMat;
	scalarPropS=Table[Prop[s,i],{i,scalarMass}]//ListToMat;

(*Lorentz structures that appear.*)
	A=4*t*u/s;
	
(*Group structures that are used*)
	CS=Contract[scalarPropS . \[Lambda]3Tensor[[1,2]], gTensor[[3]],{{1,6}}]//OrderArray[#,1,2,4,3]&;
	CT=Contract[scalarPropT . \[Lambda]3Tensor[[1,3]], gTensor[[2]],{{1,6}}]//OrderArray[#,1,4,2,3]&;
	CU=Contract[scalarPropU . \[Lambda]3Tensor[[2,3]], gTensor[[1]],{{1,6}}]//OrderArray[#,4,1,2,3]&;

(* === assembling the full result === *)
	resTot=-4*Total[CS*(u*CT+t*CU),-1];
	resTot+=-4 s*Total[CT CU,-1];
	
	If[Mod[kinFlip,2]==1,resTot=resTot/.{t->t1,u->u1}/.{t1->u,u1->t};];
	Return[resTot]
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
	
	If[bNormalizeWithDOF==False,dof=1];
	
	Return[dof];
]


ExtractOutOfEqElement[Channel_?StringQ][particleListAll_, particleListOutOfEq_, particleMasses_]:=
Block[{
	OutOfEqParticles,MatrixElements,Elements,CollElements,
	VectorMass,FermionMass,ScalarMass,
	symmetries,deltaF
},
(*incomingParticle is the particle associated with the momentum p_1*)
(*deltaFparticle is the particles whos deltaF contributions we want*)
(*
	Essentially this is generates the elements going into
	Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]
*)

(*First we generate all matrix elements*)
(*	OutOfEqParticles=Complement[Table[i,{i,1,Length[particleList]}],LightParticles];*)
(*	OutOfEqParticles=Table[i,{i,1,Length[particleListOutOfEq]}];*)

(*Divide the incoming particle by its degree's of freedom*)
	MatrixElements=Table[
		1/degreeOfFreedom[a]
		CreateMatrixElement[Channel][a,b,c,d,particleMasses],
		{a,particleListOutOfEq},
		{b,particleListAll},{c,particleListAll},{d,particleListAll}]//SparseArray;
	(*This is a list of all non-zero matrix elements*)
	Elements=MatrixElements["NonzeroPositions"];

(*If there are no matching matrix elements we return 0*)
If[Length[Elements]==0,
	Return[{}]; 
,
	MatrixElements=Map[Extract[MatrixElements, #] &, Elements];

(*We now select all elements where deltaFparticle is amongst the scattered particles*)
(*deltaF is here a list of 1 and 0s*)
(*Gives all light-particle the label specified by the first lightparticle*)
(*	If[LightParticles!={},
		deltaF=Elements/. x_?NumericQ /;MemberQ[ LightParticles,x ] -> LightParticles[[1]]; 
	,
		deltaF=Elements;
	];*)
	
(*Now we create the full list of distinct collision elements*)
	CollElements=Table[{MatrixElements[[i]],Elements[[i]]},{i,1,Length[Elements]}];

(*We now add all elements with the same deltaF list*)
	symmetries=Gather[CollElements,#1[[2]]==#2[[2]]&];
	CollElements=Map[{Total[#[[;;,1]], -1], #[[1, 2]]} &,symmetries];
	
	Return[CollElements]
]	
];


SimplifyConjugates[expr_, assumptions_List] :=
  expr //. Conjugate[x_] :> FullSimplify[
     Conjugate[x],
     Assumptions -> assumptions
     ];


ExtractOutOfEqElement[particleListAll_,particleListOutOfEq_,ParticleMasses_]:=
	Block[{
		collisionElements, collisionElementsTotal
	},
	(*incomingParticle is the particle associated with the momentum p_1*)
	(*deltaFparticle is the particles whos deltaF contributions we want*)
	(*Essentially this is generates the elements going into Sum_deltaFparticle \[Delta]C[incomingParticle,deltaFparticle]*deltaF[deltaFparticle]*)
	
	(*particleList is the complete list of particles*)
	(*First we extract the result for all subprocesses*)
	collisionElements = {
		"F1F2toF3F4",
		"F1F2toV1V2", "F1V1toF1V1", "V1V2toV3V4", "S1S2toS3S4", 
		"S1S2toF1F2", "F1S1toF1S1", "F1S1toF1V1",
		"F1F2toS1V1", "S1S2toV1V2", "S1V1toS1V1", "S1S2toS3V1"
	};
	
	collisionElementsTotal = 
	  Map[(If[bVerbose,Print[#]];
	   ExtractOutOfEqElement[#][particleListAll, particleListOutOfEq, ParticleMasses])&, 
	   collisionElements
	  ]//Flatten[#,1]&;
	
	Prop /: Conjugate[Prop[x__]] := Prop[x];
	collisionElementsTotal=SimplifyConjugates[
			collisionElementsTotal,
			VarAsum
		]//Expand;	
	UpValues[Prop] = DeleteCases[UpValues[Prop], HoldPattern[(Conjugate[Prop[_]]) :> _]];
	
	collisionElementsTotal =
	  Which[
	    bTruncateAtLeadingLog && !bTagLeadingLog,
	      TruncateAtLeadingLogarithm[collisionElementsTotal] /. {logEnhanced -> 1, powerEnhanced -> 1},
	    bTruncateAtLeadingLog && bTagLeadingLog,
	      TruncateAtLeadingLogarithm[collisionElementsTotal],
	    bTagLeadingLog,
	      TagLeadingLogarithm[collisionElementsTotal],
	    True,
	      collisionElementsTotal
	  ];	
	
	collisionElementsTotal=collisionElementsTotal/.Prop[x_,y_]->1/(x-y);
	
	Return[collisionElementsTotal//Expand];

];


(* ::Section:: *)
(*Generating matrix elements*)


(* function not in use any more *)
ExtractLightParticles[outOfEqParticleList_,lightParticleList_,particleListAll_]:=
Block[{
	position,nonEqParticles,lightParticles,
	missingParticles
},
(*
Extracting the out-of-eq particles
*)
	missingParticles={};
	(*This list includes the light particles which are added below*)
	particleListAll=Join[outOfEqParticleList,lightParticleList];
		
(*Adding all remaining particles as light*)
	Table[
		position=Position[particleListAll, _?(# == StringTake[particle, 1] &)][[;;,1]];
		nonEqParticles=Table[particleListAll[[i]][[1]],{i,position}]//Flatten[#]&;
		lightParticles={Complement[RepToIndices[PrintFieldRepPositions[particle]],nonEqParticles],StringTake[particle, 1]};
		AppendTo[missingParticles,lightParticles],
		{particle,{"Vector","Fermion","Scalar"}}
	];
	
	If[
		Length[DeleteCases[missingParticles, {{},a__}]]>0,
		Message[WallGoMatrix::missingParticles,
			If[missingParticles[[1]][[1]]=={},"None", missingParticles[[1]]],
			If[missingParticles[[2]][[1]]=={},"None", missingParticles[[2]]],
			If[missingParticles[[3]][[1]]=={},"None", missingParticles[[3]]]
		];
		Abort[];
	];
	If[
		Length[DeleteCases[missingParticles, {{},a__}]]>1,
		Message[WallGoMatrix::failmsg,
			"Multiple species declared as light in "<>ToString[particleListAll[[missingParticles]]]<>". "<>
			"Only one species allowed. Solution: make additional species explicit in particleList."];
		Abort[];
	];
]


GenerateMatrixElements[MatrixElements_,Cij_,particleListAll_,particleListOutOfEq_,particleMasses_]:=
Block[{Elements},

	Do[
		If[item[[1]] === {},
			Message[WallGoMatrix::failmsg,
				"Found empty particle species "<>ToString[item]<>" in "<>ToString[particleListAll]];
				Abort[];
		],
	    {item, particleListAll}
	];

	MatrixElements=ExtractOutOfEqElement[particleListAll,particleListOutOfEq,particleMasses];
		
(*Extract various C^{ij} components*)	
	Cij=ConstantArray[0,{Length[particleListOutOfEq],Length[particleListOutOfEq]}];
	Do[
		Elements=Extract[MatrixElements,Position[MatrixElements[[;;,2]],{i,___}]];
		Do[
			Cij[[i,j]]=Extract[Elements,#]&/@Union[
						Position[Elements[[;;,2]],{_,j,__}],
						Position[Elements[[;;,2]],{_,_,j,_}],
						Position[Elements[[;;,2]],{_,__,j}]];,
			{j,Length[particleListOutOfEq]}
		],
		{i,Length[particleListOutOfEq]}
	];

];


TagLeadingLogarithm::usage =
"TagLeadingLogarithm[MatrixElements] classifies a list of \
2-to-2 scattering matrix elements according to their enhanced structures in the \
high-energy limit, retaining only the leading logarithmic contributions.

 \[Bullet] Input: a list of elements of the form {amplitudeExpression, {i,j,k,l}},
   where the second entry encodes external state indices.

 \[Bullet] The classification follows Algorithm~1 of the accompanying publication [25xx:xxxx], applied when the option TruncateAtLeadingLog -> True is set in ExportMatrixElements.

 \[Bullet] The algorithm distinguishes three cases:
   - \!\(\*StyleBox[\"Forward channel\", FontSlant -> \"Italic\"]\): when particles a \[Congruent] c and b \[Congruent] d.
       \[Bullet] If also a \[Congruent] d (fully identical particles):
         no power enhancement;
         logarithmic enhancement from 1/t^2 and 1/u^2;
         terms \[Proportional] 1/t or 1/u are finite.
       \[Bullet] Otherwise:
         power enhancement from 1/u^2;
         logarithmic enhancement from 1/t^2 or 1/u.
   - \!\(\*StyleBox[\"Crossed channel\", FontSlant -> \"Italic\"]\): when particles a \[Congruent] d and b \[Congruent] c.
       power enhancement from 1/t^2;
       logarithmic enhancement from 1/u^2 or 1/t.
   - \!\(\*StyleBox[\"Generic case\", FontSlant -> \"Italic\"]\):
       power enhancement from 1/t^2 and 1/u^2;
       logarithmic enhancement from 1/t or 1/u.

 \[Bullet] Output: a tagged list of the form \
   {{leadingLogContribution1, {i,j,k,l}}, ...}, \
   with \!\(\*
StyleBox[\"powerEnhanced\",\nFontSlant->\"Italic\"]\) for power enhanced structures and \!\(\*
StyleBox[\"logEnhanced\",\nFontSlant->\"Italic\"]\) for logarithmic enhanced terms.\
";


TagLeadingLogarithm[MatrixElements_]:=Module[
	{
		MatrixElementsF,S,T,U
	},

	MatrixElementsF=MatrixElements/.Flatten[Map[{
			Prop[#[[1]],msq_]->#[[2]]^(-1)*Prop[#[[1]],msq],
			#[[1]]->#[[1]]*#[[2]]
		}&,{{s,S},{t,T},{u,U}}]];
	
	MatrixElementsF=Map[{
		If[#[[2]][[1]]==#[[2]][[3]]&&#[[2]][[2]]==#[[2]][[4]],
			If[#[[2]][[1]]==#[[2]][[4]],
				#[[1]]/.{
					1/U^2->logEnhanced/U^2,
					1/T^2->logEnhanced/T^2
					},
				#[[1]]/.{
					1/U^2->powerEnhanced/U^2,
					1/T^2->logEnhanced/T^2,
					1/U->logEnhanced/U
					}
			],
			If[#[[2]][[1]]==#[[2]][[4]]&&#[[2]][[2]]==#[[2]][[3]],
				#[[1]]/.{
					1/T^2->powerEnhanced/T^2,
					1/U^2->logEnhanced/U^2,
					1/T->logEnhanced/T
					},
				#[[1]]/.{
					1/U^2->powerEnhanced/U^2,
					1/T^2->powerEnhanced/T^2,
					1/T/U->logEnhanced/T/U,
					1/U->logEnhanced/U,
					1/T->logEnhanced/T
					}
			]
		],
		#[[2]]}&,
			MatrixElementsF
		];
		
	If[Max[Exponent[#, {logEnhanced,powerEnhanced}, Max]&/@MatrixElementsF[[All, 1]]] > 1,
		Message[WallGoMatrix::failmsg,
			"Found mixed log divergences some channels were not correctly leading log-filtered."];
		Abort[];
	];
		
	MatrixElementsF=Collect[MatrixElementsF,{S,T,U},Simplify]/.Thread[{T,U,S}->1];
	Return[MatrixElementsF];
]


TruncateAtLeadingLogarithm[MatrixElements_]:=Module[
	{
		MatrixElementsF
	},

	MatrixElementsF=TagLeadingLogarithm[MatrixElements];
		
	MatrixElementsF=Map[{
			+ powerEnhanced*SeriesCoefficient[#[[1]],{powerEnhanced,0,1}]
			+ logEnhanced*SeriesCoefficient[#[[1]],{logEnhanced,0,1}]
			,
			#[[2]]}&,
			MatrixElementsF
		];
		
	MatrixElementsF=DeleteCases[MatrixElementsF, {0,{a__}}];
	Return[MatrixElementsF];
]


TruncateAtLeadingLogarithmOld[MatrixElements_]:=Module[{MatrixElementsF,U,S,T},

	MatrixElementsF=MatrixElements/.Flatten[Map[{
			Power[+#[[1]] + msq_, n_?Negative]->#[[2]]^(-n)*Power[+#[[1]] + msq, n],
			Power[-#[[1]] + msq_, n_?Negative]->#[[2]]^(-n)*Power[-#[[1]] + msq, n]
		}&,{{s,S},{t,T},{u,U}}]];
	
	MatrixElementsF=Map[{
		Plus@@Table[
		+SeriesCoefficient[#[[1]]/.{T->xLarge*T,U->xLarge*U},{xLarge,Infinity,-i}]
		,
		{i,If[#[[2]][[3]]==#[[2]][[4]],2,1],2}],
		#[[2]]}&,
		MatrixElementsF];

	MatrixElementsF=Expand[MatrixElementsF]/.{S*T->0,S*U->0,T*U->0};
	MatrixElementsF=Collect[MatrixElementsF,{S,T,U},Simplify]/.Thread[{T,U,S}->1];
	MatrixElementsF=DeleteCases[MatrixElementsF, {0,{a__}}];
	
	Return[MatrixElementsF];
]
