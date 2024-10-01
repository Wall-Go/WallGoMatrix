(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: WallGoMatrix                                                          	*)

(*
      	This software is covered by the GNU General Public License 3.
       	Copyright (C) 2024-2024 Andreas Ekstedt
       	Copyright (C) 2024-2024 Oliver Gould
       	Copyright (C) 2024-2024 Joonas Hirvonen
       	Copyright (C) 2024-2024 Benoit Laurent
       	Copyright (C) 2024-2024 Lauri Niemi
       	Copyright (C) 2024-2024 Philipp Schicho
       	Copyright (C) 2024-2024 Jorinde van de Vis
*)

(* :Summary:	WallgoMatirx is an algorithm that constructs
				the Matrix Elements at Leading Log and partially NLL for generic models.	*)	

(* ------------------------------------------------------------------------ *)


BeginPackage["WallGoMatrix`"]


Unprotect@Definition;
Definition[x_Symbol] /; StringMatchQ[Context[x], "Package`" ~~ ___] :=
    StringReplace[ToString@FullDefinition[x],
        (WordCharacter .. ~~ DigitCharacter ... ~~ "`") .. ~~ s_ :> s
    ];
Protect@Definition;


(*
	Welcome banner: All credit for this part to GroupMath
*)
TexFor[text_]:=Style[text,{GrayLevel[0.3]}]
result={};
AppendTo[result,Row[{
	TexFor["GOGOGOGOGOGOGOGOGOGOGOGO "],
	TexFor["WallGoMatrix"],
	TexFor[" GOGOGOGOGOGOGOGOGOGOGGOGOGO"]}]];
AppendTo[result,Row[{"Version: "//TexFor,"0.02 (30-09-2024)"//TexFor}]];
AppendTo[result,Row[{"Authors: "//TexFor,
    Hyperlink["Andreas Ekstedt","https://inspirehep.net/authors/1799400"],", "//TexFor,
    Hyperlink["Oliver Gould","https://inspirehep.net/authors/1606373"],", "//TexFor,
    Hyperlink["Joonas Hirvonen","https://inspirehep.net/authors/1866901"],", "//TexFor,
    Hyperlink["Benoit Laurent","https://inspirehep.net/authors/1808372"],", "//TexFor,
    Hyperlink["Lauri Niemi","https://inspirehep.net/authors/1764675"],", "//TexFor,
    Hyperlink["Philipp Schicho","https://inspirehep.net/authors/1639147"],", "//TexFor,
    Hyperlink["Jorinde van de Vis","https://inspirehep.net/authors/1589979"]//TexFor
    }
  ]];
AppendTo[result,Row[{"Reference: "//TexFor,
        Hyperlink["2410.xxxxx [hep-ph]","https://arxiv.org/abs/2410.xxxxx"]//TexFor}]];
AppendTo[result,Row[{"Repository link: "//TexFor,
	Hyperlink[Mouseover[TexFor["github.com/Wall-Go/WallGoMatrix"],Style["github.com/Wall-Go/WallGoMatrix",Bold]],
	"https://github.com/DR-algo/DRalgo"]}]];
AppendTo[result,Row[{"Existing model files: "//TexFor,
	Hyperlink[Mouseover[TexFor["WallGoMatrix/examples"],Style["WallGoMatrix/examples",Bold]],
	"https://github.com/Wall-Go/WallGoMatrix/tree/main/examples"]}]];
WallGoMatrixLoad=Image[CompressedData["1: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"],"Byte",ColorSpace->"RGB",Interleaving->True];

AppendTo[result,Style["GOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGGOGOGOGOGOGOGOGO",{GrayLevel[0.3]}]];
Print[Grid[{{Image[WallGoMatrixLoad,ImageSize->140],Row[result,"\n",BaseStyle->(FontFamily->"Consolas")]}},Alignment->{Left,Center}]]
(**)



(*List of public functions*)

CreateParticle::usage=
"CreateOutOfEq[{{\!\(\*SubscriptBox[\(r\), \(1\)]\),\!\(\*SubscriptBox[\(m\), \(1\)]\)},...,{\!\(\*SubscriptBox[\(r\), \(n\)]\),\!\(\*SubscriptBox[\(m\), \(n\)]\)}},\"R\"] is a function that groups n particles into one representation for
Fermions (R=F),
Vector Bosons (R=V), and
Scalars (R=S)."
SymmetryBreaking::usage=
"Classify different scalar, fermion, and vector representations into respective particles using VEV-induced masses.
As an output particle i are given as {\!\(\*SubscriptBox[\(r\), \(i\)]\),\!\(\*SubscriptBox[\(m\), \(i\)]\)} where \!\(\*SubscriptBox[\(r\), \(i\)]\) is the label of the representation and \!\(\*SubscriptBox[\(m\), \(i\)]\) is the label of the particle mass in that representation"
ExportMatrixElements::usage=
"ExportMatrixElements[file,particleList,UserMasses,UserCouplings,ParticleName,ParticleMasses,OptionsPattern[]]

generates all possible matrix elements with the external particles specified in particleList.
The format can be specified by the option Format, with currently supported options: X=txt,json,hdf5,all,none ,
the last two options exports the result in all possible formats, and in none, respectively. A list of matrix elements is returned by the function regardless of the choosen format.";


AllocateTensors::usage="\
Creates gauge generators";
GradQuartic::usage="Creates Quartic tensors";
GradCubic::usage="Creates Cubic tensors";
GradTadpole::usage="Creates Tadpole tensors";
GradSextic::usage="Creates dim 6 tensors";
GradMass::usage="Creates Mass tensors";
CreateInvariant::usage="Creates an invariant";
CreateInvariantYukawa::usage="Creates Yukawa Tensor";
GradYukawa::usage="Creates Yukawa tensor";
GradMassFermion::usage="Creates Fermion Invariants";
CreateInvariantFermion::usage="Creates Fermion Invariants";


PerformDRhard::usage="\
Performs the dimensional reduction from hard to soft";


PrintGaugeRepPositions::usage="Prints the indices of Gauge reps";
PrintScalarRepPositions::usage="Prints the indices of Scalar reps";
PrintFermionRepPositions::usage="Prints the indices of Fermion reps";


ImportMatrixElements::usage="";


SymmetryBreakingGauge::usage="";
SymmetryBreakingFermion::usage="";
SymmetryBreakingScalar::usage="";


DefineDim6::usage="Defines a dimension 6 operator";
CompareInvariants::usage="\
Finds relations between couplings by calculating basis-invariant tensors";


$DRalgoDirectory=DirectoryName[$InputFileName];


(*
	Functions from groupmath are used to create the model.
*)
If[Global`$LoadGroupMath,
	Get["GroupMath`"];
	Print["GroupMath is an independent package, and is not part of DRalgo"];
	Print["Please Cite GroupMath: Comput.Phys.Commun. 267 (2021) 108085 \[Bullet] e-Print: 2011.01764 [hep-th]"];
];
Print["WallGoMatrix is powered by the DRalgo ModelCreation."]
Print["Please Cite DRalgo: Comput.Phys.Commun. 288 (2023) 108725 \[Bullet] e-Print: 2205.08815 [hep-ph]"];


(*
	Verbose=True removes progress messages.
	Mode=2 calculates everything, Mode=1 only calculates LO masses and couplings
	Mode=0 only calculates LO masses
*)
Options[ImportModelDRalgo]={
	Verbose -> False,Mode->2,
	Dim6->False,
	Normalization4D->False,
	AutoRG->True}


Begin["`Private`"]


(*
	Loads all functions.
*)
Get[FileNameJoin[{$DRalgoDirectory,"modelCreation.m"}]];(*Loads DRalgo model creation*)
(*Get[FileNameJoin[{$DRalgoDirectory,"HardToSoft.m"}]];(*Loads Hard -> Soft functions*)*)
Get[FileNameJoin[{$DRalgoDirectory,"matrixelements.m"}]];(*Loads matrix element creation*)


(* ::Section:: *)
(*Initialization*)


(*
	Defines internal tensors from the loaded model. Also creates help-tensors used for
	intermediate calculations.
*)
ImportModelDRalgo[GroupI_,gvvvI_,gvffI_,gvssI_,\[Lambda]1I_,\[Lambda]3I_,\[Lambda]4I_,\[Mu]ijI_,\[Mu]IJFI_,\[Mu]IJFCI_,YsffI_,YsffCI_, OptionsPattern[]]:=Module[
{
	GroupP=GroupI,gvvvP=gvvvI,gvffP=gvffI,gvssP=gvssI,\[Lambda]1IP=\[Lambda]1I,\[Lambda]3P=\[Lambda]3I,\[Lambda]4P=\[Lambda]4I,
	\[Mu]ijP=\[Mu]ijI,\[Mu]IJFP=\[Mu]IJFI,\[Mu]IJFCP=\[Mu]IJFCI,YsffP=YsffI,YsffCP=YsffCI},


If[ Global`$LoadGroupMath,
	If[!GroupMathCleared && !ValueQ[Global`$GroupMathMultipleModels],
		Remove["GroupMath`*"];
		GroupMathCleared=True;
	];
];
	\[Mu]ij=\[Mu]ijP//SparseArray//SimplifySparse;
	gvvv=gvvvP//SparseArray//SimplifySparse;
	gvss=gvssP//SparseArray//SimplifySparse;
	\[Lambda]4=\[Lambda]4P//SparseArray//SimplifySparse;
	\[Lambda]3=\[Lambda]3P//SparseArray//SimplifySparse;
	Ysff=YsffP//SparseArray//SimplifySparse;
	YsffC=YsffCP//SparseArray//SimplifySparse;
	gvff=gvffP//SparseArray//SimplifySparse;
	\[Mu]IJF=\[Mu]IJFP//SparseArray//SimplifySparse;
	\[Mu]IJFC=\[Mu]IJFCP//SparseArray//SimplifySparse;
	\[Lambda]1=\[Lambda]1IP//SparseArray//SimplifySparse;
	ns=Length[gvss[[1]]];
	nv=Length[gvvv];
	nf=Length[gvff[[1]]];
	\[Lambda]6=EmptyArray[{ns,ns,ns,ns,ns,ns}];

(*Options*)
	verbose=OptionValue[Verbose];
	mode=OptionValue[Mode]; (*If 2 everthing is calculated. And if 2 only 1-loop contributions are calculated*)
	NFMat=IdentityMatrix[nf]//SparseArray; (*This matrix is only relevant if the user wants an arbitrary number of fermion families*)
	normalization4D=OptionValue[Normalization4D];
	If[OptionValue[AutoRG],rgFac=1,rgFac=0]; (*Determines whether the RG-running is automatically incorporated in the soft masses and tadpoles*)
(*End of Options*)

	CT=False; (*Checks if counter-terms have already been calculated*)
	DefineGroup[GroupP]; (*Names Debye masses*)
	GroupDR=GroupP; (*For saving purposes*)
	CreateHelpTensors[] (*Creates recurring tensors*)
];


(*
	Defines a \[Phi]^6 operator
*)
DefineDim6[\[Lambda]6I_]:=Module[{\[Lambda]6P=\[Lambda]6I},
	If[mode>=3,
		\[Lambda]6=\[Lambda]6P//SparseArray;
	,
		Message[DRalgo::ImplementationFail, "Please set mode=3 to use this feature"];
		Abort[]
	];
];


(*
	Takes the user-defined group and names debye masses.
*)
DefineGroup[GroupI_]:=Module[{GroupP=GroupI},
(*The Debye mass matrix is \[Mu]abDef*)
(*At tree-level this matrix is 0. But after DR, thermal masses are named accoording to \[Mu]abDef.*)
	\[Mu]abDef=CreateDebyeMasses[GroupP];
	VecMassDefined=True;
];


(*
	Enables the user to add an arbitrary number of fermions.
	This works by creating a diagonal matrix. Each element NFMatP[[i,i]] corresponds to how many times the
	fermion labled by ii appears.
*)
DefineNF[NFMatP_]:=Module[{NFMatI=NFMatP},
	Do[NFMat[[i[[2]],i[[2]]]]*=i[[1]],{i,NFMatI}];(*Each element is multiplied with nF_i*)
];


(* ::Section:: *)
(*Help functions*)


OutputFormatDR[x_]:=ToExpression[StringReplace[ToString[StandardForm[x]],"WallGoMatrix`Private`"->""]];


(*
	Prints constants that appear in the effective theory.
*)
PrintConstants[]:=Module[{},
ToExpression[StringReplace[ToString[StandardForm[{Lb->(Log[\[Mu]^2/T^2]-2 Log[4 \[Pi]]+2EulerGamma),Lf->(Log[\[Mu]^2/T^2]-2 Log[4 \[Pi]]+2EulerGamma+4 Log[2])}]],"WallGoMatrix`Private`"->""]]
];


(*
	Replaces the Glaisher constant by c, which is sometimes used in the litterature. See for example hep-ph: 9508379
*)
PrintGenericBasis[]:=Module[{},
ToExpression[StringReplace[ToString[StandardForm[{Log[Glaisher]->-1/12 (Lb+2cplus-EulerGamma)}]],"WallGoMatrix`Private`"->""]]
];


(*
	Checks if two variables are identical up to a numerical factor
	For example if a=x^2 y and b=y, then there's no relation, and the function returns nothing.
	While if a=x^2 and b=5 x^2 the function returns 5
*)
CompExp2[a0_,b0_]:=Module[{a=a0,b=b0},
Comps=Solve[a v[1]- b ==0,v[1]]//Simplify//DeleteDuplicates//Select[#,UnsameQ[#,{}]&]&;
Comps/. {v[x_]->a_}:>a
]


(*
	Old function: Should maybe be removed.
	This function finds linear dependencies between the variables in list.
	Only numerical factors are included.
	The result is given in Mat. Where Mat[[i,j]]=0 if list[[i]]!=N list[[j]] where N is a number.
	Otherwise Mat[[i,j]]=N if list[[i]]=Nlist[[j]]
*)
OverallFac[list_,varMat_]:=Module[{L=list,Mat=varMat},

	For[i=1,i<Length[list],i++,
		For[j=i+1,j<Length[list]+1,j++,
			help=CompExp2[list[[i]],list[[j]]];
			TempIf=If[Length[help]>0,{NumericQ/@help}[[1,1]]&&Equal@@NumericQ/@help,False];
			If[TempIf,a=help[[1]],a=0];
				Mat[[i,j]]=a;
			];
		];
	Mat
]


(*
	Old function:Should maybe be removed.
	This function finds linear dependences between the variables in listvar
	Only relations between two variables are considered.
	This function is only used for the soft->supersoft step where the faster, and more general,
	RelationsBVariables3 has problems due to inverse powers of masses.
*)
RelationsBVariables[list_,listVar_]:=Module[{L=list,LV=listVar},
	LVTemp=LV;
	If[L[[1]]==0,
		L=Delete[L,1];
	];

	Mat=ConstantArray[0,{Length[L],Length[L]}];(*Linear-dependency matrix*)
	Mat1=OverallFac[L,Mat];
	TempVar=ConstantArray[0,{Length[LV]}];
	For[i=1,i<Length[LV],i++,
		For[j=i+1,j<Length[LV]+1,j++,
			IfTemp=Mat1[[1;;i,j]]//DeleteDuplicates//DeleteCases[#,0,Infinity]&; (*Removes cases when the numerical factor is 0 or infinity*)
			If[Length[IfTemp]==1&&Mat1[[i,j]]!=0,
				LVTemp[[j]]=Mat1[[i,j]]LV[[i]];(*Creates a list with all linear relations*)
			];
		];
	];
	LVTemp
]


(*
	Creates a list of all possible couplings that can appear in
	1-loop matching relations
	Note that this does not include couplings in debye/scalar masses
	at 1-loop level
*)
CreateBasisVanDeVis[]:=Module[{},
(* This module received founding from the fish *)

	FermionFamVar=Normal[NFMat]//Variables;
	ScalVar=\[Lambda]4//Normal//Variables;

	GaugeVarPre=Join[Normal[gvvv]//Variables,Normal[gvff]//Variables,Normal[gvss]//Variables]//DeleteDuplicates; (*Includes possible non-numeric charges*)
	GaugeCharge=Complement[GaugeVarPre,GaugeCouplingNames];(*Possible non-numeric gauge charges*)

	GaugeVarHelp=Table[i*j,{i,GaugeCouplingNames},{j,GaugeCharge}];
	GaugeVar=Join[GaugeCouplingNames,GaugeVarHelp];

	YukVar=Normal[Ysff]//Variables;
	AuxVar={1,Lb,Lf};
	AuxVar=Join[AuxVar,FermionFamVar]//DeleteDuplicates;
	varHelp=Join[ScalVar,GaugeVar,YukVar];
	t2=Table[i*j,{i,varHelp},{j,varHelp}]//Flatten[#]&//DeleteDuplicates;
	t3G=Table[i*j*k,{i,ScalVar},{j,GaugeVar},{k,GaugeVar}]//Flatten[#]&//DeleteDuplicates;
	t3F=Table[i*j*k,{i,ScalVar},{j,YukVar},{k,YukVar}]//Flatten[#]&//DeleteDuplicates;
	t4=Table[i*j*k*l,{i,GaugeVar},{j,GaugeVar},{k,GaugeVar},{l,GaugeVar}]//Flatten[#]&//DeleteDuplicates;
	t4F=Table[i*j*k*l,{i,YukVar},{j,YukVar},{k,YukVar},{l,YukVar}]//Flatten[#]&//DeleteDuplicates;
	t4FG=Table[i*j*k*l,{i,YukVar},{j,YukVar},{k,GaugeVar},{l,GaugeVar}]//Flatten[#]&//DeleteDuplicates;
	basPre=Join[varHelp,t2,t3G,t3F,t4,t4F,t4FG];

	basDR=Table[i*j*k,{i,basPre},{j,AuxVar},{k,{1,T,T^2}}]//Flatten[#]&//DeleteDuplicates;
]


{basDR};


(*
	This function finds linear dependences between the variables in list.
	This works by treating each element in list as a vector in the space spanned by basDR.
	All vectors are then rowreduced, and a minimal set of basis vectors (in list) are found.
*)
RelationsBVariables3[list_]:=Module[{L=list},
(*Creates a vector-basis*)

(*One could say that v3 and v2 are identical. With v3 being almost twice as identical as v2*)
	If[L[[1]]==0&&Length[L]>1,
		Lp=Delete[L,1];
	,
		Lp=L;
	];

	varHelp=Lp//Variables;
	varFix=#->0&/@varHelp; (*Trick to ensure that vectors are expanded properly*)

(*
Expands all elements in Lp in terms of basDR.
*)
	setVecs=Table[Coefficient[i,basDR],{i,Lp}]/.varFix;  
(*Delete columns with only 0s*)
	setVecs=Transpose[DeleteCases[Transpose[setVecs], {0 ..}, Infinity]];

(*Finds independent basis*)
	rr = setVecs // Transpose // RowReduce;
	rr=DeleteCases[rr, {0 ..}, Infinity];
(*Basis elements*)
	basisElements = Flatten[FirstPosition[#, 1, Nothing] & /@ rr];
	varBasis=Table[ \[Lambda]VL[a],{a,1,Length[basisElements]}];

(*Puts everything together*)
	LVTemp=LVTemp=ConstantArray[0,Length[Lp]];
	Do[LVTemp[[i]]=rr[[;;,i]] . varBasis,{i,1,Length[LVTemp]}];

	Return[LVTemp];
]


myPrint[args__,{style__}]:=Print[Row[{args},BaseStyle->{style}]]


(*
	Old function: Should maybe be removed.
	This function finds linear dependencies between the variables in list.
	Only numerical factors are included.
	The result is given in Mat. Where Mat[[i,j]]=0 if list[[i]]!=N list[[j]] where N is a number.
	Otherwise Mat[[i,j]]=N if list[[i]]=Nlist[[j]]
*)
OverallFac2[list_,varMat_]:=Module[{L=list,Mat=varMat},

	TotVar=L//Variables;
	DoneList=ConstantArray[0,1];

	For[i=1,i<Length[list],i++,
		Var=list[[i]]//Variables;
		varCompliment=Complement[TotVar,Var];
		For[j=i+1,j<Length[list]+1,j++,
			If[!MemberQ[DoneList, j],
				If[!CheckVariables[list[[j]],varCompliment],
					h=CompExp3[list[[i]],list[[j]]];
					Mat[[i,j]]=h;
					If[h!=0,AppendTo[DoneList,j]];
				];

			];
		];
	];
	Mat
]


(*
	This function finds linear dependencies between the variables in list.
	Only numerical factors are included.
	The result is given in Mat. Where Mat[[i,j]]=0 if list[[i]]!=N list[[j]] where N is a number.
	Otherwise Mat[[i,j]]=N if list[[i]]=Nlist[[j]]
*)
	RelationsBVariables2[list_,listVar_]:=Module[{L=list,LV=listVar},
	LVTemp=LV//FullSimplify;
	Mat=ConstantArray[0,{Length[Delete[L,1]],Length[Delete[L,1]]}];(*Linear-dependency matrix*)
	Mat1=OverallFac2[Delete[L,1],Mat];
	TempVar=ConstantArray[0,{Length[LV]}];
	For[i=1,i<Length[LV],i++,
		For[j=i+1,j<Length[LV]+1,j++,
			IfTemp=Mat1[[1;;i,j]]//DeleteDuplicates//DeleteCases[#,0,Infinity]&; (*Removes cases when the numerical factor is 0 or infinity*)
			If[Length[IfTemp]==1&&Mat1[[i,j]]!=0,
				LVTemp[[j]]=Mat1[[i,j]]LV[[i]](*Creates a list with all linear relations*)
			];
		];
	];
	LVTemp
]


CheckVariables[a_,Vars_]:=MemberQ[Boole@(MemberQ[a, #, {0, -1}, Heads -> True]&/@Vars),1, {0, -1}, Heads -> True];


(*
	Checks if two variables are identical up to a numerical factor
	For example if a=x^2 y and b=y, then there's no relation, and the function returns nothing.
	 While if a=x^2 and b=5 x^2 the function returns 5.
*)
CompExp3[a0_,b0_]:=Module[{a=a0,b=b0},
	Temp=Simplify[b/a];
	If[NumericQ[Temp],Return[Temp],Return[0]];
]


(* ::Section:: *)
(*Functions for loading and printing*)


(*
	Converts an array to a saveable form
*)
ConvertToSaveFormat[tens_SparseArray]:=Module[{},
	Return[{tens["NonzeroPositions"],tens["NonzeroValues"],Dimensions[tens]}]
];


(*
	Converts imported data, as defined by ConvertToSaveFormat, to a sparse array
*)
ConvertToSparse[arr_]:=Module[{},
	Return[SparseArray[arr[[1]]->arr[[2]],arr[[3]]]]
];


(*
	Loads the position of all scalar,gauge, and fermion representations.
*)
LoadRepPositions[repPos_]:=Module[{repPosP=repPos},
	ScalarVariablesIndices=repPosP[[1]];
	GaugeIndices=repPosP[[2]];
	FermionVariablesIndices=repPosP[[3]];
];


(*
	Loads the names of gauge couplings.
*)
LoadCouplingNames[couplingNamesI_]:=Module[{couplingNamesP=couplingNamesI},
	GaugeCouplingNames=couplingNamesP;
];


(*
	Saves a model by converting all coupling-tensors to a list.
*)
SaveModelDRalgo[modelInfo_,fileName_]:=Module[{modelInfoP=modelInfo},

	PosScalar=PrintScalarRepPositions[];
	PosVector=PrintGaugeRepPositions[];
	PosFermion=PrintFermionRepPositions[];
	PosReps={PosScalar,PosVector,PosFermion};
	tensP={modelInfoP,PosReps,GaugeCouplingNames,GroupDR,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJFC,\[Mu]IJF,Ysff,YsffC};
	SaveFile={tensP[[1]],tensP[[2]],tensP[[3]],tensP[[4]]}; (*The fourth element is the group*)
	tensP=Delete[Delete[Delete[Delete[tensP,1],1],1],1];

	Do[
		AppendTo[SaveFile,ConvertToSaveFormat[i]];
	,{i,tensP}];
	
	Export[fileName,SaveFile];
];


(*
	Loads tensors that are saved by SaveModelDRalgo
*)
LoadModelDRalgo[fileName_]:=Module[{},
	arrImp=ReadList[fileName];
	InfoText=arrImp[[1]];(*The first element is the info*)
	arrImp=Delete[arrImp,1];

	LoadRepPositions[arrImp[[1]]];(*The Second element is the repPositions*)
	arrImp=Delete[arrImp,1];

	LoadCouplingNames[arrImp[[1]]];(*The Third element is the gauge-coupling names*)
	arrImp=Delete[arrImp,1];

	ImportFile={arrImp[[1]]};(*The fourth element is the group*)
	arrImp=Delete[arrImp,1];

	Do[
		AppendTo[ImportFile,ConvertToSparse[i]];
	,{i,arrImp}];
(*Prints the info text*)
	Print[Grid[{{Row[InfoText,"\n",BaseStyle->(FontFamily->"Consolas")]}},Alignment->{Left,Center}]];
(**)
	Return[ImportFile]
];


(* ::Title:: *)
(*Matrix elements*)


(* ::Section:: *)
(*Export functions for different formats*)


(* ::Subsubsection:: *)
(*json matrix elements functions*)


makeJsonMatrixElements::usage="makeJsonMatrixElements[particles,parameters,results] converts a list of particle names {'Phi',...}, a list of particle parameters {g,...}, and a list of matrix elements results in the form {M[0,0,0,0]->g^4 s/t,...} to a JSON object in a standard format.";
makeJsonMatrixElements[particles_,parameters_,resultsI_]:=Module[
{
	particlesJson,matrixElementsJson,toString,getRelevantParameters,
	replaceSpecials,results=resultsI
},

	toString[arg_]:=If[StringQ[arg],arg,ToString[arg,InputForm]];
	replaceSpecials[arg_]:=StringReplace[arg,{"Pi"->"_pi","sReplace"->"_s","tReplace"->"_t","uReplace"->"_u"}];
	getRelevantParameters[arg_]:=Select[parameters,Not[FreeQ[arg,#]]&];
	particlesJson=Table[<|"index"->i-1,"name"->toString[particles[[i]]]|>,{i,1,Length[particles]}];
	results=results/.{s->sReplace,t->tReplace,u->uReplace};
	
	matrixElementsJson=Map[<|
		"externalParticles"->#[[1]]/.M[a__]->List[a],
		"parameters"->Map[toString,getRelevantParameters[#[[2]]]],
		"expression"->replaceSpecials[toString[PrintNonPrivate[#[[2]]]]]|>&,
		results];
		
	Return[<|"particles"->particlesJson,"matrixElements"->matrixElementsJson|>]
];


testJsonMatrixElements::usage="testJsonMatrixElements[json] tests if a JSON object is of the expected form for exporting matrix elements.";
testJsonMatrixElements[json_]:=Module[{testBool,returnString,nParticles,expectedForm},
testBool=True;
returnString="Json object matches expected schema";
(* checking head *)
If[Head[json]!=Association,
returnString="Not Association";testBool=False];
(* checking dimensions *)
If[Dimensions[json]!={2},
returnString="Dimensions not {2}";testBool=False];
(* checking top level keys *)
If[Keys[json]!={"particles","matrixElements"},
returnString="Top level keys not {'particles','matrixElements'}";testBool=False];
(* checking lower level keys *)
If[Keys[json["particles"][[1]]]!={"index","name"},
returnString="'particles' keys not {'index','name'}";testBool=False];
If[Keys[json["matrixElements"][[1]]]!={"externalParticles","parameters","expression"},
returnString="'matrixElements' keys not {'externalParticles','parameters','expressions'}";testBool=False];
(* returning results *)
{testBool, returnString}
]


splitJsonMatrixElements::usage="splitJsonMatrixElements[json] splits a JSON object containing matrix elements into a list {particleNames,parameters,results}.";
splitJsonMatrixElements[json_]:=Module[
{
	particles,matrixElements,particleIndices,particleNames,
	matrixElementIndices,matrixElementParameters,matrixElementExpressions,
	parameters,expressions,results
}
,
particles=json["particles"];
particleIndices=Map[#["index"]&,json["particles"]];
particleNames=Map[#["name"]&,json["particles"]];
matrixElements=json["matrixElements"];
matrixElementIndices=Map[#["externalParticles"]&,json["matrixElements"]];
matrixElementParameters=Map[#["parameters"]&,json["matrixElements"]];
matrixElementExpressions=Map[#["expression"]&,json["matrixElements"]];
parameters=Map[ToExpression,DeleteDuplicates[Flatten[matrixElementParameters]]];
expressions=Map[ToExpression,StringReplace[matrixElementExpressions,{RegularExpression["(\\W)_s"]->"$1s",RegularExpression["(\\W)_t"]->"$1t",RegularExpression["(\\W)_u"]->"$1u"}]];
results=Thread[matrixElementIndices->expressions];
results=Map[M[#[[1]]/.List->Sequence]->#[[2]]&,results];
{particleNames,parameters,PrintNonPrivate[results]}
];



ExportTo["json"][MatrixElement_,ParticleName_,UserCouplings_,file_]:=Block[{toExportJson},
	
(*Formatting the matrix elements*)
	toExportJson=makeJsonMatrixElements[ParticleName,UserCouplings,MatrixElement];
(*Exporting the result*)
	exportJsonMatrixElements[StringJoin[file,".json"],toExportJson];	
	Print["Results have been exported to: ", StringJoin[file,".json"]];	
]


(* reading JSON matrix elements *)
importJSONMatrixElements::usage="importJSONMatrixElements[file] imports a JSON file of matrix elements into a JSON object.";
importJSONMatrixElements[file_]:=Import[file,"RawJSON"];
(* export JSONMatrixElements *)
exportJsonMatrixElements::usage="exportJsonMatrixElements[file,jsonMatrixElements] exports a JSON object of matrix elements into a JSON file.";
exportJsonMatrixElements[file_,jsonMatrixElements_]:=Module[{test},
	If[Not[StringQ[file]],Print["File must be a string"];Return[]];
	If[StringTake[file,-5]!=".json",Print["File must end in .json"];Return[]];
	Export[file,jsonMatrixElements]
];


(* ::Subsubsection:: *)
(*hdf5 matrix elements functions*)


ExportTo["hdf5"][Cij_,OutOfEqParticles_,ParticleName_,UserCouplings_,file_]:=Block[{ExportH5,writeData,CijName,CijExport,ParticleInfo,CouplingInfo},
	
(*Metadata*)
	ParticleInfo=Table[{ToString[OutOfEqParticles[[i]]-1],ParticleName[[i]]},{i,Length[OutOfEqParticles]}];
	AppendTo[ParticleInfo,{ToString[Length[OutOfEqParticles]],"LightParticle"}];
	CouplingInfo=Table[{ToString[UserCouplings[[i]]],ToString[Symbol["c"][i-1]]},{i,1,Length[UserCouplings]}];
	
	CijExport=Cij;
	Do[CijExport[[i,j]]=Table[MatrixElemToC@k,{k,Cij[[i,j]]}];,
		{i,OutOfEqParticles},{j,OutOfEqParticles}];
(*In the hdf5 file we separate them into Cij components*)
	ExportH5=Reap[Do[
		CijName=StringJoin["MatrixElements",ParticleName[[i]],ParticleName[[j]]];
		Sow[
			writeData=Table[{ToString[FortranForm[PrintNonPrivate[a[[1]]]]],ToString[FortranForm[PrintNonPrivate[a[[2]]]]]},{a,CijExport[[i,j]]}];
			If[Length[CijExport[[i,j]]]==0,writeData=""];
			CijName -> {"Data" -> writeData}
			];
		,{i,OutOfEqParticles},{j,OutOfEqParticles}]];
	
	ExportH5=Flatten[ExportH5[[2]][[1]]];

(*Adding metadata*)
	AppendTo[ExportH5,"ParticleInfo"->{"Data"->ParticleInfo}];
	AppendTo[ExportH5,"CouplingInfo"->{"Data"->CouplingInfo}];
	
(*Exporting the reult*)
	Export[StringJoin[file,".hdf5"],ExportH5];
	Print["Results have been exported to: ", StringJoin[file,".hdf5"]];	
]


(* ::Subsubsection:: *)
(*txt matrix elements functions*)


ExportTo["txt"][MatrixElements_,OutOfEqParticles_,ParticleName_,UserCouplings_,file_]:=Block[
{
	ParticleInfo,CouplingInfo,ExportTXT,matrixElementsTXT,replaceSpecials,toString,
	sReplace,tReplace,uReplace
},

	(*Creating some metadata*)
		ParticleInfo=Table[{ToString[OutOfEqParticles[[i]]-1],ParticleName[[i]]},{i,Length[OutOfEqParticles]}];
		AppendTo[ParticleInfo,{ToString[Length[OutOfEqParticles]],"LightParticle"}];
		CouplingInfo=Table[{ToString[UserCouplings[[i]]],ToString[Symbol["c"][i-1]]},{i,1,Length[UserCouplings]}];
		
		ExportTXT=MatrixElements/.{s->sReplace,t->tReplace,u->uReplace};
		
		toString[arg_]:=If[StringQ[arg],arg,ToString[arg,InputForm]];
		replaceSpecials[arg_]:=StringReplace[arg,{"Pi"->"_pi","sReplace"->"_s","tReplace"->"_t","uReplace"->"_u"}];
	
		matrixElementsTXT=Map[
		toString[PrintNonPrivate[#[[1]]]]<>" -> "<>replaceSpecials[toString[PrintNonPrivate[#[[2]]]]]&,
		ExportTXT];
		
	(*Adding metadata to the matrix elements*)
		PrependTo[matrixElementsTXT,ParticleInfo];
		PrependTo[matrixElementsTXT,CouplingInfo];	
	
(*Exporting*)
	Export[StringJoin[file,".txt"],matrixElementsTXT];
	Print["Results have been exported to: ", StringJoin[file,".txt"]];		
]


(* ::Section:: *)
(*Exporting the results*)


PrintNonPrivate[PrivateExpression_]:=ToExpression[StringReplace[ToString[StandardForm[PrivateExpression]],"WallGoMatrix`Private`"->""]];
ReplaceMandelStam[Expression_]:=StringReplace[ToString[Expression],{"s"->"_s","t"->"_t","u"->"_u"}];


Options[ExportMatrixElements]={
	Replacements->{},
	NormalizeWithDOF->True,
	Format->"none"};


ExportMatrixElements[file_,particleList_,UserMasses_,UserCouplings_,ParticleName_,ParticleMasses_,OptionsPattern[]]:=
Block[
{
	ParticleMassesI=ParticleMasses,ExportTXT,ExportH5,
	Cij,ParticleInfo,LightParticles,particleListFull,
	CouplingInfo,MatrixElements,OutOfEqParticles,RepMasses,RepCouplings,
	FormatOptions,userFormat,MatrixElementsList,userParameters
},

(*Specifies whether the first particle should be normalized by the number of degrees of freedom*)
	normalizeDOF = OptionValue[NormalizeWithDOF];

(*Splits ParticleList into out-of-eq and light particles*)
	ExtractLightParticles[particleList,OutOfEqParticles,particleListFull,LightParticles];

(*Creates an assumption rule for simplifying Conjugate[....] terms*)
	VarAsum=#>0&/@Variables@Normal@{Ysff,gvss,gvff,gvvv,\[Lambda]4,\[Lambda]3,ParticleMasses,s,t,u}; (*All variables are assumed to be real*)
	
(*Allocates one element for each species mass to avoid errors*)	
	If[ParticleMasses[[1]]=={},ParticleMassesI[[1]]={msq}];
	If[ParticleMasses[[2]]=={},ParticleMassesI[[2]]={msq}];
	If[ParticleMasses[[3]]=={},ParticleMassesI[[3]]={msq}];

(*Extracting all matrix elements*)	
	GenerateMatrixElements[MatrixElements,Cij,particleListFull,LightParticles,ParticleMassesI,OutOfEqParticles];
	MatrixElementsList=Table[MatrixElemToC@i//.OptionValue[Replacements],{i,MatrixElements}]; (*Creates a replacement list and shifts the indices to start at 0.*)

(*Exporting the matrix elements to the choosen format*)
	FormatOptions = {"txt", "json", "hdf5", "all", "none"};
	userFormat = OptionValue[Format];
	
	(* Convert userFormat to a list if it's a single value *)
	userFormatsList = If[ListQ[userFormat], userFormat, {userFormat}];
	
	(* Replace "all" with all formats if it is present *)
	If[MemberQ[userFormatsList, "all"], userFormatsList = {"txt", "json", "hdf5"}];
	
	(* Check if all formats in the list are valid *)
	If[AllTrue[userFormatsList, MemberQ[FormatOptions, #] &],
	   (* Iterate over each format and perform the export *)
	   Do[
	     Switch[fmt,
	       "txt", ExportTo["txt"][MatrixElementsList, OutOfEqParticles, ParticleName, UserCouplings, file],
	       "hdf5", ExportTo["hdf5"][Cij, OutOfEqParticles, ParticleName, UserCouplings, file],
	       "json",
	       userParameters = Flatten[Join[UserCouplings, ParticleMasses]] // DeleteDuplicates;
	       ExportTo["json"][MatrixElementsList, ParticleName, userParameters, file]
	     ],
	     {fmt, userFormatsList}
	   ],
	   Print["The currently allowed formats are: txt, hdf5, and json"];
	   Print["Please choose a valid format"];
	];
	
	Return[PrintNonPrivate[MatrixElementsList]]
];


MatrixElemToC[MatrixElem_]:=Block[{Ind},
	Ind=MatrixElem[[2]];
	
	Return[M[Ind[[1]]-1,Ind[[2]]-1,Ind[[3]]-1,Ind[[4]]-1]->MatrixElem[[1]]]
]


MatrixElemFromC[MatrixElem_]:=Block[{Ind},
	Ind=MatrixElem[[2]];
	
	Return[M[Ind[[1]]+1,Ind[[2]]+1,Ind[[3]]+1,Ind[[4]]+1]->MatrixElem[[1]]]
]


ImportMatrixElements[file_]:=Module[
{
	jsonObject,particleNames,parameters,results
},
	(* Only implemented for JSON *)
	If[StringTake[file,-5]!=".json",Print["File must end in .json"];Return[]];

	(* Importing into JSON object *)
	jsonObject=importJSONMatrixElements[file];

	(* Splitting JSON object *)
	{particleNames,parameters,results}=splitJsonMatrixElements[jsonObject];

	(* Returning results *)
	Return[{particleNames,parameters,results}]
];


(* ::Section:: *)
(*Private constants*)


{ZLij,GvvssTSS,\[Lambda]3DSS,\[Mu]ijSSLO,\[Mu]ijSSNLO,\[Mu]ijSNLOSS,\[Lambda]3DSS,\[Lambda]KVecTSS,\[Lambda]3CTot,\[Lambda]3CSSS,ZSij,\[Mu]ijSSNLO2,Ggvvv};


{TadPoleS,ContriTadPoleSoftToHard,GgvvvSS,\[Lambda]4SMod,\[Lambda]4Tot,IdentMatPre};


{\[Beta]gvff,Zgvff,\[Mu]ijEP,gvvvEP,gvssEP,\[Lambda]4EP,\[Lambda]3EP,nsEP,nvEP,CT,HelpSolveEffectiveHardM};


{aS3D,ZijS,aV3D,ZabT,ZabL,\[Lambda]3D,GvvssL,GvvssT,\[Lambda]AA,\[Lambda]3CS,\[Mu]SijNLO,GvvsL};(*DimRed Results*)


{\[CapitalLambda]\[Lambda],\[CapitalLambda]g,Hg,Habij,HabIJF,HabIJFC,Ysij,YsijC,YTemp,YTempC,Yhelp,YhelpC};(*Private Variables*)


{\[Gamma]ij,\[Beta]mij,\[Beta]\[Lambda]ijkl,Z\[Lambda]ijkl,\[Gamma]ab,\[Beta]vvss,Zgvvss,\[Beta]gvvv,Zgvvv,\[Gamma]IJF,\[Beta]Ysij,ZYsij,\[Beta]YsijC,ZYsijC};(*CounterTerms*)


{\[Lambda]3DS,\[Lambda]KVecT,\[Lambda]KVec,\[Lambda]AAS,IdentMat,\[Mu]ijVNLO,\[Mu]ijSNLO,\[Lambda]3CSRed,\[Lambda]1IP,NFMat,NSMat};(*Private Variables*)


{nsH,nSl,\[Lambda]K,\[Lambda]4S,\[Lambda]4K,\[Lambda]x,\[Lambda]y,gAvss,gvssL,\[Mu]ijL,\[Mu]ijLight};


{HabijL,HabijVL,HabijA,HabijVA,\[Lambda]3Cx,\[Lambda]3Cy,\[Lambda]3CLight,\[Lambda]3CHeavy,\[Mu]IJF,\[Mu]IJFC,\[Mu]VabNLO,\[Mu]abDef,GroupMathCleared}


End[]
EndPackage[]
