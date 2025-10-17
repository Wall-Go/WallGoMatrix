(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: WallGoMatrix *)

(*
		This software is covered by the GNU General Public License 3.
		Copyright (C) 2024-2025 Andreas Ekstedt
		Copyright (C) 2024-2025 Oliver Gould
		Copyright (C) 2024-2025 Joonas Hirvonen
		Copyright (C) 2024-2025 Benoit Laurent
		Copyright (C) 2024-2025 Lauri Niemi
		Copyright (C) 2024-2025 Philipp Schicho
		Copyright (C) 2024-2025 Jorinde van de Vis
*)

(* :Summary:	WallGoMatrix constructs matrix elements for generic models at
				leading logarithmic (LL) and partially next-to-leading logarithmic (NLL) orders.
				It automates the computation of these elements by contracting coupling tensors
				with the appropriate kinematic factors.
*)	

(* ------------------------------------------------------------------------ *)


BeginPackage["WallGo`WallGoMatrix`"]


$WallGoMatrixOutput::usage =
"$WallGoMatrixOutputFlag contains std output.";


$WallGoMatrixOutput = $Output;
If[$ScriptCommandLine=={},Null,
	$Output=OpenWrite["/dev/null"];
];


If[! ValueQ[$GroupMathMultipleModels],
	$GroupMathMultipleModels = False];
If[! ValueQ[$LoadGroupMath],
	$LoadGroupMath = True];
If[! ValueQ[$InstallGroupMath],
	$InstallGroupMath = False];


(*
	Welcome banner
*)
authors[1] = "\!\(\*TemplateBox[{RowBox[{\"Andreas\", \" \", \"Ekstedt\"}], {URL[\"https://inspirehep.net/authors/1799400\"], None}, \"https://inspirehep.net/authors/1799400\", \"HyperlinkActionRecycled\", {\"HyperlinkActive\"}, BaseStyle -> {\"Hyperlink\"}, HyperlinkAction -> \"Recycled\"},\"HyperlinkTemplate\"]\)";
authors[2] = "\!\(\*TemplateBox[{RowBox[{\"Oliver\", \" \", \"Gould\"}], {URL[\"https://inspirehep.net/authors/1606373\"], None}, \"https://inspirehep.net/authors/1606373\", \"HyperlinkActionRecycled\", {\"HyperlinkActive\"}, BaseStyle -> {\"Hyperlink\"}, HyperlinkAction -> \"Recycled\"},\"HyperlinkTemplate\"]\)";
authors[3] = "\!\(\*TemplateBox[{RowBox[{\"Joonas\", \" \", \"Hirvonen\"}], {URL[\"https://inspirehep.net/authors/1866901\"], None}, \"https://inspirehep.net/authors/1866901\", \"HyperlinkActionRecycled\", {\"HyperlinkActive\"}, BaseStyle -> {\"Hyperlink\"}, HyperlinkAction -> \"Recycled\"},\"HyperlinkTemplate\"]\)";
authors[4] = "\!\(\*TemplateBox[{RowBox[{\"Benoit\", \" \", \"Laurent\"}], {URL[\"https://inspirehep.net/authors/1808372\"], None}, \"https://inspirehep.net/authors/1808372\", \"HyperlinkActionRecycled\", {\"HyperlinkActive\"}, BaseStyle -> {\"Hyperlink\"}, HyperlinkAction -> \"Recycled\"},\"HyperlinkTemplate\"]\)";
authors[5] = "\n\!\(\*TemplateBox[{RowBox[{\"Lauri\", \" \", \"Niemi\"}], {URL[\"https://inspirehep.net/authors/1764675\"], None}, \"https://inspirehep.net/authors/1764675\", \"HyperlinkActionRecycled\", {\"HyperlinkActive\"}, BaseStyle -> {\"Hyperlink\"}, HyperlinkAction -> \"Recycled\"},\"HyperlinkTemplate\"]\)";
authors[6] = "\!\(\*TemplateBox[{RowBox[{\"Philipp\", \" \", \"Schicho\"}], {URL[\"https://inspirehep.net/authors/1639147\"], None}, \"https://inspirehep.net/authors/1639147\", \"HyperlinkActionRecycled\", {\"HyperlinkActive\"}, BaseStyle -> {\"Hyperlink\"}, HyperlinkAction -> \"Recycled\"},\"HyperlinkTemplate\"]\)";
authors[7] = "\!\(\*TemplateBox[{RowBox[{\"Jorinde\", \" \", \"van de Vis\"}], {URL[\"https://inspirehep.net/authors/1589979\"], None}, \"https://inspirehep.net/authors/1589979\", \"HyperlinkActionRecycled\", {\"HyperlinkActive\"}, BaseStyle -> {\"Hyperlink\"}, HyperlinkAction -> \"Recycled\"},\"HyperlinkTemplate\"]\)";


(* Helper to split a long string at word boundaries near 80 characters *)
wrapText[text_String, maxLen_: 80] := Module[
  {words, lines = {}, current = ""},
  words = StringSplit[text, " "];
  Do[
    If[StringLength[current <> " " <> word] <= maxLen,
      current = If[current === "", word, current <> " " <> word],
      AppendTo[lines, current];
      current = word;
    ],
    {word, words}
  ];
  AppendTo[lines, current];
  StringRiffle[lines, "\n"]
];

pacletInfo = PacletFind["WallGo/WallGoMatrix"];
WallGoMatrixLoadImage=Image[CompressedData["1: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"],"Byte",ColorSpace->"RGB",Interleaving->True];
If[pacletInfo =!= {},
	infoTable = Normal@KeyTake[
		      pacletInfo[[1]]["PacletInfo"],
		      {"Name", "Version", "Description", "WolframVersion", "Creator", "URL"}
		      ];
	,
	infoTable = {
		"Name"->"WallGo/WallGoMatrix",
		"Description"->"Computes 2-to-2 scattering matrix amplitudes for arbitrary quantum field theories.",
		"WolframVersion"->"13.+",
		"Creator"->StringRiffle[Table[authors[i], {i, 7}], ", "],
		"URL"->Hyperlink[Mouseover["github.com/Wall-Go/WallGoMatrix",Style["github.com/Wall-Go/WallGoMatrix"]],"https://github.com/Wall-Go/WallGoMatrix"]
	};
];

AppendTo[
  infoTable,
  "Reference" -> Row[{
    Hyperlink["JHEP 04 (2025) 101", "https://doi.org/10.1007/JHEP04(2025)101"],
    " \[Bullet] e-Print: ",
    Hyperlink["2411.04970 [hep-ph]", "https://arxiv.org/abs/2411.04970"]
  }]
];
AppendTo[infoTable,"Model files"->
	Hyperlink[Mouseover["WallGoMatrix/examples",Style["WallGoMatrix/examples"]],
	"https://github.com/Wall-Go/WallGoMatrix/tree/main/examples"]];


Print@Panel@Row[{
  Graphics[WallGoMatrixLoadImage, ImageSize -> 100],
  Spacer[20],
  Column[{
      Row[{
      Style["GOGOGOGOGOGOGOGOGOGOGOGO ", Bold],
      Style["WallGoMatrix", Bold, Larger],
      Style[" GOGOGOGOGOGOGOGOGOGOGGOGOGO", Bold]
    }],
    TableForm[
      infoTable /. {
        Rule["Description", desc_] :> {"Description:", wrapText[desc]},
        Rule["Creator", desc_] :> {"Creator:", wrapText[desc]},
        Rule[k_, v_] :> {k <> ":", v}
      },
      TableAlignments -> Left
    ],
    Style["GOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGGOGOGOGOGOGOG",Bold]
  }, Spacings -> 1]
}];


(*List of public functions*)
ImportModel::usage = 
"ImportModel[Group, gvvv, gvff, gvss, \[Lambda]1, \[Lambda]3, \[Lambda]4, \[Mu]ij, \[Mu]IJ, \[Mu]IJC, Ysff, YsffC, Verbose -> False, Mode -> ModeNumber] \
loads the model and creates helper tensors.

Arguments:
  - Group: The symmetry group of the model.
  - gvvv: Structure constants.
  - gvff: Vector-fermion trilinear couplings.
  - gvss: Vector-scalar trilinear couplings.
  - \[Lambda]1: Scalar tadpole couplings.
  - \[Lambda]3: Cubic couplings.
  - \[Lambda]4: Quartic couplings.
  - \[Mu]ij: Scalar mass matrix.
  - \[Mu]IJ: Fermion mass matrices.
  - \[Mu]IJC: Extended fermion mass matrices.
  - Ysff: Yukawa couplings.
  - YsffC: Conjugate Yukawa couplings.

Options:
  - Verbose: If True, enables verbose output (default is False).
  - Mode: Specifies the operating mode using ModeNumber.

This function organizes the parameters necessary for model loading and tensor generation.";


LatexStyle[text_] := StyleBox[text, FontFamily -> "Times New Roman", FontWeight -> "Regular", FontSlant -> "Italic"];


SymmetryBreaking::usage = 
"SymmetryBreaking classifies different scalar, fermion, and vector representations into their respective particles based on VEV-induced masses.

Usage:
SymmetryBreaking[\!\(\*
StyleBox[\"vev\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)];

Arguments:
- \!\(\*
StyleBox[\"vev\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): \
List of vacuum expectation values

Output:
- Each particle i is represented as {\
\!\(\*StyleBox[\"ri\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"mi\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)}, \
where:
    - \!\(\*StyleBox[\"ri\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): \
Label for the particle's representation.
    - \!\(\*StyleBox[\"mi\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): \
Label for the particle's mass within that representation.

Options:
- VevDependentCouplings (Boolean): If True, allows for VEV-dependent couplings (default is False).

This function allows the classification of particles by type and mass using the representation-breaking mechanism provided by VEVs.\
";


CreateParticle::usage = 
"Groups \
\!\(\*StyleBox[\"n\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\) \
particles into a single representation.

Usage:
CreateParticle[{\
\!\(\*StyleBox[\"r1\", \"TI\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)\
,...,\
\!\(\*StyleBox[\"rn\", \"TI\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)}}, \
\"\!\(\*StyleBox[\"R\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)\", \
\!\(\*StyleBox[\"M\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\), \
\"\!\(\*StyleBox[\"particleName\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)\"\
]
CreateParticle[{{\!\(\*
StyleBox[\"r1\", \"TI\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"m1\", \"TI\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)},...,{\!\(\*
StyleBox[\"rn\", \"TI\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"mn\", \"TI\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)}}, \"\!\(\*
StyleBox[\"R\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)\", \!\(\*
StyleBox[\"M\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\), \"\!\(\*
StyleBox[\"particleName\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)\"]

Arguments:
- {{\!\(\*
StyleBox[\"r1\", \"TI\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"m1\", \"TI\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)},...,{\!\(\*
StyleBox[\"rn\", \"TI\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"mn\", \"TI\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)}}:
    List of particle indices (obtainable via PrintFieldRepPositions[]), where each pair {\!\(\*
StyleBox[\"ri\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"mi\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)} represents,
    - \!\(\*
StyleBox[\"ri\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): Label for the particle's representation,
    - \!\(\*
StyleBox[\"mi\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): Label for the particle's mass within that representation.
- \"\!\(\*
StyleBox[\"R\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)\": Specifies the type of particle representation:
    - \"F\" for Fermions
    - \"V\" for Vector Bosons
    - \"S\" for Scalars
- \!\(\*
StyleBox[\"M\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): The mass parameter of the grouped particle representation.
- \"\!\(\*
StyleBox[\"particleName\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)\": A string specifying the name of the particle.

Return:
- Returns a list formatted as {{\!\(\*
StyleBox[\"fieldIndices\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)}, \!\(\*
StyleBox[\"R\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\) = {V, F, S}, \!\(\*
StyleBox[\"M\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"particleName\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)}, where:
    - \!\(\*
StyleBox[\"fieldIndices\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): Indices of the grouped particles.
    - \!\(\*
StyleBox[\"R\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): Type of particle representation (Vector Boson, Fermion, Scalar).
    - \!\(\*
StyleBox[\"M\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): The mass of the grouped particle.
    - \!\(\*
StyleBox[\"particleName\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): The name of the grouped particle.

This function is useful for grouping particles into one collective representation based on their properties and type.";


TruncateAtLeadingLogarithm::usage =
"TruncateAtLeadingLogarithm[MatrixElements] classifies and truncates a list of \
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

 \[Bullet] Power-enhanced structures are tagged by \!\(\*
StyleBox[\"powerEnhanced\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\".\",\nFontWeight->\"Plain\",\nFontSlant->\"Italic\"]\)
   Logarithmically-enhanced structures are tagged by \!\(\*
StyleBox[
StyleBox[
RowBox[{\"l\", 
StyleBox[\"ogEnhanced\",\nFontWeight->\"Plain\"]}]],\nFontSlant->\"Italic\"]\).

 \[Bullet] Output: a truncated list of the form
   {{leadingLogContribution1, {i,j,k,l}}, ...},
   containing only terms contributing to leading logarithmic behaviour.\
";


ExportMatrixElements::usage = 
"Generates all possible matrix elements with the external particles specified in particleList, and exports the results to file.

Usage:
ExportMatrixElements[\!\(\*
StyleBox[\"fileName\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"particleList\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"lightParticleList\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\), OptionsPattern[]]

Arguments:
- \!\(\*
StyleBox[\"fileName\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): The name of the file to export the matrix elements to.
- \!\(\*
StyleBox[\"particleList\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): A list specifying the external particles for matrix element generation, formatted as{{\!\(\*
StyleBox[\"fieldIndices\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)},\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"R\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\) = {\"V\", \"F\", \"S\"}, \!\(\*
StyleBox[\"M\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\), \"\!\(\*
StyleBox[\"particleName\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)\"}, where:
    - \!\(\*
StyleBox[\"fieldIndices\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): Indices of the fields involved.
    - \!\(\*
StyleBox[\"R\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): Type of particle representation (Vector Boson, Fermion, Scalar).
    - \!\(\*
StyleBox[\"M\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): Mass parameter of the particle.
    - \!\(\*
StyleBox[\"particleName\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): The name of the particle.
- \!\(\*
StyleBox[\"lightParticleList\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): A list of light particles involved in the calculations, formatted similarly to \!\(\*
StyleBox[\"particleList\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\).

Options:
- Format: Specifies the export format. Supported values include:
    - \"txt\": Plain text format
    - \"json\": JSON format
    - \"hdf5\": HDF5 format
    - \"all\": Exports in all supported formats
    - \"none\": Does not export to a file (default)
- Verbose (Boolean): If True, lists channels that are currently being computed (default is False).
- TruncateAtLeadingLog (Boolean): If True, truncates the matrix elements at the leading logarithmic order (default is True).
- TagLeadingLog (Boolean): tag leading-log terms in the output being either power enhanced (powerEnhanced) or logarithmically enhanced (logEnhanced).
- NormalizeWithDOF (Boolean): If True, matrix elements are normalized by the number of degrees of freedom of the incoming particle at index 1 (default is True).
- Replacements (List): list of replacement rules applied to the generated elements.

Explanation:
- Choosing \"all\" exports the result in all possible formats, while choosing \"none\" skips exporting.

Return:
- A list of generated matrix elements is returned by the function, regardless of the chosen export format.\
";


ImportMatrixElements::usage = 
"Imports matrix elements from file.

Usage:
ImportMatrixElements[StyleBox[\"fileName\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\), OptionsPattern[]]
\!\(\*


Arguments:
- \!\(\*
StyleBox[\"fileName\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): The name of the file to export the matrix elements to.
- \!\(\*

Return:
StyleBox[\"particles\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): A list of particle names and their associated indices.
StyleBox[\"modelParameters\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): A list of model parameters.
StyleBox[\"results\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): The matrix elements.
";


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


PrintFieldRepPositions::usage = 
"PrintFieldRepPositions[\"\!\(\*
StyleBox[\"Field\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)\"] prints the indices of field representations for the given input \"\!\(\*
StyleBox[\"Field\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)\", \
where \"\!\(\*
StyleBox[\"Field\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)\" can be one of the following representations:
  - Gauge: Represents gauge fields.
  - Fermion: Represents fermionic fields.
  - Scalar: Represents scalar fields.\
";


WallGoMatrix::failmsg =
"Fatal Error! WallGoMatrix has encountered a critical issue and must abort the computation. \
Details of the problem: `1`.";

WallGoMatrix::failWarning =
"Warning! WallGoMatrix has encountered a non-fatal issue. \
Please review the following problem: `1`.";

WallGoMatrix::missingParticles = 
 "WallGoMatrix encountered a critical issue and must terminate the computation.\n\n\
Issue: Missing particle declarations.\n\n\
The following particles need to be declared as input particles in ExportMatrixElements[]:\n\
- Vector-Type particles: `1`\n\
- Fermion-Type particles: `2`\n\
- Scalar-Type particles: `3`\n\n\
Please ensure all required particles are specified before proceeding.";


DefineDim6::usage="Defines a dimension 6 operator";
CompareInvariants::usage="\
Finds relations between couplings by calculating basis-invariant tensors";


$WallGoMatrixDirectory=DirectoryName[$InputFileName];


DownloadPackage[url_, targetName_]:=Module[{zipPath, targetFolder, targetDir, targetPath},
  
  (* Define the URL and file paths *)
  targetFolder=FileBaseName[url];
  targetDir = FileNameJoin[{$UserBaseDirectory, "Applications"}];
  zipPath = FileNameJoin[{targetDir, targetFolder<>".zip"}];
  targetPath = FileNameJoin[{targetDir, targetName}];
  
  (* Download the ZIP file *)
  If[! FileExistsQ[zipPath],
    URLSave[url, zipPath];
    Print["Downloaded "<>Last[FileNameSplit[url]]]
  ];
  
  (* Unzip the file *)
  If[DirectoryQ[targetPath],
  Print[targetPath<>" already exists. Check your "<>targetName<>" installation."];
  Abort[];
  ];
  ExtractArchive[zipPath, targetDir];
  Print["Unpacked and installed "<>FileBaseName[url]<>" in Applications folder."];
  
  (* Clean up *)
  DeleteFile[zipPath];
]


(*
	Functions from GroupMath are used to create the model.
*)
If[$LoadGroupMath,
	Unprotect[BlockDiagonalMatrix];
	If[$InstallGroupMath,
		If[
			Quiet[Check[Needs["GroupMath`"], True]],
			DownloadPackage["https://renatofonseca.net/groupmath/ProgramVersions/GroupMath-1.1.2.zip","GroupMath"];
			Print[Style["GroupMath installed.","Text", Red, Bold]];
		]
	];
	Check[
		Needs["GroupMath`"];,
		Message[Get::noopen,
			"GroupMath` at "<>ToString[$UserBaseDirectory]<>"/Applications.\n"<>
			"Set WallGoMatrix`$InstallGroupMath=True for automatic installation of GroupMath"];
		Abort[];
	];
	
	Print["GroupMath is an independent package, and is not part of WallGoMatrix."];
	Print["Please cite GroupMath: Comput.Phys.Commun. 267 (2021) 108085 \[Bullet] e-Print: \!\(\*TemplateBox[{RowBox[{\"2011.01764\", \" \", \"[\", RowBox[{\"hep\", \"-\", \"th\"}], \"]\"}], {URL[\"https://arxiv.org/abs/2011.01764\"], None}, \"https://arxiv.org/abs/2011.01764\", \"HyperlinkActionRecycled\", {\"HyperlinkActive\"}, BaseStyle -> {\"Hyperlink\"}, HyperlinkAction -> \"Recycled\"},\n\"HyperlinkTemplate\"]\)"];

]


Print["WallGoMatrix is powered by the DRalgo ModelCreation."];
Print["Please cite \!\(\*TemplateBox[{\"DRalgo\", {URL[\"https://github.com/DR-algo/DRalgo\"], None}, \"https://github.com/DR-algo/DRalgo\", \"HyperlinkActionRecycled\", {\"HyperlinkActive\"}, BaseStyle -> {\"Hyperlink\"}, HyperlinkAction -> \"Recycled\"},\n\"HyperlinkTemplate\"]\): Comput.Phys.Commun. 288 (2023) 108725 \[Bullet] e-Print: \!\(\*TemplateBox[{RowBox[{\"2205.08815\", \" \", \"[\", RowBox[{\"hep\", \"-\", \"ph\"}], \"]\"}], {URL[\"https://arxiv.org/abs/2205.08815\"], None}, \"https://arxiv.org/abs/2205.08815\", \"HyperlinkActionRecycled\", {\"HyperlinkActive\"}, BaseStyle -> {\"Hyperlink\"}, HyperlinkAction -> \"Recycled\"},\n\"HyperlinkTemplate\"]\)"];


(*
	Verbose=True removes progress messages.
	Mode=2 calculates everything,
	Mode=1 only calculates LO masses and couplings
	Mode=0 only calculates LO masses
*)
Options[ImportModel]={
	Verbose->False,
	Mode->2,
	Dim6->False}


Begin["`Private`"]


(*
	Loads all functions.
*)
(*Loads DRalgo model creation*)
Get[FileNameJoin[{$WallGoMatrixDirectory,"modelCreation.m"}]];
(*Loads matrix element creation*)
Get[FileNameJoin[{$WallGoMatrixDirectory,"matrixElements.m"}]];


(* ::Section:: *)
(*Model: Initialization*)


(*
	Defines internal tensors from the loaded model. Also creates help-tensors used for
	intermediate calculations.
*)
ImportModel[GroupI_,gvvvI_,gvffI_,gvssI_,\[Lambda]1I_,\[Lambda]3I_,\[Lambda]4I_,\[Mu]ijI_,\[Mu]IJFI_,\[Mu]IJFCI_,YsffI_,YsffCI_, OptionsPattern[]]:=
Module[
{
	GroupP=GroupI,
	gvvvP=gvvvI,gvffP=gvffI,gvssP=gvssI,\[Lambda]1IP=\[Lambda]1I,\[Lambda]3P=\[Lambda]3I,\[Lambda]4P=\[Lambda]4I,
	\[Mu]ijP=\[Mu]ijI,\[Mu]IJFP=\[Mu]IJFI,\[Mu]IJFCP=\[Mu]IJFCI,YsffP=YsffI,YsffCP=YsffCI
},

If[$LoadGroupMath,
	If[!GroupMathCleared && !ValueQ[$GroupMathMultipleModels],
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
(*End of Options*)

	CreateHelpTensors[] (*Creates recurring tensors*)
];


(*
	Defines a \[Phi]^6 operator
*)
DefineDim6[\[Lambda]6I_]:=Module[{\[Lambda]6P=\[Lambda]6I},
	If[mode>=3,
		\[Lambda]6=\[Lambda]6P//SparseArray;
	,
		Message[WallGoMatrix::failmsg, "Please set mode=3 to use this feature"];
		Abort[]
	];
];


(* ::Section:: *)
(*Model: Functions for loading and printing*)


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
SaveModel[modelInfo_,fileName_]:=Module[{modelInfoP=modelInfo},

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
	Loads tensors that are saved by SaveModel
*)
LoadModel[fileName_]:=Module[{},
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


(* ::Section:: *)
(*MatrixElements: Export functions for different formats*)


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
	particles,matrixElements,particleIndices,particleNames,particleRules,
	matrixElementIndices,matrixElementParameters,matrixElementExpressions,
	parameters,expressions,results
}
,
particles=json["particles"];
particleIndices=Map[#["index"]&,json["particles"]];
particleNames=Map[#["name"]&,json["particles"]];
particleRules=Map[{#["name"]->#["index"]}&,json["particles"]];
matrixElements=json["matrixElements"];
matrixElementIndices=Map[#["externalParticles"]&,json["matrixElements"]];
matrixElementParameters=Map[#["parameters"]&,json["matrixElements"]];
matrixElementExpressions=Map[#["expression"]&,json["matrixElements"]];
parameters=Map[ToExpression,DeleteDuplicates[Flatten[matrixElementParameters]]];
expressions=Map[ToExpression,StringReplace[matrixElementExpressions,{RegularExpression["(\\W)_s"]->"$1s",RegularExpression["(\\W)_t"]->"$1t",RegularExpression["(\\W)_u"]->"$1u"}]];
results=Thread[matrixElementIndices->expressions];
results=Map[M[#[[1]]/.List->Sequence]->#[[2]]&,results];
{particleRules,parameters,PrintNonPrivate[results]}
];



ExportTo["json"][MatrixElement_,ParticleName_,UserCouplings_,file_]:=Block[
{
	toExportJson,
	particleName=ParticleName
},
	
(*Formatting the matrix elements*)
	toExportJson=makeJsonMatrixElements[particleName,UserCouplings,MatrixElement];
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


ExportTo["hdf5"][Cij_,OutOfEqParticleList_,ParticleName_,UserCouplings_,file_]:=
Block[{
	ExportH5,writeData,CijName,CijExport,
	OutOfEqParticles,
	ParticleInfo,CouplingInfo
},
	OutOfEqParticles=Range[Length[OutOfEqParticleList]];
	
(*Metadata*)
	ParticleInfo=Table[{ToString[OutOfEqParticles[[i]]-1],ParticleName[[i]]},{i,Length[OutOfEqParticles]}];
	CouplingInfo=Table[{ToString[UserCouplings[[i]]],ToString[Symbol["c"][i-1]]},{i,1,Length[UserCouplings]}];
	
	CijExport=Cij;
	Do[
		CijExport[[i,j]]=Table[MatrixElemToC@k,{k,Cij[[i,j]]}];,
		{i,OutOfEqParticles},{j,OutOfEqParticles}
	];

(*In the hdf5 file we separate them into Cij components*)
	ExportH5=Reap[Do[
		CijName=StringJoin["MatrixElements",ParticleName[[i]],ParticleName[[j]]];
		Sow[
			writeData=Table[{
				ToString[FortranForm[PrintNonPrivate[a[[1]]]]],
				ToString[FortranForm[PrintNonPrivate[a[[2]]]]]
				},{a,CijExport[[i,j]]}];
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


ExportTo["txt"][MatrixElements_,particleListAll_,UserCouplings_,file_]:=
Block[{
	ParticleInfo,CouplingInfo,
	particleNames,
	ExportTXT,matrixElementsTXT,replaceSpecials,toString,
	sReplace,tReplace,uReplace
},
	particleNames=particleListAll[[All,4]];

	(*Creating some metadata*)
		ParticleInfo=Table[{ToString[i-1],particleNames[[i]]},{i,Length[particleListAll]}];
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
(*MatrixElements: Exporting the results*)


PrintNonPrivate[PrivateExpression_]:=ToExpression[StringReplace[
	ToString[StandardForm[PrivateExpression]],
	{
		"WallGoMatrix`Private`"->"",
		"WallGo`"->"",
		"Global`"->""
	}]];
ReplaceMandelStam[Expression_]:=StringReplace[
	ToString[Expression],{"s"->"_s","t"->"_t","u"->"_u"}
	];


Options[GenerateMatrixElementsData]={
	Replacements->{},
	NormalizeWithDOF->True,
	TruncateAtLeadingLog->True,
	TagLeadingLog->False,
	Verbose->False};


Options[ExportMatrixElements]={
	Options[GenerateMatrixElementsData],
	Format->"none"}//Flatten;


extractParticleMasses[particleList_, type_, length_] := 
  Module[{resultVector = ConstantArray[Null, length]},
    Do[
      resultVector[[pos]] = entry[[3]],
      {entry, Select[particleList, #[[2]] === type &]},
      {pos, entry[[1]]}
    ];
(*    If[MemberQ[resultVector, Null],
    Message[WallGoMatrix::failmsg,
			"Missing particle mass declarations."];
		Abort[];
	];*)
	
    Return[resultVector];
  ];


GenerateMatrixElementsData::usage = 
"Computes all 2\[RightArrow]2 scattering matrix elements for the given sets of \
out-of-equilibrium and light particles. It returns an association containing \
the symbolic matrix elements, coupling tensors, particle information, and \
mass parameters, which can be exported later using ExportMatrixElements.

Usage:
GenerateMatrixElementsData[\!\(\*
StyleBox[\"particleList\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"lightParticleList\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\), OptionsPattern[]]

Arguments:
- \!\(\*
StyleBox[\"particleList\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): A list specifying the external particles for matrix element generation, formatted as{{\!\(\*
StyleBox[\"fieldIndices\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)},\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"R\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\) = {\"V\", \"F\", \"S\"}, \!\(\*
StyleBox[\"M\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\), \"\!\(\*
StyleBox[\"particleName\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\)\"}, where:
    - \!\(\*
StyleBox[\"fieldIndices\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): Indices of the fields involved.
    - \!\(\*
StyleBox[\"R\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): Type of particle representation (Vector Boson, Fermion, Scalar).
    - \!\(\*
StyleBox[\"M\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): Mass parameter of the particle.
    - \!\(\*
StyleBox[\"particleName\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): The name of the particle.
- \!\(\*
StyleBox[\"lightParticleList\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\): A list of light particles involved in the calculations, formatted similarly to \!\(\*
StyleBox[\"particleList\",\nFontFamily->\"Times New Roman\",\nFontWeight->\"Regular\",\nFontSlant->\"Italic\"]\).

The returned association includes the following keys:
  \"Cij\" \[RightArrow] coupling tensor representation
  \"MatrixElementsList\" \[RightArrow] list of computed matrix elements (after replacements)
  \"particleListAll\" \[RightArrow] combined particle list
  \"particleNames\" \[RightArrow] particle names
  \"userCouplings\" \[RightArrow] coupling symbols appearing in the matrix elements
  \"particleMasses\" \[RightArrow] particle masses grouped by spin type

Options:
- Verbose (Boolean, default \[RightArrow] False):
	If True, lists channels that are currently being computed (default is False).
- TruncateAtLeadingLog (Boolean, default \[RightArrow] True):
	If True, truncates the matrix elements at the leading logarithmic order (default is True).
- TagLeadingLog (Boolean, default \[RightArrow] False):
	tag leading-log terms in the output being either power enhanced (powerEnhanced) or logarithmically enhanced (logEnhanced).
- NormalizeWithDOF (Boolean, default \[RightArrow] True):
	If True, matrix elements are normalized by the number of degrees of freedom of the incoming particle at index 1 (default is True).
- Replacements (List, default \[RightArrow] { }):
	list of replacement rules applied to the generated elements.

Example:
  data = GenerateMatrixElementsData[outOfEqParticles, lightParticles, 
           NormalizeWithDOF \[RightArrow] True, Verbose \[RightArrow] True];
  data[\"MatrixElementsList\"]
";


GenerateMatrixElementsData[outOfEqParticleList_, lightParticleList_, OptionsPattern[]] :=
Block[
  {
    particleMasses,
    userCouplings,
    particleNames,
    particlesAll, particlesOutOfEq,
    particleListAll,
    Cij, MatrixElements,
    MatrixElementsList,
    VarAsum,
    RepMasses, RepCouplings
  },
  
  (* Read options *)
  bNormalizeWithDOF = OptionValue[NormalizeWithDOF];
  bTruncateAtLeadingLog = OptionValue[TruncateAtLeadingLog];
  bTagLeadingLog = OptionValue[TagLeadingLog];
  bVerbose = OptionValue[Verbose];
  
  (* Split particles into out-of-eq and light *)
  ExtractLightParticles[outOfEqParticleList, lightParticleList, particleListAll];
  
  (* Separate names and particles *)
  particleNames = particleListAll[[All, 4]];
  particlesOutOfEq = outOfEqParticleList[[All, 1 ;; 2]];
  particlesAll = particleListAll[[All, 1 ;; 2]];
  
  (* Collect user couplings *)
  userCouplings = Variables@Normal@{Ysff, gvss, gvff, gvvv, \[Lambda]4, \[Lambda]3} // DeleteDuplicates;
  
  (* Construct particle mass vectors *)
  particleMasses = {
    extractParticleMasses[particleListAll, "V", Length[gvff]],
    extractParticleMasses[particleListAll, "F", Length[gvff[[1]]]],
    extractParticleMasses[particleListAll, "S", Length[gvss[[1]]]]
  };
  
  (* Real variable assumptions *)
  VarAsum = {# > 0} & /@ 
    Variables@Normal@{Ysff, gvss, gvff, gvvv, \[Lambda]4, \[Lambda]3, particleMasses, s, t, u, fsign} //
    Flatten;
  
  (* Ensure particle masses non-empty *)
  Table[
    If[particleMasses[[i]] == {}, particleMasses[[i]] = {msq}],
    {i, 1, 3}
  ];
  
  (* Generate the matrix elements *)
  GenerateMatrixElements[MatrixElements, Cij, particlesAll, particlesOutOfEq, particleMasses];
  
  (* Create replacement list and shift indices *)
  MatrixElementsList = Table[MatrixElemToC@i //. OptionValue[Replacements], {i, MatrixElements}];
  MatrixElementsList = DeleteCases[MatrixElementsList, a_ -> 0];
  
  <|
    "Cij" -> Cij,
    "MatrixElementsList" -> MatrixElementsList,
    "particleListAll" -> particleListAll,
    "particleNames" -> particleNames,
    "userCouplings" -> userCouplings,
    "particleMasses" -> particleMasses
  |>
];


ExportMatrixElements[file_: None, outOfEqParticleList_, lightParticleList_, OptionsPattern[]] :=
Block[
  {
    data,
    userFormat, userFormatsList,
    FormatOptions,
    privFile = file,
    commandLineArgs, userParameters
  },
  
  (* Restore normal output *)
  $Output = $WallGoMatrixOutput;
  commandLineArgs = $ScriptCommandLine;
  
  (* Allow overriding file name from command line *)
  Do[
    If[commandLineArgs[[i]] == "--outputFile", privFile = commandLineArgs[[i + 1]]],
    {i, 2, Length[commandLineArgs] - 1}
  ];
  
  (* Always generate matrix elements *)
  data = GenerateMatrixElementsData[
    outOfEqParticleList, lightParticleList,
    NormalizeWithDOF   -> OptionValue[NormalizeWithDOF],
    TruncateAtLeadingLog -> OptionValue[TruncateAtLeadingLog],
    TagLeadingLog      -> OptionValue[TagLeadingLog],
    Verbose            -> OptionValue[Verbose],
    Replacements       -> OptionValue[Replacements]
  ];
  
  (* Export only if file name provided *)
  If[privFile =!= None && StringQ[privFile] && StringLength[privFile] > 0,
    
    FormatOptions = {"txt", "json", "hdf5", "all", "none"};
    userFormat = OptionValue[Format];
    userFormatsList = If[ListQ[userFormat], userFormat, {userFormat}];
    If[MemberQ[userFormatsList, "all"], userFormatsList = {"txt", "json", "hdf5"}];
    
    If[
      AllTrue[userFormatsList, MemberQ[FormatOptions, #] &],
      Do[
        Switch[fmt,
          "txt",
            ExportTo["txt"][
              data["MatrixElementsList"],
              data["particleListAll"],
              data["userCouplings"],
              privFile],
          "hdf5",
            ExportTo["hdf5"][
              data["Cij"],
              outOfEqParticleList,
              data["particleNames"],
              data["userCouplings"],
              privFile],
          "json",
            userParameters = Flatten[Join[data["userCouplings"], data["particleMasses"]]] // DeleteDuplicates;
            ExportTo["json"][
              data["MatrixElementsList"],
              data["particleNames"],
              userParameters,
              privFile]
        ],
        {fmt, userFormatsList}
      ],
      Message[WallGoMatrix::failWarning, 
        "Allowed export formats are: txt, hdf5, and json."]
    ];
  ,
    If[OptionValue[Verbose],
      Print["No file name provided \[LongDash] skipping export, returning generated data only."]
    ];
  ];
  
  Return[PrintNonPrivate[data["MatrixElementsList"]]]
];


MatrixElemToC[MatrixElem_]:=Block[{Ind},
	Ind=MatrixElem[[2]];
	
	Return[M[Ind[[1]]-1,Ind[[2]]-1,Ind[[3]]-1,Ind[[4]]-1]->MatrixElem[[1]]]
]


MatrixElemFromC[MatrixElem_]:=Block[{Ind},
	Ind=MatrixElem[[2]];
	
	Return[M[Ind[[1]]+1,Ind[[2]]+1,Ind[[3]]+1,Ind[[4]]+1]->MatrixElem[[1]]]
]


(* ::Section:: *)
(*MatrixElements: Importing the results*)


ImportMatrixElements[file_]:=Module[
{
	jsonObject,particles,parameters,results
},
	(* Only implemented for JSON *)
	If[StringTake[file,-5]!=".json",Print["File must end in .json"];Return[]];

	(* Importing into JSON object *)
	jsonObject=importJSONMatrixElements[file];

	(* Splitting JSON object *)
	{particles,parameters,results}=splitJsonMatrixElements[jsonObject];

	(* Returning results *)
	Return[{particles,parameters,results}]
];


(* ::Section:: *)
(*Private constants*)


{CT}


{\[CapitalLambda]\[Lambda],\[CapitalLambda]g,Hg,Habij,HabIJF,HabIJFC,Ysij,YsijC,YTemp,YTempC,Yhelp,YhelpC};(*Private Variables*)


{\[Lambda]1IP};(*Private Variables*)


{\[Mu]IJF,\[Mu]IJFC,GroupMathCleared}


End[]
EndPackage[]
