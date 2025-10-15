(* ::Package:: *)

Quit[];


(* ::Section:: *)
(*FeynRules evaluation*)


$FeynRulesPath=DirectoryName[FindFile["FeynRules`"]];
<<FeynRules`


(* load model file *)
model="n-scalars";
modelFile =FileNameJoin[{ NotebookDirectory[],model<>".fr"}];
LoadModel[modelFile]


L


(* choosing output location *)
outputFiles=FileNameJoin[NotebookDirectory[]<>model];
(* writing FeynArts output *)
WriteFeynArtsOutput[L,Output->outputFiles,CouplingRename->False]
(* deleting empty directory made by FeynRules *)
DeleteDirectory[outputFiles]


SetDirectory[NotebookDirectory[]]
Directory[]
RunProcess[{"bash", "../replaceFeynCalc.sh", "n-scalars.mod"}]
RunProcess[{"bash", "../replaceFeynCalc.sh", "n-scalars.gen"}]
