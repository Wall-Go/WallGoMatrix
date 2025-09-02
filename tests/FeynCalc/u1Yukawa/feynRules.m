(* ::Package:: *)

Quit[];


(* ::Section:: *)
(*FeynRules evaluation*)


$FeynRulesPath=DirectoryName[FindFile["FeynRules`"]];
<<FeynRules`


(* load model file *)
model="u1-higgs-yukawa";
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
RunProcess[{"bash", "../replaceFeynCalc.sh", "u1-higgs-yukawa.mod"}]
RunProcess[{"bash", "../replaceFeynCalc.sh", "u1-higgs-yukawa.gen"}]



