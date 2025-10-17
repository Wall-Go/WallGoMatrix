(* ::Package:: *)

$FeynRulesPath=DirectoryName[FindFile["FeynRules`"]];
<<FeynRules`


(* load model file *)
model="yukawa";
modelFile =FileNameJoin[{ NotebookDirectory[],model<>".fr"}];
LoadModel[modelFile]


(* choosing output location *)
outputFiles=FileNameJoin[NotebookDirectory[]<>model];
(* writing FeynArts output *)
WriteFeynArtsOutput[L,Output->outputFiles,CouplingRename->False]
(* deleting empty directory made by FeynRules *)
DeleteDirectory[outputFiles]



