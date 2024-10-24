(* ::Package:: *)

(* Builds a Paclet *)
pacletDirectory=FileNameJoin[{NotebookDirectory[],"src/"}];
pacletFilename=CreatePacletArchive[pacletDirectory];
