(* ::Package:: *)

(* Builds a Paclet *)


(* Define the temporary directory *)
tempDir = CreateDirectory[FileNameJoin[{NotebookDirectory[], "build", "TempPaclet"}]];

(* Copy Kernel directory and PacletInfo.m to the temporary directory *)
CopyDirectory[FileNameJoin[{NotebookDirectory[], "Kernel"}], FileNameJoin[{tempDir, "Kernel"}]];
CopyFile[FileNameJoin[{NotebookDirectory[], "PacletInfo.m"}], FileNameJoin[{tempDir, "PacletInfo.m"}]];

(* Set the temporary directory as the paclet directory *)
pacletFilename=CreatePacletArchive[tempDir];
DeleteDirectory[tempDir, DeleteContents -> True];
Print["Temporary directory deleted."];
