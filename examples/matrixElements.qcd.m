(* ::Package:: *)

Quit[];


SetDirectory[NotebookDirectory[]];
(*Put this if you want to create multiple model-files with the same kernel*)
$GroupMathMultipleModels=True;
$LoadGroupMath=True;
<<../DRalgo/DRalgo.m
<<../src/matrixElements.m


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


(* ::Title:: *)
(*A model with 6 quarks and 1 gluon*)


(* ::Subtitle:: *)
(*UserInput*)


(*
In DRalgo fermions are Weyl.
So to create one Dirac we need
one left-handed and
one right-handed fermoon
*)


(*Below
rep 1-6 are quarks,
rep 7 is a gluon
*)
Rep1=CreateOutOfEq[{1,2},"F"];
Rep2=CreateOutOfEq[{3,4},"F"];
Rep3=CreateOutOfEq[{5,6},"F"];
Rep4=CreateOutOfEq[{7,8},"F"];
Rep5=CreateOutOfEq[{9,10},"F"];
Rep6=CreateOutOfEq[{11,12},"F"];
RepGluon=CreateOutOfEq[{1},"V"];
(*check*)
(*Rep2=CreateOutOfEq[{3,4,...,12},"F"];*)


ParticleList={Rep1,RepGluon,Rep2,Rep3,Rep4,Rep5,Rep6};
(*
These particles do not have out-of-eq contributions
*)
LightParticles={3,4,5,6,7};


(*Defining various masses and couplings*)


VectorMass=Table[mg2,{i,1,Length[gvff]}];
FermionMass=Table[mq2,{i,1,Length[gvff[[1]]]}];
(*
up to the user to make sure that the same order is given in the python code
*)
UserMasses={mq2,mg2};
UserCouplings={gs};


SetDirectory[NotebookDirectory[]];
ParticleName={"Top","Gluon"};
MatrixElements=ExportMatrixElements["matrixElements.qcd",ParticleList,LightParticles,UserMasses,UserCouplings,ParticleName];


MatrixElements


(*comparison with https://arxiv.org/pdf/hep-ph/0302165.pdf*)
(*q1q2->q1q2*)
(*5 light quarks*)
M[0,2,0,2]/.MatrixElements/.{c[0]->1}(*/.{-u->-t}/.{t*u->-s*t}*)
5*2/(2*CA)*(8*dF^2*CF^2/dA((s^2+u^2)/ttsq))/.{ttsq->(-t+msq[1])^2}/.{Nf->3,CA->3,CF->4/3,dF->3,dA->8}
%-%%//Simplify
(*q1q1->gg*)
M[0,0,1,1]/.MatrixElements/.{c[0]->1}(*/.{-u->-t}/.{t*u->-s*t}*)
1/(2*2*CA)*(8*dF*CF^2(u/tt+t/uu)-8*dF*CF*CA((t^2+u^2)/s^2))/.{tt->(-t+msq[0])^2/t,uu->(-u+msq[0])^2/u}/.{Nf->3,CA->3,CF->4/3,dF->3}
%-%%//Simplify
(*q1 g->q1 g*)
M[0,1,0,1]/.MatrixElements/.{c[0]->1}(*/.{-u->-t}/.{t*u->-s*t}*)//Simplify
1/(2*CA)*(-8*dF*CF^2(u/s +s/uu)+8*dF*CF*CA((s^2+u^2)/ttsq))/.{ttsq->(-t+msq[1])^2,tt->(-t+msq[1])^2/t,uu->(-u+msq[0])^2/u}/.{Nf->3,CA->3,CF->4/3,dF->3}//Simplify
%-%%//Simplify


(*crosscheck*)
(*tt->tt*)
M[0,0,0,0]/.MatrixElements/.{c[0]->1}(*/.{-u->-t}/.{t*u->-s*t}*)
(*result from AMY paper*)
(*initial sates and finial states are either particle or anti-particle*)
+8*dF^2*CF^2/dA((s^2+u^2)/ttsq+(s^2+t^2)/uusq)+16*dF*CF(CF-CA/2)(s^2/(t*u));
(*mixed initial and finial states*)
+8*dF^2*CF^2/dA((s^2+u^2)/ttsq+(u^2+t^2)/s^2)+16*dF*CF(CF-CA/2)(u^2/(t*s));
1/(2*CA)*(xx*%%+%)/.{ttsq->(-t+msq[1])^2,uusq->(-u+msq[1])^2}/.{Nf->3,CA->3,CF->4/3,dF->3,dA->8}
bb %-%%%%//FullSimplify
(*result only agrees by assuming factor 1/2 for same particles in final state*)
%/.{xx->xx,bb->xx}/.{xx->1}


(*g g->g g*)
M[1,1,1,1]/.MatrixElements/.{c[0]->1}//Simplify
(*Final states are the same -> divide out 1/2 symmetry factor*)
1/(2*2(CA^2-1))(16*dA*CA^2(3-s*u/ttsq-s*t/uusq-t*u/s^2));
%/.{s*t/uusq->1/4-1/4(s-t)^2/(-u+msq[1])^2,s*u/ttsq->1/4-1/4(s-u)^2/(-t+msq[1])^2}/.{Nf->3,dA->8,CA->3,CF->4/3,dF->3}//Simplify
%%%-%//Simplify


Import["MatrixElem.hdf5"]


Import["MatrixElem.hdf5","MatrixElementsTopTop"]


Import["MatrixElem.hdf5","CouplingInfo"]


Import["MatrixElem.hdf5","ParticleInfo"]
