#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(**)


s[a_]:=1-Exp[-a^2]


k=Sqrt[100];
nmax=30;
shift=(*3*) Pi/(4k);
a=nmax*2 Pi/k+shift; (* this is the bounday between naive integration and approximation*)


(* ::Text:: *)
(*Try to compute this numerically*)


FourierTransform[s[r]/r^2/.r->Sqrt[r1^2+r2^2],{r1,r2},{k1,k2}]/.k1->Sqrt[k^2-k2^2]//FullSimplify
analytic=%//N//ScientificForm[#,10]&


(* ::Input:: *)
(**)


(* ::Section:: *)
(**)


(* ::Text:: *)
(*Here the integration is done without any enhancement*)


Install["/home/tomoki/Saturation-Model/saturation-ver3/math/levin"]


Plot[s[r]BesselJ[0,r k]/r,{r,0, a}]
NIntegrate[s[r]BesselJ[0,r k]/r,{r,0, a}]


table={};
min=0;
max=shift;
LevinReset[];
Do[max+=Pi/k;
(*Print[Plot[s[r]BesselJ[0,r k]/r,{r,min,max}]];*)
val= NIntegrate[s[r]BesselJ[0,r k]/r,{r,min,max}];
LevinAdd[val];
AppendTo[table ,{i,val}];
min=max;
,{i,0,nmax-1}]
sum1=Map[Plus@@(table[[1;;(#[[1]]+1),2]])&,table]



sum1[[nmax]]
LevinAccel[nmax-1,0]
LevinAccel[nmax-3,2]
LevinAccel[nmax-6,5]
LevinAccel[nmax-8,7]
analytic


nmax
LevinAccel[nmax-6-14,5]
LevinSum[nmax-14]


(* ::Text:: *)
(*Now, after 15 sectors, Levin 5 point acceleration produce correct result!!!*)


(* ::Text:: *)
(*Note if integral was done in 2Pi/k, somehow result is bad...*)
(*Also, shift=3Pi/4 is fine.*)


ltable=Table[{i-1,LevinAccel[nmax-6-(nmax-i),5]},{i,6,nmax}];
table;
ListPlot[{%%,%},PlotLegends->{"Levin","Direct" }]


ltable=Table[{i-1,LevinAccel[nmax-7-(nmax-i),6]//N[Round[#,10^(Round[Log10[Abs[#]],1]-4)]]&(*//ScientificForm[#,4]&*),table[[i,2]]//N[Round[#,10^(Round[Log10[Abs[#]],1]-4)]]&(*//ScientificForm[#,4]&*)},{i,10,30,3}]


Export["/home/tomoki/Saturation-Model/saturation-ver3/levin.dat", ltable, "Table"]

