#!/usr/bin/env wolframscript
(* ::Package:: *)

points= Import["./f2-plots.m"];

pointsl=Map[Log[10,#]&,points,{2}];

plot=ListPlot[pointsl,PlotStyle->PointSize[Small], ImageSize->Full];

Export["PlotF2.jpg",plot];
