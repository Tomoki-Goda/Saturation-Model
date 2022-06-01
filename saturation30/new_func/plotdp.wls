#!/usr/bin/env wolframscript
(* ::Package:: *)

points= Import["./dp-plots.m"];

pointsl=Map[Log[10,#]&,points,{2}];

plot=ListPlot[pointsl,PlotStyle->PointSize[Small], ImageSize->Full];

Export["PlotDP.jpg",plot];
