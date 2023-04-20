# Saturation-Model

it is a folder for the project on the GBW/ BGK models.
Written by Tomoki Goda
Most functions are in "Function".

To use one needs to write "control.h", in which control macros are defined. 
default values of them can be found in the file "control-default.h"

Makefile requires environmental variable DIR set to where the "control.h" is.

Examples can be found in Run-kt.

Generally, the following are necessary

CUBA (T. Hahn)
Gnu Scientific Library (GSL)
ROOT::Minuit2

KaTie (A.van Hameren,2016) is required for the contents of "dijet", and gluplot if one wants to run 

$ gnuplot plotting.txt.

Grid-Interpolation is just a copy of KS gluon.
