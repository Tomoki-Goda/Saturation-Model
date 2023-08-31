# Saturation-Model
<p> Saturation-Model version-3 </p>
<p>
it is a folder for the project on the GBW/ BGK models.
Written by Tomoki Goda. 
Most functions are in "saturation-ver3/Function".
</p>

<p>
To use one needs to write "control.h", in which control macros are defined. 
default values of them can be found in the file "control-default.h"
</p>
<p>
Makefile requires environmental variable DIR set to where the "control.h" is.
</p>
<p>
Examples can be found in Run-kt.
</p>
<p>
Generally, the following are necessary<br>
CUBA (T. Hahn)<br>
Gnu Scientific Library (GSL)<br>
ROOT::Minuit2<br>
</p>

<p>
KaTie (A.van Hameren,2016) is required for the contents of "dijet", and gluplot if one wants to run <br>
<code>$ gnuplot plotting.txt</code>.
</p>
<p>
Grid-Interpolation is just a copy of KS gluon.
</p>

Folders may contain INSTRUCTION files to tell you how to use them.
