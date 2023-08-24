#! /usr/bin/env bash

DATADIR=./dijet-HERA-plot

gnuplot -p -persist <<EOF
set terminal epslatex standalone size 2.75in,2.5in  level3
set logscale xy
set xrange [5:50]
set yrange [1e-4:1e+2]
set format y "%T"


set rmargin 0.5
set tmargin 0.5
set lmargin 2.5
set bmargin 1.5

unset label
set object rectangle from 5.1,1.1e-4 to 15,5e-4 fillcolor rgbcolor '#FFFFFF' fs solid noborder front
set output '${DATADIR}/plots/${5}figure16-1.tex'
set label "\\\\small \$5.5<Q^2<8\\\\;\\\\mathrm{GeV}^2\$" at 5.5,2e-4 front
plot '${DATADIR}/${1}/fig161.hst' using 1:4:3 with filledcurves lt rgb "blue" fs transparent pattern 5 t '${3}', \
	 '${DATADIR}/${1}/fig161.hst' using 1:3 with lines lt rgb "blue" t '', \
	 '${DATADIR}/${1}/fig161.hst' using 1:4 with lines lt rgb "blue" t '', \
	 '${DATADIR}/${2}/fig161.hst' using 1:4:3 with filledcurves lt rgb "red" fs transparent pattern 4 t '${4}', \
	 '${DATADIR}/${2}/fig161.hst' using 1:3 with lines lt rgb "red" t '', \
	 '${DATADIR}/${2}/fig161.hst' using 1:4 with lines lt rgb "red" t '', \
	 '${DATADIR}/${1}/fig161.hst' using 1:2 with lines lt rgb "blue" lw 1.5 t '', \
	 '${DATADIR}/${2}/fig161.hst' using 1:2 with lines lt rgb "red" lw 1.5 t '', \
     'dijet-HERA/subtable-1.txt' using (\$1+\$2)/2:6:(\$6*(1+0.01*\$4)):(\$6*(1+0.01*\$5)) with yerrorbars t 'Data' \
		lt rgb 'green' 
	 
unset label
set output '${DATADIR}/plots/${5}figure16-2.tex'
set label "\\\\small\$8<Q^2<11\\\\;\\\\mathrm{GeV}^2\$" at 5.5,2e-4 front
plot '${DATADIR}/${1}/fig162.hst' using 1:4:3 with filledcurves lt rgb "blue" fs transparent pattern 5 t '${3}', \
	 '${DATADIR}/${1}/fig162.hst' using 1:3 with lines lt rgb "blue" t '', \
	 '${DATADIR}/${1}/fig162.hst' using 1:4 with lines lt rgb "blue" t '', \
	 '${DATADIR}/${2}/fig162.hst' using 1:4:3 with filledcurves lt rgb "red" fs transparent pattern 4 t  '${4}', \
	 '${DATADIR}/${2}/fig162.hst' using 1:3 with lines lt rgb "red" t '', \
	 '${DATADIR}/${2}/fig162.hst' using 1:4 with lines lt rgb "red" t '', \
	 '${DATADIR}/${1}/fig162.hst' using 1:2 with lines lt rgb "blue" lw 1.5 t '', \
	 '${DATADIR}/${2}/fig162.hst' using 1:2 with lines lt rgb "red" lw 1.5 t '', \
     'dijet-HERA/subtable-2.txt' using (\$1+\$2)/2:6:(\$6*(1+0.01*\$4)):(\$6*(1+0.01*\$5)) with yerrorbars t 'Data' \
		lt rgb 'green' 
	 

unset label
set output '${DATADIR}/plots/${5}figure16-3.tex'
set label "\\\\small\$11<Q^2<16\\\\;\\\\mathrm{GeV}^2\$" at 5.5,2e-4 front
plot '${DATADIR}/${1}/fig163.hst' using 1:4:3 with filledcurves lt rgb "blue" fs transparent pattern 5 t '${3}', \
	 '${DATADIR}/${1}/fig163.hst' using 1:3 with lines lt rgb "blue" t '', \
	 '${DATADIR}/${1}/fig163.hst' using 1:4 with lines lt rgb "blue" t '', \
	 '${DATADIR}/${2}/fig163.hst' using 1:4:3 with filledcurves lt rgb "red" fs transparent pattern 4 t '${4}', \
	 '${DATADIR}/${2}/fig163.hst' using 1:3 with lines lt rgb "red" t '', \
	 '${DATADIR}/${2}/fig163.hst' using 1:4 with lines lt rgb "red" t '', \
	 '${DATADIR}/${1}/fig163.hst' using 1:2 with lines lt rgb "blue" lw 1.5 t '', \
	 '${DATADIR}/${2}/fig163.hst' using 1:2 with lines lt rgb "red" lw 1.5 t '', \
     'dijet-HERA/subtable-3.txt' using (\$1+\$2)/2:6:(\$6*(1+0.01*\$4)):(\$6*(1+0.01*\$5)) with yerrorbars t 'Data' \
		lt rgb 'green' 
	 
unset label
set output '${DATADIR}/plots/${5}figure16-4.tex'
set label "\\\\small\$16<Q^2<22\\\\;\\\\mathrm{GeV}^2\$" at 5.5,2e-4 front
plot '${DATADIR}/${1}/fig164.hst' using 1:4:3 with filledcurves lt rgb "blue" fs transparent pattern 5 t '${3}', \
	 '${DATADIR}/${1}/fig164.hst' using 1:3 with lines lt rgb "blue" t '', \
	 '${DATADIR}/${1}/fig164.hst' using 1:4 with lines lt rgb "blue" t '', \
	 '${DATADIR}/${2}/fig164.hst' using 1:4:3 with filledcurves lt rgb "red" fs transparent pattern 4 t  '${4}', \
	 '${DATADIR}/${2}/fig164.hst' using 1:3 with lines lt rgb "red" t '', \
	 '${DATADIR}/${2}/fig164.hst' using 1:4 with lines lt rgb "red" t '', \
	 '${DATADIR}/${1}/fig164.hst' using 1:2 with lines lt rgb "blue" lw 1.5 t '', \
	 '${DATADIR}/${2}/fig164.hst' using 1:2 with lines lt rgb "red" lw 1.5 t '', \
     'dijet-HERA/subtable-4.txt' using (\$1+\$2)/2:6:(\$6*(1+0.01*\$4)):(\$6*(1+0.01*\$5)) with yerrorbars t 'Data' \
		lt rgb 'green' 

unset label
set output '${DATADIR}/plots/${5}figure16-5.tex'
set label "\\\\small\$22<Q^2<30\\\\;\\\\mathrm{GeV}^2\$" at 5.5,2e-4 front
plot '${DATADIR}/${1}/fig165.hst' using 1:4:3 with filledcurves lt rgb "blue" fs transparent pattern 5 t '${3}', \
	 '${DATADIR}/${1}/fig165.hst' using 1:3 with lines lt rgb "blue" t '', \
	 '${DATADIR}/${1}/fig165.hst' using 1:4 with lines lt rgb "blue" t '', \
	 '${DATADIR}/${2}/fig165.hst' using 1:4:3 with filledcurves lt rgb "red" fs transparent pattern 4 t '${4}', \
	 '${DATADIR}/${2}/fig165.hst' using 1:3 with lines lt rgb "red" t '', \
	 '${DATADIR}/${2}/fig165.hst' using 1:4 with lines lt rgb "red" t '', \
	 '${DATADIR}/${1}/fig165.hst' using 1:2 with lines lt rgb "blue" lw 1.5 t '', \
	 '${DATADIR}/${2}/fig165.hst' using 1:2 with lines lt rgb "red" lw 1.5 t '', \
     'dijet-HERA/subtable-5.txt' using (\$1+\$2)/2:6:(\$6*(1+0.01*\$4)):(\$6*(1+0.01*\$5)) with yerrorbars t 'Data' \
		lt rgb 'green' 

unset label
set output '${DATADIR}/plots/${5}figure16-6.tex'
set label "\\\\small\$30<Q^2<42\\\\;\\\\mathrm{GeV}^2\$" at 5.5,2e-4 front
plot '${DATADIR}/${1}/fig166.hst' using 1:4:3 with filledcurves lt rgb "blue" fs transparent pattern 5 t '${3}', \
	 '${DATADIR}/${1}/fig166.hst' using 1:3 with lines lt rgb "blue" t '', \
	 '${DATADIR}/${1}/fig166.hst' using 1:4 with lines lt rgb "blue" t '', \
	 '${DATADIR}/${2}/fig166.hst' using 1:4:3 with filledcurves lt rgb "red" fs transparent pattern 4 t  '${4}', \
	 '${DATADIR}/${2}/fig166.hst' using 1:3 with lines lt rgb "red" t '', \
	 '${DATADIR}/${2}/fig166.hst' using 1:4 with lines lt rgb "red" t '', \
	 '${DATADIR}/${1}/fig166.hst' using 1:2 with lines lt rgb "blue" lw 1.5 t '', \
	 '${DATADIR}/${2}/fig166.hst' using 1:2 with lines lt rgb "red" lw 1.5 t '', \
     'dijet-HERA/subtable-6.txt' using (\$1+\$2)/2:6:(\$6*(1+0.01*\$4)):(\$6*(1+0.01*\$5)) with yerrorbars t 'Data' \
		lt rgb 'green' 


unset label
set output '${DATADIR}/plots/${5}figure16-7.tex'
set label "\\\\small\$42<Q^2<60\\\\;\\\\mathrm{GeV}^2\$" at 5.5,2e-4 front
plot '${DATADIR}/${1}/fig167.hst' using 1:4:3 with filledcurves lt rgb "blue" fs transparent pattern 5 t '${3}', \
	 '${DATADIR}/${1}/fig167.hst' using 1:3 with lines lt rgb "blue" t '', \
	 '${DATADIR}/${1}/fig167.hst' using 1:4 with lines lt rgb "blue" t '', \
	 '${DATADIR}/${2}/fig167.hst' using 1:4:3 with filledcurves lt rgb "red" fs transparent pattern 4 t '${4}', \
	 '${DATADIR}/${2}/fig167.hst' using 1:3 with lines lt rgb "red" t '', \
	 '${DATADIR}/${2}/fig167.hst' using 1:4 with lines lt rgb "red" t '', \
	 '${DATADIR}/${1}/fig167.hst' using 1:2 with lines lt rgb "blue" lw 1.5 t '', \
	 '${DATADIR}/${2}/fig167.hst' using 1:2 with lines lt rgb "red" lw 1.5 t '', \
     'dijet-HERA/subtable-7.txt' using (\$1+\$2)/2:6:(\$6*(1+0.01*\$4)):(\$6*(1+0.01*\$5)) with yerrorbars t 'Data' \
		lt rgb 'green' 
	 
unset label
set output '${DATADIR}/plots/${5}figure16-8.tex'
set label "\\\\small\$60<Q^2<80\\\\;\\\\mathrm{GeV}^2\$" at 5.5,2e-4 front
plot '${DATADIR}/${1}/fig168.hst' using 1:4:3 with filledcurves lt rgb "blue" fs transparent pattern 5 t '${3}', \
	 '${DATADIR}/${1}/fig168.hst' using 1:3 with lines lt rgb "blue" t '', \
	 '${DATADIR}/${1}/fig168.hst' using 1:4 with lines lt rgb "blue" t '', \
	 '${DATADIR}/${2}/fig168.hst' using 1:4:3 with filledcurves lt rgb "red" fs transparent pattern 4 t '${4}', \
	 '${DATADIR}/${2}/fig168.hst' using 1:3 with lines lt rgb "red" t '', \
	 '${DATADIR}/${2}/fig168.hst' using 1:4 with lines lt rgb "red" t '', \
	 '${DATADIR}/${1}/fig168.hst' using 1:2 with lines lt rgb "blue" lw 1.5 t '', \
	 '${DATADIR}/${2}/fig168.hst' using 1:2 with lines lt rgb "red" lw 1.5 t '', \
     'dijet-HERA/subtable-8.txt' using (\$1+\$2)/2:6:(\$6*(1+0.01*\$4)):(\$6*(1+0.01*\$5)) with yerrorbars t 'Data' \
		lt rgb 'green' 
EOF

for i in ./dijet-HERA-plot/plots/*.tex; do pdflatex ${i} --output-directory=./dijet-HERA-plot; done

