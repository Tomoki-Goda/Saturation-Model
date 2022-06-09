#! /usr/bin/env bash


gcc Auto-Plot-DP.c -o Auto-Plot -lm
#gcc makelist.c -o makelist
gcc Auto-Plot-f2.c -o Auto-F2 -lm

dir="NewRun"
##########################################################

#./makelist 7 ./Archive0106/GBS/*massless*/result.txt ./Archive0106/GBS/massless.txt

######################################################################
./Auto-Plot 50 5 ./${dir}/*/Mass*

./Auto-Plot 50 3 ./${dir}/*/Mass*

#./Auto-Plot 500 5 ./${dir}/*/Mass*

#./Auto-Plot 500 3 ./${dir}/*/Mass*

#####################################################################
#./Auto-F2 5 ./${dir}/*/Mass*

#./Auto-F2 3 ./${dir}/*/Mass*


