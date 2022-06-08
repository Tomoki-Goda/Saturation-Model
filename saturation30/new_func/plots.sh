#! /usr/bin/env bash


gcc Auto-Plot-DP.c -o Auto-Plot -lm
#gcc makelist.c -o makelist
gcc Auto-Plot-f2.c -o Auto-F2 -lm

dir="NewRun"
##########################################################


#./makelist 7 ./Archive0106/GBS/*massless*/result.txt ./Archive0106/GBS/massless.txt
#./makelist 7 ./Archive0106/GBS/*massive*/result.txt ./Archive0106/GBS/massive.txt


#./makelist 7 ./Archive0106/GBSrfix/*massless*/result.txt ./Archive0106/GBSrfix/massless.txt
#./makelist 7 ./Archive0106/GBSrfix/*massive*/result.txt ./Archive0106/GBSrfix/massive.txt


#./makelist 5 ./Archive0106/GBSPert/*massless*/result.txt ./Archive0106/GBSPert/massless.txt
#./makelist 5 ./Archive0106/GBSPert/*massive*/result.txt ./Archive0106/GBSPert/massive.txt


#./makelist 5 ./Archive0106/GBSPertrfix/*massless*/result.txt ./Archive0106/GBSPertrfix/massless.txt
#./makelist 5 ./Archive0106/GBSPertrfix/*massive*/result.txt ./Archive0106/GBSPertrfix/massive.txt

######################################################################
./Auto-Plot 100 5 ./${dir}/*/Mass*

./Auto-Plot 100 3 ./${dir}/*/Mass*

./Auto-Plot 500 5 ./${dir}/*/Mass*

./Auto-Plot 500 3 ./${dir}/*/Mass*

#####################################################################
./Auto-F2 5 ./${dir}/*/Mass*

./Auto-F2 3 ./${dir}/*/Mass*

#./Auto-F2 5 ./${dir}/*/*LCB

#./Auto-F2 3 ./${dir}/*/*LCB
