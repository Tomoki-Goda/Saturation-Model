#! /usr/bin/env bash

gcc Automation_Run.c -o Auto_Run
gcc Auto-Plot-DP.c -o Auto-Plot -lm
gcc makelist.c -o makelist
##########################################################
./Auto_Run ./GBWIntegrated/*LCB
./makelist 3 ./GBWIntegrated/*massless*/result.txt ./GBWIntegrated/massless.txt
./makelist 3 ./GBWIntegrated/*massive*/result.txt ./GBWIntegrated/massive.txt
./Auto-Plot 500 3 ./GBWIntegrated/*LCB
./Auto-Plot 50 3 ./GBWIntegrated/*LCB
./Auto-Plot 500 5 ./GBWIntegrated/*LCB
./Auto-Plot 50 5 ./GBWIntegrated/*LCB

./Auto_Run ./GBSPert/*LCB 
./makelist 5 ./GBSPert/*massless*/result.txt ./GBSPert/massless.txt
./makelist 5 ./GBSPert/*massive*/result.txt ./GBSPert/massive.txt
./Auto-Plot 500 3 ./GBSPert/*LCB
./Auto-Plot 50 3 ./GBSPert/*LCB
./Auto-Plot 500 5 ./GBSPert/*LCB
./Auto-Plot 50 5 ./GBSPert/*LCB


./Auto_Run ./GBSPertrfix/*LCB 
./makelist 5 ./GBSPertrfix/*massless*/result.txt ./GBSPertrfix/massless.txt
./makelist 5 ./GBSPertrfix/*massive*/result.txt ./GBSPertrfix/massive.txt
./Auto-Plot 500 3 ./GBSPertrfix/*LCB
./Auto-Plot 50 3 ./GBSPertrfix/*LCB
./Auto-Plot 500 5 ./GBSPertrfix/*LCB
./Auto-Plot 50 5 ./GBSPertrfix/*LCB

./Auto_Run ./GBS/*LCB 
./makelist 7 ./GBS/*massless*/result.txt ./GBS/massless.txt
./makelist 7 ./GBS/*massive*/result.txt ./GBS/massive.txt
./Auto-Plot 500 3 ./GBS/*LCB
./Auto-Plot 50 3 ./GBS/*LCB
./Auto-Plot 500 5 ./GBS/*LCB
./Auto-Plot 50 5 ./GBS/*LCB


./Auto_Run ./GBSrfix/*LCB 
./makelist 7 ./GBSrfix/*massless*/result.txt ./GBSrfix/massless.txt
./makelist 7 ./GBSrfix/*massive*/result.txt ./GBSrfix/massive.txt
./Auto-Plot 500 3 ./GBSrfix/*LCB
./Auto-Plot 50 3 ./GBSrfix/*LCB
./Auto-Plot 500 5 ./GBSrfix/*LCB
./Auto-Plot 50 5 ./GBSrfix/*LCB





######################################################################


