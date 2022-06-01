#! /usr/bin/env bash

gcc Automation_Run.c -o Auto_Run
gcc Auto-Plot-DP.c -o Auto-Plot
##########################################################
./Auto_Run ./GBS/*LCB 
./makelist 7 ./GBS/*massless*/result.txt ./GBS/massless.txt
./makelist 7 ./GBS/*massive*/result.txt ./GBS/massive.txt

./Auto_Run ./GBSrfix/*LCB 
./makelist 7 ./GBSrfix/*massless*/result.txt ./GBSrfix/massless.txt
./makelist 7 ./GBSrfix/*massive*/result.txt ./GBSrfix/massive.txt

./Auto_Run ./GBSPert/*LCB 
./makelist 7 ./GBSPert/*massless*/result.txt ./GBSPert/massless.txt
./makelist 7 ./GBSPert/*massive*/result.txt ./GBSPert/massive.txt

./Auto_Run ./GBSPertrfix/*LCB 
./makelist 7 ./GBSPert/*massless*/result.txt ./GBSPert/massless.txt
./makelist 7 ./GBSPert/*massive*/result.txt ./GBSPert/massive.txt

######################################################################
./Auto-Plot ./GBS/*LCB
./Auto-Plot ./GBSrfix/*LCB
./Auto-Plot ./GBSPert/*LCB
./Auto-Plot ./GBS/*LCB

