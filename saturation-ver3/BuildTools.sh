#! /usr/bin/env bash

DIR="./Utilities"
#DIR=.

gcc ${DIR}/Auto_control.c -o Auto-Control 
gcc ${DIR}/Append-control.c -o Append


gcc ${DIR}/Write_TeX.c -o Write_TeX
gcc ${DIR}/make-table.c -o make-table
