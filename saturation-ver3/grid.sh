#! /usr/bin/env bash

export DIR

WW=/home/tomoki/Saturation-Model/saturation-ver3/RunWW
PREV=/home/tomoki/Saturation-Model/saturation-ver3/Run2808
RUN3=/home/tomoki/Saturation-Model/saturation-ver3/Run3

DIR=${WW}/GBW
make ww -silent
${DIR}/ww -i ${RUN3}/GBW/result.txt -o ${DIR}/GBWr.dat

DIR=${WW}/BGK
make ww -silent
${DIR}/ww -i ${RUN3}/BGK/result.txt -o ${DIR}/BGKr.dat

DIR=${WW}/GBWS
make ww -silent
${DIR}/ww -i ${RUN3}/GBW/result.txt -o ${DIR}/GBWSr.dat

DIR=${WW}/BGKS
make ww -silent
${DIR}/ww -i ${RUN3}/BGK/result.txt -o ${DIR}/BGKSr.dat

###########################################


DIR=${WW}/fixa-bjorx
make ww -silent
${DIR}/ww -i ${RUN3}/fixa-bjorx/result.txt -o ${DIR}/GBWkt.dat

DIR=${WW}/fixa-bjorx-BGK
make ww -silent
${DIR}/ww -i ${RUN3}/fixa-bjorx-BGK/result.txt -o ${DIR}/BGKkt.dat

DIR=${WW}/fixa-bjorx-Sud
make ww -silent
${DIR}/ww -i ${RUN3}/fixa-bjorx/result.txt -o ${DIR}/GBWSkt.dat

DIR=${WW}/fixa-bjorx-BGK-Sud
make ww -silent
${DIR}/ww -i ${RUN3}/fixa-bjorx-BGK/result.txt -o ${DIR}/BGKSkt.dat

