#!/usr/bin/env bash


#DIR="/home/tomoki/Saturation-Model/saturation30/new_func/Models"

DIR="./Models"

tar -czvf "${DIR}-previous.tar.gz" ${DIR}
rm -r ${DIR}/*G*
################# GBW  #################################

NDIR="${DIR}/GBW"
mkdir ${NDIR}

mkdir ${NDIR}/Massless
${DIR}/../Auto-Control -dir ${NDIR}/Massless -lmass 0.0 -qup 50 500 -sudakov 0 -rfix 0 -model 0

mkdir ${NDIR}/Massive
${DIR}/../Auto-Control -dir ${NDIR}/Massive -lmass 0.0196 -qup 50 500 -sudakov 0 -rfix 0 -model 0

#mkdir ${NDIR}/Masslessrfix
#${DIR}/../Auto-Control -dir ${NDIR}/Masslessrfix -lmass 0.0 -qup 50 500 -sudakov 0 -rfix 1 -model 0

#mkdir ${NDIR}/Massiverfix
#${DIR}/../Auto-Control -dir ${NDIR}/Massiverfix -lmass 0.0196 -qup 50 500 -sudakov 0 -rfix 1 -model 0
######################################################################

################BGK#######################################
NDIR="${DIR}/BGK"
mkdir ${NDIR}

mkdir ${NDIR}/Massless
${DIR}/../Auto-Control -dir ${NDIR}/Massless -lmass 0.0 -qup 50 500 -sudakov 0 -rfix 0 -model 1

mkdir ${NDIR}/Massive
${DIR}/../Auto-Control -dir ${NDIR}/Massive -lmass 0.0196 -qup 50 500 -sudakov 0 -rfix 0 -model 1

mkdir ${NDIR}/Masslessrfix
${DIR}/../Auto-Control -dir ${NDIR}/Masslessrfix -lmass 0.0 -qup 50 500 -sudakov 0 -rfix 1 -model 1

mkdir ${NDIR}/Massiverfix
${DIR}/../Auto-Control -dir ${NDIR}/Massiverfix -lmass 0.0196 -qup 50 500 -sudakov 0 -rfix 1 -model 1
################Variants of BGK c mass=1.3 and new star prescription######################
NDIR="${DIR}/BGK-GWIA"

cp -r ${DIR}/BGK ${NDIR}

NDIR="${DIR}/BGK-13"

cp -r ${DIR}/BGK ${NDIR}

NDIR="${DIR}/BGK-GWIA-13"

cp -r ${DIR}/BGK ${NDIR}

################GBWS#######################################
NDIR="${DIR}/GBWS"
mkdir ${NDIR}

mkdir ${NDIR}/Massless
${DIR}/../Auto-Control -dir ${NDIR}/Massless -lmass 0.0 -qup 50 500 -sudakov 1 -rfix 0 -model 22

mkdir ${NDIR}/Massive
${DIR}/../Auto-Control -dir ${NDIR}/Massive -lmass 0.0196 -qup 50 500 -sudakov 1 -rfix 0 -model 22

mkdir ${NDIR}/Masslessrfix
${DIR}/../Auto-Control -dir ${NDIR}/Masslessrfix -lmass 0.0 -qup 50 500 -sudakov 1 -rfix 1 -model 22

mkdir ${NDIR}/Massiverfix
${DIR}/../Auto-Control -dir ${NDIR}/Massiverfix -lmass 0.0196 -qup 50 500 -sudakov 1 -rfix 1 -model 22
################Variants of GBWS c mass=1.3 and new star prescription######################
NDIR="${DIR}/GBWS-GWIA"

cp -r ${DIR}/GBWS ${NDIR}

NDIR="${DIR}/GBWS-13"

cp -r ${DIR}/GBWS ${NDIR}

NDIR="${DIR}/GBWS-GWIA-13"

cp -r ${DIR}/GBWS ${NDIR}

################BGKS#######################################
NDIR="${DIR}/BGKS"
mkdir ${NDIR}

mkdir ${NDIR}/Massless
${DIR}/../Auto-Control -dir ${NDIR}/Massless -lmass 0.0 -qup 50 500 -sudakov 1 -rfix 0 -model 3

mkdir ${NDIR}/Massive
${DIR}/../Auto-Control -dir ${NDIR}/Massive -lmass 0.0196 -qup 50 500 -sudakov 1 -rfix 0 -model 3

mkdir ${NDIR}/Masslessrfix
${DIR}/../Auto-Control -dir ${NDIR}/Masslessrfix -lmass 0.0 -qup 50 500 -sudakov 1 -rfix 1 -model 3

mkdir ${NDIR}/Massiverfix
${DIR}/../Auto-Control -dir ${NDIR}/Massiverfix -lmass 0.0196 -qup 50 500 -sudakov 1 -rfix 1 -model 3
################Variants of BGKS c mass=1.3 and new star prescription######################
NDIR="${DIR}/BGKS-GWIA"
cp -r ${DIR}/BGKS ${NDIR}

NDIR="${DIR}/BGKS-13"
cp -r ${DIR}/BGKS ${NDIR}

NDIR="${DIR}/BGKS-GWIA-13"
cp -r ${DIR}/BGKS ${NDIR}

NDIR="${DIR}/BGKS-GWIA-13-fixC"
cp -r ${DIR}/BGKS ${NDIR}

${DIR}/../Append "#define MASS_C2 1.69" ${DIR}/*13*/M*/M*

${DIR}/../Append "#define STAR 1" ${DIR}/*GWIA*/M*/M*

${DIR}/../Append "#define INDEPENDENT_C 0" ${DIR}/*fixC*/M*/M*

#${DIR}/../Append "#define PRINT_PROGRESS 1" ${DIR}/*G*/Massless/*500*


#for i in ./Models/*G*; do ./Auto_Run ${i}/M*/M*; done


