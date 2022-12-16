#! /usr/bin/env bash

DIR=../Run

GBW=${DIR}/GBW/Mass0.0-Qup650-Model0-Sud0
GBWS=${DIR}/GBWS-Fix-S/Mass0.0-Qup650-Model22-Sud1
BGK=${DIR}/BGK/Mass0.0-Qup650-Model1-Sud0
BGKS=${DIR}/BGKS-Fix-S/Mass0.0-Qup650-Model3-Sud1
SAVE=../Document/Plots
PLOT=./Plotting

${PLOT}/visualize-dipole.py \
	-s ${SAVE}/dipole \
	${GBW} \
	${GBWS} \
	${BGK} \
	${BGKS}

#${PLOT}/visualize-slope.py \
#	-s ${SAVE}/slope \
#	${GBW} \
#	${GBWS} \
#	${BGK} \
#	${BGKS}

#${PLOT}/visualize-F2.py \
#	-s ${SAVE}/F2 \
#	${GBW} \
#	${GBWS} \
#	${BGK} \
#	${BGKS}

#${PLOT}/visualize-FL.py \
#	-s ${SAVE}/FL \
#	${GBW} \
#	${GBWS} \
#	${BGK} \
#	${BGKS}

${PLOT}/visualize-gluon.py \
	-s ${SAVE}/gluon \
	${GBW} \
	${GBWS} \
	${BGK} \
	${BGKS}

${PLOT}/visualize-gluon-x.py \
	-s ${SAVE}/gluon-x \
	${GBW} \
	${GBWS} \
	${BGK} \
	${BGKS}

${PLOT}/visualize-ww-gluon.py \
	-s ${SAVE}/ww-gluon \
	${GBW} \
	${GBWS} \
	${BGK} \
	${BGKS}

#${PLOT}/visualize-critical.py \
#	-s ${SAVE}/critical \
#	${GBW} \
#	${GBWS} \
#	${BGK} \
#	${BGKS}
#./visualize-ww-critical.py \
#	-s ${SAVE}/ww-critical-GBW \
#	${GBW} \
#	${GBWS} \
#	${BGK} \
#	${BGKS}
	
