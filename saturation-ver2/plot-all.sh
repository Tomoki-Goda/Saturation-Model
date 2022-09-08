#! /usr/bin/env bash

DIR=../Run2808

GBW=${DIR}/GBW/Mass0.0-Qup650-Model0-Sud0
GBWS=${DIR}/GBWS-Fix-S/Mass0.0-Qup650-Model22-Sud1
BGK=${DIR}/BGK/Mass0.0-Qup650-Model1-Sud0
BGKS=${DIR}/BGKS-Fix-S/Mass0.0-Qup650-Model3-Sud1
SAVE=/media/tomoki/TOMOKI-USB/Saturation-Model/Saturation-Notes/Run2808

./visualize-dipole.py \
	-s ${SAVE}/dipole-GBW \
	${GBW} \
	${GBWS}
./visualize-slope.py \
	-s ${SAVE}/slope-GBW \
	${GBW} \
	${GBWS}

./visualize-F2.py \
	-s ${SAVE}/F2-GBW \
	${GBW} \
	${GBWS}

./visualize-FL.py \
	-s ${SAVE}/FL-GBW \
	${GBW} \
	${GBWS}

./visualize-gluon.py \
	-s ${SAVE}/gluon-GBW \
	${GBW} \
	${GBWS}
./visualize-gluon-x.py \
	-s ${SAVE}/gluon-x-GBW \
	${GBW} \
	${GBWS}
./visualize-ww-gluon.py \
	-s ${SAVE}/ww-gluon-GBW \
	${GBW} \
	${GBWS}
./visualize-critical.py \
	-s ${SAVE}/critical-GBW \
	${GBW} \
	${GBWS}
	
./visualize-dipole.py \
	-s ${SAVE}/dipole-BGK \
	${BGK} \
	${BGKS}

./visualize-slope.py \
	-s ${SAVE}/slope-BGK \
	${BGK} \
	${BGKS}

./visualize-F2.py \
	-s ${SAVE}/F2-BGK \
	${BGK} \
	${BGKS}

./visualize-FL.py \
	-s ${SAVE}/FL-BGK \
	${BGK} \
	${BGKS}

./visualize-gluon.py \
	-s ${SAVE}/gluon-BGK \
	${BGK} \
	${BGKS}
./visualize-gluon-x.py \
	-s ${SAVE}/gluon-x-BGK \
	${BGK} \
	${BGKS}
./visualize-ww-gluon.py \
	-s ${SAVE}/ww-gluon-BGK \
	${BGK} \
	${BGKS}
./visualize-critical.py \
	-s ${SAVE}/critical-BGK \
	${BGK} \
	${BGKS}
