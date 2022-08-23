#! /usr/bin/env bash

GBW=../Run1708/GBW/Mass0.0-Qup650-Model0-Sud0
GBWS=../Run1708/GBWS-Fix-S/Mass0.0-Qup650-Model22-Sud1
BGK=../Run1708/BGK/Mass0.0-Qup650-Model1-Sud0
BGKS=../Run1708/BGKS-Fix-S/Mass0.0-Qup650-Model3-Sud1
SAVE=/media/tomoki/TOMOKI-USB/Saturation-Model/Saturation-Notes/Run1708

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
./visualize-gluon.py \
	-s ${SAVE}/gluon-GBW \
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
./visualize-gluon.py \
	-s ${SAVE}/gluon-BGK \
	${BGK} \
	${BGKS}
./visualize-critical.py \
	-s ${SAVE}/critical-BGK \
	${BGK} \
	${BGKS}
