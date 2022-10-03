#! /usr/bin/env bash


#./Utilities/plotf2.py -s ../Run1708/GBWS-Fix-S/Mass0.0-Qup650-Model22-Sud1/grid-plot ../Run1708/GBW/Mass0.0-Qup650-Model0-Sud0/fcn.txt ../Run1708/GBWS-Fix-S/Mass0.0-Qup650-Model22-Sud1/fcn.txt ../Run1708/GBWS-Fix-S/Mass0.0-Qup650-Model22-Sud1/data.txt

#./Utilities/plotf2.py -s ../Run1708/BGKS-Fix-S/Mass0.0-Qup650-Model3-Sud1/grid-plot ../Run1708/BGK/Mass0.0-Qup650-Model1-Sud0/fcn.txt ../Run1708/BGKS-Fix-S/Mass0.0-Qup650-Model3-Sud1/fcn.txt ../Run1708/BGKS-Fix-S/Mass0.0-Qup650-Model3-Sud1/data.txt

#SAVEDIR=..
RESDIR=../Run
SAVEDIR=../Document/Plots

./Utilities/plot-heavy-f2.py -s ${SAVEDIR}/F2-data-GBW ${RESDIR}/GBW/Mass0.0-Qup650-Model0-Sud0/F2c ${RESDIR}/GBWS-Fix-S/Mass0.0-Qup650-Model22-Sud1/F2c

./Utilities/plot-heavy-f2.py -s ${SAVEDIR}/F2-data-BGK ${RESDIR}/BGK/Mass0.0-Qup650-Model1-Sud0/F2c ${RESDIR}/BGKS-Fix-S/Mass0.0-Qup650-Model3-Sud1/F2c


