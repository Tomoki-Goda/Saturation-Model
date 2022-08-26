#! /usr/bin/env bash

dir=../Run1108


./Write_TeX \
	${dir}/BGK/Mass0.0196-Qup100*/result.txt \
	${dir}/BGKS-Fix-S/Mass0.0196-Qup100*/result.txt \
	${dir}/BGK/Mass0.0-Qup100*/result.txt \
	${dir}/BGKS-Fix-S/Mass0.0-Qup100*/result.txt \
	${dir}/BGKS-Fix-S-100.tex

./Write_TeX \
	${dir}/BGK/Mass0.0196-Qup650*/result.txt \
	${dir}/BGKS-Fix-S/Mass0.0196-Qup650*/result.txt \
	${dir}/BGK/Mass0.0-Qup650*/result.txt \
	${dir}/BGKS-Fix-S/Mass0.0-Qup650*/result.txt \
	${dir}/BGKS-Fix-S-650.tex


#for i in 100 650
#do
#	for j in 0.0 0.0196
#	do
#		./Write_TeX \
#		 ../Run2207/GBW/Mass${j}-Qup${i}-*/result.txt \
#		 ../Run2207/GBWS-Fix-S/Mass${j}-Qup${i}-*/result.txt \
#		 ../Run2207/GBWS-Fix-C/Mass${j}-Qup${i}-*/result.txt \
#		 ../Run2207/GBWS-Fix-rmax/Mass${j}-Qup${i}-*/result.txt \
#		 ../Run2207/GBWS/Mass${j}-Qup${i}-*/result.txt \
#		 ../Saturation-Model2207/GBWS-Fix-S-${j}-${i}.tex
#	done
#done

./Write_TeX \
	${dir}/GBW/Mass0.0-Qup100-*/result.txt \
	${dir}/GBWS-Fix-S/Mass0.0-Qup100-*/result.txt \
	${dir}/GBW/Mass0.0196-Qup100-*/result.txt \
	${dir}/GBWS-Fix-S/Mass0.0196-Qup100-*/result.txt \
	${dir}/GBWS-Fix-S-100.tex
./Write_TeX \
	${dir}/GBW/Mass0.0-Qup650-*/result.txt \
	${dir}/GBWS-Fix-S/Mass0.0-Qup650-*/result.txt \
	${dir}/GBW/Mass0.0196-Qup650-*/result.txt \
	${dir}/GBWS-Fix-S/Mass0.0196-Qup650-*/result.txt \
	${dir}/GBWS-Fix-S-650.tex
