#! /usr/bin/env bash

dir=../Run2808
savedir=/media/tomoki/TOMOKI-USB/Saturation-Model/Saturation-Notes/Run2808

./make-table \
	-v "sigma_0 A_g lambda_g C mu02 chisq/dof" \
	${dir}/BGK/Mass0.0196-Qup650*/result.txt \
	${dir}/BGKS-Fix-S/Mass0.0196-Qup650*/result.txt \
	${dir}/BGK/Mass0.0-Qup650*/result.txt \
	${dir}/BGKS-Fix-S/Mass0.0-Qup650*/result.txt \
	-o ${savedir}/BGKS-Fix-S.tex


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

./make-table \
	-v "sigma_0 x_0 lambda chisq/dof" \
	${dir}/GBW/Mass0.0196-Qup650-*/result.txt \
	${dir}/GBWS-Fix-S/Mass0.0196-Qup650-*/result.txt \
	${dir}/GBW/Mass0.0-Qup650-*/result.txt \
	${dir}/GBWS-Fix-S/Mass0.0-Qup650-*/result.txt \
	-o ${savedir}/GBWS-Fix-S.tex

./make-table \
	-v "chisq/dof" \
	${dir}/GBW/Mass0.0-Qup5-*/result.txt \
	${dir}/GBWS-Fix-S/Mass0.0-Qup5-*/result.txt \
	${dir}/GBW/Mass0.0-Qup25-*/result.txt \
	${dir}/GBWS-Fix-S/Mass0.0-Qup25-*/result.txt \
	${dir}/GBW/Mass0.0-Qup50-*/result.txt \
	${dir}/GBWS-Fix-S/Mass0.0-Qup50-*/result.txt \
	${dir}/GBW/Mass0.0-Qup100-*/result.txt \
	${dir}/GBWS-Fix-S/Mass0.0-Qup100-*/result.txt \
	${dir}/GBW/Mass0.0-Qup650-*/result.txt \
	${dir}/GBWS-Fix-S/Mass0.0-Qup650-*/result.txt \
	-o ${savedir}/GBWS-Q2.tex
	
./make-table \
	-v "chisq/dof" \
	${dir}/BGK/Mass0.0-Qup5-*/result.txt \
	${dir}/BGKS-Fix-S/Mass0.0-Qup5-*/result.txt \
	${dir}/BGK/Mass0.0-Qup25-*/result.txt \
	${dir}/BGKS-Fix-S/Mass0.0-Qup25-*/result.txt \
	${dir}/BGK/Mass0.0-Qup50-*/result.txt \
	${dir}/BGKS-Fix-S/Mass0.0-Qup50-*/result.txt \
	${dir}/BGK/Mass0.0-Qup100-*/result.txt \
	${dir}/BGKS-Fix-S/Mass0.0-Qup100-*/result.txt \
	${dir}/BGK/Mass0.0-Qup650-*/result.txt \
	${dir}/BGKS-Fix-S/Mass0.0-Qup650-*/result.txt \
	-o ${savedir}/BGKS-Q2.tex

