#! /usr/bin/env bash

./table "sigma_0, x_0, lambda,chisq/dof" \
	"GBW-r,GBW-r-Massive, GBW, GBW-Thresh, GBW-Run, GBW-Run-Thresh" \
	./Run2/GBW/result.txt \
	./Run2/GBW-Massive/result.txt \
	./Run2/fixa-bjorx/result.txt \
	./Run2/fixa-bjorx2/result.txt \
	./Run2/runa-bjorx-4/result.txt \
	./Run2/runa-bjorx2-4/result.txt \
	./kt-formula-report1/GBW-parameters.tex
./table "sigma_0, A_g, lambda_g,C1,mu102, chisq/dof" \
	"BGK-r, BGK,BGK-Thresh" \
	./Run2/BGK/result.txt \
	./Run2/fixa-bjorx-BGK/result.txt \
	./Run2/fixa-bjorx-BGK2/result.txt \
	./kt-formula-report1/BGK-parameters.tex
