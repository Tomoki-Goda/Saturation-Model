#! /usr/bin/env bash

./table "sigma_0, x_0, lambda,chisq/dof" \
	"$r$-GBW,$r$-GBW-massive,$k_t$-GBW,rc-$k_t$-GBW" \
	./Run2/GBW/result.txt \
	./Run2/GBW-Massive/result.txt \
	./Run2/fixa-bjorx/result.txt \
	./Run2/runa-bjorx-4/result.txt \
	./kt-formula-report1/GBW-parameters.tex
./table "sigma_0, A_g, lambda_g,C1,mu102, chisq/dof" \
	"$r$-BGK, $k_t$-BGK" \
	./Run2/BGK/result.txt \
	./Run2/fixa-bjorx-BGK/result.txt \
	./kt-formula-report1/BGK-parameters.tex
