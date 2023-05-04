#! /usr/bin/env bash

gnuplot plotting.txt

for i in plots/*.tex 
do
	pdflatex --output-directory='../Paper/plots' ${i}
done
