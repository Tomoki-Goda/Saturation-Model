#! /usr/bin/env bash 

rm ./dijet-HERA-plot/plots/*
rm ./dijet-HERA-plot/*.pdf

./plotting-scale.sh BGKkt BGKSkt "BGK \$k_t\$" "BGK \$k_t\$ + S" BGKSkt-ren

./plotting.sh  GBWr BGKr GBWkt BGKkt  "GBW \$r\$" "BGK \$r\$"  "GBW \$k_t\$" "BGK \$k_t\$" rkt

./plotting3.sh BGKr BGKkt BGKSkt  "BGK \$r\$"  "BGK \$k_t\$" "BGK \$k_t\$ +S" BGKrktS

for i in ./dijet-HERA-plot/plots/*.tex; do pdflatex ${i} --output-directory=./dijet-HERA-plot; done
rm ./dijet-HERA-plot/*log ./dijet-HERA-plot/*aux
