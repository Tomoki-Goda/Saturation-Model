#! /usr/bin/env bash 


./plotting-scale.sh BGKkt BGKSkt "BGK \$k_t\$" "BGK \$k_t\$ + S" BGKSkt

./plotting.sh  GBWr BGKr GBWkt BGKkt  "GBW \$r\$" "BGK \$r\$"  "GBW \$k_t\$" "BGK \$k_t\$" rkt

./plotting3.sh BGKr BGKkt BGKSkt  "BGK \$r\$"  "BGK \$k_t\$" "BGK \$k_t\$ +S" BGKrktS

for i in ./dijet-HERA-plot/plots/*.tex; do pdflatex ${i} --output-directory=./dijet-HERA-plot; done
