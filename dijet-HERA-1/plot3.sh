#! /usr/bin/env bash
./plotting.txt GBWkt "GBW \$k_t\$" BGKkt "BGK \$k_t\$" BGKGBW 

./plotting.txt BGKr "BGK \$r\$" BGKkt "BGK \$k_t\$" BGK 

#./plotting.txt BGKSkt "BGK+S " BGKkt "BGK" BGKS 

./plotting.txt GBWr "GBW \$r\$" GBWkt "GBW \$k_t\$" GBW 

#./plotting.txt GBWSkt "GBW+S " GBWkt "GBW" GBWS 

#./plotting.txt BGKSkt "BGK+S " KS "KS" KS 

#./plotting.txt BGKSkt "BGK+S " rcBK "rcBK" rcBK 
./plotting.txt GBWkt "GBW \$k_t\$" rcGBWkt "rcGBW \$k_t\$"  rcGBWkt

./plotting.txt BGKkt "BGK \$k_t\$" BGKSkt "BGK+S \$k_t\$" BGKkt

./plotting.txt GBWkt "GBW \$k_t\$" GBWSkt "GBW+S \$k_t\$" GBWkt

##############################################################

./plotting3.txt GBWkt "GBW \$k_t\$" rcGBWkt "rcGBW \$k_t\$" rcGBWSkt "rcGBW+S \$k_t\$" rcGBW

./plotting3.txt BGKr "BGK \$r\$" BGKkt "BGK \$k_t\$" BGKSkt "BGK+S \$k_t\$" BGK3

./plotting3.txt GBWr "GBW \$r\$" GBWkt "GBW \$k_t\$" GBWSkt "GBW+S \$k_t\$" GBW3

./plotting3.txt KS "KS " rcBK "rcBK" BGKSkt "BGK+S \$k_t\$" rcBKKS
