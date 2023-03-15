#! /usr/bin/env bash
#./plot-F2.py -l 'rc-$k_t$:$k_t$:$r$' -c 'r:magenta:b' -p '- -- -.' ./Run3/runa-bjorx-4/data.txt ./Run3/runa-bjorx-4/F2.txt ./Run3/fixa-bjorx/F2.txt ./Run3/GBW/F2.txt
#./plot-hist.py -l 'rc-$k_t$:$k_t$:$r$' -c 'r:magenta:b' -p '- -- -.' ./Run3/runa-bjorx-4/data.txt ./Run3/fixa-bjorx/data.txt ./Run3/GBW/data.txt

#./plot-F2.py -l '$k_t$:$r$' -c 'r:b' -p '- -.' ./Run3/fixa-bjorx-BGK/data.txt ./Run3/fixa-bjorx-BGK/F2.txt ./Run3/BGK/F2.txt
#./plot-hist.py -l '$k_t$:$r$' -c 'r:b' -p '- -.' ./Run3/fixa-bjorx-BGK/data.txt ./Run3/BGK/data.txt

./plot-F2.py -L -l 'rc-$k_t$:$k_t$:$r$' -c 'r:magenta:b' -p '- -- -.' ./Run3/runa-bjorx-4/data.txt ./Run3/runa-bjorx-4/F2.txt ./Run3/fixa-bjorx/F2.txt ./Run3/GBW/F2.txt

./plot-F2.py -L -l '$k_t$:$r$' -c 'r:b' -p '- -.' ./Run3/fixa-bjorx-BGK/data.txt ./Run3/fixa-bjorx-BGK/F2.txt ./Run3/BGK/F2.txt
