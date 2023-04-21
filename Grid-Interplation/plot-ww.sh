#! /usr/bin/env bash

./KSgluon-example ../gluon-grids/GBWr.dat 2 ./GBWr.txt 17
./KSgluon-example ../gluon-grids/GBWkt.dat 2 ./GBWkt.txt 17
./KSgluon-example ../gluon-grids/BGKr.dat 2 ./BGKr.txt 17
./KSgluon-example ../gluon-grids/BGKkt.dat 2 ./BGKkt.txt 17

for i in 17 67
do
./KSgluon-example ../gluon-grids/GBWSr.dat 3 ./GBWSr-${i}.txt ${i}
./KSgluon-example ../gluon-grids/GBWSkt.dat 3 ./GBWSkt-${i}.txt ${i}
./KSgluon-example ../gluon-grids/BGKSr.dat 3 ./BGKSr-${i}.txt ${i}
./KSgluon-example ../gluon-grids/BGKSkt.dat 3 ./BGKSkt-${i}.txt ${i}
done

./KSgluon-example /usr/local/share/tmdlib/KS-WeizWill-2017/KS-WeizWill-2017_0000.dat 2 ./KS.txt 17
