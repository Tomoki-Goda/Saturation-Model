#! /usr/bin/bash 

rm "./*.o"

for i in "test" "testcheb" "graph" "main" "simps" "bessels" "simps2d"
do rm ${i}
done
