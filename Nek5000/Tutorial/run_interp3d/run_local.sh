#!/bin/bash 



# create a mesh for interpolation first

rm ZSTAT/a*
rm ZSTAT/b*

cp ../run_post3d/duct.ma2 .
cp ../run_post3d/duct.re2 .

cp ../run_post3d/a* ZSTAT/
cp ../run_post3d/b* ZSTAT/

casename=duct
rm  -f $casename.sch
echo $casename > SESSION.NAME
echo $PWD/ >> SESSION.NAME

rm log_*

mpirun -n 4 ./nek5000a | tee log_01a.txt

rm -rf output_a
mkdir output_a
mv ZSTAT/L* ZSTAT/U* ZSTAT/V* ZSTAT/W* output_a
mv history2.txt output_a

mpirun -n 4 ./nek5000b | tee log_01b.txt

rm -rf output_b
mkdir output_b
mv ZSTAT/L* ZSTAT/U* ZSTAT/V* ZSTAT/W* output_b
mv history2.txt output_b
