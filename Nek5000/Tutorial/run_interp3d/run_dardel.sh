#!/bin/bash 

# for dardel: 
# salloc -N 1 -t 02:00:00 -A 2021-29 -p main
# ml swap PrgEnv-cray/8.1.0 PrgEnv-gnu/8.1.0


# create a mesh for interpolation first
cp ../run_tutorial/duct.ma2 .
cp ../run_tutorial/duct.re2 .
cp ../run_post3d/a* ZSTAT/
cp ../run_post3d/b* ZSTAT/

casename=duct
rm  -f $casename.sch
echo $casename > SESSION.NAME
echo $PWD/ >> SESSION.NAME

rm log_*

srun -n 128 ./nek5000a | tee log_01a.txt

rm -rf output_a
mkdir output_a
mv ZSTAT/L* ZSTAT/U* ZSTAT/V* ZSTAT/W* output_a
mv history2.txt output_a

srun -n 128 ./nek5000b | tee log_01b.txt

rm -rf output_b
mkdir output_b
mv ZSTAT/L* ZSTAT/U* ZSTAT/V* ZSTAT/W* output_b
mv history2.txt output_b
