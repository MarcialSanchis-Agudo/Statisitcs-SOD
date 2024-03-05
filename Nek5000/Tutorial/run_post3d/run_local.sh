#!/bin/bash 

casename=duct
rm  -f $casename.sch
echo $casename > SESSION.NAME
echo $PWD/ >> SESSION.NAME

mpirun -n 4 ./nek5000 | tee log_01.txt
