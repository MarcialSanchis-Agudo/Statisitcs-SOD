#!/bin/bash 

# salloc -N 1 -t 02:00:00 -A 2021-29 -p main

casename=duct
rm  -f $casename.sch
echo $casename > SESSION.NAME
echo $PWD/ >> SESSION.NAME

srun -n 4 ./nek5000 | tee log_01.txt
