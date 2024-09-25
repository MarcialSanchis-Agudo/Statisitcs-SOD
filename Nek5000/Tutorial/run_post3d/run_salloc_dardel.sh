#!/bin/bash 

# salloc -N 1 -t 00:10:00 -A naiss2023-3-13 -p main

casename=duct
rm  -f $casename.sch
echo $casename > SESSION.NAME
echo $PWD/ >> SESSION.NAME

srun -n 128 ./nek5000 | tee log_01.txt
