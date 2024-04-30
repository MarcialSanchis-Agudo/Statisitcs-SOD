#!/bin/bash 


casename=duct
rm  -f $casename.sch
echo $casename > SESSION.NAME
echo $PWD/ >> SESSION.NAME

srun -n 256 ./nek5000 | tee log_01.txt

folder=output_01
mkdir $folder
mv lap* la2* s??duct0.f* t??duct0.f* pts* duct0.f* log_* $folder


