#!/bin/bash 

if [ $1 = clean ]
then
    rm rs6duct0.f* duct0.f* s??duct0.f000* t??duct0.f000* lap* la2* pts* log*
else
    casename=duct
    rm  -f $casename.sch
    echo $casename > SESSION.NAME
    echo $PWD/ >> SESSION.NAME

    mpirun -n 4 ./nek5000 | tee log_09.txt

    #mkdir output_run_04
    #mv lap* la2* s??duct0.f* t??duct0.f* pts* duct0.f* log_* output_run_04
fi

