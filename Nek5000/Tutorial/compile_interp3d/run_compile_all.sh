#!/bin/bash

./compile_script_dardel.sh --clean

cp duct.usr_a duct.usr

./compile_script_dardel.sh  --compile

mv nek5000 ../run_interp3d/nek5000a

cp duct.usr_b duct.usr

./compile_script_dardel.sh --compile

mv nek5000 ../run_interp3d/nek5000b
