#!/bin/bash -l

cd build
./compile_script.sh --clean

./compile_script.sh --compile

cd ..
cp build/nek5000 .

./run_local.sh

