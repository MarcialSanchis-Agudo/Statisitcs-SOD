#!/usr/bin/env bash 
if [ -z "$1" -o -z "$2" -o -z "$3" -o -z "$4" ]; then  
	echo "--||Usage: ./run.sh <CASE> <UTAU> <NU> <PORDER>"  
	exit 
fi 
SOD2D_SRC_DIR=/home/mech/sanchis/sanchis/sod2d_gitlab_save
	
../../gmsh $1.geo -order $4 -o $1.msh -0  
pvpython python_getCoord.py $1 
pvpython python_checkYp.py $1 $2 $3 1 
rm -f *Coord.dat *.coord chan.box 

python3 ../gmsh2sod2d.py channel -p 2 -r $4
mpirun -np 1 ../tool_meshConversorPar channel.dat
