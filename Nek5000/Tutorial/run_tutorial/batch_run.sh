#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A naiss2023-3-13

# The name of the script is myjob
#SBATCH -J Stat18

# The partition
#SBATCH -p main

# 10 hours wall-clock time will be given to this job
#SBATCH -t 15:00:00

# Number of nodes
#SBATCH --nodes=2
#SBATCh --ntasks-per-node=128
#SBATCH --mail-user=walesiak@kth.se
#SBATCH --mail-type=ALL



# Run the executable named myexe
# and write the output into my_output_file
# ml PDC
# ml PrgEnv-gnu
# Compiling

# cd build
# ./compile_script.sh --clean

# ./compile_script.sh --compile

# Change folder and copy exec
# cd ..
# cp build/nek5000 .

# Run
./run_dardel.sh
