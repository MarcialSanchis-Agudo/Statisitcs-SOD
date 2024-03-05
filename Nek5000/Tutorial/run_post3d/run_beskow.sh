#!/bin/bash 
# The -l above is required to get the full environment with modules

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH --exclusive
#SBATCH -A 2020-3-5

# The name of the script is myjob
#SBATCH -J run

# Wall-clock time will be given to this job
#SBATCH -t 04:59:00

# Number of nodes
#SBATCH -N 16
# Number of MPI processes per node (the following is actually the default)
#SBATCH --ntasks-per-node=32
# Number of MPI processes.
#SBATCH -n 512

#SBATCH -e error_file.e
#SBATCH -o output_file.o

casename=duct
rm  -f $casename.sch
echo $casename > SESSION.NAME
echo $PWD/ >> SESSION.NAME

##module add openmpi/intel

# Run the executable named myexe 
# and write the output into my_output_file
srun -n 512 ./nek5000 > obstacle.output
