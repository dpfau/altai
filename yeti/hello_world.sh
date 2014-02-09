#!/bin/sh
# hello_world.sh
# Test for submitting jobs on GPU-enabled nodes

# Torque directives
#PBS -N hello_world
#PBS -W group_list=yetistats
#PBS -l walltime=00:05:00,mem=400mb,other=gpu
#PBS -M dbp2112@columbia.edu
#PBS -m abe
#PBS -V

#set output and error directories
#PBS -o localhost:/vega/stats/users/dbp2112/ahrens/results/
#PBS -e localhost:/vega/stats/users/dbp2112/ahrens/results/

echo "hello world"

#End of script