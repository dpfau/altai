#!/bin/sh
# hello_world.sh
# Test script to do some sanity checks

# Torque directives
#PBS -N hello_world
#PBS -W group_list=yetistats
#PBS -l walltime=1:00:00:00,mem=8192mb,other=gpu
#PBS -M dbp2112@columbia.edu
#PBS -m abe
#PBS -V

#set output and error directories
#PBS -o localhost:/vega/stats/users/dbp2112/ahrens/results/
#PBS -e localhost:/vega/stats/users/dbp2112/ahrens/results/

echo $RNG
if [[ -z "$RNG" ]]; then
	echo "Hello World"
else
	echo "Goodbye World"
fi

#End of script
