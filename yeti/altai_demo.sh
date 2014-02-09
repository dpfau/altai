#!/bin/sh
# altai_demo.sh
# Torque script to run GPU-optimized Altai pipeline on spontaneous zebrafish data.

# Torque directives
#PBS -N altai_demo
#PBS -W group_list=yetistats
#PBS -l walltime=24:00:00,mem=8gb,other=gpu
#PBS -M dbp2112@columbia.edu
#PBS -m abe
#PBS -V

#set output and error directories
#PBS -o localhost:/vega/stats/users/dbp2112/ahrens/results/
#PBS -e localhost:/vega/stats/users/dbp2112/ahrens/results/

matlab-R2012b -nosplash -nodisplay -nodesktop -r 'demo' > ../../results/matoutfile