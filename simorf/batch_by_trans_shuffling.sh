#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBTTCH --job-name="hw"

cd /projects/b1080/hy/simorf/simorf/

s=$1   # species
olf=$2 # list/golub/xaa 
o=$3   # output/golub/xaa

python batch_by_trans_shuffling.py -s human -olf $olf -sn 2000 -o $o

