#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBTTCH --job-name="hw"

cd /projects/b1080/hy/simorf/simorf/

s=$1   # species
olf=$2 # list/golub/xaa 
sn=$3
o=$4  # output/golub/xaa

python batch_by_trans_shuffling.py -s $s -olf $olf -sn $sn -o $o

