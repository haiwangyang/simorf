#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBTTCH --job-name="hw"

cd /projects/b1080/hy/simorf/simorf/

s=$1   # species
olf=$2 # list/golub/xaa 
o=$3   # output/golub/xaa

python batch_by_random_int.py -s human -olf $olf -sn 1000 -tp 100 -tn 1000 -o $o

