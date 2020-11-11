#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=90000
#SBATCH -t 48:00:00
#SBTTCH --job-name="hw"

cd /projects/b1080/hy/simorf/simorf/

s=$1   # species
olf=$2 # list/golub/xaa 
sn=$3
tp=$4
tn=$5
o=$6   # output/golub/xaa

python batch_by_random_int.py -s $s -olf $olf -sn $sn -tp $tp -tn $tn -o $o

