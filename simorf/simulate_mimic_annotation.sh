#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 48:00:00
#SBTTCH --job-name="hw"

cd /projects/b1080/hy/simorf/simorf/

s=$1   # species
sc=$2
tp=$3
tn=$4

python simulate_mimic_annotation.py -s $s -sc $sc -tp $tp -tn $tn

