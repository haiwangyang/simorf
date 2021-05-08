cd /projects/b1080/hy/simorf/simorf/
for sc in `echo ATG CTG GTG TTG ACG`; do
	sbatch simulate_mimic_annotation.sh human $sc 100 1000
done

