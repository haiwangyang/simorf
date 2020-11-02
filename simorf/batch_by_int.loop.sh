for f in `cat list/golub/file.list`; do
    sbatch batch_by_int.sh human list/golub/$f output/golub/$f
done


