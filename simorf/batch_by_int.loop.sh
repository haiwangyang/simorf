folder=human
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_int.sh human list/$folder/$f output/$folder/$f
done


