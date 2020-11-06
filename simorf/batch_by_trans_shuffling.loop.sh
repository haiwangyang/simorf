folder=human_shu
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_trans_shuffling.sh human list/$folder/$f output/$folder/$f
done


