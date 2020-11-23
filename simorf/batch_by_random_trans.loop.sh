sbatch batch_by_random_trans.sh human list/golub.orf_id.use.list 2000 output/golub_trans.transFDRs.sim2000.v2.txt


folder=human_trans
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_random_trans.sh human list/$folder/$f 2000 output/$folder/$f
done

