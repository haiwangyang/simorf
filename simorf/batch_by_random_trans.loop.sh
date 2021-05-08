sbatch batch_by_random_trans.sh human list/golub.orf_id.use.list 2000 output/golub_trans.transFDRs.sim2000.v2.txt


folder=human_trans
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_random_trans.sh human list/$folder/$f 2000 output/$folder/$f
done

folder=mouse_trans
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_random_trans.sh mouse list/$folder/$f 2000 output/$folder/$f
done

folder=zebrafish_trans
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_random_trans.sh zebrafish list/$folder/$f 2000 output/$folder/$f
done

folder=worm_trans
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_random_trans.sh worm list/$folder/$f 2000 output/$folder/$f
done

folder=yeast_trans
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_random_trans.sh yeast list/$folder/$f 2000 output/$folder/$f
done


folder=arabidopsis_trans
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_random_trans.sh arabidopsis list/$folder/$f 2000 output/$folder/$f
done

