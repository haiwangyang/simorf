folder=human_trans
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_trans_shuffling.sh human list/$folder/$f 2000 output/$folder/$f
done

folder=mouse_trans
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_trans_shuffling.sh mouse list/$folder/$f 2000 output/$folder/$f
done

folder=zebrafish_trans
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_trans_shuffling.sh zebrafish list/$folder/$f 2000 output/$folder/$f
done

folder=worm_trans
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_trans_shuffling.sh worm list/$folder/$f 2000 output/$folder/$f
done

folder=arabidopsis_trans
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_trans_shuffling.sh arabidopsis list/$folder/$f 2000 output/$folder/$f
done

folder=yeast_trans
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_trans_shuffling.sh yeast list/$folder/$f 2000 output/$folder/$f
done


