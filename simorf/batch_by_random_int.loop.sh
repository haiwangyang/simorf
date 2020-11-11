folder=human
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_random_int.sh human list/$folder/$f 2000 100 1000 output/$folder/$f
done

folder=mouse_int
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_random_int.sh mouse list/$folder/$f 2000 200 1000 output/$folder/$f
done

folder=zebrafish_int
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_random_int.sh zebrafish list/$folder/$f 2000 3000 1000 output/$folder/$f
done

folder=worm_int
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_random_int.sh worm list/$folder/$f 2000 0 1000 output/$folder/$f
done

folder=yeast_int
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_random_int.sh yeast list/$folder/$f 2000 0 1 output/$folder/$f
done

folder=arabidopsis_int
for f in `cat list/$folder/file.list`; do
    sbatch batch_by_random_int.sh arabidopsis list/$folder/$f 2000 0 1000 output/$folder/$f
done
