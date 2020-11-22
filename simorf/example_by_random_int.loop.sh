folder=human_int
for f in `cat list/$folder/file.list`; do
    sbatch example_by_random_int.sh human list/$folder/$f 1 100 1000 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done

folder=mouse_int
for f in `cat list/$folder/file.list`; do
    sbatch example_by_random_int.sh mouse list/$folder/$f 1 200 1000 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done

folder=zebrafish_int
for f in `cat list/$folder/file.list`; do
    sbatch example_by_random_int.sh zebrafish list/$folder/$f 1 3000 1000 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done


folder=worm_int
for f in `cat list/$folder/file.list`; do
    sbatch example_by_random_int.sh worm list/$folder/$f 1 0 1000 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done

folder=yeast_int
for f in `cat list/$folder/file.list`; do
    sbatch example_by_random_int.sh yeast list/$folder/$f 1 0 1 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done

folder=arabidopsis_int
for f in `cat list/$folder/file.list`; do
    sbatch example_by_random_int.sh arabidopsis list/$folder/$f 1 0 1000 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done
