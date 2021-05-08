folder=human_trans
for f in `cat list/$folder/file.list`; do
    sbatch example_by_random_trans.sh human list/$folder/$f 1 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done

