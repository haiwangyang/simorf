folder=human_int
for f in `cat list/$folder/file.list`; do
    sbatch example_by_random_int.sh human list/$folder/$f 1 100 1000 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done

folder=mouse_int
for f in `cat list/$folder/file.list`; do
    sbatch example_by_random_int.sh mouse list/$folder/$f 1 200 1000 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done
