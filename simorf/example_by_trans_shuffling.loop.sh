folder=human_shu
for f in `cat list/$folder/file.list`; do
    sbatch example_by_trans_shuffling.sh human list/$folder/$f 1 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done


folder=mouse_shu
for f in `cat list/$folder/file.list`; do
    sbatch example_by_trans_shuffling.sh mouse list/$folder/$f 1 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done


