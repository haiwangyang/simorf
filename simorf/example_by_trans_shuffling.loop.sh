folder=human_trans
for f in `cat list/$folder/file.list`; do
    sbatch example_by_trans_shuffling.sh human list/$folder/$f 1 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done


folder=mouse_shu
for f in `cat list/$folder/file.list`; do
    sbatch example_by_trans_shuffling.sh mouse list/$folder/$f 1 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done

folder=zebrafish_shu
for f in `cat list/$folder/file.list`; do
    sbatch example_by_trans_shuffling.sh zebrafish list/$folder/$f 1 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done

folder=worm_shu
for f in `cat list/$folder/file.list`; do
    sbatch example_by_trans_shuffling.sh worm list/$folder/$f 1 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done

folder=arabidopsis_shu
for f in `cat list/$folder/file.list`; do
    sbatch example_by_trans_shuffling.sh arabidopsis list/$folder/$f 1 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done

folder=yeast_shu
for f in `cat list/$folder/file.list`; do
    sbatch example_by_trans_shuffling.sh yeast list/$folder/$f 1 output/$folder/$f.example.fa output/$folder/$f.example.pep_len
done

