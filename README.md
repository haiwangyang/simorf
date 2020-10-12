# simorf
The tool is meant to conduct simulation to test if an given Open Reading Frame (ORF) is significantly longer or shorter than expectation

# usage
## 0. prepare annotation
put species fasta and genePred annotation in ./annotation folder
cd /projects/b1080/hy/simorf/simorf/annotation/
for s in `echo human mouse zebrafish worm yeast arabidopsis`; do ln -s /projects/b1080/hy/ribo/ribo/$s/annotation/genome.fa $s.fa; done

ln -s /projects/b1080/hy/ribo/ribo/human/annotation/annotation.gencodeplus.reduced.genePred human.annotation.genePred
ln -s /projects/b1080/hy/ribo/ribo/mouse/annotation/annotation.genePred mouse.annotation.genePred
ln -s /projects/b1080/hy/ribo/ribo/zebrafish/annotation/annotation.genePred zebrafish.annotation.genePred
ln -s /projects/b1080/hy/ribo/ribo/worm/annotation/annotation.added.genePred worm.annotation.genePred
ln -s /projects/b1080/hy/ribo/ribo/yeast/annotation/annotation.added.genePred yeast.annotation.genePred
ln -s /projects/b1080/hy/ribo/ribo/arabidopsis/annotation/annotation.added.genePred arabidopsis.annotation.genePred

## 1. prepare intergenic and intron region files for different species
python prepare_intergenic_region_files.py -s human

## 2. prepare exon structure file
python prepare_exon_structure_file.py -s human -ti ENST00000370418.7

## 3. run three simulations

### (a) by intergenic region
python simulate_by_intergenic_region.py -s human -oi "ENST00000370418.7:chr10:-|3|1863:63:90|uORF|ATG|253|1630|63|90" -esf data/exon_structure/human/ENST00000370418.7.exon_structure -cs 252 -ce 1629 -sn 100

### (b) by intron
python simulate_by_intron.py -s human -oi "ENST00000370418.7:chr10:-|3|1863:63:90|uORF|ATG|253|1630|63|90" -esf example/ENST00000370418.7.exon_structure -cs 252 -ce 1629 -sn 100

### (c) by transcript shuffling
python simulate_by_transcript_shuffling.py -s human -oi "ENST00000528722.5:chr11:+|2|954:7:88|uORF|CTG|566|824|7|88" -ts AAGTGGCTGGAGTCGGGAGGCTGGAAAGAGACTCCGAGAAAGTACCAGCGGAAGGCGGCCGCCGCTACGGCGATTCGCAGGGAGTAGCAGACGAAGACGGTGGCCGCCGCACTAGCCACCACGTGTGGAGGATAAACGGTCTACACGGCCATTCCGGCGCCGAGTCTAGGGAAAGAGTTAGCGACGACGGGGAAAGAAAATGTGAAGAGAGCGACCGCCGCTCCAGGGTCGCTGCAGGAAGCCTAAGTGCAGACGCCGGCTTCTCCCGCAGTGACTTGAGAAGGGTTCCAGGAAGAAGGTGCATTTTGGCAGCATACATGATGCAGTACGAGCTGGAGATGTAAAGCAGCTTTCAGAAATAGTGGTACGTGGAGCCAGCATTAATGAACTTGATGTTCTCCATAAGTTTACCCCTTTACATTGGGCAGCACATTCTGGAAGTTTGGAGTGTCTTCATTGGCTGCTCTGGCATGGAGCTGATATCACACACGTAACAACGAGAGGTTGGACAGCATCTCACATAGCTGCAATCAGGGGTCAGGATGCTTGTGTACAGGCTCTTATAATGAATGGAGCAAATCTGACAGCCCAGGATGACCGGGGATGCACTCCTTTACATCTTGCTGCAACTCATGGACATTCTTTCACTTTACAAATAATGCTCCGAAGTGGAGTGGATCCCAGTGTGACTGATAAGAGAGAATGGAGACCTGTGCATTATGCAGCTTTTCATGGGCGGCTTGGCTGCTTGCAACTTCTTGTTAAATGGGGTTGTAGCATAGAAGATGTGGACTACAATGGAAACCTTCCAGAACCTCCTTAGATCCCTGTGGAGCCTCTGATTCCTTGGCTTGGCATGCTGGACTTACATAATTTGGCACCTACCTACCTGCATGGCCTCATCTTTCTATATCTTTACACTACTGTCCAGGATCATTTTGTTTTTTTCTGAAG -cs 565 -ce 823 -sn 10000

