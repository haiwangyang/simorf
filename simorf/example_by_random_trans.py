"""

A script to batch test ORFs by int random fetching method

Example script:
python -i example_by_random_trans.py -s human -olf list/golub/xaa -sn 1 -o1 output/golub/xaa.example.fa -o2 output/golub/xaa.example.pep_len

Parameter description:
s  = species
olf = orf_list_file
sn = simulation_num
"""

__version__ = "1.0.1"
__email__ = "haiwang.yang@northwestern.edu"


from timeit import default_timer as timer
from common_functions import Orf, Transcript2, get_elements, get_A2B, object_to_pickle, pickle_to_object
from simulate_by_random_trans import example_trans, fetch_exon_structure
import subprocess
import argparse
import random
import os
import re
import pickle
from Bio.Seq import Seq
from pyfaidx import Fasta


if __name__ == "__main__":
    # fetching parameters
    parser = argparse.ArgumentParser(description='please provide information')
    parser.add_argument('-s', '--species', type=str) # human
    parser.add_argument('-olf', '--orf_list_file', type=str) # orf list file path
    parser.add_argument('-sn', '--simulation_num', type=int) # simulation num
    parser.add_argument('-tp', '--top_pos', type=int) # pick the position of longest intergenic/intron to boost speed
    parser.add_argument('-tn', '--top_num', type=int) # pick the number of longest intergenic/intron to boost speed
    parser.add_argument('-o1', '--output1', type=str) # example orf cds seq
    parser.add_argument('-o2', '--output2', type=str) # example orf pep len

    args = parser.parse_args()
    species = args.species
    orf_list_file = args.orf_list_file
    simulation_num = args.simulation_num
    output1 = args.output1
    output2 = args.output2

    s_fasta = Fasta("./data/trans/" + species + ".shuffled_concatenated_transcript.fa")
    s_ids = list(s_fasta.keys())
    dct_trans_GT_AG = pickle_to_object("/projects/b1080/hy/simorf/simorf/data/trans/" + species + ".dct_trans_GT_AG.pickle")

    cds_range_file = "annotation/" + species + ".transcript.cds_range"
    transcript_id_to_transcript_len = get_A2B(cds_range_file, 1, 3)
    transcript_id_to_canonical_start = get_A2B(cds_range_file, 1, 4)
    transcript_id_to_UTR3_len = get_A2B(cds_range_file, 1, 6)

    with open(output1, "w") as w1, open(output2, "w") as w2:
        for orf_id in get_elements(orf_list_file):
            print(orf_id)
            orf = Orf(species, orf_id)
            transcript_id = orf.transcript_id

            exon_structure_file = "data/exon_structure/" + species + "/" + transcript_id + ".exon_structure"
            if not os.path.exists(exon_structure_file):
                #os.system("python prepare_exon_structure_file.py -s " + species + " -ti " + transcript_id)
                subprocess.run(["python", "prepare_exon_structure_file.py", "-s", species, "-ti", transcript_id])
            dct_exon_structure = fetch_exon_structure(exon_structure_file)
            exon_structure = dct_exon_structure[orf.transcript_id]
           

            orf_type = orf.orf_type
            if orf_type != "noncoding":
                transcript_len = int(transcript_id_to_transcript_len[transcript_id])
                canonical_start = int(transcript_id_to_canonical_start[transcript_id])
                UTR3_len = int(transcript_id_to_UTR3_len[transcript_id])
                canonical_end = transcript_len - UTR3_len
                lst_sim_all_orf_cds, lst_sim_main_orf_cds, lst_sim_uorf_cds, lst_sim_overlapping_uorf_cds, lst_sim_all_orf_pep_len, lst_sim_main_orf_pep_len, lst_sim_uorf_pep_len, lst_sim_overlapping_uorf_pep_len  = example_trans(species, orf_id, exon_structure_file, canonical_start, canonical_end,  simulation_num, s_fasta, s_ids, dct_trans_GT_AG)
            else:
                transcript_len = int(transcript_id_to_transcript_len[transcript_id])
                lst_sim_all_orf_cds, lst_sim_main_orf_cds, lst_sim_uorf_cds, lst_sim_overlapping_uorf_cds, lst_sim_all_orf_pep_len, lst_sim_main_orf_pep_len, lst_sim_uorf_pep_len, lst_sim_overlapping_uorf_pep_len = example_trans(species, orf_id, exon_structure_file, 0, 0,  simulation_num, s_fasta, s_ids, dct_trans_GT_AG)
            
            i = 0
            for _ in lst_sim_all_orf_cds:
                w1.write(">" + orf_id + "__all" + str(i)  + "\n" + _ + "\n")
                i += 1

            i = 0
            for _ in lst_sim_main_orf_cds:
                w1.write(">" + orf_id + "__main" + str(i)  + "\n" + _ + "\n")
                i += 1

            i = 0
            for _ in lst_sim_uorf_cds:
                w1.write(">" + orf_id + "__uORF" + str(i)  + "\n" + _ + "\n")
                i += 1

            i = 0
            for _ in lst_sim_overlapping_uorf_cds:
                w1.write(">" + orf_id + "__ouORF" + str(i)  + "\n" + _ + "\n")
                i += 1
 
            w2.write(orf_id + "\tall"  + "\t" + ",".join([str(_) for _ in lst_sim_all_orf_pep_len]) + "\n")
            w2.write(orf_id + "\tmain"  + "\t" + ",".join([str(_) for _ in lst_sim_main_orf_pep_len]) + "\n")
            w2.write(orf_id + "\tuORF"  + "\t" + ",".join([str(_) for _ in lst_sim_uorf_pep_len]) + "\n")
            w2.write(orf_id + "\touORF"  + "\t" + ",".join([str(_) for _ in lst_sim_overlapping_uorf_pep_len]) + "\n")


            
