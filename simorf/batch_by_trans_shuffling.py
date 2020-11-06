"""

A script to batch test ORFs by int random fetching method

Example script:
python -i batch_by_trans_shuffling.py -s human -olf list/golub/xaa -sn 1000 -o output/golub/xaa

Parameter description:
s  = species
olf = orf_list_file
sn = simulation_num
"""

__version__ = "1.0.1"
__email__ = "haiwang.yang@northwestern.edu"


from timeit import default_timer as timer
from common_functions import Orf, Transcript2, get_elements, get_A2B, object_to_pickle, pickle_to_object
from simulate_by_trans_shuffling import simulation_shu
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
    parser.add_argument('-o', '--output', type=str)

    args = parser.parse_args()
    species = args.species
    orf_list_file = args.orf_list_file
    simulation_num = args.simulation_num
    output = args.output

    transcript_fasta = Fasta("annotation/" + species + ".transcript.fa")
    cds_range_file = "annotation/" + species + ".transcript.cds_range"
    transcript_id_to_transcript_len = get_A2B(cds_range_file, 1, 3)
    transcript_id_to_canonical_start = get_A2B(cds_range_file, 1, 4)
    transcript_id_to_UTR3_len = get_A2B(cds_range_file, 1, 6)

    with open(output, "w") as w:
        w.write("\t".join(["orf_id", "orf_type", "orf_pep_len", "FDR_S_all", "FDR_L_all", "FDR_S_main", "FDR_L_main", "FDR_S_uORF", "FDR_L_uORF", "FDR_S_ouORF", "FDR_L_ouORF", "shorter_all", "longer_all", "total_all", "shorter_main", "longer_main", "total_main", "shorter_uORF", "longer_uORF", "total_uORF", "shorter_ouORF", "longer_ouORF", "total_ouORF"]) + "\n")
        for orf_id in get_elements(orf_list_file):
            print(orf_id)
            orf = Orf(species, orf_id)
            transcript_id = orf.transcript_id
            transcript_seq = str(transcript_fasta[transcript_id])

            orf_type = orf.orf_type
            if orf_type != "noncoding":
                transcript_len = int(transcript_id_to_transcript_len[transcript_id])
                canonical_start = int(transcript_id_to_canonical_start[transcript_id])
                UTR3_len = int(transcript_id_to_UTR3_len[transcript_id])
                canonical_end = transcript_len - UTR3_len
                pep_len, shorter_all, longer_all, total_all, shorter_main, longer_main, total_main, shorter_uorf, longer_uorf, total_uorf, shorter_ouorf, longer_ouorf, total_ouorf = simulation_shu(species, orf_id, transcript_seq, canonical_start, canonical_end, simulation_num)
            else:
                canonical_start, canonical_end = 0, 0
                pep_len, shorter_all, longer_all, total_all, shorter_main, longer_main, total_main, shorter_uorf, longer_uorf, total_uorf, shorter_ouorf, longer_ouorf, total_ouorf = simulation_shu(species, orf_id, transcript_seq, canonical_start, canonical_end, simulation_num)

            
            fdr_S_all = round(shorter_all/total_all, 16)
            fdr_L_all = round(longer_all/total_all, 16)

            fdr_S_main = round(shorter_main/total_main, 16)
            fdr_L_main = round(longer_main/total_main, 16)

            fdr_S_uorf = round(shorter_uorf/total_uorf, 16)
            fdr_L_uorf = round(longer_uorf/total_uorf, 16)

            fdr_S_ouorf = round(shorter_ouorf/total_ouorf, 16)
            fdr_L_ouorf = round(longer_ouorf/total_ouorf, 16)

            w.write("\t".join([str(_) for _ in [orf_id, orf_type, orf.pep_len, fdr_S_all, fdr_L_all, fdr_S_main, fdr_L_main, fdr_S_uorf, fdr_L_uorf, fdr_S_ouorf, fdr_L_ouorf, shorter_all, longer_all, total_all, shorter_main, longer_main, total_main, shorter_uorf, longer_uorf, total_uorf, shorter_ouorf, longer_ouorf, total_ouorf]]) + "\n")
             

