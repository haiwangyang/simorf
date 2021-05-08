""" 
A script to generate exon structure file for certain species and transcript

Example script:
python prepare_whole_trans.py -s human -p 1000000

Parameter description:
s  = species
"""

__version__ = "1.0.1"
__email__ = "haiwang.yang@northwestern.edu"

import re
import argparse
from Bio.Seq import Seq
import random
from os import path
from pyfaidx import Fasta
from common_functions import get_revcom_dna, object_to_pickle, pickle_to_object


def get_concatenated_transcript_seq(species):
    fasta = Fasta("./annotation/" + species + ".transcript.fa")
    concatenated_transcript_seq = ''
    for transcript_id in fasta.keys():
        transcript_seq = str(fasta[transcript_id])
        concatenated_transcript_seq += transcript_seq
    return concatenated_transcript_seq

def get_shuffled_seq(seq):
    lst = list(seq)
    random.shuffle(lst)
    return "".join(lst)

def get_shuffled_seq2(seq):
    lst = list(seq)
    n = len(lst)
    return "".join([lst[random.randrange(n)] for _ in range(n)])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='please provide information')
    parser.add_argument('-s', '--species', type=str) # human
    parser.add_argument('-p', '--partition', type=int) # human

    args = parser.parse_args()
    species = args.species
    partition = args.partition

    concatenated_transcript_seq = get_concatenated_transcript_seq(species)
    shuffled_concatenated_transcript_seq = get_shuffled_seq(concatenated_transcript_seq)
    i = 0
    with open("data/trans/" + species + ".shuffled_concatenated_transcript.fa", "w") as w:
        while len(shuffled_concatenated_transcript_seq) > 0:
            piece, shuffled_concatenated_transcript_seq = shuffled_concatenated_transcript_seq[:partition], shuffled_concatenated_transcript_seq[partition:]
            w.write(">shuffled_concatenated_transcript_" +  str(i) + "\n" + piece + "\n")
            i += 1


    s_fasta = Fasta("data/trans/" + species + ".shuffled_concatenated_transcript.fa")

    #import re
    dct_trans_GT_AG = dict()
    dct_trans_GT_AG['GT'] = dict()
    dct_trans_GT_AG['AG'] = dict()
    
    for s_id in s_fasta.keys():
        dct_trans_GT_AG[s_id] = dict()
        dct_trans_GT_AG[s_id]['GT'] = []
        dct_trans_GT_AG[s_id]['AG'] = []
    
    for s_id in s_fasta.keys():
        s_seq = str(s_fasta[s_id])
        for m in re.finditer(r"GT", s_seq):
            dct_trans_GT_AG[s_id]['GT'].append(m.start())
    
        for m in re.finditer(r"AG", s_seq):
            dct_trans_GT_AG[s_id]['AG'].append(m.start())
    
    object_to_pickle(dct_trans_GT_AG, "data/trans/" + species + ".dct_trans_GT_AG.pickle")

