"""

A script to concatenate intergenic regions e.g., yeast

Example script:
python -i batch_by_random_int.py -s human -olf list/golub/xaa -sn 1000 -tp 100 -tn 1000 -o output/golub/xaa

Parameter description:
s  = species
olf = orf_list_file
sn = simulation_num
"""

__version__ = "1.0.1"
__email__ = "haiwang.yang@northwestern.edu"


from timeit import default_timer as timer
from common_functions import Orf, Transcript2, get_elements, get_A2B, object_to_pickle, pickle_to_object
from simulate_by_random_int import simulation_int, generate_merged_int, fetch_exon_structure
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
    parser.add_argument('-f', '--fasta', type=str) # intergenic fasta: /projects/b1080/hy/simorf/simorf/data/intergenic/yeast/yeast.annotation.intergenic.fa
    
    args = parser.parse_args()
    species = args.species
    fasta = args.fasta

    intergenic_fasta = Fasta(fasta)
    concatenated_intergenic = ""
    for name in intergenic_fasta.keys():
        concatenated_intergenic += str(intergenic_fasta[name])
   
    intergenic_id = "concatenated"
    with open("/projects/b1080/hy/simorf/simorf/data/intergenic/" + species + ".annotation.intergenic.fa", "w") as w:
        w.write(">" + intergenic_id + "\n" + concatenated_intergenic + "\n")
    with open("/projects/b1080/hy/simorf/simorf/data/intron/" + species + ".annotation.intron.fa", "w") as w:
        w.write(">" + intergenic_id + "\n" + concatenated_intergenic + "\n")

    dct_intergenic_GT_AG = {}
    dct_intergenic_GT_AG[intergenic_id] = dict()
    dct_intergenic_GT_AG[intergenic_id]['GT'] = []
    dct_intergenic_GT_AG[intergenic_id]['AG'] = []
    for m in re.finditer(r"GT", concatenated_intergenic):
        dct_intergenic_GT_AG[intergenic_id]['GT'].append(m.start())

    for m in re.finditer(r"AG", concatenated_intergenic):
        dct_intergenic_GT_AG[intergenic_id]['AG'].append(m.start())

    object_to_pickle(dct_intergenic_GT_AG, "/projects/b1080/hy/simorf/simorf/data/intergenic/" + species + ".dct_intergenic_GT_AG.pickle")
    object_to_pickle(dct_intergenic_GT_AG, "/projects/b1080/hy/simorf/simorf/data/intron/" + species + ".dct_intron_GT_AG.pickle")
