from simulate_by_random_int import simulation_int
from pyfaidx import Fasta
from common_functions import Orf, pickle_to_object, get_A2B
import sys

species = "human"
top_pos = 100
top_num = 1000
#orf_id = "ENST00000451962.5:chr7:-|12|687:88:247|noncoding|ATG|0|0|88|247"
orf_id = sys.argv[1]

transcript_id = Orf(species, orf_id).transcript_id
transcript_len = int(get_A2B("annotation/" + species + ".transcript.cds_range", 1, 3)[transcript_id])
canonical_start = int(get_A2B("annotation/" + species + ".transcript.cds_range", 1, 4)[transcript_id])
UTR3_len = int(get_A2B("annotation/" + species + ".transcript.cds_range", 1, 6)[transcript_id])
canonical_end = transcript_len - UTR3_len

int_fasta = Fasta("./data/int/" + species + ".annotation.int.from" + str(top_pos) + "top" + str(top_num) + ".fa")
int_ids = list(int_fasta.keys())
dct_int_GT_AG = pickle_to_object("/projects/b1080/hy/simorf/simorf/data/int/" + species + ".dct_int_GT_AG.from" + str(top_pos) + "top" + str(top_num) + ".pickle")

exon_structure_file = "data/exon_structure/" + species + "/" + transcript_id + ".exon_structure"

print("simulation_num", "rep", "fdr_L_all", "fdr_L_main")
for simulation_num in sorted([10**_ for _ in list(range(5))] + [1000*_ for _ in list(range(2,10))]):
    for rep in range(10):
        pep_len, shorter_all, longer_all, total_all, shorter_main, longer_main, total_main, shorter_uorf, longer_uorf, total_uorf, shorter_ouorf, longer_ouorf, total_ouorf = simulation_int(species, orf_id, exon_structure_file, canonical_start, canonical_end,  simulation_num, int_fasta, int_ids, dct_int_GT_AG)

        fdr_L_all = round(longer_all/total_all, 16)
        fdr_L_main = round(longer_main/total_main, 16)
        print(simulation_num, rep, fdr_L_all, fdr_L_main)

