"""

A script to test if more simulation number could generate more accurate FDR by transcript shuffling method

Example script:
python test_simulation_depth_by_transcript_shuffling.py -s human -oi "ENST00000621318.4:chr7:-|39|3084:647:1052|uORF|ATG|1250|2882|647|1052"

Parameter description:
s  = species
oi  = orf_id

"""

__version__ = "1.0.1"
__email__ = "haiwang.yang@northwestern.edu"

from simulate_by_transcript_shuffling import simulation_shu, Orf
from timeit import default_timer as timer
from common_functions import get_A2B
from pyfaidx import Fasta
import argparse

def main():
    # fetching parameters
    parser = argparse.ArgumentParser(description='please provide information')
    parser.add_argument('-s', '--species', type=str) # human
    parser.add_argument('-oi', '--orf_id', type=str) # ENST00000528722.5:chr11:+|2|954:7:88|uORF|CTG|566|824|7|88
    args = parser.parse_args()
    species = args.species
    orf_id = args.orf_id
    
    transcript_id = Orf(species, orf_id).transcript_id
    transcript_len = int(get_A2B("annotation/" + species + ".transcript.cds_range", 1, 3)[transcript_id])
    canonical_start = int(get_A2B("annotation/" + species + ".transcript.cds_range", 1, 4)[transcript_id])
    UTR3_len = int(get_A2B("annotation/" + species + ".transcript.cds_range", 1, 6)[transcript_id])
    canonical_end = transcript_len - UTR3_len
    fasta = Fasta("annotation/" + species + ".transcript.fa")
    transcript_seq = str(fasta[transcript_id])
    for i in list(range(10,100,1)) + list(range(100,1000,10)) + list(range(1000,10001,100)):
        #pep_len, shorter0, longer0, total0, shorter1, longer1, total1, simulations0, simulations1 = simulation_shu("human", "ENST00000528722.5:chr11:+|2|954:7:88|uORF|CTG|566|824|7|88", "AAGTGGCTGGAGTCGGGAGGCTGGAAAGAGACTCCGAGAAAGTACCAGCGGAAGGCGGCCGCCGCTACGGCGATTCGCAGGGAGTAGCAGACGAAGACGGTGGCCGCCGCACTAGCCACCACGTGTGGAGGATAAACGGTCTACACGGCCATTCCGGCGCCGAGTCTAGGGAAAGAGTTAGCGACGACGGGGAAAGAAAATGTGAAGAGAGCGACCGCCGCTCCAGGGTCGCTGCAGGAAGCCTAAGTGCAGACGCCGGCTTCTCCCGCAGTGACTTGAGAAGGGTTCCAGGAAGAAGGTGCATTTTGGCAGCATACATGATGCAGTACGAGCTGGAGATGTAAAGCAGCTTTCAGAAATAGTGGTACGTGGAGCCAGCATTAATGAACTTGATGTTCTCCATAAGTTTACCCCTTTACATTGGGCAGCACATTCTGGAAGTTTGGAGTGTCTTCATTGGCTGCTCTGGCATGGAGCTGATATCACACACGTAACAACGAGAGGTTGGACAGCATCTCACATAGCTGCAATCAGGGGTCAGGATGCTTGTGTACAGGCTCTTATAATGAATGGAGCAAATCTGACAGCCCAGGATGACCGGGGATGCACTCCTTTACATCTTGCTGCAACTCATGGACATTCTTTCACTTTACAAATAATGCTCCGAAGTGGAGTGGATCCCAGTGTGACTGATAAGAGAGAATGGAGACCTGTGCATTATGCAGCTTTTCATGGGCGGCTTGGCTGCTTGCAACTTCTTGTTAAATGGGGTTGTAGCATAGAAGATGTGGACTACAATGGAAACCTTCCAGAACCTCCTTAGATCCCTGTGGAGCCTCTGATTCCTTGGCTTGGCATGCTGGACTTACATAATTTGGCACCTACCTACCTGCATGGCCTCATCTTTCTATATCTTTACACTACTGTCCAGGATCATTTTGTTTTTTTCTGAAG", 565, 823, i)
        pep_len, shorter0, longer0, total0, shorter1, longer1, total1, simulations0, simulations1 = simulation_shu(species, orf_id, transcript_seq, canonical_start, canonical_end, i)
        fdr_l_0 = round(longer0/total0, 3)
        fdr_l_1 = round(longer1/total1, 3)
        print(i, fdr_l_0, fdr_l_1)


if __name__ == "__main__":
    start = timer()
    main()
    end = timer()
    print("time used: ", round(end - start, 5), " sec")    
