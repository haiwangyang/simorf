"""

A script to test if a given ORF is shorter or longer than expected by transcript shuffling method

Example script:
python simulate_by_trans_shuffling.py -s human -oi "ENST00000528722.5:chr11:+|2|954:7:88|uORF|CTG|566|824|7|88" -ts AAGTGGCTGGAGTCGGGAGGCTGGAAAGAGACTCCGAGAAAGTACCAGCGGAAGGCGGCCGCCGCTACGGCGATTCGCAGGGAGTAGCAGACGAAGACGGTGGCCGCCGCACTAGCCACCACGTGTGGAGGATAAACGGTCTACACGGCCATTCCGGCGCCGAGTCTAGGGAAAGAGTTAGCGACGACGGGGAAAGAAAATGTGAAGAGAGCGACCGCCGCTCCAGGGTCGCTGCAGGAAGCCTAAGTGCAGACGCCGGCTTCTCCCGCAGTGACTTGAGAAGGGTTCCAGGAAGAAGGTGCATTTTGGCAGCATACATGATGCAGTACGAGCTGGAGATGTAAAGCAGCTTTCAGAAATAGTGGTACGTGGAGCCAGCATTAATGAACTTGATGTTCTCCATAAGTTTACCCCTTTACATTGGGCAGCACATTCTGGAAGTTTGGAGTGTCTTCATTGGCTGCTCTGGCATGGAGCTGATATCACACACGTAACAACGAGAGGTTGGACAGCATCTCACATAGCTGCAATCAGGGGTCAGGATGCTTGTGTACAGGCTCTTATAATGAATGGAGCAAATCTGACAGCCCAGGATGACCGGGGATGCACTCCTTTACATCTTGCTGCAACTCATGGACATTCTTTCACTTTACAAATAATGCTCCGAAGTGGAGTGGATCCCAGTGTGACTGATAAGAGAGAATGGAGACCTGTGCATTATGCAGCTTTTCATGGGCGGCTTGGCTGCTTGCAACTTCTTGTTAAATGGGGTTGTAGCATAGAAGATGTGGACTACAATGGAAACCTTCCAGAACCTCCTTAGATCCCTGTGGAGCCTCTGATTCCTTGGCTTGGCATGCTGGACTTACATAATTTGGCACCTACCTACCTGCATGGCCTCATCTTTCTATATCTTTACACTACTGTCCAGGATCATTTTGTTTTTTTCTGAAG -cs 565 -ce 823 -sn 10000

Parameter description:
s  = species
oi  = orf_id
ts = transcript_seq
cs = canonical_start
ce = canonical_end
sn = simulation_num

"""

__version__ = "1.0.1"
__email__ = "haiwang.yang@northwestern.edu"

from timeit import default_timer as timer
from collections import Counter
import argparse
import random
import re

def findall(string, substring):
    return [m.start() for m in re.finditer(substring, string)]

def obtain_all_orfs_in_string_with_perticular_start_codon(string, start_codon):
    lst = []
    start_positions = findall(string, start_codon)

    TAA_positions = findall(string, "TAA")
    TGA_positions = findall(string, "TGA")
    TAG_positions = findall(string, "TAG")

    start_positions = sorted(start_positions)
    end_positions = sorted(TAA_positions + TGA_positions + TAG_positions)
    for i in start_positions:
        for j in [j for j in end_positions if j > i]:
            d = j - i
            if d % 3 == 0:
                tag = ",".join([str(i + 1), str(j + 4)])
                lst.append([i, j + 3])
                break
    return(lst)

def get_shorter_longer_fdrs(len_pep, lst_sim):
    shorter_than_count = 0
    for i in lst_sim:
        if len_pep >= i:  # >?
            shorter_than_count += 1
        else:
            break

    longer_than_count = 0
    for i in lst_sim[::-1]:
        if i >= len_pep: # >?
            longer_than_count += 1
        else:
            break
    return(shorter_than_count, longer_than_count, len(lst_sim))




class Transcript:
    def __init__(self, species, transcript_id, transcript_seq):
        self.species = species
        self.transcript_id = transcript_id
        self.transcript_seq = transcript_seq

    def get_shuffled_seq(self):
        lst = list(self.transcript_seq)
        random.shuffle(lst)
        return "".join(lst)        

    def get_shuffled_seq2(self):
        lst = list(self.transcript_seq)
        n = len(lst)
        return "".join([lst[random.randrange(n)] for _ in range(n)])
        
class Orf:
    def __init__(self, species, orf_id):
        self.species = species
        self.orf_id = orf_id
        self.parse_orf_id()

    def parse_orf_id(self):
        temp = self.orf_id.split("|")
        self.transcript_id, self.chrom, self.strand = temp[0].split(":")
        self.transcript_len, self.orf_start, self.orf_end = [int(_) for _ in temp[2].split(":")]
        self.orf_type = temp[3]
        self.start_codon = temp[4]
        self.pep_len = (self.orf_end - self.orf_start - 3) // 3 

def simulation_shu(species, orf_id, transcript_seq, canonical_start, canonical_end, simulation_num):
    orf = Orf(species, orf_id)
    transcript =Transcript(species, orf.transcript_id, transcript_seq)

    # distributions of three backgrounds (longest main orfs; all uorfs; all overlapping uorfs)
    lst_sim_all_orf_pep_len = []          # all ORFs in simulated transcript
    lst_sim_main_orf_pep_len = []         # main longest ORFs in simulated transcript
    lst_sim_uorf_pep_len = []             # ORFs in simulated 5'UTR (based on longest ATG-ORFs in simulated transcript)
    lst_sim_overlapping_uorf_pep_len = [] # span original canonical start but in simulated transcript

    for sim in range(simulation_num):
        shuffled_transcript_seq = transcript.get_shuffled_seq2()
        lst0 = sorted(obtain_all_orfs_in_string_with_perticular_start_codon(shuffled_transcript_seq, orf.start_codon), key=lambda x: x[0] - x[1])
        if len(lst0) > 0:
            for _ in lst0:
                lst_sim_all_orf_pep_len.append((_[1] - _[0] - 3) // 3)

            # longest
            lst_sim_main_orf_pep_len.append((lst0[0][1] - lst0[0][0] - 3) // 3)
        else:
            lst_sim_all_orf_pep_len.append(0)
            lst_sim_main_orf_pep_len.append(0)

        if orf.start_codon != "ATG": # if ORF have a different start_codon than ATG
             lst = sorted(obtain_all_orfs_in_string_with_perticular_start_codon(shuffled_transcript_seq, "ATG"), key=lambda x: x[0] - x[1])
        else:
            lst = lst0

        if len(lst) > 0:
            sim_main_orf_start, sim_main_orf_end = lst[0][0], lst[0][1]
            sim_main_orf_pep_len = (sim_main_orf_end - sim_main_orf_start - 3) // 3
            lst_sim_main_orf_pep_len.append(sim_main_orf_pep_len)

            sim_5UTR_seq = shuffled_transcript_seq[:sim_main_orf_start]
            sim_main_orf_seq = shuffled_transcript_seq[sim_main_orf_start:sim_main_orf_end]
            sim_3UTR_seq = shuffled_transcript_seq[sim_main_orf_end:]

            if_overlapping_uorf = False
            for s, e in lst:
                if s < canonical_start < e:
                    if_overlapping_uorf = True
                    lst_sim_overlapping_uorf_pep_len.append((e - s - 3)//3)
            if not if_overlapping_uorf:
                lst_sim_overlapping_uorf_pep_len.append(0)

            lst2 = obtain_all_orfs_in_string_with_perticular_start_codon(sim_5UTR_seq, orf.start_codon)
            if len(lst2) > 0:
                for s, e in lst2:
                    lst_sim_uorf_pep_len.append((e - s - 3)//3)
            else:
                lst_sim_uorf_pep_len.append(0)
        else: # if no orf were found in simulation
            lst_sim_main_orf_pep_len.append(0)
            lst_sim_overlapping_uorf_pep_len.append(0)
            lst_sim_uorf_pep_len.append(0)

    lst_sim_all_orf_pep_len.sort()
    lst_sim_main_orf_pep_len.sort()
    lst_sim_uorf_pep_len.sort()
    lst_sim_overlapping_uorf_pep_len.sort()

    shorter_all, longer_all, total_all = get_shorter_longer_fdrs(orf.pep_len, lst_sim_all_orf_pep_len)
    shorter_main, longer_main, total_main = get_shorter_longer_fdrs(orf.pep_len, lst_sim_main_orf_pep_len)
    shorter_uorf, longer_uorf, total_uorf = get_shorter_longer_fdrs(orf.pep_len, lst_sim_uorf_pep_len)
    shorter_ouorf, longer_ouorf, total_ouorf = get_shorter_longer_fdrs(orf.pep_len, lst_sim_overlapping_uorf_pep_len)

    #print("lst_sim_all_orf_pep_len:", lst_sim_all_orf_pep_len)
    #print("lst_sim_main_orf_pep_len:", lst_sim_main_orf_pep_len)
    #print("lst_sim_uorf_pep_len:", lst_sim_uorf_pep_len)
    #print("lst_sim_overlapping_uorf_pep_len", lst_sim_overlapping_uorf_pep_len)

    return orf.pep_len, shorter_all, longer_all, total_all, shorter_main, longer_main, total_main, shorter_uorf, longer_uorf, total_uorf, shorter_ouorf, longer_ouorf, total_ouorf


###################
### main script ###
###################
def main():
    # fetching parameters
    parser = argparse.ArgumentParser(description='please provide information')
    parser.add_argument('-s', '--species', type=str) # human
    parser.add_argument('-oi', '--orf_id', type=str) # ENST00000528722.5:chr11:+|2|954:7:88|uORF|CTG|566|824|7|88
    parser.add_argument('-ts', '--transcript_seq', type=str) # transcfipt_seq (intron-less)
    parser.add_argument('-cs', '--canonical_start', type=int) # canonical_start
    parser.add_argument('-ce', '--canonical_end', type=int) # canonical_end
    parser.add_argument('-sn', '--simulation_num', type=int) # simulation_num
    
    args = parser.parse_args()
    species = args.species
    orf_id = args.orf_id
    simulation_num = args.simulation_num
    transcript_seq = args.transcript_seq
    canonical_start = args.canonical_start
    canonical_end = args.canonical_end
    
    pep_len, shorter_all, longer_all, total_all, shorter_main, longer_main, total_main, shorter_uorf, longer_uorf, total_uorf, shorter_ouorf, longer_ouorf, total_ouorf = simulation_shu(species, orf_id, transcript_seq, canonical_start, canonical_end, simulation_num)
    print(pep_len, shorter_all, longer_all, total_all, shorter_main, longer_main, total_main, shorter_uorf, longer_uorf, total_uorf
, shorter_ouorf, longer_ouorf, total_ouorf)

if __name__ == "__main__":
    start = timer()
    main()
    end = timer()
    print("time used: ", round(end - start, 5), " sec")
