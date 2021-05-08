"""

A script to test if a given ORF is shorter or longer than expected by intergenic region fetching method

Example script:
python simulate_by_random_UTR3.py -s human -oi "ENST00000370418.7:chr10:-|3|1863:63:90|uORF|ATG|253|1630|63|90" -esf data/exon_structure/human/ENST00000370418.7.exon_structure -cs 252 -ce 1629 -sn 1000 -tp 1 -tn 1000

Parameter description:
s  = species
oi  = orf_id
ts = transcript_seq
cs = canonical_start
ce = canonical_end
sn = simulation_num
esf= exon structure file path

"""

__version__ = "1.0.1"
__email__ = "haiwang.yang@northwestern.edu"

from timeit import default_timer as timer
import argparse
import random
import os
import re
import pickle
from Bio.Seq import Seq
from pyfaidx import Fasta

def get_revcom_dna(dna):
    return(str(Seq(dna).reverse_complement()))

def object_to_pickle(obj, pkl_path):
    with open(pkl_path, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

def pickle_to_object(pkl_path):
    with open(pkl_path, "rb") as input:
        return(pickle.load(input))

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




class Transcript2:
    def __init__(self, species, transcript_id, exon_structure, canonical_start, canonical_end):
        self.species = species
        self.transcript_id = transcript_id
        self.exon_structure = exon_structure
        # below will slow down the script, is just for testing purpose
        #self.transcript_seq = self.get_transcript_seq()
        #self.cds_seq = self.get_cds_seq()
  
    def get_shuffled_seq(self):
        lst = list(self.transcript_seq)
        random.shuffle(lst)
        return "".join(lst)

    def get_transcript_seq(self):
        transcript_seq = ''
        for _, iexon_seq, __ in self.exon_structure:
            transcript_seq += iexon_seq
        return transcript_seq

    def get_cds_seq(self):
        return self.transcript_seq[canonical_start:canonical_end]

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

def fetch_exon_structure(exon_structure_file):
    dct = {}
    with open(exon_structure_file, "r") as r:
        for line in r.readlines():
            transcript_id, GT, iexon_seq, AG = line.rstrip().split("\t")
            if not transcript_id in dct:
                dct[transcript_id] = []
            dct[transcript_id].append([GT, iexon_seq, AG])
    return dct

def get_int_simulation(exons, int_ids, int_fasta, dct_int_GT_AG):
    if exons == []:
        pass
    else:
        whole_SIM_seq = ""
        for left, useq, right in exons: # for each exon
            L = len(useq)
            HL = L//2
            if left == "AG" and right == "GT":
                lh_useq =  useq[:HL]   # AG-left-half
                lh_L = HL
                rh_useq =  useq[HL:]   # right-half-GT
                rh_L = L - HL

                rGTloc, rAGloc = "", ""
                while rGTloc == "" or rAGloc == "":
                    rid = random.choice(int_ids) # random trans id
                    riseq = str(int_fasta[rid])
                    riL = len(riseq)

                    GTlocs = dct_int_GT_AG[rid]['GT']
                    AGlocs = dct_int_GT_AG[rid]['AG']
                    if len(GTlocs) > 0 and len(AGlocs) > 0:
                        rGTloc = random.choice(GTlocs)      # random GT loc
                        rAGloc = random.choice(AGlocs)      # random AG loc
                        if rGTloc < rh_L: # rGTloc is not enough
                            rGTloc = ""
                        if lh_L > (riL - rAGloc - 2): # if rAGloc is not enough
                            rAGloc = ""

                GTseq = riseq[(rGTloc-rh_L):rGTloc]
                AGseq = riseq[(rAGloc + 2):(rAGloc + lh_L + 2)]
                SIM_iexonseq = AGseq + GTseq
            elif left == "AG" and right == "-":
                rAGloc = ""
                while rAGloc == "":
                    rid = random.choice(int_ids) # random trans id
                    riseq = str(int_fasta[rid])
                    riL = len(riseq)

                    AGlocs = dct_int_GT_AG[rid]['AG']
                    if len(AGlocs) > 0:
                        rAGloc = random.choice(AGlocs)      # random AG loc
                        if L > (riL - rAGloc - 2): # if rAGloc is not enough
                            rAGloc = ""

                SIM_iexonseq = riseq[(rAGloc + 2):(rAGloc + L + 2)]
            elif left == "-" and right == "GT":
                rGTloc = ""
                while rGTloc == "":
                    rid = random.choice(int_ids) # random trans id
                    riseq = str(int_fasta[rid])
                    riL = len(riseq)

                    GTlocs = dct_int_GT_AG[rid]['GT']
                    if len(GTlocs) > 0:
                        rGTloc = random.choice(GTlocs)      # random GT loc
                        if rGTloc < L: # rGTloc is not enough
                            rGTloc = ""

                SIM_iexonseq = riseq[(rGTloc-L):rGTloc]
            else: # non-splicing exon use random int
                ranloc = ""
                while ranloc == "":
                    rid = random.choice(int_ids)
                    riseq = str(int_fasta[rid])
                    riL = len(riseq)

                    ranloc = random.choice(list(range(riL)))
                    if L > (riL - ranloc): # if intergnic is not long enough
                        ranloc = ""

                SIM_iexonseq = riseq[(ranloc):(ranloc + L)]

            whole_SIM_seq += SIM_iexonseq # concatenate simulated iexon to whole
    return whole_SIM_seq

def simulation_int(species, orf_id, exon_structure_file, canonical_start, canonical_end,  simulation_num, int_fasta, int_ids, dct_int_GT_AG):
    orf = Orf(species, orf_id)
    dct_exon_structure = fetch_exon_structure(exon_structure_file)
    exon_structure = dct_exon_structure[orf.transcript_id]

    transcript =Transcript2(species, orf.transcript_id, exon_structure, canonical_start, canonical_end)
    # distributions of three backgrounds (longest main orfs; all uorfs; all overlapping uorfs)

    lst_sim_all_orf_pep_len = []          # all ORFs in simulated transcript
    lst_sim_main_orf_pep_len = []         # main longest ORFs in simulated transcript
    lst_sim_uorf_pep_len = []             # ORFs in simulated 5'UTR (based on longest ATG-ORFs in simulated transcript)
    lst_sim_overlapping_uorf_pep_len = [] # span original canonical start but in simulated transcript

    for sim in range(simulation_num):
        simulated_transcript_seq = get_int_simulation(exon_structure, int_ids, int_fasta, dct_int_GT_AG)
        lst0 = sorted(obtain_all_orfs_in_string_with_perticular_start_codon(simulated_transcript_seq, orf.start_codon), key=lambda x: x[0] - x[1])
        if len(lst0) > 0:
            # all
            for _ in lst0:
                lst_sim_all_orf_pep_len.append((_[1] - _[0] - 3) // 3)

            # longest
            lst_sim_main_orf_pep_len.append((lst0[0][1] - lst0[0][0] - 3) // 3)
        else:
            lst_sim_all_orf_pep_len.append(0)
            lst_sim_main_orf_pep_len.append(0)

        if orf.start_codon != "ATG": # if ORF have a different start_codon than ATG
            lst = sorted(obtain_all_orfs_in_string_with_perticular_start_codon(simulated_transcript_seq, "ATG"), key=lambda x: x[0] - x[1])
        else:
            lst = lst0

        if len(lst) > 0:
            sim_main_orf_start, sim_main_orf_end = lst[0][0], lst[0][1]
            #sim_main_orf_pep_len = (sim_main_orf_end - sim_main_orf_start - 3) // 3
            #lst_sim_main_orf_pep_len.append(sim_main_orf_pep_len)

            sim_5UTR_seq = simulated_transcript_seq[:sim_main_orf_start]
            sim_main_orf_seq = simulated_transcript_seq[sim_main_orf_start:sim_main_orf_end]
            sim_3UTR_seq = simulated_transcript_seq[sim_main_orf_end:]

            if_overlapping_uorf = False
            for s, e in lst0:
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

    return orf.pep_len, shorter_all, longer_all, total_all, shorter_main, longer_main, total_main, shorter_uorf, longer_uorf, total_uorf, shorter_ouorf, longer_ouorf, total_ouorf
    #return orf.pep_len, shorter0, longer0, total0, shorter1, longer1, total1, shorter2, longer2, total2, ",".join([str(_) for _ in lst_sim_main_orf_pep_len]), ",".join([str(_) for _ in lst_sim_uorf_pep_len]), ",".join([str(_) for _ in lst_sim_overlapping_uorf_pep_len])

def example_int(species, orf_id, exon_structure_file, canonical_start, canonical_end,  simulation_num, int_fasta, int_ids, dct_int_GT_AG):
    orf = Orf(species, orf_id)
    dct_exon_structure = fetch_exon_structure(exon_structure_file)
    exon_structure = dct_exon_structure[orf.transcript_id]

    transcript =Transcript2(species, orf.transcript_id, exon_structure, canonical_start, canonical_end)
    # distributions of three backgrounds (longest main orfs; all uorfs; all overlapping uorfs)

    lst_sim_all_orf_pep_len = []          # all ORFs in simulated transcript
    lst_sim_main_orf_pep_len = []         # main longest ORFs in simulated transcript
    lst_sim_uorf_pep_len = []             # ORFs in simulated 5'UTR (based on longest ATG-ORFs in simulated transcript)
    lst_sim_overlapping_uorf_pep_len = [] # span original canonical start but in simulated transcript

    lst_sim_all_orf_cds = []          # all ORFs in simulated transcript
    lst_sim_main_orf_cds = []         # main longest ORFs in simulated transcript
    lst_sim_uorf_cds = []             # ORFs in simulated 5'UTR (based on longest ATG-ORFs in simulated transcript)
    lst_sim_overlapping_uorf_cds = [] # span original canonical start but in simulated transcript

    for sim in range(simulation_num):
        simulated_transcript_seq = get_int_simulation(exon_structure, int_ids, int_fasta, dct_int_GT_AG)
        lst0 = sorted(obtain_all_orfs_in_string_with_perticular_start_codon(simulated_transcript_seq, orf.start_codon), key=lambda x: x[0] - x[1])
        if len(lst0) > 0:
            # all
            for _ in lst0:
                lst_sim_all_orf_pep_len.append((_[1] - _[0] - 3) // 3)
                lst_sim_all_orf_cds.append(simulated_transcript_seq[_[0]:_[1]])

            # longest
            lst_sim_main_orf_pep_len.append((lst0[0][1] - lst0[0][0] - 3) // 3)
            lst_sim_main_orf_cds.append(simulated_transcript_seq[lst0[0][0]:lst0[0][1]])
        else:
            lst_sim_all_orf_pep_len.append(0)
            lst_sim_main_orf_pep_len.append(0)

        if orf.start_codon != "ATG": # if ORF have a different start_codon than ATG
            lst = sorted(obtain_all_orfs_in_string_with_perticular_start_codon(simulated_transcript_seq, "ATG"), key=lambda x: x[0] - x[1])
        else:
            lst = lst0

        if len(lst) > 0:
            sim_main_orf_start, sim_main_orf_end = lst[0][0], lst[0][1]
            #sim_main_orf_pep_len = (sim_main_orf_end - sim_main_orf_start - 3) // 3
            #lst_sim_main_orf_pep_len.append(sim_main_orf_pep_len)
            

            sim_5UTR_seq = simulated_transcript_seq[:sim_main_orf_start]
            sim_main_orf_seq = simulated_transcript_seq[sim_main_orf_start:sim_main_orf_end]
            sim_3UTR_seq = simulated_transcript_seq[sim_main_orf_end:]

            if_overlapping_uorf = False
            for s, e in lst0:
                if s < canonical_start < e:
                    if_overlapping_uorf = True
                    lst_sim_overlapping_uorf_pep_len.append((e - s - 3)//3)
                    lst_sim_overlapping_uorf_cds.append(simulated_transcript_seq[s:e])
            if not if_overlapping_uorf:
                lst_sim_overlapping_uorf_pep_len.append(0)

            lst2 = obtain_all_orfs_in_string_with_perticular_start_codon(sim_5UTR_seq, orf.start_codon)
            if len(lst2) > 0:
                for s, e in lst2:
                    lst_sim_uorf_pep_len.append((e - s - 3)//3)
                    lst_sim_uorf_cds.append(sim_5UTR_seq[s:e])
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

    return lst_sim_all_orf_cds, lst_sim_main_orf_cds, lst_sim_uorf_cds, lst_sim_overlapping_uorf_cds, lst_sim_all_orf_pep_len, lst_sim_main_orf_pep_len, lst_sim_uorf_pep_len, lst_sim_overlapping_uorf_pep_len

def get_GT_AG_dct(fasta):
    dct_GT_AG = {}
    for _ in fasta.keys():
        dct_GT_AG[_] = dict()
        dct_GT_AG[_]['GT'] = []
        dct_GT_AG[_]['AG'] = []
        for m in re.finditer(r"GT", str(fasta[_])):
            dct_GT_AG[_]['GT'].append(m.start())

        for m in re.finditer(r"AG", str(fasta[_])):
            dct_GT_AG[_]['AG'].append(m.start())
    return dct_GT_AG


def generate_merged_UTR3(species, top_pos, top_num):
    """ generate a small set of intergenic and intron mixed data to boost simulation speed """
    half_top_num = top_num//2
    UTR3_fasta = Fasta("./data/UTR3/" + species + ".annotation.UTR3.fa")
    UTR3_ids = list(UTR3_fasta.keys())

    longest_UTR3_ids = sorted(UTR3_ids, key=lambda x: -len(UTR3_fasta[x]))[top_pos:(top_pos+half_top_num)]

    with open("./data/UTR3/" + species + ".annotation.UTR3.from" + str(top_pos) + "top" + str(top_num) + ".fa", "w") as w:
        for UTR3_id in longest_UTR3_ids:
            w.write(">" + UTR3_id + "\n" + str(UTR3_fasta[UTR3_id]) + "\n")

    dct_UTR3_GT_AG = get_GT_AG_dct(UTR3_fasta)
    object_to_pickle(dct_UTR3_GT_AG, "/projects/b1080/hy/simorf/simorf/data/UTR3/" + species + ".dct_UTR3_GT_AG.pickle")

    dct_longest_UTR3_GT_AG = {}
    for _ in longest_UTR3_ids:
        dct_longest_UTR3_GT_AG[_] = {}
        dct_longest_UTR3_GT_AG[_]["GT"] = dct_UTR3_GT_AG[_]["GT"]
        dct_longest_UTR3_GT_AG[_]["AG"] = dct_UTR3_GT_AG[_]["AG"]
    object_to_pickle(dct_longest_UTR3_GT_AG, "/projects/b1080/hy/simorf/simorf/data/UTR3/" + species + ".dct_UTR3_GT_AG.from" + str(top_pos) + "top" + str(top_num) + ".pickle")

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
    parser.add_argument('-esf', '--exon_structure_file', type=str) # exon structure file with context
    parser.add_argument('-tp', '--top_pos', type=int) # pick the position of longest intergenic/intron to boost speed
    parser.add_argument('-tn', '--top_num', type=int) # pick the number of longest intergenic/intron to boost speed

    args = parser.parse_args()
    species = args.species
    orf_id = args.orf_id
    simulation_num = args.simulation_num
    transcript_seq = args.transcript_seq
    canonical_start = args.canonical_start
    canonical_end = args.canonical_end
    exon_structure_file = args.exon_structure_file
    top_pos = args.top_pos
    top_num = args.top_num

    orf = Orf(species, orf_id)
    
    
    dct_exon_structure = fetch_exon_structure(exon_structure_file)
    exon_structure = dct_exon_structure[orf.transcript_id]

    if not os.path.exists("/projects/b1080/hy/simorf/simorf/data/UTR3/" + species + ".dct_int_GT_AG.from" + str(top_pos) + "top" + str(top_num) + ".pickle"):    
        generate_merged_UTR3(species, top_pos, top_num) 
    
    #int_fasta = Fasta("./data/int/" + species + ".annotation.int.from" + str(top_pos) + "top" + str(top_num) + ".fa")
    #int_ids = list(int_fasta.keys())

    #dct_int_GT_AG = pickle_to_object("/projects/b1080/hy/simorf/simorf/data/int/" + species + ".dct_int_GT_AG.from" + str(top_pos) + "top" + str(top_num) + ".pickle")

    #transcript =Transcript2(species, orf.transcript_id, exon_structure, canonical_start, canonical_end)


    #pep_len, shorter_all, longer_all, total_all, shorter_main, longer_main, total_main, shorter_uorf, longer_uorf, total_uorf, shorter_ouorf, longer_ouorf, total_ouorf = simulation_int(species, orf_id, exon_structure_file, canonical_start, canonical_end,  simulation_num, int_fasta, int_ids, dct_int_GT_AG)    
    #print(pep_len, shorter_all, longer_all, total_all, shorter_main, longer_main, total_main, shorter_uorf, longer_uorf, total_uorf, sh#orter_ouorf, longer_ouorf, total_ouorf)
    #
if __name__ == "__main__":
    start = timer()
    main()
    end = timer()
    print("time used: ", round(end - start, 5), " sec") 
    
    
    
    
    
