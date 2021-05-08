"""

A script to test if a given ORF is shorter or longer than expected by intergenic region fetching method

Example script:
python -i simulate_mimic_annotation.py -s human -sc ATG -tp 100 -tn 10000 -nc intron
python -i simulate_mimic_annotation.py -s human -sc ATG -tp 100 -tn 1000 -nc intergenic
python -i simulate_mimic_annotation.py -s human -sc ATG -tp 1 -tn 1000 -nc UTR3

Parameter description:
s  = species
sc = stop_codon
tp = top_pos
tn = top_num

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

def cds2pep(cds):
    return(str(Seq(cds).translate()))

def get_revcom_dna(dna):
    return(str(Seq(dna).reverse_complement()))

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


###################
### main script ###
###################
def main():
    # fetching parameters
    parser = argparse.ArgumentParser(description='please provide information')
    parser.add_argument('-s', '--species', type=str) # human
    parser.add_argument('-sc', '--start_codon', type=str) # ATG
    parser.add_argument('-esf', '--exon_structure_file', type=str) # exon structure file with context
    parser.add_argument('-tp', '--top_pos', type=int) # pick the position of longest intergenic/intron to boost speed
    parser.add_argument('-tn', '--top_num', type=int) # pick the number of longest intergenic/intron to boost speed
    parser.add_argument('-nc', '--noncoding', type=str) # noncoding type: intron or intergenic

    args = parser.parse_args()
    species = args.species
    start_codon = args.start_codon
    top_pos = args.top_pos
    top_num = args.top_num
    noncoding = args.noncoding

    if noncoding == "intergenic":
        folder = "int"
        noncoding_tag = "intg_"
    elif noncoding == "intron":
        folder = "int"
        noncoding_tag = "ENST"
    elif noncoding == "UTR3":
        folder = "UTR3"
        noncoding_tag = "ENST"

    genome_fasta = Fasta("annotation/" + species + ".fa")
    int_fasta = Fasta("./data/" + folder + "/" + species + ".annotation." + folder + ".from" + str(top_pos) + "top" + str(top_num) + ".fa")
    #int_fasta = Fasta("./data/" + noncoding + "/" + species + ".annotation." + noncoding + ".fa")
    int_ids = [_ for _ in list(int_fasta.keys()) if noncoding_tag in _ if "chr" in _ if "-" in _]
    with open("mim/" + species + "." + start_codon + ".mim_" + noncoding + "_orf.gtf", "w") as w0, open("mim/" + species + "." + start_codon + ".mim_" + noncoding + "_orf.cds.fa", "w") as w1, open("mim/" + species + "." + start_codon + ".mim_" + noncoding + "_orf.pep.fa", "w") as w2:
        for rid in int_ids:
            print(rid)
            #rid = random.choice(int_ids) # random trans id
            if noncoding == "intergenic":
                _, __ = rid.split("::")
                chrom, other = __.split(":")
                s0, e0 = other.split("(")[0].split("-")
                strand = other.split("(")[1].replace(")", "")
            else:
                transcript_id, chrom, s0, e0, strand_raw = rid.split("__")
                strand = strand_raw.split("(")[0]
            s0, e0 = int(s0), int(e0)
            #print(chrom, s0, e0)
            riseq  = str(int_fasta[rid])
            riL = len(riseq)
            lst = sorted(obtain_all_orfs_in_string_with_perticular_start_codon(riseq, start_codon), key=lambda x: x[0] - x[1])
            for s1, e1 in lst:
                orf_cds = riseq[s1:e1]
                orf_pep = cds2pep(orf_cds)
                if strand == "+":
                    #testseq = str(genome_fasta[chrom])[(s0 + s1):(s0 + e1)]
                    start = str(s0 + s1 + 1)
                    end = str(s0 + e1)
                    MIM_id = 'MIM' + '_' + chrom + '_' + strand + '_' + start + '_' + end
                    attribute = 'gene_id "' + MIM_id + '"; transcript_id "' + MIM_id + '"; gene_name "' + MIM_id + '"'
                    w0.write("\t".join([chrom, "ME", "CDS", start, end, ".", strand, "0", attribute]) + "\n")
                else:
                    #testseq = str(genome_fasta[chrom])[(e0 - e1):(e0 - s1)]
                    start = str(e0 - e1)
                    end = str(e0 - s1)
                    MIM_id = 'MIM' + '_' + chrom + '_' + strand + '_' + start + '_' + end
                    attribute = 'gene_id "' + MIM_id + '"; transcript_id "' + MIM_id + '"; gene_name "' + MIM_id + '"'
                    w0.write("\t".join([chrom, "ME", "CDS", start, end, ".", strand, "0", attribute]) + "\n")
                w1.write(">" +  MIM_id + "\n" + orf_cds + "\n")
                w2.write(">" +  MIM_id + "\n" + orf_pep + "\n")

if __name__ == "__main__":
    start = timer()
    main()
    end = timer()
    print("time used: ", round(end - start, 5), " sec") 
        
        
        
        
        
