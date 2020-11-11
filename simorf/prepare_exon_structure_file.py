""" 
A script to generate exon structure file for certain species and transcript

Example script:
python prepare_exon_structure_file.py -s human -ti ENST00000370418.7
python prepare_exon_structure_file.py -s mouse -ti ENSMUST00000183805
python prepare_exon_structure_file.py -s zebrafish -ti ENSDART00000177727
python prepare_exon_structure_file.py -s worm -ti K09C6.2e.2
python prepare_exon_structure_file.py -s arabidopsis -ti AT1G02640.1
python prepare_exon_structure_file.py -s yeast -ti YGR177C

Parameter description:
s  = species
ti = transcript_id
"""

__version__ = "1.0.1"
__email__ = "haiwang.yang@northwestern.edu"

import argparse
from Bio.Seq import Seq
from os import path
from pyfaidx import Fasta
from common_functions import get_revcom_dna, object_to_pickle, pickle_to_object


def generate_transcript_fa(species):
    fasta = Fasta("./annotation/" + species + ".fa")
    transcript_id_to_strand = {}
    transcript_id_to_iexons = {}
    with open("./annotation/" + species + ".annotation.genePred", "r") as r:
        for line in r.readlines():
            transcript_id, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds = line.rstrip().split("\t")
            transcript_id_to_strand[transcript_id] = strand
            if not transcript_id in transcript_id_to_iexons:
                transcript_id_to_iexons[transcript_id] = []
            exonmap = list(zip(exonStarts.rstrip(",").split(","), exonEnds.rstrip(",").split(",")))
            for s, e in exonmap:
                s, e = int(s), int(e)
                transcript_id_to_iexons[transcript_id].append(str(fasta[chrom][s:e]).upper())
    
    with open("./annotation/" + species + ".transcript.fa", "w") as w:
        for transcript_id in transcript_id_to_iexons.keys():
            if transcript_id_to_strand[transcript_id] == "+":
                transcript_seq = ''.join(transcript_id_to_iexons[transcript_id])
            else:
                transcript_seq = ''
                for iexon in transcript_id_to_iexons[transcript_id][::-1]:
                    transcript_seq += get_revcom_dna(iexon)
            w.write(">" + transcript_id + "\n" + transcript_seq + "\n")

def generate_transcript_cds_range(species):
    fasta = Fasta("./annotation/" + species + ".fa")
    with open("./annotation/" + species + ".annotation.genePred", "r") as r, open("./annotation/" + species + ".transcript.cds_range", "w") as w:
        for line in r.readlines():
            transcript_id, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds = line.rstrip().split("\t")
            if cdsStart != cdsEnd:
                cdsStart, cdsEnd = int(cdsStart), int(cdsEnd)
                exonmap = list(zip(exonStarts.rstrip(",").split(","), exonEnds.rstrip(",").split(",")))
                left_utr_size, right_utr_size, transcript_size = 0, 0, 0
                for s, e in exonmap:
                    s, e = int(s), int(e)
                    
                    """ test purpose """
                    #iexon_seq = str(fasta[chrom])[s:e]
                    #_ = 0
                    #for i in range(s, e):
                    #    if i < cdsStart:
                    #        print(str(i) +  iexon_seq[_] + "l ", end="")
                    #    elif cdsEnd <= i:
                    #        print(str(i) +  iexon_seq[_] + "r ", end="")
                    #    else:
                    #        print(str(i) + iexon_seq[_] + "c ", end="")
                    #    _ += 1

                    left_utr_size += sum([1 for i in range(s, e) if i < cdsStart])
                    right_utr_size += sum([1 for i in range(s, e) if cdsEnd <= i])
                    transcript_size += sum([1 for i in range(s, e)])
                    cds_size = transcript_size - left_utr_size - right_utr_size
                if strand == "+":
                    w.write("\t".join([transcript_id, "coding", str(transcript_size), str(left_utr_size), str(cds_size), str(right_utr_size)]) + "\n")
                else:
                    w.write("\t".join([transcript_id, "coding", str(transcript_size), str(right_utr_size), str(cds_size), str(left_utr_size)]) + "\n")
            else: # for non-coding transcripts
                cdsStart, cdsEnd = int(cdsStart), int(cdsEnd)
                exonmap = list(zip(exonStarts.rstrip(",").split(","), exonEnds.rstrip(",").split(",")))
                left_utr_size, right_utr_size, transcript_size = 0, 0, 0
                for s, e in exonmap:
                    s, e = int(s), int(e)
                    transcript_size += sum([1 for i in range(s, e)])
                w.write("\t".join([transcript_id, "noncoding", str(transcript_size), str(0), str(0), str(0)]) + "\n")

def generate_dct_exon_pickle(species):
    dct_exon = dict()
    print("fetching exon length for " + species)

    with open("./annotation/" + species + ".annotation.genePred", "r") as r:
        for line in r.readlines():
            transcript_id, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds = line.rstrip().split("\t")
            #if transcript_id == "ENST00000182527.3":
            if True:
                txStart = int(txStart)
                txEnd = int(txEnd)
                cdsStart = int(cdsStart)
                cdsEnd = int(cdsEnd)
                exonCount = int(exonCount)
                if cdsStart == cdsEnd:
                    coding = False
                else:
                    coding = True

                dct_exon[transcript_id] = []
                if True:
                    dct_exon[transcript_id] = []
                    left_range = [(txStart, cdsStart)]
                    right_range = [(cdsEnd, txEnd)]
                    exonmap = list(zip(exonStarts.rstrip(",").split(","), exonEnds.rstrip(",").split(",")))
                    for exon in exonmap:
                        eL = int(exon[1]) - int(exon[0])
                        if strand == "+":
                            dct_exon[transcript_id].append(eL)
                        elif strand == "-":
                            dct_exon[transcript_id].insert(0, eL)

    object_to_pickle(dct_exon, "./annotation/" + species + ".dct_exon.pickle")

def generate_dct_exon_structure_pickle(dct_exon):
    dct_exon_structure = dict()
    fasta = Fasta("./annotation/" + species + ".transcript.fa")
    for transcript_id, exons_raw in dct_exon.items():
        exons = [_ for _ in exons_raw if not _ == -1]
        dct_exon_structure[transcript_id] = []
        if exons == []:
            pass # no exon pass
        else:
            transcript_seq = str(fasta[transcript_id])
            if len(exons) == 1: # one exon - 2 scenarios
                L = exons[0]
                if len(exons_raw) == 1: # no splicing
                    dct_exon_structure[transcript_id].append(["-", transcript_seq[:L].upper(), "-"])
                else:                       # yes splicing
                    dct_exon_structure[transcript_id].append(["-", transcript_seq[:L].upper(), "GT"])
            else:                  # multiple exon - 3 scenarios
                for i, L in enumerate(exons):
                    this_seq, transcript_seq = transcript_seq[:L], transcript_seq[L:]
                    if i == 0:
                        dct_exon_structure[transcript_id].append(["-", this_seq.upper(), "GT"])
                    elif i == (len(exons) - 1):
                        dct_exon_structure[transcript_id].append(["AG", this_seq.upper(), "-"])
                    else:
                        dct_exon_structure[transcript_id].append(["AG", this_seq.upper(), "GT"])
    object_to_pickle(dct_exon_structure, "./annotation/" + species + ".dct_exon_structure.pickle")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='please provide information')
    parser.add_argument('-s', '--species', type=str) # human
    parser.add_argument('-ti', '--transcript_id', type=str) # python ENST00000370418.7
    args = parser.parse_args()
    species = args.species
    transcript_id = args.transcript_id

    if not path.exists("./annotation/" + species + ".transcript.fa"):
        #print("generating transcript fasta for " + species)
        generate_transcript_fa(species)

    if not path.exists("./annotation/" + species + ".transcript.cds_range"):
        #print("generating transcript cds range for " + species)
        generate_transcript_cds_range(species)
    #print("transcript cds range for " + species + " is ready")
    
    if not path.exists("./annotation/" + species + ".dct_exon.pickle"):
        #print("generating exon pickle for " + species)
        generate_dct_exon_pickle(species)
    #print("loading exon pickle for " + species)
    dct_exon = pickle_to_object("./annotation/" + species + ".dct_exon.pickle")

    if not path.exists("./annotation/" + species + ".dct_exon_structure.pickle"):
        #print("generating exon structure pickle for " + species)
        generate_dct_exon_structure_pickle(dct_exon)
    #print("loading exon structure pickle for " + species)
    dct_exon_structure = pickle_to_object("./annotation/" + species + ".dct_exon_structure.pickle")

    print("writing exon structure file for " + species + " " + transcript_id)
    with open("./data/exon_structure/" + species + "/" + transcript_id + ".exon_structure", "w") as w:
        for GT, iexon_seq, AG in dct_exon_structure[transcript_id]:
            w.write(transcript_id + "\t" + GT + "\t" + iexon_seq + "\t" + AG + "\n")
