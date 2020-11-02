__version__ = "1.0.1"
__email__ = "haiwang.yang@northwestern.edu"

import re
import pickle
from Bio.Seq import Seq


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

def get_elements(filename):
    """ get elements from file, each line is an element """
    elements = list()
    with open(filename, 'r') as f:
        lines = f.readlines()
        return [e.rstrip() for e in lines]

def get_A2B(filename, colKey, colValue):
    """ get dict of A2B from tab-separated table """
    dct = dict()
    with open(filename, 'r') as f:
        for line in f.readlines():
            temp = line.rstrip().split("\t")
            dct[temp[colKey - 1]] = temp[colValue - 1]
    return(dct)
