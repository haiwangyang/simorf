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
        self.orf_start -= 1 # convert to 0 based
        self.orf_end -= 1   # convert to 0 based
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

def get_A2Last(filename, colKey):
    """ get dict of A2Last from tab-separated table """
    dct = dict()
    with open(filename, 'r') as f:
        for line in f.readlines():
            temp = line.rstrip().split("\t")
            dct[temp[colKey - 1]] = temp[-1]
    return(dct)

def get_intersection_union_jaccard(exonmap1, exonmap2):
    """ get jaccard between exon map
        1.3 and 2.4
        jaccard = 2/4 = 0.5
    """
    union_sum = 0
    intersection_sum = 0

    dct1 = dict()
    for se in exonmap1:
        s, e = se
        for i in range(int(s), int(e) + 1):
            if not i in dct1.keys():
                dct1[i] = 0
            dct1[i] += 1

    dct2 = dict()
    for se in exonmap2:
        s, e = se
        for i in range(int(s), int(e) + 1):
            if not i in dct2.keys():
                dct2[i] = 0
            dct2[i] += 1

    st = set()
    for ii in [dct1.keys(), dct2.keys()]:
        for i in ii:
            st.add(i)

    union_sum = len(st)
    for i in st:
        if i in dct1.keys() and i in dct2.keys():
            intersection_sum += 1

    j = intersection_sum / union_sum
    return(intersection_sum, union_sum, j)
