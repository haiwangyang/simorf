__version__ = "1.0.1"
__email__ = "haiwang.yang@northwestern.edu"

import re
import pickle
from Bio.Seq import Seq

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
