from simulate_by_intergenic_region import simulation_int0, Orf
from pyfaidx import Fasta
from common_functions import pickle_to_object, get_A2B
from timeit import default_timer as timer

species = "human"
#orf_id = "ENST00000370418.7:chr10:-|3|1863:63:90|uORF|ATG|253|1630|63|90"
orf_id = "ENST00000528722.5:chr11:+|2|954:7:88|uORF|CTG|566|824|7|88"
transcript_id = Orf(species, orf_id).transcript_id
transcript_len = int(get_A2B("annotation/" + species + ".transcript.cds_range", 1, 3)[transcript_id])
canonical_start = int(get_A2B("annotation/" + species + ".transcript.cds_range", 1, 4)[transcript_id])
UTR3_len = int(get_A2B("annotation/" + species + ".transcript.cds_range", 1, 6)[transcript_id])
canonical_end = transcript_len - UTR3_len

intergenic_fasta = Fasta("./data/intergenic/" + species + ".annotation.intergenic.fa")
intergenic_ids = list(intergenic_fasta.keys())
dct_intergenic_GT_AG = pickle_to_object("/projects/b1080/hy/ribo/ribo/uORF_paper/annotation/" + species + ".dct_intergenic_GT_AG.pickle")

for i in list(range(10,100,10)) + list(range(100,1000,100)) + list(range(1000,10001,1000)):
    start = timer()
    pep_len, shorter0, longer0, total0, shorter1, longer1, total1, simulations0, simulations1 = simulation_int0(species, orf_id, "data/exon_structure/human/" + transcript_id + ".exon_structure", canonical_start, canonical_end, i, intergenic_fasta, intergenic_ids, dct_intergenic_GT_AG)
    fdr_l_0 = round(longer0/total0, 3)
    fdr_l_1 = round(longer1/total1, 3)
    end = timer()
    print(i, round(end - start, 5))

