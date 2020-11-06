"""

A script to classify orf based on gtf annotation gene_type and look up table
and add basic orf and transcript/gene info


Example script:
python add_basic_info_to_orf_list.py -s human -olf list/human.orf_id.list

Parameter description:
s  = species
olf  = orf list file
"""

__version__ = "1.0.1"
__email__ = "haiwang.yang@northwestern.edu"


from timeit import default_timer as timer
from common_functions import Orf, Transcript2, get_elements, get_A2B, object_to_pickle, pickle_to_object
import argparse
import os
import re

def fetch_id_in_gtf_attribute(string, tag_id):
    if tag_id in string:
        return(string.rstrip().split(tag_id + " \"")[1].split("\"")[0])
    else:
        return("")

if __name__ == "__main__":
    # fetching parameters
    parser = argparse.ArgumentParser(description='please provide information')
    parser.add_argument('-s', '--species', type=str) # FDR tables
    parser.add_argument('-olf', '--orf_list_file', type=str) # orf list file

    args = parser.parse_args()
    species = args.species
    orf_list_file = args.orf_list_file

    gene_type_to_gene_classification = get_A2B("annotation/all_gene_type_to_classfication.txt", 1, 2)
    cds_range_file = "annotation/" + species + ".transcript.cds_range"
    transcript_id_to_UTR5_len = get_A2B(cds_range_file, 1, 4)
    transcript_id_to_CDS_len = get_A2B(cds_range_file, 1, 5)
    transcript_id_to_UTR3_len = get_A2B(cds_range_file, 1, 6)

    transcript_id_to_gene_id = {}
    gene_id_to_gene_classification = {}
    gene_id_to_gene_name = {}
    with open("annotation/" + species + ".annotation.gtf", "r") as r:
        for line in r.readlines():
            if not line.startswith("#"):
                chrom, source, feature, start, end, score, strand, frame, attribute = line.rstrip().split("\t")
                if feature == "gene":
                    gene_id = fetch_id_in_gtf_attribute(attribute, "gene_id")
                    gene_name = fetch_id_in_gtf_attribute(attribute, "gene_name")
                    gene_id_to_gene_name[gene_id] = gene_name
                    if "gene_biotype" in attribute:
                        gene_type = fetch_id_in_gtf_attribute(attribute, "gene_biotype")
                    elif "gene_type" in attribute:
                        gene_type = fetch_id_in_gtf_attribute(attribute, "gene_type")
                    else:
                        gene_type = "NA"

                    gene_id_to_gene_classification[gene_id] = [gene_type, gene_type_to_gene_classification.get(gene_type, "UNKNOWN")]

                if feature == "exon":
                    gene_id = fetch_id_in_gtf_attribute(attribute, "gene_id")
                    transcript_id = fetch_id_in_gtf_attribute(attribute, "transcript_id")
                    transcript_id_to_gene_id[transcript_id] = gene_id
                    
    orf_list = get_elements(orf_list_file)
    with open(orf_list_file.replace(".list", ".with_info.list"), "w") as w:
        w.write("\t".join(["orf_id", "orf_type", "orf_classification", "orf_start_codon", "orf_pep_len", "transcript_id", "transcript_len", "UTR5_len", "canonical_len", "UTR3_len", "gene_id", "gene_name", "gene_type", "gene_classification"]) + "\n") 
        for orf_id in orf_list:
            orf = Orf(species, orf_id)
            transcript_id = orf.transcript_id
            UTR5_len = transcript_id_to_UTR5_len[transcript_id]
            CDS_len = transcript_id_to_CDS_len[transcript_id]
            UTR3_len = transcript_id_to_UTR3_len[transcript_id]
            gene_id = transcript_id_to_gene_id[transcript_id]
            gene_name = gene_id_to_gene_name[gene_id]
            gene_type, gene_classification = gene_id_to_gene_classification[gene_id]
            orf_classification = "canonical"
            if "uORF" in orf.orf_type:
                orf_classification = orf.orf_type
            elif gene_classification == "lncRNA":
                orf_classification = "lncRNA"
            elif gene_classification == "pseudogene":
                orf_classification = "pseudogene"
            lst_to_print = [orf_id, orf.orf_type, orf_classification, orf.start_codon, orf.pep_len, transcript_id, orf.transcript_len, UTR5_len, CDS_len, UTR3_len, gene_id, gene_name, gene_type, gene_classification]
            w.write("\t".join([str(_) for _ in lst_to_print]) + "\n")
