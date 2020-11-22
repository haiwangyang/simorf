"""

A script to classify human orf based on gtf annotation gene_type and look up table
and add basic orf and transcript/gene info


Example script:
python add_info_to_human_orf_list.py -s human -olf list/human.orf_id.list -iff output/human.intFDRs.sim2000.txt -sff output/human.shuFDRs.sim2000.txt -o output/human.FDRs.with_info.txt

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
    parser.add_argument('-iff', '--int_fdr_file', type=str) # int fdr files output/human.intFDRs.sim1000.merged.txt
    parser.add_argument('-sff', '--shu_fdr_file', type=str) # shu fdr files output/human.shuFDRs.sim1000.merged.txt
    parser.add_argument('-o', '--output', type=str) # output

    args = parser.parse_args()
    species = args.species
    orf_list_file = args.orf_list_file
    int_fdr_file = args.int_fdr_file
    shu_fdr_file = args.shu_fdr_file
    output = args.output

    gene_type_to_gene_classification = get_A2B("annotation/all_gene_type_to_classfication.txt", 1, 2)
    cds_range_file = "annotation/" + species + ".transcript.cds_range"
    transcript_id_to_UTR5_len = get_A2B(cds_range_file, 1, 4)
    transcript_id_to_CDS_len = get_A2B(cds_range_file, 1, 5)
    transcript_id_to_UTR3_len = get_A2B(cds_range_file, 1, 6)

    orf_id_to_int_FDR_L_all = get_A2B(int_fdr_file,  1, 5)
    orf_id_to_int_FDR_L_main = get_A2B(int_fdr_file,  1, 7)
    orf_id_to_int_FDR_L_uORF = get_A2B(int_fdr_file,  1, 9)
    orf_id_to_int_FDR_L_ouORF = get_A2B(int_fdr_file,  1, 11)

    orf_id_to_shu_FDR_L_all = get_A2B(shu_fdr_file,  1, 5)
    orf_id_to_shu_FDR_L_main = get_A2B(shu_fdr_file,  1, 7)
    orf_id_to_shu_FDR_L_uORF = get_A2B(shu_fdr_file,  1, 9)
    orf_id_to_shu_FDR_L_ouORF = get_A2B(shu_fdr_file,  1, 11)

    orf_id_to_phylocsf = get_A2B("features/human.phylocsf.txt", 1, 2)
    orf_id_to_optimizedcodon = get_A2B("features/human.optimizedcodon.txt", 1, 4)
    orf_id_to_r_hydrophobic = get_A2B("features/human.orf.aa_feature", 1, 3)
    orf_id_to_r_amphipathic = get_A2B("features/human.orf.aa_feature", 1, 4)
    orf_id_to_r_polar = get_A2B("features/human.orf.aa_feature", 1, 5)
    orf_id_to_r_charged = get_A2B("features/human.orf.aa_feature", 1, 6)
    orf_id_to_kozak_score = get_A2B("features/human.orf.kozak_feature", 1, 2)
    orf_id_to_orf_upstream_MFE = get_A2B("features/human.orf.MFE_feature", 1, 2)
    orf_id_to_orf_MFE = get_A2B("features/human.orf.MFE_feature", 1, 3)
    expression_columns = ["A549_GSE101760_SRX3028093","Muscle_GSE103308_SRX3146514","HeLa_GSE105082_SRX3292091","Fibroblast_GSE42509_SRX207945","Brain_GSE51424_SRX690810","U2OS_GSE56924_SRX522481","HCT116_GSE58207_SRX569214","HEK293T_GSE59095_SRX646183","Kidney_GSE59820_SRX663300","ESC_GSE62247_SRX730865","HeLa_GSE63591_SRX767133","Breast_GSE69923_SRX1059911","Erythroid_GSE89183_SRX2268404","Huh7_GSE94454_SRX2536403","tau"]
    dct_expression = {}
    for _, colname in enumerate(expression_columns):
        dct_expression[colname] = get_A2B("features/human.orf.expression_feature", 1, _ + 2)

    orf_id_to_westernblot_intensity = get_A2B("features/human.orf.westernblot_feature", 1, 2)

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
    with open(output, "w") as w:
        w.write("\t".join(["orf_id", "orf_type", "orf_classification", "orf_upstream_len", "orf_start_codon", "orf_pep_len", "westernblot_intensity", "phylocsf", "optimizedcodon", "r_hydrophobic", "r_amphipathic", "r_polar", "r_charged", "kozak_score", "orf_upstream_MFE", "orf_MFE", "int_FDR_L_all", "int_FDR_L_main", "int_FDR_L_uORF", "int_FDR_L_ouORF", "shu_FDR_L_all", "shu_FDR_L_main", "shu_FDR_L_uORF", "shu_FDR_L_ouORF", "transcript_id", "transcript_len", "UTR5_len", "canonical_len", "UTR3_len", "gene_id", "gene_name", "gene_type", "gene_classification"] + expression_columns) + "\n") 
        for orf_id in orf_list:
            orf = Orf(species, orf_id)
            int_FDR_L_all = orf_id_to_int_FDR_L_all[orf_id]
            int_FDR_L_main = orf_id_to_int_FDR_L_main[orf_id]
            int_FDR_L_uORF = orf_id_to_int_FDR_L_uORF[orf_id]
            int_FDR_L_ouORF = orf_id_to_int_FDR_L_ouORF[orf_id]

            shu_FDR_L_all = orf_id_to_shu_FDR_L_all[orf_id]
            shu_FDR_L_main = orf_id_to_shu_FDR_L_main[orf_id]
            shu_FDR_L_uORF = orf_id_to_shu_FDR_L_uORF[orf_id]
            shu_FDR_L_ouORF = orf_id_to_shu_FDR_L_ouORF[orf_id]


            westernblot_intensity = orf_id_to_westernblot_intensity.get(orf_id, "NA")
            phylocsf = orf_id_to_phylocsf.get(orf_id, "-10")
            optimizedcodon = orf_id_to_optimizedcodon.get(orf_id, "0")
            r_hydrophobic = orf_id_to_r_hydrophobic[orf_id]
            r_amphipathic = orf_id_to_r_amphipathic[orf_id]
            r_polar = orf_id_to_r_polar[orf_id]
            r_charged = orf_id_to_r_charged[orf_id]

            kozak_score = orf_id_to_kozak_score[orf_id]
            orf_upstream_MFE = orf_id_to_orf_upstream_MFE[orf_id]
            orf_MFE = orf_id_to_orf_MFE[orf_id]

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
            lst_to_print = [orf_id, orf.orf_type, orf_classification, orf.orf_start, orf.start_codon, orf.pep_len, westernblot_intensity, phylocsf, optimizedcodon, r_hydrophobic, r_amphipathic, r_polar, r_charged, kozak_score, orf_upstream_MFE, orf_MFE, int_FDR_L_all, int_FDR_L_main, int_FDR_L_uORF, int_FDR_L_ouORF, shu_FDR_L_all, shu_FDR_L_main, shu_FDR_L_uORF, shu_FDR_L_ouORF, transcript_id, orf.transcript_len, UTR5_len, CDS_len, UTR3_len, gene_id, gene_name, gene_type, gene_classification]
            w.write("\t".join([str(_) for _ in lst_to_print]) + "\t")
            w.write("\t".join([dct_expression[_][orf_id] for _ in expression_columns]) + "\n")
