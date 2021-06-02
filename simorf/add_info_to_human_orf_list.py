"""

A script to classify human orf based on gtf annotation gene_type and look up table
and add basic orf and transcript/gene info


Example script:
python add_info_to_human_orf_list.py -s human -olf list/human.orf_id.list.all -iff output/human.intFDRs.sim2000.latest.txt -sff output/human.transFDRs.sim2000.latest.txt -o output/human.FDRs.with_info.txt

python add_info_to_human_orf_list.py -s human -olf list/human.orf_id.list.all -iff output/human.intFDRs.sim2000.extra.txt -sff output/human.transFDRs.sim2000.extra.txt -o output/human.FDRs.with_info.extra.txt

s=human
python add_info_to_human_orf_list.py -s $s -olf list/${s}_int/test_old -iff output/${s}_int/test_old -sff output/${s}_trans/test_old -o output/${s}.FDRs.with_info.test_old.txt

python -i add_info_to_human_orf_list.py -s $s -olf list/${s}_int/test_replace -iff output/${s}_int/test_replace -sff output/${s}_trans/test_replace -o output/${s}.FDRs.with_info.test_replace.txt


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

    orf_id_to_int_FDR_L_all = get_A2B(int_fdr_file,  1, 13)
    orf_id_to_int_FDR_L_all_T = get_A2B(int_fdr_file,  1, 14)
    orf_id_to_int_FDR_L_main = get_A2B(int_fdr_file,  1, 16)
    orf_id_to_int_FDR_L_main_T = get_A2B(int_fdr_file,  1, 17)
    orf_id_to_int_FDR_L_uORF = get_A2B(int_fdr_file,  1, 19)
    orf_id_to_int_FDR_L_uORF_T = get_A2B(int_fdr_file,  1, 20)
    orf_id_to_int_FDR_L_ouORF = get_A2B(int_fdr_file,  1, 22)
    orf_id_to_int_FDR_L_ouORF_T = get_A2B(int_fdr_file,  1, 23)

    orf_id_to_shu_FDR_L_all = get_A2B(shu_fdr_file,  1, 13)
    orf_id_to_shu_FDR_L_all_T = get_A2B(shu_fdr_file,  1, 14)
    orf_id_to_shu_FDR_L_main = get_A2B(shu_fdr_file,  1, 16)
    orf_id_to_shu_FDR_L_main_T = get_A2B(shu_fdr_file,  1, 17)
    orf_id_to_shu_FDR_L_uORF = get_A2B(shu_fdr_file,  1, 19)
    orf_id_to_shu_FDR_L_uORF_T = get_A2B(shu_fdr_file,  1, 20)
    orf_id_to_shu_FDR_L_ouORF = get_A2B(shu_fdr_file,  1, 22)
    orf_id_to_shu_FDR_L_ouORF_T = get_A2B(shu_fdr_file,  1, 23)

    orf_id_to_chrom = get_A2B("features/human.orf.genePred", 1, 2)
    orf_id_to_strand = get_A2B("features/human.orf.genePred", 1, 3)
    orf_id_to_orf_start = get_A2B("features/human.orf.genePred", 1, 6)
    orf_id_to_orf_end = get_A2B("features/human.orf.genePred", 1, 7)
    orf_id_to_GSE_tissue_exist = get_A2B("features/human.orf.GSE_tissue_exist", 1, 41)

    orf_id_to_phylocsf = get_A2B("features/human.phylocsf.txt", 1, 2) 
    orf_id_to_optimizedcodon = get_A2B("features/human.optimizedcodon.txt", 1, 4)
    orf_id_to_r_hydrophobic = get_A2B("features/human.orf.aa_feature", 1, 3)
    orf_id_to_r_amphipathic = get_A2B("features/human.orf.aa_feature", 1, 4)
    orf_id_to_r_polar = get_A2B("features/human.orf.aa_feature", 1, 5)
    orf_id_to_r_charged = get_A2B("features/human.orf.aa_feature", 1, 6)
    orf_id_to_kozak_score = get_A2B("features/human.orf.kozak_feature", 1, 2)
    orf_id_to_MFE_AUG_upstream = get_A2B("features/human.orf.MFE_feature", 1, 2) # U
    orf_id_to_MFE_AUG_downstream = get_A2B("features/human.orf.MFE_feature", 1, 3) # U
    orf_id_to_var_syn = get_A2B("features/human.orf.variance_count_per_aa", 1, 2)
    orf_id_to_var_mis = get_A2B("features/human.orf.variance_count_per_aa", 1, 3)
    orf_id_to_var_cid = get_A2B("features/human.orf.variance_count_per_aa", 1, 4)
    orf_id_to_var_spa = get_A2B("features/human.orf.variance_count_per_aa", 1, 5)
    orf_id_to_var_spd = get_A2B("features/human.orf.variance_count_per_aa", 1, 6)
    orf_id_to_var_did = get_A2B("features/human.orf.variance_count_per_aa", 1, 7)
    orf_id_to_var_sog = get_A2B("features/human.orf.variance_count_per_aa", 1, 8)
    orf_id_to_var_sal = get_A2B("features/human.orf.variance_count_per_aa", 1, 9)
    orf_id_to_var_fsh = get_A2B("features/human.orf.variance_count_per_aa", 1, 10)

    expression_columns = ["EMTAB7247h_Brain","EMTAB7247h_Liver","EMTAB7247h_Testis","GSE101760_A549","GSE103308_Muscle","GSE105082_HeLa","GSE125218_HEK293T","GSE125218_HeLa","GSE125218_K562","GSE129061_HepG2","GSE129061_K562","GSE131650_iPSC","GSE143263_A375","GSE143263_Bcell","GSE143263_HCT116","GSE42509_Fibroblast","GSE51424_BrainTumor","GSE56924_U2OS","GSE58207_HCT116","GSE59095_HEK293T","GSE59817_MCF10A","GSE59820_Kidney","GSE61742_Lymphoblastoid","GSE62247_ESC","GSE63591_HeLa","GSE64962_Fibroblast","GSE65885_Fibroblast","GSE65885_MCF10AERSrc","GSE67902_RPE1","GSE69923_MCF7","GSE69923_T47D","GSE70211_HEK293T","GSE71763_CHL1","GSE73136_HEK293","GSE78961_ESC","GSE78961_Neuron","GSE80156_HEK293","GSE89183_Erythroid","GSE94454_Huh7","GSE94460_HEK293","GSE_all"]
    dct_expression = {}
    for _, colname in enumerate(expression_columns):
        dct_expression[colname] = get_A2B("features/human.orf.RPKM_feature", 1, _ + 3)

    orf_id_to_verified = get_A2B("features/human.orf.verified", 1, 2)
    orf_id_to_westernblot_intensity = get_A2B("features/human.orf.westernblot_feature", 1, 2)
    orf_id_to_ucsc_matched_stop = get_A2B("features/human.orf.ucsc_matched", 1, 2)
    orf_id_to_ucsc_matched_stopstart = get_A2B("features/human.orf.ucsc_matched", 1, 3)
    orf_id_to_ucsc_matched_start = get_A2B("features/human.orf.ucsc_matched", 1, 4)
    orf_id_to_gencode_matched_start = get_A2B("features/human.orf.ucsc_matched", 1, 5)
    orf_id_to_gencode_matched_stop = get_A2B("features/human.orf.ucsc_matched", 1, 6)

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
        w.write("\t".join(["orf_id", "chrom", "strand", "orf_start", "orf_end", "GSE_tissue_exist", "if_GSE_tissue_multi_exist","if_candidate", "verified", "if_gencode_matched_start", "if_gencode_matched_stop", "ucsc_matched_start", "ucsc_matched_stop", "ucsc_matched_stopstart", "orf_type", "orf_classification", "orf_upstream_len", "orf_start_codon", "orf_pep_len", "westernblot_intensity", "phylocsf", "optimizedcodon", "r_hydrophobic", "r_amphipathic", "r_polar", "r_charged", "kozak_score", "MFE_AUG_upstream", "MFE_AUG_downstream", "int_FDR_L_all", "int_FDR_L_main", "int_FDR_L_uORF", "int_FDR_L_ouORF", "shu_FDR_L_all", "shu_FDR_L_main", "shu_FDR_L_uORF", "shu_FDR_L_ouORF", "max_FDR_L_all", "transcript_id", "transcript_len", "UTR5_len", "canonical_len", "UTR3_len", "gene_id", "gene_name", "gene_type", "gene_classification", "synonymous_variant", "missense_variant", "conservative_inframe_deletion", "splice_acceptor_variant", "splice_donor_variant", "disruptive_inframe_deletion", "stop_gained", "start_lost", "frameshift_variant"] + expression_columns) + "\n")

        for orf_id in orf_list:
            orf = Orf(species, orf_id)
            chrom = orf_id_to_chrom[orf_id]
            strand = orf_id_to_strand[orf_id]
            orf_start = orf_id_to_orf_start[orf_id]
            orf_end = orf_id_to_orf_end[orf_id]
            GSE_tissue_exist = orf_id_to_GSE_tissue_exist[orf_id]
            if_GSE_tissue_multi_exist = 0
            if int(GSE_tissue_exist) > 1:
                if_GSE_tissue_multi_exist = 1
            
            if float(orf_id_to_int_FDR_L_all[orf_id]) > 0:
                int_FDR_L_all = float(orf_id_to_int_FDR_L_all[orf_id]) / float(orf_id_to_int_FDR_L_all_T[orf_id])
            else:
                int_FDR_L_all = 1 / float(orf_id_to_int_FDR_L_all_T[orf_id])

            if float(orf_id_to_int_FDR_L_main[orf_id]) > 0:
                int_FDR_L_main = float(orf_id_to_int_FDR_L_main[orf_id]) / float(orf_id_to_int_FDR_L_main_T[orf_id])
            else:
                int_FDR_L_main = 1 / float(orf_id_to_int_FDR_L_main_T[orf_id])

            if float(orf_id_to_int_FDR_L_uORF[orf_id]) > 0:
                int_FDR_L_uORF = float(orf_id_to_int_FDR_L_uORF[orf_id]) / float(orf_id_to_int_FDR_L_uORF_T[orf_id])
            else:
                int_FDR_L_uORF = 1 / float(orf_id_to_int_FDR_L_uORF_T[orf_id])

            if float(orf_id_to_int_FDR_L_ouORF[orf_id]) > 0:
                int_FDR_L_ouORF = float(orf_id_to_int_FDR_L_ouORF[orf_id]) / float(orf_id_to_int_FDR_L_ouORF_T[orf_id])
            else:
                int_FDR_L_ouORF = 1 / float(orf_id_to_int_FDR_L_ouORF_T[orf_id])

            if float(orf_id_to_shu_FDR_L_all[orf_id]) > 0:
                shu_FDR_L_all = float(orf_id_to_shu_FDR_L_all[orf_id]) / float(orf_id_to_shu_FDR_L_all_T[orf_id])
            else:
                shu_FDR_L_all = 1 / float(orf_id_to_shu_FDR_L_all_T[orf_id])

            if float(orf_id_to_shu_FDR_L_main[orf_id]) > 0:
                shu_FDR_L_main = float(orf_id_to_shu_FDR_L_main[orf_id]) / float(orf_id_to_shu_FDR_L_main_T[orf_id])
            else:
                shu_FDR_L_main = 1 / float(orf_id_to_shu_FDR_L_main_T[orf_id])

            if float(orf_id_to_shu_FDR_L_uORF[orf_id]) > 0:
                shu_FDR_L_uORF = float(orf_id_to_shu_FDR_L_uORF[orf_id]) / float(orf_id_to_shu_FDR_L_uORF_T[orf_id])
            else:
                shu_FDR_L_uORF = 1 / float(orf_id_to_shu_FDR_L_uORF_T[orf_id])

            if float(orf_id_to_shu_FDR_L_ouORF[orf_id]) > 0:
                shu_FDR_L_ouORF = float(orf_id_to_shu_FDR_L_ouORF[orf_id]) / float(orf_id_to_shu_FDR_L_ouORF_T[orf_id])
            else:
                shu_FDR_L_ouORF = 1 / float(orf_id_to_shu_FDR_L_ouORF_T[orf_id])

            max_FDR_L_all = max(int_FDR_L_all, shu_FDR_L_all)

            verified = orf_id_to_verified.get(orf_id, "-")
            ucsc_matched_stop = orf_id_to_ucsc_matched_stop[orf_id]
            ucsc_matched_stopstart = orf_id_to_ucsc_matched_stopstart[orf_id]
            ucsc_matched_start = orf_id_to_ucsc_matched_start[orf_id]
            gencode_matched_start = orf_id_to_gencode_matched_start[orf_id]
            gencode_matched_stop = orf_id_to_gencode_matched_stop[orf_id]

            westernblot_intensity = orf_id_to_westernblot_intensity.get(orf_id, "NA")
            phylocsf = orf_id_to_phylocsf.get(orf_id, "-10")
            if phylocsf == "":
                phylocsf = "-10"
            optimizedcodon = orf_id_to_optimizedcodon.get(orf_id, "0")
            r_hydrophobic = orf_id_to_r_hydrophobic[orf_id]
            r_amphipathic = orf_id_to_r_amphipathic[orf_id]
            r_polar = orf_id_to_r_polar[orf_id]
            r_charged = orf_id_to_r_charged[orf_id]

            kozak_score = orf_id_to_kozak_score[orf_id]
            MFE_AUG_upstream = orf_id_to_MFE_AUG_upstream[orf_id]
            MFE_AUG_downstream = orf_id_to_MFE_AUG_downstream[orf_id]

            var_syn = orf_id_to_var_syn.get(orf_id, 0)
            var_mis = orf_id_to_var_mis.get(orf_id, 0)
            var_cid = orf_id_to_var_cid.get(orf_id, 0)
            var_spa = orf_id_to_var_spa.get(orf_id, 0)
            var_spd = orf_id_to_var_spd.get(orf_id, 0)
            var_did = orf_id_to_var_did.get(orf_id, 0)
            var_sog = orf_id_to_var_sog.get(orf_id, 0)
            var_sal = orf_id_to_var_sal.get(orf_id, 0)
            var_fsh = orf_id_to_var_fsh.get(orf_id, 0)

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
            if_candidate = 0
            if verified != "-" or ucsc_matched_stopstart != "-":
                if_candidate = 1
            lst_to_print = [orf_id, chrom, strand, orf_start, orf_end, GSE_tissue_exist, if_GSE_tissue_multi_exist, if_candidate, verified, gencode_matched_start, gencode_matched_stop, ucsc_matched_start, ucsc_matched_stop, ucsc_matched_stopstart, orf.orf_type, orf_classification, orf.orf_start, orf.start_codon, orf.pep_len, westernblot_intensity, phylocsf, optimizedcodon, r_hydrophobic, r_amphipathic, r_polar, r_charged, kozak_score, MFE_AUG_upstream, MFE_AUG_downstream, int_FDR_L_all, int_FDR_L_main, int_FDR_L_uORF, int_FDR_L_ouORF, shu_FDR_L_all, shu_FDR_L_main, shu_FDR_L_uORF, shu_FDR_L_ouORF, max_FDR_L_all, transcript_id, orf.transcript_len, UTR5_len, CDS_len, UTR3_len, gene_id, gene_name, gene_type, gene_classification, var_syn, var_mis, var_cid, var_spa, var_spd, var_did, var_sog, var_sal, var_fsh]
            w.write("\t".join([str(_) for _ in lst_to_print]) + "\t")
            w.write("\t".join([dct_expression[_][orf_id] for _ in expression_columns]) + "\n")
