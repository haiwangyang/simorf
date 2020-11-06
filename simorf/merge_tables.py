"""

A script to merge FDR tables of different batch

Example script:
python merge_tables.py -t output/human.FDRs.sim1000.batchA.txt,output/human.FDRs.sim1000.batchB.txt -o output/human.FDRs.sim1000.merged.txt

Parameter description:
t  = FDR tables
o  = output
"""

__version__ = "1.0.1"
__email__ = "haiwang.yang@northwestern.edu"


from timeit import default_timer as timer
from common_functions import Orf, Transcript2, get_elements, get_A2B, object_to_pickle, pickle_to_object
import argparse
import os
import re


if __name__ == "__main__":
    # fetching parameters
    parser = argparse.ArgumentParser(description='please provide information')
    parser.add_argument('-t', '--tables', type=str) # FDR tables
    parser.add_argument('-o', '--output', type=str) # output

    args = parser.parse_args()
    tables = args.tables
    output = args.output

    print("collecting data from table")
    dct = {} 
    for table in tables.split(","):
        print("\t", table)
        with open(table, "r") as r:
            for line in r.readlines()[1:]:
                elements = line.rstrip().split("\t")
                orf_id = elements[0]
                if not orf_id in dct:
                   dct[orf_id] = []
                dct[orf_id].append(elements[1:])

    print("writing merged table") 
    colnames = ["orf_id","orf_type","orf_pep_len","FDR_S_all","FDR_L_all","FDR_S_main","FDR_L_main","FDR_S_uORF","FDR_L_uORF","FDR_S_ouORF","FDR_L_ouORF","shorter_all","longer_all","total_all","shorter_main","longer_main","total_main","shorter_uORF","longer_uORF","total_uORF","shorter_ouORF","longer_ouORF","total_ouORF"]
    with open(output, "w") as w:
        w.write("\t".join(colnames) + "\n")
        
        for orf_id in dct.keys():
            dct1_sum = {}
            dct1_sum["orf_id"] = orf_id
            for items in dct[orf_id]:
                dct1 = dict(zip(colnames[1:], items))
                for k in dct1.keys():
                    if "orf" in k:
                        dct1_sum[k] = dct1[k]
                    elif "shorter" in k or "longer" in k or "total" in k:
                        if not k in dct1_sum:
                            dct1_sum[k] = 0
                        dct1_sum[k] += int(dct1[k])

            dct1_sum["FDR_S_all"] = round(dct1_sum["shorter_all"] /  dct1_sum["total_all"], 16)
            dct1_sum["FDR_L_all"] = round(dct1_sum["longer_all"] /  dct1_sum["total_all"], 16)
            dct1_sum["FDR_S_main"] = round(dct1_sum["shorter_main"] /  dct1_sum["total_main"], 16)
            dct1_sum["FDR_L_main"] = round(dct1_sum["longer_main"] /  dct1_sum["total_main"], 16)
            dct1_sum["FDR_S_uORF"] = round(dct1_sum["shorter_uORF"] /  dct1_sum["total_uORF"], 16)
            dct1_sum["FDR_L_uORF"] = round(dct1_sum["longer_uORF"] /  dct1_sum["total_uORF"], 16)
            dct1_sum["FDR_S_ouORF"] = round(dct1_sum["shorter_ouORF"] /  dct1_sum["total_ouORF"], 16)
            dct1_sum["FDR_L_ouORF"] = round(dct1_sum["longer_ouORF"] /  dct1_sum["total_ouORF"], 16)

            w.write("\t".join([str(dct1_sum[_]) for _ in colnames]) + "\n")


            
                
