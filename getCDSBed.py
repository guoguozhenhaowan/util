#!/usr/bin/env python3

import os
import sys
import gzip
import argparse
import pandas as pd
from pyfaidx import Fasta

def get_cds_region(cds_start, cds_end, exon_start_l, exon_end_l):
    """get cds region of transcript"""
    lower_index = -1
    upper_index = -1

    for i in range(len(exon_start_l)):
        if lower_index == -1 and  exon_end_l[i] >= cds_start:
            lower_index = i
        if upper_index == -1 and exon_end_l[i] >= cds_end:
            upper_index = i
    
    lower_reg = 0;
    upper_reg = 0;
    for i in range(len(exon_start_l)):
        if i < lower_index and i < upper_index:
            lower_reg += exon_end_l[i] - exon_start_l[i]
            upper_reg += exon_end_l[i] - exon_start_l[i]
            continue
        if i == lower_index:
            lower_reg += cds_start - exon_start_l[i]
        if i < upper_index:
            upper_reg += exon_end_l[i] - exon_start_l[i]
        if i == upper_index:
            upper_reg += cds_end - exon_start_l[i]
            break

    return lower_reg, upper_reg


def main():
    """main function to accept commandline arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gen", help="refGene table", required=True)
    parser.add_argument("-c", "--chr", help="onlyget Chrom[default: True]", action="store_false", required=False)
    parser.add_argument("-o", "--out", help="output bed path", required=False)
    args = parser.parse_args()
    # parse args
    ref_gene_with_version = args.gen
    out_p = os.path.join(os.getcwd(), "reg.bed")
    if args.out is not None:
        out_p = args.out
    tmp_str = "#Transcript\tStart\tEnd\tStrand\tChrom\tTransVer\n"
    # parse refGene table (refGene official table with transcript version column added
    ref_gene_df = pd.read_table(ref_gene_with_version, sep='\t', header=0)
    valid_chr = ["chr" + str(i) for i in range(1, 23)]
    valid_chr.extend(["chrX", "chrY", "chrM"])
    for index, row in ref_gene_df.iterrows():
        cds_start = int(row['cdsStart'])
        cds_end = int(row['cdsEnd'])
        exon_start_l = [int(k) for k in row['exonStarts'].rstrip(",").split(",")]
        exon_end_l = [int(k) for k in row['exonEnds'].rstrip(",").split(",")]
        strand = row['strand']
        trs_name = row['name']
        trs_ver = row['version']
        chrom = row['chrom']
        if args.chr and chrom not in valid_chr:
            continue
        start, end = get_cds_region(cds_start, cds_end, exon_start_l, exon_end_l)
        rec_str = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(trs_name, start, end, strand, chrom, trs_ver)
        tmp_str += rec_str
    
    with open(out_p, "w") as fw:
        fw.write(tmp_str)

if __name__ == "__main__":
    main()

