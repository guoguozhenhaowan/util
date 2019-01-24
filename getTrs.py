#!/usr/bin/env python3

import os
import sys
import gzip
import argparse
import pandas as pd
from pyfaidx import Fasta

def get_trs_seq(ref_fa, chrom, exon_start_l, exon_end_l):
    """get seq of transcript"""
    cds_seq = ''
    for i in range(len(exon_start_l)):
        cds_seq += str(ref_fa[chrom][exon_start_l[i]:exon_end_l[i]])
    return cds_seq 

def main():
    """main function to accept commandline arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--ref", help="refGeonme fa", required=True)
    parser.add_argument("-g", "--gen", help="refGene table", required=True)
    parser.add_argument("-c", "--chr", help="onlyget Chrom[default: True]", action="store_false", required=False)
    parser.add_argument("-o", "--out", help="output fa path", required=False)
    args = parser.parse_args()
    # parse args
    ref_genome_fa = args.ref
    ref_fa = Fasta(ref_genome_fa)
    ref_gene_with_version = args.gen
    out_p = os.path.join(os.getcwd(), "refgene.fa")
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
        cds_seq = get_trs_seq(ref_fa, chrom, exon_start_l, exon_end_l)
        rec_str = ">{0} {1}\n{2}\n".format(trs_name, trs_ver, cds_seq)
        tmp_str += rec_str
    
    with open(out_p, "w") as fw:
        fw.write(tmp_str)

if __name__ == "__main__":
    main()
