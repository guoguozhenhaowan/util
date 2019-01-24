#!/usr/bin/env python3

import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import os
import argparse
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "wljlinksly"

def get_gb(acc_no, out_dir):
    """get genbank record"""
    out_file = os.path.join(out_dir, "{0}.gb".format(acc_no))
    if os.access(out_file, os.F_OK):
        return 0
    handle = Entrez.efetch(db="nucleotide", id=acc_no, rettype="gb", retmode="text")
    with open(out_file, "w") as fw:
        fw.write(handle.read())
    handle.close()

def main():
    parser = argparse.ArgumentParser(usage="get genbank record from NCBI")
    parser.add_argument("-i", "--input", help="input transcript id list", required=True)
    parser.add_argument("-o", "--output", help="output parent dir", required=False)
    args = parser.parse_args()
    tr_id_f = args.input
    out_dir = os.getcwd()
    if args.output is not None:
        out_dir = args.output
        if not os.access(out_dir, os.F_OK):
            os.makedirs(out_dir)
    with open(tr_id_f, "r") as fr:
        for line in fr:
            tr_id = line.strip()
            get_gb(tr_id, out_dir)

if __name__ == "__main__":
    main()
