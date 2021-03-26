#!/usr/bin/env python3.7

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import sys

sys.path.append('..')

import os
import glob
import pandas as pd

# sanger analysis for LDLR wt data
# Info:
# Row A to F -> Should include tile 2 region of the wildtype LDLR sequence.
# Because, it is only the Tile 2 portion, both reads should cover the full nt sequence.
# Row G&H -> Should include full wildtype LDLR sequence. Due to the length of full LDLR nt sequence,
# the reads will most likely not cover the entire ORF length. So, we are just looking for mismatches
# in the nt sequence that is covered from the M13F and M13R primers.


def main(sanger_data, full_gene_fasta, tile2_fasta, output_dir):
    """

    """
    tile2_rows = "ABCDEF"
    full_gene_rows = "GH"
    # go through all the seq files
    df = pd.DataFrame([])
    for f in glob.glob(f"{sanger_data}/*.seq"):
        basename = os.path.basename(f)
        # if it is in A-F
        well = basename[0]
        sample_name = basename.split(".")[0]
        if well in tile2_rows:
            # blast to get overall summary
            cmd = f'blastn -db {tile2_fasta} -out test.txt -outfmt "10 qseqid sseqid qlen slen bitscore mismatch" -query ' + f
            os.system(cmd)
            # blast to get each alignment
            cmd = f'blastn -db {tile2_fasta} -out {output_dir}/{sample_name}.txt -outfmt "0 qseqid sseqid qlen slen ' \
                  f'bitscore mismatch" -query {f}'
            os.system(cmd)
        else:
            # blast to get overall summary
            cmd = f'blastn -db {full_gene_fasta} -out test.txt -outfmt "10 qseqid sseqid qlen slen bitscore mismatch" -query ' + f
            os.system(cmd)
            # blast to get each alignment
            cmd = f'blastn -db {full_gene_fasta} -out {output_dir}/{sample_name}.txt -outfmt "0 qseqid sseqid qlen ' \
                  f'slen bitscore mismatch" -query {f}'
            # print(cmd)
            os.system(cmd)

        if os.path.getsize("./test.txt") == 0: continue

        tmp = pd.read_csv("./test.txt", sep=",", header=None)
        tmp["f_name"] = f
        if df.empty:
            df = tmp
        else:
            df = df.append(tmp)
    df.columns = ["Query seq-ID", "Subject seq-ID", "Query length", "Subject length", "Bitscore", "n_mismatch",
                  "file_name"]
    output_file = os.path.join(output_dir, "sanger_output_raw.csv")
    df.to_csv(output_file, index=False)

if __name__ == '__main__':
    sanger_data = "/Users/roujia/Desktop/210323_LDLR_sanger/sanger_seq_1852542/"
    output_dir = "/Users/roujia/Desktop/210323_LDLR_sanger/"

    # make fasta db for tile 2 and full gene
    # from fasta file make blast db
    fasta_LDLR = "/Users/roujia/Desktop/210323_LDLR_sanger/LDLR_ref.fasta"
    cmd = f"makeblastdb -in {fasta_LDLR} -parse_seqids -dbtype nucl"
    os.system(cmd)

    fasta_tile2 = "/Users/roujia/Desktop/210323_LDLR_sanger/tile2_ref.fasta"
    cmd = f"makeblastdb -in {fasta_tile2} -parse_seqids -dbtype nucl"
    os.system(cmd)

    main(sanger_data, fasta_LDLR, fasta_tile2, output_dir)