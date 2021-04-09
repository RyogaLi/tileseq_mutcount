#!/usr/bin/env python3.7

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca
import re
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


def main(sanger_data, fasta_dir, output_dir):
    """

    """
    # LDLR_ref = f"{fasta_dir}/LDLR_ref.fasta"
    # tile2_ref = f"{fasta_dir}/tile2_ref.fasta"
    # tile2_rows = "ABCDEF"
    tile3_ref = f"{fasta_dir}/tile3_ref.fasta"
    tile4_ref = f"{fasta_dir}/tile4_ref.fasta"

    # go through all the seq files
    t3 = re.compile("(A[0-9]+$)|(B[1-5]$)")
    df = pd.DataFrame([])
    for f in glob.glob(f"{sanger_data}/*.seq"):
        basename = os.path.basename(f)
        # if it is in A-F
        sample_name = basename.split(".")[0]
        well = sample_name.split("_")[0]
        match = re.match(t3, well)
        if match:
            # blast to get overall summary
            cmd = f'blastn -db {tile3_ref} -out test.txt -outfmt "10 qseqid sseqid qlen slen bitscore mismatch" -query ' + f
            os.system(cmd)
            # blast to get each alignment
            cmd = f'blastn -db {tile3_ref} -out {output_dir}/{sample_name}.txt -outfmt "0 qseqid sseqid qlen slen ' \
                  f'bitscore mismatch" -query {f}'
            os.system(cmd)
        else:
            print(well)
            # blast to get overall summary
            cmd = f'blastn -db {tile4_ref} -out test.txt -outfmt "10 qseqid sseqid qlen slen bitscore mismatch" ' \
                  f'-query {f}'
            os.system(cmd)
            # blast to get each alignment
            cmd = f'blastn -db {tile4_ref} -out {output_dir}/{sample_name}.txt -outfmt "0 qseqid sseqid qlen ' \
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
    output_file = os.path.join(output_dir, "sanger_output_raw_t3_t4.csv")
    df.to_csv(output_file, index=False)

if __name__ == '__main__':
    sanger_data = "/Users/roujia/Desktop/210323_LDLR_sanger/sanger_seq_1857750/"
    output_dir = "/Users/roujia/Desktop/210323_LDLR_sanger/"
    fasta_dir = "/Users/roujia/Desktop/210323_LDLR_sanger/fasta/"
    # make fasta db for tile 2 and full gene
    # from fasta file make blast db
    fasta_LDLR = "/Users/roujia/Desktop/210323_LDLR_sanger/fasta/LDLR_ref.fasta"
    cmd = f"makeblastdb -in {fasta_LDLR} -parse_seqids -dbtype nucl"
    os.system(cmd)

    fasta_tile2 = "/Users/roujia/Desktop/210323_LDLR_sanger/fasta/tile2_ref.fasta"
    cmd = f"makeblastdb -in {fasta_tile2} -parse_seqids -dbtype nucl"
    os.system(cmd)

    fasta_tile3 = "/Users/roujia/Desktop/210323_LDLR_sanger/fasta/tile3_ref.fasta"
    cmd = f"makeblastdb -in {fasta_tile3} -parse_seqids -dbtype nucl"
    os.system(cmd)

    fasta_tile4 = "/Users/roujia/Desktop/210323_LDLR_sanger/fasta/tile4_ref.fasta"
    cmd = f"makeblastdb -in {fasta_tile4} -parse_seqids -dbtype nucl"
    os.system(cmd)

    main(sanger_data, fasta_dir, output_dir)