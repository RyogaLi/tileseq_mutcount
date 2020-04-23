#!/usr/bin/env python3.6

# THis script is used to align fatq files to reference sequence
# Input: list of fastq files
# Output: samfiles

import os

def make_ref(name, ref_seq, ref_path, bowtie_build):
    """
    given the reference sequence and ref_path
    make fasta file and build from fasta
    """
    ref_fasta = os.path.join(ref_path, name+".fasta")
    with open(ref_fasta, "w") as fasta:
        fasta.write(">"+name+"\n")
        fasta.write(ref_seq+"\n")
    build_cmd =f"{bowtie_build}--quiet -f {ref_fasta} {os.path.join(ref_path, name)}"
    os.system(build_cmd)

    return os.path.join(ref_path, name)

def align_main(ref, r1, r2, sam_path, bowtie, shfile):
    """
    ref: reference fatsa file
    r1: input fastq file - R1
    r2: input fastq file - R2
    sam_path: path to sam files
    bowtie2: path to bowtie2
    shfile: write the command to this file
    return shfile
    """
    log_file = os.path.join(sam_path, os.path.basename(shfile).replace(".sh", ".log"))

    r1_sam_file = os.path.join(sam_path, os.path.basename(r1).replace(".fastq.gz", ".sam"))
    r2_sam_file = os.path.join(sam_path, os.path.basename(r2).replace(".fastq.gz", ".sam"))

    r1_cmd = f"{bowtie}--no-head --norc --no-sq --rdg 12,1 --rfg 12,1 --local -x {ref} -U {r1} -S {r1_sam_file}"
    r2_cmd = f"{bowtie}--no-head --nofw --no-sq --rdg 12,1 --rfg 12,1 --local -x {ref} -U {r2} -S {r2_sam_file}"

    with open(shfile, "w") as f:
        f.write(r1_cmd + "\n")
        f.write(r2_cmd + "\n")
    os.system(f"chmod 755 {shfile}")
    return r1_sam_file, r2_sam_file, log_file
