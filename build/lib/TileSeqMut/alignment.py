#!/usr/bin/env python3.7

# THis script is used to align fatq files to reference sequence
# Input: list of fastq files
# Output: samfiles

import os

def make_ref(name, ref_seq, ref_path):
    """
    given the reference sequence and ref_path
    make fasta file and build from fasta
    """
    if name == "phix":
        ref_fasta = os.path.join(ref_path, "phix.fasta")
    else:
        ref_fasta = os.path.join(ref_path, name+".fasta")
        with open(ref_fasta, "w") as fasta:
            fasta.write(">"+name+"\n")
            fasta.write(ref_seq+"\n")
    build_cmd =f"bowtie2-build --quiet -f {ref_fasta} {os.path.join(ref_path, name)}"
    os.system(build_cmd)

    return os.path.join(ref_path, name)

def align_main(ref, r1, r2, sam_path, shfile, rc=False, header=""):
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

    if not rc: # do not need to check for reverse complement
        r1_cmd = f"bowtie2 --no-head --norc --no-sq --rdg 12,1 --rfg 12,1 --local -x {ref} -U {r1} -S {r1_sam_file}"
        r2_cmd = f"bowtie2 --no-head --nofw --no-sq --rdg 12,1 --rfg 12,1 --local -x {ref} -U {r2} -S {r2_sam_file}"
    else:
        r1_cmd = f"bowtie2 --no-head --no-sq --rdg 12,1 --rfg 12,1 --local -x {ref} -U {r1} -S {r1_sam_file}"
        r2_cmd = f"bowtie2 --no-head --no-sq --rdg 12,1 --rfg 12,1 --local -x {ref} -U {r2} -S {r2_sam_file}"

    with open(shfile, "w") as f:
        # header is only used on galen
        if header != "":
            f.write(header)
        f.write(r1_cmd + "\n")
        f.write(r2_cmd + "\n")
    os.system(f"chmod 755 {shfile}")
    return r1_sam_file, r2_sam_file, log_file
