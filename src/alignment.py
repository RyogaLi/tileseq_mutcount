#!/usr/bin/python3.6

# THis script is used to align fatq files to reference sequence
# Input: list of fastq files 
# Output: samfiles

import os
import settings

class Alignment(object):

    def __init__(self, fasta, fastq, output_sam, log):
        self._fasta = fasta # path to fasta reference files 
        self._fastq = fastq # path to fastq reference files
        self._output_sam = output_sam # path to save the sam output files
        self._log = log # alignemnt log (saved in the same folder as sam files 

    def _align(self):

        # align with bowtie 
        # path to bowtie and bowtie build can be found in settings 
    
        # check if bowtie index file were built 
        basename = os.path.basename(self._fasta)
        bt = basename.replace(".fasta", "")

        # align fastq 
        
        # direct output sam file to output sam folder
        pass


def align_main(ref, fastq, output_sam, shfiles):
    """
    ref: reference fatsa file
    fastq: input fastq file
    output_sam: path to sam files
    shfiles: path to folder to save files to be submitted to the cluster
    """

    pass


if __name__ == "__main__":

    # for testing this script

