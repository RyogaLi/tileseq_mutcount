#!/usr/bin/env python3.6

import pandas as pd
import os
import math
import glob
from Bio import SeqIO

class SeqParser(object):

    def __init__(self, ref_dir, mut2_func, gene):
        # read from mut2func file
        self._gene = gene
        self._tiles = mut2_func

        fasta = os.path.join(ref_dir, gene+".fasta")
        fasta_coding = os.path.join(ref_dir, gene+"_coding.fasta")
        self._seq = list(SeqIO.parse(fasta, "fasta"))[0].seq
        self._seq_coding = list(SeqIO.parse(fasta_coding, "fasta"))[0].seq
        self._start_pos = self._seq.find(self._seq_coding)

    def _map_loc(self):
        """
        read mut2func file, from the file map sample names and tile info
        """
        pass

    def _get_samples(self):
        """
        from the tile doc, get samples and corresponding tiles
        * the first 5 columns are fixed
        """
        # get samples
        samples = self._tiles.iloc[:,5:]
        tiles = self._tiles.tile
        conditions = samples.columns

        samples.index = tiles

        #print(samples.to_dict("index"))
        return samples


    def _create_fasta(self):
        """
        based on the tile information, get ref sequence and make fasta file
        self._tiles: file contains all the file information of that gene
            region, start, stop, len
        """

        # convert start and stop pos to nt pos
        # include 30bp from upstream and 30bp from downstream
        self._tiles["nt_start"] = self._tiles.start *3 -3 +self._start_pos-30 # inclusive
        self._tiles["nt_end"] = self._tiles.end * 3+self._start_pos+30


        # select corresponding sequence and write to fasta
        for index, row in self._tiles.iterrows():
            f_name = "tile_"+str(row.tile)+".fasta"
            print(f_name)
            with open(ref_dir+f_name, "w") as fa:
                fa.write(">"+self._gene+"\n")
                seq = self._seq[int(row.nt_start):int(row.nt_end)]
                row["seq"] = seq
                fa.write(str(seq))
            cmd = parameters.BOWTIE2_BUILD+ref_dir+f_name+" --quiet "+os.path.join(ref_dir,"tile_"+str(row.tile))
            os.system(cmd)


    def _make_sh_files(self, fastq_dir, sh_dir, sam_dir):
        """
        in the tiles file, user can add as many columns as they want for all the fastq files
        for each fastq file, use corresponding tile fasta reference and map them
        fastq_dir: directory that contains fastq files
        sh_dir: directory that contains sh files (for sge submission)
        sam_dir: direcotry for output sam/bam files
        """
        samples = self._get_samples()

        for index, row in samples.iterrows():
            tile = "tile_" + str(index)
            tile_ref = ref_dir + tile
            for i in set(row.values.tolist()):
                if math.isnan(i): continue
                # get R1 from fastq dir
                r1 = glob.glob(fastq_dir+str(int(i))+"*_R1_*")

                # get R2 from fastq dir
                r2 = glob.glob(fastq_dir+str(int(i))+"*_R2_*")

if __name__ == "__main__":
    ref_dir = "./MTHFR_3del/bowtieIndex/"
    tiles = "./MTHFR_tiles.csv"
    fastq_dir = "../MTHFR/190510_MTHFR_3del_strain/WT/"
    sh_dir = "./MTHFR_3del/shfiles/"
    sam_dir = "/MTHFR_3del/samfiles/"

    tiles = pd.read_csv(tiles, sep="\t")
    sp = SeqParser(ref_dir,tiles, "MTHFR")
    #sp._get_samples()
    #sp._create_fasta()
    sp._make_sh_files(fastq_dir, sh_dir, sam_dir)
