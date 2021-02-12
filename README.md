## TileSeq mutation count package

This package is made to parse input sequecning files (fastq) with user provided parameter.csv file.
Output of this pipeline is mutation counts for each pair of fastq files.

## Dependencies

`python 3.7/3.8 (tested mainly under py3.7)`

`R 3.4.4+`

`Bowtie2 Bowtie2-build`

## Installation 
Please use conda to set up the environment before installing the package: 

`conda install -n <env_name>`

You will also need the script `csv2json.R` which can be installed via installing [tileseqMave](https://github.com/jweiletileseqMave). Make sure `csv2json.R` can be found in `$PATH`

The alpha version is available by running:

`python -m pip install TileSeqMut`

To update to the newest stable release:

`python -m pip install TileSeqMut==0.4.201`

### Execution
---

After installation, you can run the package: 

```
tileseq_mut -p ~/path/to/paramSheet.csv -o ~/path/to/output_folder -f ~/path/to/fastq_file_folder/ -name
 name_of_the_run 
```

**Examples:**

``` bash
# on DC
tileseq_mut -p $HOME/dev/tilseq_mutcount/190506_param_MTHFR.csv -o $HOME/dev/tilseq_mutcount/output/ -f $HOME
/tileseq_data/WT/ -name MTHFR_test

# on BC2
tileseq_mut -p $HOME/dev/tilseq_mutcount/190506_param_MTHFR.csv -o $HOME/dev/tilseq_mutcount/output/ -f $HOME
/tileseq_data/WT/ -name MTHFR_test -env BC2
```
* This command will analyze fastq files in the folder: `~/tileseq_data/WT/` and make a time stamped output folder
 with the prefix: `MTHFR_test` in `$HOME/dev/tilseq_mutcount/output/` (Using all default parameters, see below)


**Parameters**

* Run `tileseq_mut --help`

``` bash
(py37) [rli@dc06 DC_jobs]$ tileseq_mut -h
usage: tileseq_mut [-h] [-f FASTQ] -o OUTPUT -p PARAM -n NAME
                   [--skip_alignment] [-r1 R1] [-r2 R2] [-log LOG_LEVEL]
                   [-env ENVIRONMENT] [-at AT] [-mt MT] [-c C] [-b BASE]
                   [-test] [-rc] [-override]

TileSeq mutation counts

optional arguments:
  -h, --help            show this help message and exit
  -f FASTQ, --fastq FASTQ
                        Path to all fastq files you want to analyze
  -o OUTPUT, --output OUTPUT
                        Output folder
  -p PARAM, --param PARAM
                        csv paramter file
  -n NAME, --name NAME  Name for this run
  --skip_alignment      skip alignment for this analysis, ONLY submit jobs for
                        counting mutations in existing output folder
  -r1 R1                r1 SAM file
  -r2 R2                r2 SAM file
  -log LOG_LEVEL, --log_level LOG_LEVEL
                        set log level: debug, info, warning, error, critical.
                        (default = info)
  -env ENVIRONMENT, --environment ENVIRONMENT
                        The cluster used to run this script (default = DC)
  -at AT                Alignment time (default = 8h)
  -mt MT                Mutation call time (default = 48h)
  -c C                  Number of cores to use for mutation counting
  -b BASE, --base BASE  ASCII code base
  -test                 Turn on test mode
  -rc                   Turn on rc mode, both direction of the reads will be
                        aligned to the reference. Variant calling will be
                        performed on all the reads that are aligned, regardless of their direction (BE
                        CAREFUL!)
  -override, --sr_Override
                        Provide this argument when there is only one replicate

```

**Start the run**

* Once the run starts, it will first submit alignment jobs to the cluster and keep tracking of all the submitted
 alignment jobs. Once all the jobs are finished, the pipeline will submit another batch of jobs for mutation calling.
 
 * if you want to skip alignment and only do mutation calls for existing sam files you can run the following command:
 
* **Example of skipping alignment:**

```
tileseq_mut -p $HOME/dev/tilseq_mutcount/190506_param_MTHFR.csv -o /home/rothlab1/rli/dev/tilseq_mutcount/output
/190506_MTHFR_WT_2020-01-29-17-07-04/ --skip_alignment
```

### Input files
---

`/path/to/fastq/` - Full path to input fastq files

`parameters.csv` - CSV file contains information for this run (please see example
[here](https://docs.google.com/spreadsheets/d/1tIblmIFgOApPNzWN2KUwj8BKzBiJ1pOL7R4AOUGrqvE/edit?usp=sharing)
).
This file is required to be comma-seperated and saved in csv format.


### Output files
---

One output folder is created for each run. The output folder are named with `name_time-stamp`

Within each output folder, the following files and folders will be generated:

`./main.log` - main logging file for alignment

`./args.log` - arguments for this run

`./ref/` - Reference fasta file and bowtie2 index

`./env_aln_sh/` - Bash scripts for submitting the alignment jobs

`./sam_files/` - Alignment output and log files for the raw fastq files

`./name_time-stamped_mut_count/` - Mutation counts in each sample are saved in csv files

    - `./main.log` - Main log file for mutation calling

    - `./args.log` - command line arguments
    
    - `./info.csv` - Meta information for each sample: sequencing depth, tile starts/ends and # of reads mapped outside of the targeted tile

    - `./count_sample_*.csv` - Raw mutation counts for each sample. With meta data in header. Variants are represented in hgvs format

    - `./env_mut/` - Bash scripts for summitting the mutation count jobs, also log files for each sample.

The count_sample_\*\*.csv is passed to tileseqMave for further analysis

### Alignment
---

The pipeline takes the sequence in the parameter file as reference and align the fastq files
to the whole reference sequence. This is the sequence specified by user in the parameter file.

For each pair of fastq files (R1 and R2), the pipeline submits one alignment job to the cluster. In the folder `env_sh` you can find all the scripts that were submitted to the cluster when you run `main.py`.

Alignments were done using `Bowtie2` with following parameters:

```
~/bowtie2 --no-head --norc --no-sq --local -x {ref} -U {r1} -S {r1_sam_file}
~/bowtie2 --no-head --nofw --no-sq --local -x {ref} -U {r2} -S {r2_sam_file}
```

### Mutation Calls
---

From each pair of sam files we count mutations for each sample.

We first filter out reads that did not map to reference or reads that are outside of the tile. Then pass the rest of the reads to `count_mut.py`. Please read the wiki page about how to call mutations using CIGAR string and MD:Z tag.

In order to eliminate sequencing errors. We apply a posterior probability cut-off. The posterior probability of a mutation was calculated using the Phred scores provided in SAM files.
