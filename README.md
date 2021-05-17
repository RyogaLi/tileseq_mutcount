## TileSeq mutation count package

This package is made to parse input sequecning files (fastq) with user provided parameter.csv file.
Output of this pipeline is mutation counts for each pair of fastq files.

## Dependencies

`python 3.7/3.8 (tested mainly under py3.7)`

`R 3.4.4+`

`Bowtie2 Bowtie2-build`

## Installation 
Please use conda to set up the environment before installing the package: 

* if you don't have python3.7 installed:  `conda install python==3.7`

* create a python3.7 environment:  `conda create -n py37 python=3.7`

* activate an environment: `conda activate py37`

You will also need the script `csv2json.R` which can be installed via installing [tileseqMave](https://github.com/jweiletileseqMave). Make sure `csv2json.R` can be found in `$PATH`

To install the newest stable release:

`python -m pip install TileSeqMut==0.6.3`

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

TileSeq mutation counts

optional arguments:
  -h, --help            show this help message and exit
  -f FASTQ, --fastq FASTQ
                        Path to all fastq files you want to analyze
  -o OUTPUT, --output OUTPUT
                        Output folder
  -p PARAM, --param PARAM
                        csv or json paramter file
  -n NAME, --name NAME  Name for this run (required)
  --skip_alignment      skip alignment for this analysis, ONLY submit jobs for
                        counting mutations in existing output folder. 
  -r1 R1                r1 SAM file
  -r2 R2                r2 SAM file
  -log LOG_LEVEL, --log_level LOG_LEVEL
                        set log level: debug, info, warning, error, critical.
                        (default = info)
  -env ENVIRONMENT, --environment ENVIRONMENT
                        The cluster used to run this script (default = DC)
  -at AT                Alignment time (default = 8h)
  -mt MT                Mutation call time (default = 36h)
  -c C                  Number of cores to use for mutation counting (default = 16)
  -b BASE, --base BASE  ASCII code base (default = 33)
  -rc                   Turn on rc mode, both direction of the reads will be
                        aligned to the reference. Variant calling will be
                        performed on all the reads that are aligned, regardless of their direction (BE
                        CAREFUL!)
  -override, --sr_Override
                        Provide this argument when there is only one replicate
  --posteriorQC         Turn on posterior QC mode, this requires more memory and runtime, please change the arguments
                        accordingly
  --wt_override         When no wt conditions defined in the parameter sheet, turn on this option will treat
                        EVERYTHING as WT. Phred scores will be adjusted based on the first replicate. 
                        Only use this when you also have --calibratePhredWT on.

  --calibratePhredWT    When this parameter is provided, use WT samplese (first replicate) to calibrate phred scores.
  --calibratePhredPhix  When this parameter is provided, use phix alignments to calibrate phred scores. Note that you
                        need to make sure Undetermined_**.fastq.gz files are also in the fastq directory. 

  --resubmit            For a finished run, batch resubmit failed scripts (if any). Use the *_mut_count dir as -o. 
                        Also need to provide --skip_alignment

```

**Start the run**

* Once the run starts, it will first submit alignment jobs to the cluster and keep tracking of all the submitted
 alignment jobs. Once all the jobs are finished, the pipeline will automatically submit another batch of jobs for
  mutation calling.
 
 * if you want to skip alignment and only do mutation calls for existing sam files you can run the following command:
 
* **Example of skipping alignment:**

```
tileseq_mut -p ~/dev/tilseq_mutcount/190506_param_MTHFR.csv -o /home/rothlab1/rli/dev/tilseq_mutcount/output
/190506_MTHFR_WT_2020-01-29-17-07-04/ --skip_alignment -n rerun_mut_count
```

 * if you want to resubmit mutation calls for failed samples:
 
* **Example of resubmit:**

```
tileseq_mut -p ~/dev/tilseq_mutcount/190506_param_MTHFR.csv -o /home/rothlab1/rli/dev/tilseq_mutcount/output
/190506_MTHFR_WT_2020-01-29-17-07-04//190506_MTHFR_WT_2020-01-29-17-07-04_mut_count/ --resubmit --skip_alignment -n
 resub
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

`./env_jobs/` - Bash scripts for submitting the alignment jobs

`./sam_files/` - Alignment output (SAM) and log files for bowtie2

`./name_time-stamped_mut_count/` - Mutation counts in each sample are saved in csv files
    
    - `./count_sample_*.csv` - Raw mutation counts for each sample. With meta data in header. Variants are represented in hgvs format

    - `./env_jobs/` - Bash scripts for summitting the mutation count jobs, also log files for each sample. 
    
        * `./env_jobs/*.log` - log file from the cluster
        * `./env_jobs/*.sh`  - submission script
    
    - `coverage_sample_name.csv` - File contains read counts for each position in the tile. There are
     four columns in the file: 
     pos: nt position of the tile 
     m_both: Number of variants found on both reads covering the site 
     m_r1: Number of variants found only on Read 1
     m_r2: Number of variants found only on Read 2
     passed: Number of variants passed filter
    
    - `*_R1/R2_calibrate_phred.csv` - Calibrated Phred scores (if applicable)

The count_sample_\*\*.csv is passed to tileseqMave for further analysis

### Alignment
---

The pipeline takes the sequence in the parameter file as reference and align the fastq files
to the whole reference sequence. This is the sequence specified by user in the parameter file.

For each pair of fastq files (R1 and R2), the pipeline submits one alignment job to the cluster. In the folder `env_sh` you can find all the scripts that were submitted to the cluster when you run `main.py`.

Alignments were done using `Bowtie2` with following parameters:

```
bowtie2 --no-head --norc --no-sq --rdg 12,1 --rfg 12,1 --local -x {ref} -U {r1} -S {r1_sam_file}
bowtie2 --no-head --nofw --no-sq --rdg 12,1 --rfg 12,1 --local -x {ref} -U {r2} -S {r2_sam_file}
```

If `-rc` is provided, the following parameters are used. BE CAREFUL! In this case, the reads are aligned to both fw
 and rc reference and variants are called regardless of which strand the read mapped to.

```
bowtie2 --no-head --no-sq --rdg 12,1 --rfg 12,1 --local -x {ref} -U {r1} -S {r1_sam_file}
bowtie2 --no-head --no-sq --rdg 12,1 --rfg 12,1 --local -x {ref} -U {r2} -S {r2_sam_file}
```

### Mutation Calls
---

From each pair of sam files we count mutations for each sample.

We first filter out reads that did not map to reference or reads that are outside of the tile. Then pass the rest of the reads to `count_mut.py`. Please read the wiki page about how to call mutations using CIGAR string and MD:Z tag.

In order to eliminate sequencing errors. We apply a posterior probability cut-off. The posterior probability of a mutation was calculated using the Phred scores provided in SAM files.
