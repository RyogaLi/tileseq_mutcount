## TileSeq mutation count package

This package is made to parse input sequecning files (fastq) with user provided parameter.csv file.
Output of this pipeline is mutation counts for each pair of fastq files.

## Dependencies

`python 3.7/3.8 (tested mainly under py3.7)`

`R 3.6+` and [tileseqMave](https://github.com/jweiletileseqMave)

`Bowtie2` and `Bowtie2-build` (need to be in the same folder)

## Installation 
Please use conda to set up the environment before installing the package: 

* if you don't have python3.7 installed:  `conda install python==3.7`

* create a python3.7 environment:  `conda create -n py37 python=3.7`

* activate an environment: `conda activate py37`

You will also need the scripts: `csv2json.R` and `calibratePhred.R` which can be installed via installing [tileseqMave](https://github.com/jweiletileseqMave). Make sure `csv2json.R` can be found in `$PATH`

To install the package (please visit [pypi page](https://pypi.org/project/TileSeqMut/) for the newest version number):

`python -m pip install TileSeqMut==newest_version`

### Execution
---

After installation, you can run the package (on DC): 

```
tileseq_mut -p /path/to/paramSheet.csv -o /path/to/output_folder -f ~/path/to/fastq_file_folder/ -name name_of_the_run 
```

**Examples:**

``` bash
# on DC
tileseq_mut -p $HOME/dev/tilseq_mutcount/190506_param_MTHFR.csv -o $HOME/dev/tilseq_mutcount/output/ -f $HOME
/tileseq_data/WT/ -name MTHFR_test -env DC

# on GALEN
tileseq_mut -p $HOME/dev/tilseq_mutcount/190506_param_MTHFR.csv -o $HOME/dev/tilseq_mutcount/output/ -f $HOME
/tileseq_data/WT/ -name MTHFR_test
```
* This command will analyze fastq files in the folder: `~/tileseq_data/WT/` and make a time stamped output folder
 with the prefix: `MTHFR_test` in `$HOME/dev/tilseq_mutcount/output/` (Using all default parameters, see below)


**Parameters**

* Run `tileseq_mut --help`

``` bash
(py37) [rli@dc06 DC_jobs]$ tileseq_mut -h

TileSeq mutation counts

Required arguments:
  -f FASTQ, --fastq FASTQ
                        Path to all fastq files you want to analyze
  -o OUTPUT, --output OUTPUT
                        Output folder
  -p PARAM, --param PARAM
                        csv paramter file
  -n NAME, --name NAME  Name for this run

Optional arguments
  -h, --help            show this help message and exit
  --skip_alignment      skip alignment for this analysis, ONLY submit jobs for
                        counting mutations in existing output folder
  -r1 R1                r1 SAM file
  -r2 R2                r2 SAM file
  -log LOG_LEVEL, --log_level LOG_LEVEL
                        set log level: debug, info, warning, error, critical.
                        (default = info)
  -env ENVIRONMENT, --environment ENVIRONMENT
                        The cluster used to run this script (default = GALEN)
  -at AT                Alignment time (default = 8h)
  -mt MT                Mutation call time (default = 36h)
  -mm MM                Mutation call request memory (default = 15GB)
  -c C                  Number of cores (default = 8)
  -b BASE, --base BASE  ASCII code base (default = 33)
  -test                 Turn on testing mode
  -rc                   Turn on rc mode, both direction of the reads will be
                        aligned to the reference. Variant calling will be
                        performed on all the reads that are aligned (BE
                        CAREFUL!)
  --sr_Override         Provide this argument when there is only one replicate
  --posteriorQC         Turn on posterior QC mode, this requires more memory
                        and runtime, please change the variable accordingly
  --wt_override         When no wt conditions defined in the parameter sheet,
                        turn on this option will treat EVERYTHING as wt (not
                        recommended). Phred scores will be adjusted based on
                        the first replicate
  --calibratePhredWT    When this parameter is provided, use wt to calibrate
                        phred scores
  --calibratePhredPhix  When this parameter is provided, use phix alignments
                        to calibrate phred scores
  --resubmit            For a finished run, batch resubmit failed scripts (if
                        any)
  --version VERSION     print version and exit
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
/190506_MTHFR_WT_2020-01-29-17-07-04/190506_MTHFR_WT_2020-01-29-17-07-04_mut_count/ --resubmit -n resub
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

`./env_jobs/` - Bash scripts for submitting the alignment jobs. If any of the sam file(s) didn't generate correctly, use `sbatch *.sh` to resubmit job for this sample

`./sam_files/` - Alignment output (SAM) and log files for bowtie2. Also contains alignment log file. Use this to check alignment rate.

`./name_time-stamped_mut_count/` - Mutation counts in each sample are saved in csv files
    
    - `./count_sample_*.csv` - Raw mutation counts for each sample. With meta data in header. Variants are represented in hgvs format
  
    - `./offreads_samplename.csv` - Contains information for all the reads that did not cover 90% of the tile

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

The count_sample_\*\*.csv is passed to tileseqMave for further analysis. If, for some reason, some of the count_sample_*.csv file did not generate correctly, you can resubmmit the job with the argument `--resubmit` (see example above)

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


### Other pacakges
---

The following pacakges are also installed with TileSeqMut. 

`posterior_QC`
* Process posterior probability files, generate posteriorQC output. It can only be used when the run was executed with `--posteriorQC`
* To run posterior_QC:

```
usage: posterior_QC [-h] -i INPUT -p PARAM

TileSeq mutation posterior QC, this will generate posterior_QC folder in
*_mut_count

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to folder that contains mutation counts
                        (*_mut_count)
  -p PARAM, --param PARAM
                        json paramter file
```

`random_ds`
* Random downsampling of sequencing data. Note that after downsampling, the output fastq files are saved in the same
 dir as the original fastq files. You will need to gzip them.
* To run random_ds:
```
usage: random_ds [-h] [-input INPUT] [-n N]

optional arguments:
  -h, --help    show this help message and exit
  -input INPUT  input folder contains original fastq files
  -n N          number of reads to downsample to
```

`makePRC`
* Before you run this, please install: [yogiroc](https://github.com/jweile/yogiroc)
* The script plots balanced Precision Recall Curve (PRC) given the score file (*_simple_aa.csv) and gene name. Plots
 are made using R package [yogiroc](https://github.com/jweile/yogiroc)
* It gets variants from public databases: Clinvar and gnomAD. The pathogenic/likely pathogenic variants from clinvar
 are used as pathogenic set. The benign/likely benign from clinvar and common variants from gnomAD (MAF > 1e-4) are
  used as benign set
* GRCh38 is used
* Currently only supports [VARITY](http://varity.varianteffect.org/) as a comparison
* To run makePRC:
```
usage: makePRC [-h] -s SCORES -g GENE [-c CLINVAR] [-o OUTPUT]
               [-r RANGE RANGE] [-v VARITY]

Make PRC curve using DMS scores

optional arguments:
  -h, --help            show this help message and exit
  -s SCORES, --scores SCORES
                        Input score file to make prc curve. score file ends
                        with _simple_aa.csv (required)
  -g GENE, --gene GENE  Gene symbol (required)
  -c CLINVAR, --clinvar CLINVAR
                        Path to Clinvar data, if not provided, the clinvar
                        database will be download to output dir, this might
                        take sometime.
  -o OUTPUT, --output OUTPUT
                        Output folder, if not specified, output plots will be
                        saved with score file
  -r RANGE RANGE, --range RANGE RANGE
                        Two integers to indicate the start/end of the targeted
                        region. If specified, only variants in this range will
                        be included. e.g -r 0 180 means variants in the range
                        of (0, 180] will be included for PRC curve
  -v VARITY, --varity VARITY
                        File contains hgvsp and VARITY scores (VARITY_R and
                        VARITY_ER) must be in columns.
```

`mergeRuns`
* Merge mut_counts file from two runs
* The script will ONLY add counts from samples with the same condition name, tile, time point and replicate
* The script output merged counts files in the output dir, also output a csv file contains the new file names
* To run mergeRuns
```
usage: mergeRuns [-h] [-p1 PARAMONE] [-p2 PARAMTWO] [-d1 DIR1] [-d2 DIR2] -o
                 OUTPUT [--covOverride] [--subtractWT] [-log LOG_LEVEL]

TileSeq mutation counts - Merge two runs

optional arguments:
  -h, --help            show this help message and exit
  -p1 PARAMONE, --paramOne PARAMONE
                        Path to parameter sheet (JSON) for the first run
  -p2 PARAMTWO, --paramTwo PARAMTWO
                        Path to parameter sheet (JSON) for the second run
  -d1 DIR1, --dir1 DIR1
                        Path to *_mut_count folder for the first run
  -d2 DIR2, --dir2 DIR2
                        Path to *_mut_count folder for the second run
  -o OUTPUT, --output OUTPUT
                        Output folder
  --covOverride         Ignore coverage files, please use this if there's no coverage_* file for either of the run
  --subtractWT          Subtract wt counts from non-select before merging
  -log LOG_LEVEL, --log_level LOG_LEVEL
                        set log level: debug, info, warning, error, critical.
                        (default = debug)
```
