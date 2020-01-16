## TileSeq mutation count package

This package is made to parse input sequecning files (fastq) with user provided parameter.csv file.
Output of this pipeline is mutation counts for each pair of fastq files. 

### Execution
---

To run this pipeline on GURU (SGE)


To run this pipeline on GALEN (Slurm)


### Input files
---

`/path/to/fastq/` - Full path to input fastq files 

`param.csv` - CSV file contains information for this run (please see example
[here]:https://docs.google.com/spreadsheets/d/1tIblmIFgOApPNzWN2KUwj8BKzBiJ1pOL7R4AOUGrqvE/edit?usp=sharing ).
This file is required to be comma-seperated and saved in csv format. 


### Output files
---

One ouptut folder is created for each run. The output folder are named with `project_name-time-stamp`

Within each output folder, the following files and folders:

`main.log` - main logging output 

`ref` - Reference fasta file and bowtie2 index

`ds_fastq` - Raw fastq files are downsampled to `n` reads and the downsampled fastq files are saved
in this folder

`sam_files` - Alignemnt output and log files for the raw fastq files

`ds_sam_files` - Alignment output and log files for the down-sampled fastq files


### Alignment
---

The pipeline takes the sequence in the parameter file as reference and try to align the fastq files
to the whole reference sequence. 
