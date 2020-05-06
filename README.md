# snakemake_pipelines
This is a snake pipeline for RNA-seq analysis of repetitive element expression. The Pipeline steps are as follows:

1. Initial QC - FastQC
2. Adapter/read trimming -  cutadapt
3. 2nd pass QC on trimmed data - FastQC
4. Read alignment - STAR
5. Read calling - telescope and TEtranscripts
6. Differential expression - DESeq2
7. Combine all samples into two final output tables -- one for telescope and one for TEtranscripts

## Requirements:
Raw data:
- Raw files must be in the format `<sample name>_<replicate>_<read>.<fastq file extension>`
  -	`<sample name>` - can be anything you want
  -	`<replicate>` - a single number (all replicaetes for a given sample must be sequential, i.e. 1,2,3, NOT 1,3,4,)
  -	`<read>` - R1 for read 1 and R2 for read 2
  -	`<fastq file extension>` - set in the config file. Must be gzipped fastq files (e.g., fastq.gz, fq.gz)!
- All samples must have the same number of replicates

## Folder structure:
- Raw data should be located in folder called `raw_files`
- `RNAseq.Snakefile`, `RNAseq.Snakemake.cluster.config.yaml`, and `RNAseq.Snakemake.config.yaml` must be in the working directory
- You need to create a new directory called `outputs` (sbatch outputs will go here)
- All other subfolders will be created automatically

## To run the pipeline:
1. start a new tmux session with the command `tmux new-session -s <session name>`
2. Run the following command:
```
snakemake -s RNAseq.Snakefile -j 100 --configfile RNAseq.Snakemake.config.yaml --cluster-config RNAseq.Snakemake.cluster.config.yaml --cluster "sbatch -o {cluster.output} -e {cluster.err} -p {cluster.p} -N {cluster.N} -J {cluster.jobName} -t {cluster.time} --mail-user={cluster.mail-user} --mail-type={cluster.mail-type}"
```
3. Snakemake will create two `sample_table.txt` files located in `telescope/` and `TEtranscripts/`. Add columns to these files with any extra information about each sample. These will appear as extra columns in the final `all.samples.DESeq2.tibble.tsv files`.
