"""
test
This is a snake pipeline for RNA-seq analysis of repetitive element expression. The Pipeline steps are as follows:

1. Initial QC - FastQC
2. Adapter/read trimming -  cutadapt
3. 2nd pass QC on trimmed data - FastQC
4. Read alignment - STAR
5. Read calling - telescope
6. Differential expression - DESeq2

REQUIREMENTS:
Raw data:
-Raw files must be in the format <sample name>_<replicate>_<read>.<fastq file extension>
	<sample name> - can be anything you want
	<replicate> - a single number (all replicaetes for a given sample must be sequential, i.e. 1,2,3, NOT 1,3,4,)
	<read> - R1 for read 1 and R2 for read 2
	<fastq file extension> - set in the config file. Must be gzipped fastq files (e.g., fastq.gz, fq.gz)!
-All samples must have the same number of replicates

Folder structure:
-Raw data should be located in folder called raw_files
-This file, RNAseq.Snakemake.cluster.config.yaml, and RNAseq.Snakemake.config.yaml must be in the working directory
-You need to create a new directory called output (sbatch outputs will go here)

TO RUN THE PIPELINE:
1. start a new tmux session with the command tmux new-session -s <session name>
2. Run the following command:
snakemake -s RNAseq.Snakefile -j 100 --configfile RNAseq.Snakemake.config.yaml --cluster-config RNAseq.Snakemake.cluster.config.yaml --cluster "sbatch -o {cluster.output} -e {cluster.err} -p {cluster.p} -N {cluster.N} -J {cluster.jobName} -t {cluster.time} --mail-user={cluster.mail-user} --mail-type={cluster.mail-type}"
"""


"""
TO DO:
- add functionality for single strand
- test that config file values are correct
- put pipeline and config files in a new folder in .../software
	- make the config files generic (and fail if not edited) and write-protect them 
"""

import numpy as np

#configfile: 'RNAseq_snakemake.config'
localrules: Telescope_DESeq

controlSample = config['control_sample']
replicates = list(range(1, int(config['num_replicates']) + 1))
working_dir = config['working_dir']
raw_file_ext = config['raw_file_extension']

#run test that config params are correct
#run test on file naming convention

SAMPLE_IDS = glob_wildcards(working_dir + 'raw_files/{sample}_{replicate}_{read}.' + raw_file_ext).sample

print("starting Snakemake...")
print("SAMPLE_IDS = " + str(SAMPLE_IDS))

#run test to confirm SAMPLE_IDS exist
assert len(SAMPLE_IDS) > 1, "< 1 sample found!"

rule all:
	input:
		expand(working_dir + 'FastQC/{sample}_{replicate}_{read}_fastqc.html', sample = SAMPLE_IDS, replicate = replicates, read = ['R1','R2']),
		expand(working_dir + 'FastQC_2/{sample}_{replicate}_{read}.cutadapt.q20.minlen1_fastqc.html', sample = SAMPLE_IDS, replicate = replicates, read = ['R1','R2']),
		# expand(working_dir + 'cutadapt/{sample}_{replicate}_{read}.cutadapt.q20.minlen1.' + raw_file_ext, sample = SAMPLE_IDS, replicate = replicates, read = ['R1', 'R2']),
		# expand(working_dir + 'RNAseq.STAR/RNAseq.STAR.{file}.Aligned.out.bam', file = SAMPLE_IDS),
		# expand(working_dir + 'telescope/{sample}_{replicate}-telescope_report.tsv', sample = SAMPLE_IDS, replicate = replicates),
		expand(working_dir + 'telescope/{sample}.telescope.count.table.DESeq2.tsv', sample = SAMPLE_IDS),
		expand(working_dir + 'TEtranscripts/{sample}.TEtranscripts.DESeq_gene_TE_analysis.txt', sample = SAMPLE_IDS)


rule FastQC:
	input:
		working_dir + 'raw_files/{sample}_{replicate}_{read}.' + raw_file_ext,

	output:
		working_dir + 'FastQC/{sample}_{replicate}_{read}_fastqc.html',
		working_dir + 'FastQC/{sample}_{replicate}_{read}_fastqc.zip'
	shell:
		'''
		ml fastQC/0.11.5 &&
				fastqc -o FastQC {input}
		'''
rule CutAdapt:
	input:
		read1 = working_dir + 'raw_files/{sample}_{replicate}_R1.' + raw_file_ext,
		read2 = working_dir + 'raw_files/{sample}_{replicate}_R2.' + raw_file_ext
	output:
		read1 = working_dir + 'cutadapt/{sample}_{replicate}_R1.cutadapt.q20.minlen1.' + raw_file_ext,
		read2 = working_dir + 'cutadapt/{sample}_{replicate}_R2.cutadapt.q20.minlen1.' + raw_file_ext
	shell:
		'''
		ml python/2.7.6

		cutadapt -a {config[fw_adapter]} -A {config[rev_adapter]} -q 20 --minimum-length 1 -o {output.read1} -p {output.read2} {input.read1} {input.read2} > cutadapt/{wildcards.sample}_{wildcards.replicate}.cutadapt.report.txt
		'''
rule FastQC_pass2:
	input:
	  working_dir + 'cutadapt/{sample}_{replicate}_{read}.cutadapt.q20.minlen1.' + raw_file_ext

	output:
	  working_dir + 'FastQC_2/{sample}_{replicate}_{read}.cutadapt.q20.minlen1_fastqc.html',
	  working_dir + 'FastQC_2/{sample}_{replicate}_{read}.cutadapt.q20.minlen1_fastqc.zip'
	shell:
		'''
		ml fastQC/0.11.5 &&
			fastqc -o FastQC_2 {input}
		'''

rule STAR:
	input: 
		read1 = working_dir + 'cutadapt/{sample}_{replicate}_R1.cutadapt.q20.minlen1.' + raw_file_ext,
		read2 = working_dir + 'cutadapt/{sample}_{replicate}_R2.cutadapt.q20.minlen1.' + raw_file_ext
	output:
		working_dir + 'RNAseq.STAR/RNAseq.STAR.{sample}_{replicate}.Aligned.out.bam'
	shell:
		'''
		ml star/2.6
		STAR --runThreadN 16 --genomeDir /groups/chiappinellilab/genomes/hg38/STAR.hg38.index/ \
		--sjdbGTFfile /groups/chiappinellilab/genomes/hg38/gencode.annotation/gencode.v21.primary.assembly.only.annotation.gtf \
		--sjdbOverhang 100 \
		--readFilesIn {input.read1} {input.read2} \
		--readFilesCommand zcat --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 --outFileNamePrefix RNAseq.STAR/RNAseq.STAR.{wildcards.sample}_{wildcards.replicate}.
		'''


rule TEtranscripts:
	input: 
		treatment_files = expand(working_dir + 'RNAseq.STAR/RNAseq.STAR.{{sample}}_{replicate}.Aligned.out.bam', replicate = replicates),
		control_files = expand(working_dir + 'RNAseq.STAR/RNAseq.STAR.' + controlSample + '_{replicate}.Aligned.out.bam', replicate = replicates)
	output: working_dir + 'TEtranscripts/{sample}.TEtranscripts.DESeq_gene_TE_analysis.txt'
	shell:
		'''
		ml python/2.7.6
		ml R/3.4.2
		ml gcc/8.1.0
		ml xz/5.2.3
		TEtranscripts --format BAM --mode multi --stranded reverse -t {input.treatment_files} -c {input.control_files} \
		--GTF /groups/chiappinellilab/genomes/hg38/gencode.annotation/gencode.v21.primary.assembly.only.annotation.gtf \
		--TE /groups/chiappinellilab/genomes/hg38/RepeatMasker/rmsk.hg38.sorted.BOTHfeature.gtf \
		--project TEtranscripts/{wildcards.sample}.TEtranscripts.DESeq

		'''

rule Telescope:
	input: working_dir + 'RNAseq.STAR/RNAseq.STAR.{sample}_{replicate}.Aligned.out.bam'
	output: working_dir + 'telescope/{sample}_{replicate}-telescope_report.tsv'
	shell:
		'''
		ml miniconda
		source activate telescope_env

		telescope assign {input} /groups/chiappinellilab/genomes/hg38/HERV_L1_rmsk.hg38.gtf \
		--outdir telescope/ --exp_tag {wildcards.sample}_{wildcards.replicate}

		'''

rule Telescope_DESeq:
	input:
		treatment_reports = expand(working_dir + 'telescope/{{sample}}_{replicate}-telescope_report.tsv', replicate = replicates),
		control_reports = expand(working_dir + 'telescope/' + controlSample + '_{replicate}-telescope_report.tsv', replicate = replicates),
		script = '/groups/chiappinellilab/software/tkanholm/make.TCGA.telescope.DESeq2.input.filter.baseMean.10.py',
		annotation = '/groups/chiappinellilab/genomes/hg38/HERV_L1_rmsk.hg38.gtf'
	output:
		treat_files_list = temp(working_dir + 'telescope/{sample}.treat_files.txt'),
		cntrl_files_list = temp(working_dir + 'telescope/{sample}.cntrl_files.txt'),
		cntTable = working_dir + 'telescope/{sample}.telescope.count.table.tsv',
		DESeq2_script = working_dir + 'telescope/{sample}.telescope.count.table.DESeq2.Rscript.R',
		DESeq2_output = working_dir + 'telescope/{sample}.telescope.count.table.DESeq2.tsv'
	params:
		cntrl_sample = controlSample,
		workingDir = working_dir
	shell:
		'''
		ml python
		ml R/3.4.2
		ml gcc/8.1.0

		ls {params.workingDir}telescope/*{wildcards.sample}*telescope_report.tsv > {output.treat_files_list}
		ls {params.workingDir}telescope/*{params.cntrl_sample}*telescope_report.tsv > {output.cntrl_files_list}
		
		#  make the DESeq count table and generate a DESeq script
		python {input.script} {output.treat_files_list} {output.cntrl_files_list} {input.annotation} {params.workingDir}telescope/ -o {wildcards.sample}.telescope.count.table
		# run the generated script
		Rscript {output.DESeq2_script}

		'''

rule Combine_tables:
	input: expand(working_dir + 'telescope/{sample}.telescope.count.table.DESeq2.tsv', sample = SAMPLE_IDS)
	output: working_dir + 'all.samples.telescope.DESeq2.tibble.tsv'
	shell:
		'''
		ml python
		python /groups/chiappinellilab/software/jimcdonald/make.TCGA.telescope.DESeq2.tibble.py /lustre/groups/chiappinellilab/processed_data/melbourneBrainTumors/telescope/ -n all.samples.telescope. -x .count.table.DESeq2.tsv -a /groups/chiappinellilab/genomes/hg38/HERV_L1_rmsk.hg38.gtf  -p TRUE -s /lustre/groups/chiappinellilab/processed_data/melbourneBrainTumors/telescope/sample_table.txt
		'''
