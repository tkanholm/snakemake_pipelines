"""
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
-You need to create a new directory called outputs (sbatch outputs will go here)
-All other subfolders will be created automatically

TO RUN THE PIPELINE:
1. start a new tmux session with the command tmux new-session -s <session name>
2. Run the following command:
snakemake -s RNAseq.Snakefile -j 100 --configfile RNAseq.Snakemake.config.yaml --cluster-config RNAseq.Snakemake.cluster.config.yaml --cluster "sbatch -o {cluster.output} -e {cluster.err} -p {cluster.p} -N {cluster.N} -J {cluster.jobName} -t {cluster.time} --mail-user={cluster.mail-user} --mail-type={cluster.mail-type}"
3. Snakemake will create two sample_table.txt files located in telescope/ and TEtranscripts/. Add columns to these files with
	any extra information about each sample. These will appear as extra columns in the final all.samples.DESeq2.tibble.tsv
	files.
"""


"""
TO DO:
- test that config file values are correct
- make output folder automagically
- add functionality for single strand
	- can this be done with a flag, or do we need a whole new pipeline?
"""

import numpy as np

#configfile: 'RNAseq_snakemake.config'
localrules: Telescope_DESeq, make_sample_tables, Combine_tables_telescope, Combine_tables_TEtranscripts

controlSample = config['control_sample']
replicates = list(range(1, int(config['num_replicates']) + 1))
working_dir = config['working_dir']
raw_file_ext = config['raw_file_extension']

#run test that config params are correct
#run test on file naming convention

SAMPLE_IDS = glob_wildcards(working_dir + 'raw_files/{sample}_{replicate}_{read}.' + raw_file_ext).sample
# this will grab multiple of the same SAMPLE_ID (one for each replicate/read combo). Grab just the unique ones
# convert to a set (only takes unique values)
SAMPLE_ID_set = set(SAMPLE_IDS)
# convert back to a list
SAMPLE_IDS = list(SAMPLE_ID_set)

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
		# expand(working_dir + 'telescope/{sample}.telescope.count.table.DESeq2.tsv', sample = SAMPLE_IDS),
		# expand(working_dir + 'TEtranscripts/{sample}.TEtranscripts.DESeq_gene_TE_analysis.txt', sample = SAMPLE_IDS)
		working_dir + 'all.samples.telescope.DESeq2.tibble.tsv',
		working_dir + 'all.samples.TEtranscripts.DESeq2.tibble.tsv'


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
		python {input.script} {output.cntrl_files_list} {output.treat_files_list} {input.annotation} {params.workingDir}telescope/ -o {wildcards.sample}.telescope.count.table
		# run the generated script
		Rscript {output.DESeq2_script}

		'''
rule make_sample_tables:
	input:
	output:
		telescope_table = working_dir + 'telescope/sample_table.txt',
		TEtranscripts_table = working_dir + 'TEtranscripts/sample_table.txt'
	run:
		## telescope sample table ##
		# table header
		table_str = "DESeq_output_file" + "\t" + "Sample_name" + "\n"

		for sampleName in SAMPLE_IDS:
			table_str = table_str + sampleName + ".telescope.count.table.DESeq2.tsv" "\t" + sampleName + "\n"

		print("Generated sample_table for telescope. Add columns to this file with optional sample data to be added to the final " +
			"all.samples.telescope.DESeq2.tibble.tsv output")
		print(table_str)
		sample_table_file = open(output.telescope_table, "w")
		sample_table_file.write(table_str)

		## TEtranscripts sample table ##
		# table header
		table_str = "DESeq_output_file" + "\t" + "Sample_name" + "\n"

		for sampleName in SAMPLE_IDS:
			table_str = table_str + sampleName + ".TEtranscripts.DESeq_gene_TE_analysis.txt" "\t" + sampleName + "\n"

		print("Generated sample_table for TEtranscripts. Add columns to this file with optional sample data to be added to the final " +
			"all.samples.telescope.DESeq2.tibble.tsv output")
		print(table_str)
		sample_table_file = open(output.TEtranscripts_table, "w")
		sample_table_file.write(table_str)



rule Combine_tables_telescope:
	input: 
		telescope_files = expand(working_dir + 'telescope/{sample}.telescope.count.table.DESeq2.tsv', sample = SAMPLE_IDS),
		sample_table = working_dir + 'telescope/sample_table.txt',
		script = '/groups/chiappinellilab/software/jimcdonald/make.TCGA.telescope.DESeq2.tibble.py'
	output: working_dir + 'all.samples.telescope.DESeq2.tibble.tsv'
	params:
		workingDir = working_dir
	shell:
		'''
		ml python
		python {input.script} \
		{params.workingDir}telescope/ \
		-n all.samples.telescope -x .count.table.DESeq2.tsv -a /groups/chiappinellilab/genomes/hg38/HERV_L1_rmsk.hg38.gtf  \
		-p TRUE -s {input.sample_table}
		'''
rule Combine_tables_TEtranscripts:
	input: 
		TEtranscripts_files = expand(working_dir + 'TEtranscripts/{sample}.TEtranscripts.DESeq_gene_TE_analysis.txt', sample = SAMPLE_IDS),
		sample_table = working_dir + 'TEtranscripts/sample_table.txt',
		script = '/groups/chiappinellilab/software/tkanholm/make.TCGA.telescope.DESeq2.tibble.py'
	output: working_dir + 'all.samples.TEtranscripts.DESeq2.tibble.tsv'
	params:
		workingDir = working_dir
	shell:
		'''
		ml python
		python {input.script} \
		{params.workingDir}TEtranscripts/ \
		-n all.samples.TEtranscripts -x .DESeq_gene_TE_analysis.txt   \
		-p TRUE -s {input.sample_table}
		'''
