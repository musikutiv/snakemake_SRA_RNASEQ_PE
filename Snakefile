#
configfile: "config.yaml"

import pandas as pd

sample_table = pd.read_table(config["runtable"], sep=",")
samples = list(sample_table['Sample Name'].unique())

# everything that need network connection for
localrules: prefetch, get_index_files

rule all:
    input:
        expand("{SRS}_out/{SRS}_rsem.genes.results", SRS=samples)

rule STAR:
    input:
        files1=lambda wildcards: expand("raw/{sample}_1.fastq.gz",
            sample=list(sample_table[sample_table['Sample Name']==wildcards.SRS].Run)),
        files2=lambda wildcards: expand("raw/{sample}_2.fastq.gz",
            sample=list(sample_table[sample_table['Sample Name']==wildcards.SRS].Run)),
        index="star_index/Genome"
    output:
        temp("{SRS}_out/{SRS}.Aligned.toTranscriptome.out.bam"),
        temp("{SRS}_out/{SRS}.Aligned.out.bam"),

    threads: 16

    params:
        files1=lambda wildcards: ",".join(expand("raw/{sample}_1.fastq.gz",
            sample=list(sample_table[sample_table['Sample Name']==wildcards.SRS].Run))),
        files2=lambda wildcards: ",".join(expand("raw/{sample}_2.fastq.gz",
            sample=list(sample_table[sample_table['Sample Name']==wildcards.SRS].Run))),
        quantMode="TranscriptomeSAM GeneCounts",
        outSAMtype="BAM Unsorted",
        limitBAMsortRAM="40000000000",
        readFilesCommand="gunzip -c",
        SRS=lambda wildcards: wildcards.SRS
    shell:
        """
        STAR --genomeDir star_index \
        --sjdbGTFfile genome.gtf \
        --runThreadN {threads} \
        --readFilesCommand {params.readFilesCommand} \
        --quantMode {params.quantMode} --outSAMtype {params.outSAMtype} \
        --limitBAMsortRAM {params.limitBAMsortRAM} \
        --readFilesIn {params.files1} {params.files2} \
        --outFileNamePrefix {params.SRS}_out/{params.SRS}.
        """

rule RSEM:
    input:
        "{SRS}_out/{SRS}.Aligned.toTranscriptome.out.bam",
        index="star_index/Genome"

    output:
        "{SRS}_out/{SRS}_rsem.genes.results",
        temp("{SRS}_out/{SRS}_rsem.transcript.bam"),

    threads: 16

    params:
        SRS=lambda wildcards: wildcards.SRS
    shell:
        """
        rsem-calculate-expression --bam --paired-end --strandedness none -p {threads} {params.SRS}_out/{params.SRS}.Aligned.toTranscriptome.out.bam star_index/genome {params.SRS}_out/{params.SRS}_rsem
        """

rule fastq_dump:
    input:
        "sra/{SRR}.sra"
    output:
        "raw/{SRR}_1.fastq.gz",
        "raw/{SRR}_2.fastq.gz"
    shell:
        """
        fastq-dump --split-3 --gzip {input} --outdir raw
        """

rule prefetch:
    output:
        "sra/{SRR}.sra"
    shell:
        "prefetch -O sra {wildcards.SRR}"

rule star_index:
	input:
		fasta="genome.fa",
		gtf="genome.gtf"
	output:
		"star_index/Genome",
		"star_index/genome.transcripts.fa",
		"star_index/chrNameLength.txt"
	threads: 30
	shell:
		"""
		rsem-prepare-reference --gtf {input.gtf} --star -p {threads} {input.fasta} star_index/genome
		"""

rule get_index_files:
    params:
        genomeLocation=config['genomeFTP'],
        gtfLocation=config['gtfFTP']
    output:
        "genome.fa",
        "genome.gtf"
    shell:
        """
        wget -O genome.fa.gz {params.genomeLocation}
        gunzip genome.fa.gz
        wget -O genome.gtf.gz {params.gtfLocation}
        gunzip genome.gtf.gz
        """

onsuccess:
        shell("rm *.out")
