configfile: "config.yaml"

workdir: config['work_dir']      ## Configure Working Directory

import os
import pandas as pd

## for the sake of speed, fastqc do it outside the snakemake pipe, todo
# v2 change to bowtie2 
###########################################
## paths for pipeline and/or reference data
work_dir = config["work_dir"]
species = config['species']
binsize = 50 if species != "yeast" else 1
##################################################q
## read in sample and corresponding fq files talbe
## and aggregatio talbe
SAMPLES = (
    pd.read_csv(config["samples"], sep="\t")
    .set_index("sample_id", drop=False)
    .sort_index()
)

effectiveSize = 2862010428
###############################
## read in refrence files' info
REF = pd.read_csv(config["ref_files"], sep="\t", header = None, index_col = 0)
# blacklist = REF.loc["blacklist"][1]   ## ENCODE blacklist

threads = config["threads"]
mapq_cutoff = config['MAPQ']

##  get corresponding bwa_index
def get_bwa_index():
    return REF.loc["bwa_index"][1]

def get_bowtie2_index():
    return REF.loc["bowtie2_index"][1]

def get_raw_fastq_se(wildcards):
    if config["paired-end"] == False:
        return SAMPLES.loc[wildcards.sample]["R1"].split(",")
    else:
        return ""

def get_raw_fastq_pe_R1(wildcards):
    if config["paired-end"]:
        return SAMPLES.loc[wildcards.sample]["R1"].split(",")
    else:
        return ""

def get_raw_fastq_pe_R2(wildcards):
    if config["paired-end"]:
        return SAMPLES.loc[wildcards.sample]["R2"].split(",")
    else:
        return ""

def get_renamed_fastq(wildcards):
    if config["paired-end"]:
        R1 = "renamed_fq/{}_R1.fastq.gz".format(wildcards.sample)
        R2 = "renamed_fq/{}_R2.fastq.gz".format(wildcards.sample)
        return R1 + R2
    else:
        return "{}/renamed_fq/{}.fastq.gz".format(work_dir,wildcards.sample)

def get_fastq_4trim(wildcards):
    if config["paired-end"]:
        R1 = "renamed_fq/{}_R1.fastq.gz".format(wildcards.sample),
        R2 = "renamed_fq/{}_R2.fastq.gz".format(wildcards.sample),
        return R1 + R2
    else:
        return "renamed_fq/{}.fastq.gz".format(wildcards.sample)

##################################
## get trimmed fastq files for BWA
def get_trimmed_fastq(wildcards):
    if config["paired-end"]:
        R1 = "trimmed_fq/{}_R1_val_1.fq.gz".format(wildcards.sample),
        R2 = "trimmed_fq/{}_R2_val_2.fq.gz".format(wildcards.sample),
        return R1 + R2
    else:
        return "trimmed_fq/{}_trimmed.fq.gz".format(wildcards.sample)

########################################
## get dedup bam files for meth_qc_quant
def get_dedup_bam(wildcards):
    if config["paired-end"]:
        return "dedup_bam_pe/{}_dedup.bam".format(wildcards.sample)
    else:
        return "dedup_bam_se/{}_dedup.bam".format(wildcards.sample)

def get_fastqc_stats(wildcards):
    if config["paired-end"]:
        r1_raw  = expand("fastqc_pe/{samples}_R1_fastqc.zip", samples = wildcards.sample),
        r2_raw  = expand("fastqc_pe/{samples}_R2_fastqc.zip", samples = wildcards.sample),
        r1_trim = expand("fastqc_pe/{samples}_R1_val_1_fastqc.zip", samples = wildcards.sample),
        r2_trim = expand("fastqc_pe/{samples}_R2_val_2_fastqc.zip", samples = wildcards.sample),
        return r1_raw + r2_raw + r1_trim + r2_trim
    else:
         r1_raw  = expand("fastqc_se/{samples}_fastqc.zip", samples = wildcards.sample),
         r1_trim = expand("fastqc_se/{samples}_trimmed_fastqc.zip", samples = wildcards.sample),
         return r1_raw + r1_trim

def get_dedup_bam_stats(wildcards):
    if config["paired-end"] and config["add_umi"]:
        return expand("dedup_bam_umi_pe/{samples}_dedup.bam.stats.txt", samples = wildcards.sample)
    elif config["paired-end"] and config["add_umi"] == False:
        return expand("dedup_bam_pe/{samples}_dedup.bam.stats.txt", samples = wildcards.sample)
    elif config["paired-end"] == False and config["add_umi"]:
        return expand("dedup_bam_umi_se/{samples}_dedup.bam.stats.txt", samples = wildcards.sample)
    else:
        return expand("dedup_bam_se/{samples}_dedup.bam.stats.txt", samples = wildcards.sample)

def get_rule_all_input():
    bam_out = expand("dedup_bam_pe/{sample}_dedup.bam", sample=SAMPLES["sample_id"].values.tolist())
    frag_out = expand("fragment_size/{sample}_insert_size_metrics.txt",sample=SAMPLES["sample_id"].values.tolist())
    bw_out1 = expand("bigWig/{sample}_dedup.CPM.bw", sample=SAMPLES["sample_id"].values.tolist())
    bw_out2 = expand("bigWig/{sample}_dedup.watson."+str(binsize)+".bw", sample=SAMPLES["sample_id"].values.tolist())
    bw_out3 = expand("bigWig/{sample}_dedup.crick."+str(binsize)+".bw", sample=SAMPLES["sample_id"].values.tolist())
    return bam_out+frag_out+bw_out1+bw_out2+bw_out3


rule all:
    input:
        get_rule_all_input()


## paired-end
rule merge_and_rename_fq_pe:
    input:
        R1 = get_raw_fastq_pe_R1,
        R2 = get_raw_fastq_pe_R2,
    output:
        "renamed_fq/{sample}_R1.fastq.gz",
        "renamed_fq/{sample}_R2.fastq.gz",
    conda:
        "envs/cutandrun_smk.yaml"
    shell:
        "ln -s {input.R1} {output[0]} && "
        "ln -s {input.R2} {output[1]} "

rule trim_galore_pe:
    input:
        get_fastq_4trim
    output:
        temp("trimmed_fq/{sample}_R1_val_1.fq.gz"),
        temp("trimmed_fq/{sample}_R2_val_2.fq.gz"),
        "trimmed_fq/{sample}_R1.fastq.gz_trimming_report.txt",
        "trimmed_fq/{sample}_R2.fastq.gz_trimming_report.txt",
    params:
        ## path needs to be full path: failed to recognize space
        path = work_dir + "/trimmed_fq"
    threads: threads
    log:
        "logs/{sample}_trim_galore_pe.log"
    conda:
        "envs/cutandrun_smk.yaml"
    shell:
        "(trim_galore -q {mapq_cutoff} --stringency 3 --length 20 "
        "--cores {threads} --paired -o {params.path} {input}) 2> {log}"

################
## BWA alignment
################
rule bwa_map:
    input:
        get_bwa_index(),
        get_trimmed_fastq
    output:
        temp_bam=temp("raw_bam/{sample}.bam"),
        stats="raw_bam/{sample}_bam.stats.txt"
    threads: threads
    resources: mem_mb=80000
    log:
        log_bwa="logs/{sample}_bwa_map.log",
        log_stats="logs/{sample}_stats.log"
    conda:
        "envs/cutandrun_smk.yaml"
    shell:
        """
        (bwa mem -M -K 100000000 -t {threads} {input} | \
        samtools view -Sb --threads {threads} - > {output.temp_bam}) 2> {log.log_bwa}
        samtools stats -@ {threads} {output.temp_bam} > {output.stats}
        """
##########################################
## raw bams without any filtering
## fixmate, sort, index and stats bam file
rule samtools_sort_index_stats:
    input:
        "raw_bam/{sample}.bam"
    output:
        #temp(bam = "raw_bam/{sample}_sorted.bam"), doesn't work
        #bai = "raw_bam/{sample}_sorted.bam.bai",
        #stat= "raw_bam/{sample}_sorted.bam.stats.txt"
        "raw_bam/{sample}_sorted.bam",
        "raw_bam/{sample}_sorted.bam.stats.txt"
    threads: threads
    resources: mem_mb=80000
    conda:
        "envs/cutandrun_smk.yaml"
    shell:
        ## --threads flag failed
        "(samtools fixmate -@ {threads} -m {input} - | "
        "samtools sort  -@ {threads} -o {output[0]} && "
        "samtools index -@ {threads} {output[0]} && "
        "samtools stats -@ {threads} {output[0]} > {output[1]})"

##########################################################################
## to filter out unmapped & non-uniquely mapped, not properly paired reads
## Deduplication with markup, index and stats deduplicated file
###########################################################################
# 2828 unique map, 524 allow mutli-mapping (h 3 )
rule samtools_markdup_stats_pe:
    input:
        "raw_bam/{sample}_sorted.bam"
    output:
        bam = "dedup_bam_pe/{sample}_dedup.bam",
        #bai = "dedup_bam/{sample}_dedup.bam.bai",
        stat= "dedup_bam_pe/{sample}_dedup.bam.stats.txt"
    threads: threads
    conda:
        "envs/cutandrun_smk.yaml"
    shell:
        "(samtools view -b -F 524 --threads {threads} {input} | "
        "samtools markdup -@ {threads} -r - {output.bam} && "
        "samtools index -@ {threads} {output.bam} && "
        "samtools stats -@ {threads} {output.bam} > {output.stat})"

rule generate_dedup_bw:
    input:
        "dedup_bam_pe/{sample}_dedup.bam"
    output:
        "bigWig/{sample}_dedup.CPM.bw"
    threads: threads
    params:
        eSize = effectiveSize
    conda:
        "envs/cutandrun_smk.yaml"
    shell:
        "bamCoverage --bam {input} -o {output} --binSize 1 --normalizeUsing CPM --extendReads -p {threads}"

rule insert_size:
    input:
        #"dedup_bam_pe/{sample}_dedup.bam"
        get_dedup_bam
    output:
        txt = "fragment_size/{sample}_insert_size_metrics.txt",
        hist = "fragment_size/{sample}_insert_size_histogram.pdf"
    log:
        "logs/{sample}_picard_insert_size.log"
    conda:
        "envs/cutandrun_smk.yaml"
    shell:
        "(picard CollectInsertSizeMetrics M=0.05 I={input} O={output.txt} "
        "H={output.hist}) 2> {log}"

rule plus_minus_bam:
    input:
        "dedup_bam_pe/{sample}_dedup.bam"
    output:
        plus_bam = "dedup_bam_pe/{sample}_dedup.watson.bam",
        minus_bam = "dedup_bam_pe/{sample}_dedup.crick.bam",
    threads: threads
    conda:
        "envs/cutandrun_smk.yaml"
    shell:
        """
        # Split BAM into plus and minus strand based on flag values
        samtools view -h {input} | tee >(awk '$2 == 99 || $2 == 147 || $1 ~ /^@/' | samtools view -Sb -@ {threads} - > {output.plus_bam}) | awk '$2 == 83 || $2 == 163 || $1 ~ /^@/' | samtools view -Sb -@ {threads} - > {output.minus_bam}
        
        # Index
        samtools index -@ {threads} {output.plus_bam}
        samtools index -@ {threads} {output.minus_bam}
        """

rule plus_minus_bw:
    input:
        plus_bam = "dedup_bam_pe/{sample}_dedup.watson.bam",
        minus_bam = "dedup_bam_pe/{sample}_dedup.crick.bam"
    output:
        plus_bw = "bigWig/{sample}_dedup.watson."+str(binsize)+".bw",
        minus_bw = "bigWig/{sample}_dedup.crick."+str(binsize)+".bw"
    params:
        bin_size = binsize
    threads: threads
    conda:
        "envs/cutandrun_smk.yaml"
    shell:
        """
        # Convert BAM to BigWig for plus strand
        bamCoverage -b {input.plus_bam} -o {output.plus_bw} --binSize {params.bin_size} --extendReads --numberOfProcessors {threads}
        # Convert BAM to BigWig for minus strand
        bamCoverage -b {input.minus_bam} -o {output.minus_bw} --binSize {params.bin_size} --extendReads --numberOfProcessors {threads}
        """