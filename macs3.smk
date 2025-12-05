# Snakemake workflow for CUT&RUN MACS3 peak calling
# Handles WT, POLE3KO, and POLE4_KO samples with biological replicates and technical replicates

import os
import re
import glob
from pathlib import Path

# Define broad and narrow histone marks
BROAD_MARKS = ["H3K27me3", "H3K9me3"]
NARROW_MARKS = ["H3K4me3"]

# Define input directory structure
INPUT_DIR = "dedup_bam_pe"  # Directory containing input BAM files
RESULTS_DIR = "results/macs3"  # Directory for results

# Get all sample BAM files
SAMPLE_BAMS = glob.glob(f"{INPUT_DIR}/*_dedup.bam")
SAMPLES = [os.path.basename(bam).replace("_dedup.bam", "") for bam in SAMPLE_BAMS]

# Parse sample information
def parse_sample_name(sample):
    # Expected format examples:
    # WT_1_H3K27me3_CR_2_dedup (biological sample 1, technical replicate 2)
    # POLE3KO_2_H3K4me3_CR_1_dedup (biological sample 2, technical replicate 1)
    # POLE4_KO_1_H3K9me3_CR_1_dedup (biological sample 1, technical replicate 1)
    
    components = sample.split("_")
    
    # Identify histone mark
    for component in components:
        if component in BROAD_MARKS or component in NARROW_MARKS:
            histone_mark = component
            break
    else:
        histone_mark = "unknown"
    
    # Extract condition and biological sample identifier
    if "POLE4_KO" in sample:
        # Extract biological sample number for POLE4_KO
        match = re.search(r"POLE4_KO_(\d+)", sample)
        if match:
            bio_sample = match.group(1)
            condition = f"POLE4_KO_{bio_sample}"
        else:
            condition = "POLE4_KO"
    elif "POLE3KO" in sample:
        # Extract biological sample number for POLE3KO
        match = re.search(r"POLE3KO_(\d+)", sample)
        if match:
            bio_sample = match.group(1)
            condition = f"POLE3KO_{bio_sample}"
        else:
            condition = "POLE3KO"
    else:
        # Extract biological sample number for WT
        match = re.search(r"WT_(\d+)", sample)
        if match:
            bio_sample = match.group(1)
            condition = f"WT_{bio_sample}"
        else:
            condition = "WT"
    
    # Find technical replicate number (CR_X)
    rep_match = re.search(r"_CR_(\d+)_", sample)
    if rep_match:
        replicate = rep_match.group(1)
    else:
        replicate = "1"
    
    return {
        "sample": sample,
        "condition": condition,  # Now includes biological sample ID
        "histone_mark": histone_mark,
        "replicate": replicate   # Technical replicate
    }

# Parse all samples
SAMPLE_INFO = [parse_sample_name(sample) for sample in SAMPLES]

# Group samples by condition and histone mark to identify replicates for merging
CONDITION_MARK_SAMPLES = {}
for info in SAMPLE_INFO:
    key = (info["condition"], info["histone_mark"])
    if key not in CONDITION_MARK_SAMPLES:
        CONDITION_MARK_SAMPLES[key] = []
    CONDITION_MARK_SAMPLES[key].append(info["sample"])

# Define the final target files
# Define the final target files
def get_broad_peak_files():
    files = []
    for (condition, mark), samples in CONDITION_MARK_SAMPLES.items():
        if mark in BROAD_MARKS and len(samples) > 0:
            # Individual sample peak files
            for sample in samples:
                files.append(f"{RESULTS_DIR}/{sample}_peaks.broadPeak")
            
            # Common peaks if there are technical replicates
            if len(samples) > 1:
                files.append(f"{RESULTS_DIR}/{condition}_{mark}_common_peaks.broadPeak")
    return files

def get_narrow_peak_files():
    files = []
    for (condition, mark), samples in CONDITION_MARK_SAMPLES.items():
        if mark in NARROW_MARKS and len(samples) > 0:
            # Individual sample peak files
            for sample in samples:
                files.append(f"{RESULTS_DIR}/{sample}_peaks.narrowPeak")
            
            # Common peaks if there are technical replicates
            if len(samples) > 1:
                files.append(f"{RESULTS_DIR}/{condition}_{mark}_common_peaks.narrowPeak")
    return files



# Define the target rule
rule all:
    input:
        get_broad_peak_files(),
        get_narrow_peak_files()

# MACS3 callpeak for broad histone marks - MODIFIED: corrected output filename
rule epic2_callpeak_broad:
    input:
        bam = INPUT_DIR + "/{sample}_dedup.bam"
    output:
        peaks = RESULTS_DIR + "/{sample}_peaks_ep.broadPeak"
    params:
        outdir = RESULTS_DIR,
        name = "{sample}"
    conda:
        "./envs/macs3_bedtools.yaml"
    shell:
        """
        epic2 --treatment {input.bam} \
            --output {output.peaks} \
            -gn mm10 \
            --guess-bampe
        """

# MACS3 callpeak for narrow histone marks
rule macs3_callpeak_narrow:
    input:
        bam = INPUT_DIR + "/{sample}_dedup.bam"
    output:
        peaks = RESULTS_DIR + "/{sample}_peaks.narrowPeak"
    params:
        outdir = RESULTS_DIR,
        name = "{sample}"
    conda:
        "./envs/macs3_bedtools.yaml"
    shell:
        """
        macs3 callpeak -t {input.bam} \
            --outdir {params.outdir} \
            -n {params.name} \
            -g mm \
            -f BAMPE \
            -q 0.05
        """


