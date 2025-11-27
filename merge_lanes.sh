#!/bin/bash


LANE_DIR_L5="/mnt/data2/userdata/co_worker_labs/wanglu/rawdata/cutrun/B16_H3K36me3_2510/LBFC20251677-02/20250929_LH00524_0381_A22WY3GLT4"
LANE_DIR_L2="/mnt/data2/userdata/co_worker_labs/wanglu/rawdata/cutrun/B16_H3K36me3_2510/LBFC20251677-02/20251008_LH00708_0220_A235VMGLT4"


OUTPUT_DIR="/mnt/data2/userdata/co_worker_labs/wanglu/rawdata/cutrun/B16_H3K36me3_2510/MERGED_FASTQS"


mkdir -p "$OUTPUT_DIR" || { echo "$OUTPUT_DIR"; exit 1; }



find "$LANE_DIR_L5" -maxdepth 1 -name "*.R1.fastq.gz" -print0 | while IFS= read -r -d $'\0' R1_PATH_L5; do
    

    FILENAME_L5=$(basename "$R1_PATH_L5")
    SAMPLE_PREFIX=$(echo "$FILENAME_L5" | sed -E 's/LZM0923_L5_//; s/\.R1\.fastq\.gz//')
    

    R1_PATH_L2="${LANE_DIR_L2}/LZM0923_L2_${SAMPLE_PREFIX}.R1.fastq.gz"

    if [[ -f "$R1_PATH_L2" ]]; then
        OUTPUT_R1="${OUTPUT_DIR}/${SAMPLE_PREFIX}_Merged.R1.fastq.gz"
        cat "$R1_PATH_L5" "$R1_PATH_L2" > "$OUTPUT_R1"
    else
        continue 
    fi

    R2_PATH_L5="${LANE_DIR_L5}/LZM0923_L5_${SAMPLE_PREFIX}.R2.fastq.gz"
    R2_PATH_L2="${LANE_DIR_L2}/LZM0923_L2_${SAMPLE_PREFIX}.R2.fastq.gz"
    
    if [[ -f "$R2_PATH_L5" && -f "$R2_PATH_L2" ]]; then
        OUTPUT_R2="${OUTPUT_DIR}/${SAMPLE_PREFIX}_Merged.R2.fastq.gz"
        cat "$R2_PATH_L5" "$R2_PATH_L2" > "$OUTPUT_R2"
    else
    fi

done

echo "----------------------------------------------------"
