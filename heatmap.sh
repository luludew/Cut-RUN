#!/bin/bash
# 设置 Conda 环境和参考基因组文件路径
CONDA_ENV="cutandrun_smk"
BED_FILE="/mnt/data2/userdata/co_worker_labs/wanglu/ref/ref_files_mouse/mm10.TSS.sorted.bed"
PEAK_FILE="/mnt/data2/userdata/co_worker_labs/wanglu/cutrun/B16/smk_cutrun/pipeline/callpeaks/results/peaks"

# 创建输出目录
mkdir -p ../heatmap
mkdir -p ../computeMatrix
mkdir -p ../bamCoverage_log2_results

# 定义统一的颜色方案
HEATMAP_COLOR_SCHEME='Greens'

find . -name "*_dedup.bam" | while read bwfile; do
    file=$(basename "$bwfile")
    sample=${file%%_dedup*}
    
    # 定义处理后的标准BED文件路径
    # 注意：为每种修饰生成不同后缀，以避免混淆，但实际上它们都是标准BED
    PROCESSED_PEAK_FILE_H3K9me3="${PEAK_FILE}/${sample}_peaks.bed6"
    PROCESSED_PEAK_FILE_H3K27me3="${PEAK_FILE}/${sample}_peaks.bed6" # 和H3K9me3使用相同后缀，但实际是独立文件
    # 对于H3K4me3，其BED_FILE已经是TSS Bed文件，无需处理


    # === 新增：将原始峰文件处理为标准BED6格式的函数 ===
    # 输入：原始文件路径，输出文件路径
    # 注意：此函数假定原始文件包含标题行且格式为：
    # #Chromosome   Start   End   ChIPCount   Score   Strand
    # 并且 Score 列是浮点数，需要缩放到0-1000
    convert_to_bed6() {
        local input_file="$1"
        local output_file="$2"
        
        if [[ ! -f "$input_file" ]]; then
            echo "Error: Original peak file $input_file not found. Cannot convert to BED6."
            return 1
        fi

        # 使用 awk 处理：
        # - FNR > 1: 跳过第一行 (头部)
        # - printf: 格式化输出为 BED6
        #   - $1: chrom (Chromosome)
        #   - $2: chromStart (Start) - BED 是 0-based
        #   - $3: chromEnd (End) - BED 是 1-based
        #   - "Peak_" NR-1: name (生成一个简单的Peak名称，例如 Peak_1, Peak_2...)
        #   - sprintf("%.0f", ($5 / max_score) * 1000): score (将第五列的Score缩放到0-1000，四舍五入)
        #   - $6: strand (Strand)
        # 先找到最大 score 以进行缩放（这是一个简化的方法，更严谨的需要预先计算所有文件的最大值）
        # 这里我们假设最大 score 不会超过一个合理的大值，比如 200。
        # 如果您的 Score 实际值经常远大于 200，请调整 max_score_val。
        local max_score_val=200 # 假设一个最大分数，用于将score缩放到0-1000。请根据实际数据调整！
        awk -v OFS='\t' -v max_s="$max_score_val" 'FNR > 1 {
            score_scaled = ($5 / max_s) * 1000;
            if (score_scaled > 1000) score_scaled = 1000; # 确保不超过1000
            if (score_scaled < 0) score_scaled = 0;     # 确保不低于0
            printf "%s\t%d\t%d\tPeak_%d\t%.0f\t%s\n", $1, $2, $3, FNR-1, score_scaled, $6
        }' "$input_file" > "$output_file"
        echo "Converted to BED6: $output_file"
        return 0
    }
    # ==================================================


    if [[ "$file" == *"H3K9me3"* ]]; then
        bam_file="${sample}_dedup.bam"
        
        # 处理H3K9me3的峰文件为BED6
        if ! convert_to_bed6 "${PEAK_FILE}/${sample}_peaks_ep.broadPeak" "$PROCESSED_PEAK_FILE_H3K9me3"; then
            echo "Skipping $sample due to BED6 conversion error."
            continue
        fi

        # 添加 bamCoverage 步骤
        echo "Running bamCoverage for $sample (H3K9me3)..."        
        conda run -n "$CONDA_ENV" bamCoverage \
            -b "$bam_file" \
            -o "../bamCoverage_log2_results/${sample}.bigwig" \
            --binSize 100 \
            --smoothLength 1000 \
            --normalizeUsing CPM \
            -p 8
        
        compare_bw="../bamCoverage_log2_results/${sample}.bigwig"
        
        echo "$sample contains H3K9me3 modification, using scale-regions mode..."
        # computeMatrix (scale-regions 模式)
        conda run -n "$CONDA_ENV" computeMatrix scale-regions \
            -p 60 \
            -R "$PROCESSED_PEAK_FILE_H3K9me3" \
            -S "$compare_bw" \
            --regionBodyLength 5000 \
            --upstream 5000 \
            --downstream 5000 \
            --skipZeros \
            --binSize 100 \
            -o "../computeMatrix/matrix1_${sample}_scaled.gz" \
            --outFileSortedRegions "../computeMatrix/regions1_${sample}_genes.bed" \
            --verbose
            
    # H3K4me3 部分保持不变，因为它使用BED_FILE (TSS bed文件)
    # 假设 BED_FILE 已经是标准BED格式且没有头部
    # 如果需要运行，请取消以下所有行的注释
    # elif [[ "$file" == *"H3K4me3"* ]]; then 
    #     bam_file="${sample}_dedup.bam"
        
    #     echo "Running bamCoverage for $sample (H3K4me3)..."        
    #     conda run -n "$CONDA_ENV" bamCoverage \
    #         -b "$bam_file" \
    #         -o "../bamCoverage_log2_results/${sample}.bigwig" \
    #         --binSize 100 \
    #         --smoothLength 600 \
    #         --normalizeUsing CPM \
    #         -p 8
        
    #     compare_bw="../bamCoverage_log2_results/${sample}.bigwig"
        
    #     echo "$sample contains H3K4me3 modification, using reference-point mode with TSS focus..."
    #     conda run -n "$CONDA_ENV" computeMatrix reference-point \
    #         -p 60 \
    #         -R "${BED_FILE}" \
    #         -S "$compare_bw" \
    #         --referencePoint TSS \
    #         --upstream 2000 \
    #         --downstream 2000 \
    #         --skipZeros \
    #         --binSize 50 \
    #         -o "../computeMatrix/matrix1_${sample}_centered.gz" \
    #         --outFileSortedRegions "../computeMatrix/regions1_${sample}_genes.bed" \
    #         --verbose
            
    elif [[ "$file" == *"H3K27me3"* ]]; then
        bam_file="${sample}_dedup.bam"

        # 处理H3K27me3的峰文件为BED6
        if ! convert_to_bed6 "${PEAK_FILE}/${sample}_peaks_ep.broadPeak" "$PROCESSED_PEAK_FILE_H3K27me3"; then
            echo "Skipping $sample due to BED6 conversion error."
            continue
        fi
        
        # 为 H3K27me3 添加 bamCoverage 步骤
        echo "Running bamCoverage for $sample (H3K27me3)..."        
        conda run -n "$CONDA_ENV" bamCoverage \
            -b "$bam_file" \
            -o "../bamCoverage_log2_results/${sample}.bigwig" \
            --binSize 100 \
            --smoothLength 1000 \
            --normalizeUsing CPM \
            -p 8
        
        compare_bw="../bamCoverage_log2_results/${sample}.bigwig"
        
        echo "$sample contains H3K27me3 modification, using scale-regions mode..."
        # computeMatrix (scale-regions 模式用于 H3K27me3，更适用于宽域)
        conda run -n "$CONDA_ENV" computeMatrix scale-regions \
            -p 60 \
            -R "$PROCESSED_PEAK_FILE_H3K27me3" \
            -S "$compare_bw" \
            --regionBodyLength 5000 \
            --upstream 5000 \
            --downstream 5000 \
            --skipZeros \
            --binSize 100 \
            -o "../computeMatrix/matrix1_${sample}_scaled.gz" \
            --outFileSortedRegions "../computeMatrix/regions1_${sample}_genes.bed" \
            --verbose
            
    else
        echo "$sample did not identify a specific histone modification, skipping processing."
        continue # 如果没有识别到特定的组蛋白修饰，跳过当前文件
    fi
    
    # 确定正确的矩阵文件路径和绘图参数
    # 初始化 Y 轴和 Z 轴（颜色条）范围变量
    Y_MIN_VAL=0
    Y_MAX_VAL=0
    HEATMAP_Z_MIN=0
    HEATMAP_Z_MAX=0
    X_LABEL='' # 初始化X轴标签
    REGIONS_LABEL='' # 初始化 regionsLabel
    MATRIX_FILE="" # 初始化矩阵文件路径

    if [[ "$file" == *"H3K27me3"* ]]; then
        MATRIX_FILE="../computeMatrix/matrix1_${sample}_scaled.gz"
        X_LABEL='Peak region'
        Y_MIN_VAL=0
        Y_MAX_VAL=0.35
        HEATMAP_Z_MIN=0
        HEATMAP_Z_MAX=0.55
        REGIONS_LABEL='peaks'
        
    elif [[ "$file" == *"H3K9me3"* ]]; then
        MATRIX_FILE="../computeMatrix/matrix1_${sample}_scaled.gz"
        X_LABEL='Peak region'
        Y_MIN_VAL=0 
        Y_MAX_VAL=0.45
        HEATMAP_Z_MIN=0
        HEATMAP_Z_MAX=0.55
        REGIONS_LABEL='peaks'
        
    elif [[ "$file" == *"H3K4me3"* ]]; then
        MATRIX_FILE="../computeMatrix/matrix1_${sample}_centered.gz"
        X_LABEL='Distance from TSS (bp)'
        Y_MIN_VAL=0 
        Y_MAX_VAL=4
        HEATMAP_Z_MIN=0
        HEATMAP_Z_MAX=16
        REGIONS_LABEL='genes'
    fi

    # 只有当 MATRIX_FILE 被成功设定（即匹配到H3K4me3/H3K9me3/H3K27me3）时才进行绘图
    if [[ -n "$MATRIX_FILE" ]]; then
        echo "Generating plots for $sample..."
        
        # plotHeatmap (PNG)
        conda run -n "$CONDA_ENV" plotHeatmap \
            --regionsLabel "$REGIONS_LABEL" \
            --xAxisLabel "$X_LABEL" \
            -m "$MATRIX_FILE" \
            -out "../heatmap/${sample}_Heatmap.png" \
            --yMin "$Y_MIN_VAL" --yMax "$Y_MAX_VAL" \
            --zMin "$HEATMAP_Z_MIN" --zMax "$HEATMAP_Z_MAX" \
            --colorMap "$HEATMAP_COLOR_SCHEME"
        
        # plotHeatmap (PDF)
        conda run -n "$CONDA_ENV" plotHeatmap \
            -m "$MATRIX_FILE" \
            -out "../heatmap/${sample}_Heatmap.pdf" \
            --xAxisLabel "$X_LABEL" \
            --regionsLabel "$REGIONS_LABEL" \
            --plotFileFormat pdf \
            --dpi 720 \
            --yMin "$Y_MIN_VAL" --yMax "$Y_MAX_VAL" \
            --zMin "$HEATMAP_Z_MIN" --zMax "$HEATMAP_Z_MAX" \
            --colorMap "$HEATMAP_COLOR_SCHEME"
            
        # plotProfile (PNG)
        conda run -n "$CONDA_ENV" plotProfile \
            -m "$MATRIX_FILE" \
            -out "../heatmap/${sample}_Profile.png" \
            --regionsLabel "$REGIONS_LABEL" \
            --averageType mean \
            --yMin "$Y_MIN_VAL" --yMax "$Y_MAX_VAL"
            
        # plotProfile (PDF)
        conda run -n "$CONDA_ENV" plotProfile \
            -m "$MATRIX_FILE" \
            -out "../heatmap/${sample}_Profile.pdf" \
            --regionsLabel "$REGIONS_LABEL" \
            --plotFileFormat pdf \
            --perGroup \
            --dpi 720 \
            --averageType mean \
            --yMin "$Y_MIN_VAL" --yMax "$Y_MAX_VAL"
    fi
    
    echo "Processing complete: $sample"
done

echo "All processing is complete."