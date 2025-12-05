import os
import pandas as pd
import yaml

# 配置文件
configfile: "config_k9.yaml"

# 加载样本信息
samples = pd.DataFrame(config["samples"])
comparisons = config["comparisons"]
genome = config.get("genome", "mm10")

# 输出目录结构
OUTDIR = config.get("outdir", "results")
FILEDIR = config.get("filedir", "results")

# 创建所需目录
for dir in [f"{OUTDIR}/diffbind", f"{OUTDIR}/annotation", f"{OUTDIR}/motifs", f"{OUTDIR}/logs/diffbind", 
                  f"{OUTDIR}/logs/annotation", f"{OUTDIR}/logs/motifs"]:
    os.makedirs(dir, exist_ok=True)

# 修正组织性: 定义常量和辅助函数
MARK_PEAK_CALLER_MAP = {
    "H3K4me3": "macs3",
    "H3K9me3": "epic2",
    "H3K27me3": "epic2"
}

def get_peak_file_pattern(mark, tool):
    if tool == "macs3":
        return f"{{sample}}_peaks.narrowPeak"
    elif tool == "epic2":
        return f"{{sample}}_peaks_ep.broadPeak"
    else:
        raise ValueError(f"Unknown peak caller: {tool}")

# 最终目标文件
rule all:
    input:
        # DiffBind输出
        expand(f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks.bed", comparison=comparisons.keys()),
        expand(f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks_upregulated.bed", comparison=comparisons.keys()),
        expand(f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks_downregulated.bed", comparison=comparisons.keys()),


         # 常规注释输出
        expand(f"{OUTDIR}/annotation/{{comparison}}_annotated.txt", comparison=comparisons.keys()),
        expand(f"{OUTDIR}/annotation/{{comparison}}_annotated_up.txt", comparison=comparisons.keys()),
        expand(f"{OUTDIR}/annotation/{{comparison}}_annotated_down.txt", comparison=comparisons.keys()),

        # TE注释输出
        expand(f"{OUTDIR}/annotation_TE/{{comparison}}_TE_annotated.txt", comparison=comparisons.keys()),
        expand(f"{OUTDIR}/annotation_TE/{{comparison}}_TE_annotated_up.txt", comparison=comparisons.keys()),
        expand(f"{OUTDIR}/annotation_TE/{{comparison}}_TE_annotated_down.txt", comparison=comparisons.keys()),
        # Motif分析输出
        expand(f"{OUTDIR}/motifs/up/{{comparison}}/knownResults.html", comparison=comparisons.keys()),
        expand(f"{OUTDIR}/motifs/down/{{comparison}}/knownResults.html", comparison=comparisons.keys()),
        # 通路富集输出
        expand(f"{OUTDIR}/enrichment/up/{{comparison}}_GO_enrichment.csv", comparison=comparisons.keys()),
        expand(f"{OUTDIR}/enrichment/up/{{comparison}}_reactome_enrichment.csv", comparison=comparisons.keys()),
        expand(f"{OUTDIR}/enrichment/up/{{comparison}}_GO_dotplot.pdf", comparison=comparisons.keys()),
        expand(f"{OUTDIR}/enrichment/up/{{comparison}}_reactome_dotplot.pdf", comparison=comparisons.keys()),

        expand(f"{OUTDIR}/enrichment/down/{{comparison}}_GO_enrichment.csv", comparison=comparisons.keys()),
        expand(f"{OUTDIR}/enrichment/down/{{comparison}}_reactome_enrichment.csv", comparison=comparisons.keys()),
        expand(f"{OUTDIR}/enrichment/down/{{comparison}}_GO_dotplot.pdf", comparison=comparisons.keys()),
        expand(f"{OUTDIR}/enrichment/down/{{comparison}}_reactome_dotplot.pdf", comparison=comparisons.keys())


# 生成DiffBind样本表
rule generate_samplesheet:
    output:
        f"{OUTDIR}/diffbind/samplesheet.csv"
    run:
        sample_data = []
        for idx, row in samples.iterrows():
            mark = "H3K9me3"
            peak_caller = "epic2"  # 或其他工具可以根据需要修改
            peak_file_pattern = get_peak_file_pattern(mark, peak_caller)
            sample_data.append({
                "SampleID": f"{row['sample']}",
                "Condition": row["condition"],
                "Replicate": row.get("replicate", 1),
                "bamReads": f"{FILEDIR}/dedup_bam_pe/{row['sample']}_dedup.bam",
                "Peaks": f"{FILEDIR}/callpeaks/results/{peak_caller}/{peak_file_pattern.format(sample=row['sample'])}",
                "PeakCaller": "bed"
            })
        pd.DataFrame(sample_data).to_csv(output[0], index=False)

# ===== 1. DiffBind差异分析 =====
##记得配置HOMER基因组 perl /mnt/data2/userdata/wanglu/RNAseq/RIMA/miniconda3/envs/macs3_bedtools/share/homer/.//configureHomer.pl -install mm10

rule diffbind_analysis:
    input:
        samplesheet = f"{OUTDIR}/diffbind/samplesheet.csv"
    output:
        bed = f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks.bed",
        bed_up = f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks_upregulated.bed",  # 使用OUTDIR变量
        bed_down = f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks_downregulated.bed",  # 使用OUTDIR变量
        report = f"{OUTDIR}/diffbind/{{comparison}}_report.pdf"
    params:
        comparison = "{comparison}",
        fdr = config.get("fdr_cutoff", 0.05),
        log2fc = config.get("log2fc_cutoff", 1),
        # 调整组别顺序以匹配DiffBind的对比逻辑
        group1 = lambda wildcards: comparisons[wildcards.comparison][1],  # 处理组
        group2 = lambda wildcards: comparisons[wildcards.comparison][0]   # 对照组
    log:
        f"{OUTDIR}/logs/diffbind/{{comparison}}.log"
    script:
        "scripts/run_diffbind.R"
'''
# ===== 2. HOMER注释 =====
rule annotate_peaks:
    input:
        f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks.bed"
    output:
        f"{OUTDIR}/annotation/{{comparison}}_annotated.txt"
    params:
        genome = genome,
        additional_params = config.get("annotation_params", "")
    log:
        f"{OUTDIR}/logs/annotation/{{comparison}}.log"
    shell:
        """
        annotatePeaks.pl {input} {params.genome} \
            -mask \
            {params.additional_params} \
            > {output} 2> {log}
        """
'''


# ===== 2. ChIPseeker注释 =====

rule annotate_peaks:
    input:
        peaks = f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks.bed"
    output:
        annotation = f"{OUTDIR}/annotation/{{comparison}}_annotated.txt",  # 添加统一前缀
        plot_dir = directory(f"{OUTDIR}/annotation/{{comparison}}_annotated_plots")  # 修改目录
    params:
        genome = genome,
        organism = config.get("organism", "mm10"),  # 默认为mm10，根据需要修改
        promoter_region = config.get("promoter_region", "c(-3000, 3000)"),  # 默认启动子区域为TSS上下游3kb
        additional_params = config.get("annotation_params", "")
    log:
        f"{OUTDIR}/logs/annotation/{{comparison}}.log"
    script:
        "scripts/chipseeker_annotate.R"  # R脚本路径

rule annotate_peaks_up:
    input:
        peaks = f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks_upregulated.bed"
    output:
        annotation = f"{OUTDIR}/annotation/{{comparison}}_annotated_up.txt",  # 调整为后缀式命名
        plot_dir = directory(f"{OUTDIR}/annotation/{{comparison}}_annotated_up_plots")
    params:
        genome = genome,
        organism = config.get("organism", "mm10"),
        promoter_region = config.get("promoter_region", "c(-3000, 3000)"),
        additional_params = config.get("annotation_params", "")
    log:
        f"{OUTDIR}/logs/annotation/{{comparison}}_up.log"
    script:
        "scripts/chipseeker_annotate.R"


rule annotate_peaks_down:
    input:
        peaks = f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks_downregulated.bed"
    output:
        annotation = f"{OUTDIR}/annotation/{{comparison}}_annotated_down.txt",  # 调整为后缀式命
        plot_dir = directory(f"{OUTDIR}/annotation/{{comparison}}_annotated_down_plots")
    params:
        genome = genome,
        organism = config.get("organism", "mm10"),
        promoter_region = config.get("promoter_region", "c(-3000, 3000)"),
        additional_params = config.get("annotation_params", "")
    log:
        f"{OUTDIR}/logs/annotation/{{comparison}}_down.log"
    script:
        "scripts/chipseeker_annotate.R"


# ===== 4. 转座元件(TE)注释 =====

rule annotate_TE:
    input:
        peaks = f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks.bed",
        gtf = config["repeatmasker_gtf"]  # 必须作为input声明
    output:
        annotation = f"{OUTDIR}/annotation_TE/{{comparison}}_TE_annotated.txt",
        plot_dir = directory(f"{OUTDIR}/annotation_TE/{{comparison}}_TE_annotated_plots")  # 统一添加annotated标识
    params:
        analysis_type = "all"  # 添加类型参数
    log:
        f"{OUTDIR}/logs/annotation_TE/{{comparison}}.log"
    script:
        "scripts/annotate_TE.R"

rule annotate_TE_up:
    input:
        peaks = f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks_upregulated.bed",
        gtf = config["repeatmasker_gtf"]  # 取消注释，必须作为input
    output:
        annotation = f"{OUTDIR}/annotation_TE/{{comparison}}_TE_annotated_up.txt",
        plot_dir = directory(f"{OUTDIR}/annotation_TE/{{comparison}}_TE_annotated_up_plots")
    params:
        analysis_type = "up"  # 添加类型参数
    log:
        f"{OUTDIR}/logs/annotation_TE/{{comparison}}_up.log"
    script:
        "scripts/annotate_TE.R"

rule annotate_TE_down:
    input:
        peaks = f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks_downregulated.bed",
        gtf = config["repeatmasker_gtf"]  # 取消注释，必须作为input
    output:
        annotation = f"{OUTDIR}/annotation_TE/{{comparison}}_TE_annotated_down.txt",
        plot_dir = directory(f"{OUTDIR}/annotation_TE/{{comparison}}_TE_annotated_down_plots")
    params:
        analysis_type = "down"  # 添加类型参数
    log:
        f"{OUTDIR}/logs/annotation_TE/{{comparison}}_down.log"
    script:
        "scripts/annotate_TE.R"


# ===== 3. HOMER Motif分析 =====


#perl /mnt/data2/userdata/wanglu/RNAseq/RIMA/miniconda3/envs/macs3_bedtools/share/homer/.//configureHomer.pl -install vertebrates
from snakemake.io import directory
rule find_up_motifs:
    input:
        f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks_upregulated.bed"
    output:
        html = f"{OUTDIR}/motifs/up/{{comparison}}/knownResults.html",
        outdir = directory(f"{OUTDIR}/motifs/up/{{comparison}}")
    params:
        genome = genome,
    log:
        f"{OUTDIR}/logs/motifs/up/{{comparison}}.log"
    shell:
        """
        mkdir -p {output.outdir}
        findMotifsGenome.pl {input} {params.genome} {output.outdir} -preparse > {log} 2>&1
        touch {output.html}
        """


rule find_down_motifs:
    input:
        f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks_downregulated.bed"
    output:
        html = f"{OUTDIR}/motifs/down/{{comparison}}/knownResults.html",
        outdir = directory(f"{OUTDIR}/motifs/down/{{comparison}}")
    params:
        genome = genome,
    log:
        f"{OUTDIR}/logs/motifs/down/{{comparison}}.log"
    shell:
        """
        mkdir -p {output.outdir}
        findMotifsGenome.pl {input} {params.genome} {output.outdir} -preparse > {log} 2>&1
        touch {output.html}
        """


########### 4.GO分析


rule chipseeker_clusterprofiler_up:
    input:
        peaks = f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks_upregulated.bed",
    output:
        enr_go = f"{OUTDIR}/enrichment/up/{{comparison}}_GO_enrichment.csv",
        enr_reactome = f"{OUTDIR}/enrichment/up/{{comparison}}_reactome_enrichment.csv",
        dotplot_go = f"{OUTDIR}/enrichment/up/{{comparison}}_GO_dotplot.pdf",
        dotplot_reactome = f"{OUTDIR}/enrichment/up/{{comparison}}_reactome_dotplot.pdf"
    params:
        org = config.get("organism_code", "mmu"),  # 如 mmu 或 hsa
        ont = config.get("go_ont", "BP"),         # GO 本体: BP/CC/MF
        pcut = config.get("pvalue_cutoff", 0.05),
        qcut = config.get("qvalue_cutoff", 0.2)
    log:
        f"{OUTDIR}/logs/enrichment/up/{{comparison}}.log"
    script:
        "scripts/chipseeker_clusterProfiler.R"

rule chipseeker_clusterprofiler_down:
    input:
        peaks = f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks_downregulated.bed",
    output:
        enr_go = f"{OUTDIR}/enrichment/down/{{comparison}}_GO_enrichment.csv",
        enr_reactome = f"{OUTDIR}/enrichment/down/{{comparison}}_reactome_enrichment.csv",
        dotplot_go = f"{OUTDIR}/enrichment/down/{{comparison}}_GO_dotplot.pdf",
        dotplot_reactome = f"{OUTDIR}/enrichment/down/{{comparison}}_reactome_dotplot.pdf"
    params:
        org = config.get("organism_code", "mmu"),  # 如 mmu 或 hsa
        ont = config.get("go_ont", "BP"),         # GO 本体: BP/CC/MF
        pcut = config.get("pvalue_cutoff", 0.05),
        qcut = config.get("qvalue_cutoff", 0.2)
    log:
        f"{OUTDIR}/logs/enrichment/down/{{comparison}}.log"
    script:
        "scripts/chipseeker_clusterProfiler.R"


# ===== 可选：生成报告 =====
rule generate_report:
    input:
        diffbind = expand(f"{OUTDIR}/diffbind/{{comparison}}_diffpeaks.bed", comparison=comparisons.keys()),
        annotation = expand(f"{OUTDIR}/annotation/{{comparison}}_annotated.txt", comparison=comparisons.keys()),
        motifs = expand(f"{OUTDIR}/motifs/{{comparison}}/knownResults.html", comparison=comparisons.keys())
    output:
        report = f"{OUTDIR}/final_report.html"
    params:
        project_name = config.get("project_name", "ChIP-seq Analysis")
    script:
        "scripts/generate_report.R"

