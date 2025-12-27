import glob
import os
from os.path import join
from pathlib import Path
from gitdb.util import dirname

wildcard_constraints:
    sample = r"[^/]+"

configfile: "config.yaml"

OUTDIR   = Path(config.get("outdir", ".")).resolve()
DATAROOT = Path(config["data_root"]).resolve()
SUBPATH  = Path(config.get("fastq_subpath", "basecalling2/pass"))
REF      = Path(config["reference"]).resolve()
REF_BASENAME = REF.name

def _discover_samples():
    names = []
    for d in DATAROOT.iterdir():
        sp = d / SUBPATH
        if sp.is_dir() and any(sp.rglob("*.fastq.gz")):
            names.append(d.name)
    return sorted(names)

SAMPLES = _discover_samples()

def fastq_list(wc):
    # All per-sample FASTQs, recursive just in case
    return sorted(str(p) for p in (DATAROOT / wc.sample / SUBPATH).rglob("*.fastq.gz"))

def find_dorado_inputs(wc):
    base = config["data_root"]  # e.g. "/mnt/.../UltraLong"
    s = wc.sample
    cand_dirs = [
        f"{base}/{s}/basecalling/pass/bams/",
        f"{base}/{s}/basecalling/fail",
        f"{base}/{s}/pass",
        f"{base}/{s}/fail",
        f"{base}/{s}",
    ]
    files = []
    for d in cand_dirs:
        if os.path.isdir(d):
            files += sorted(glob.glob(os.path.join(d, "*.bam")))
            files += sorted(glob.glob(os.path.join(d, "*.fastq"))) \
                  +  sorted(glob.glob(os.path.join(d, "*.fq"))) \
                  +  sorted(glob.glob(os.path.join(d, "*.fastq.gz"))) \
                  +  sorted(glob.glob(os.path.join(d, "*.fq.gz")))
    if not files:
        raise ValueError(f"No Dorado outputs found for sample {s} under {base}")
    return files
# ---------------------- targets ----------------------
rule all:
    input:
        # reference indices (written into OUTDIR/ref/)
        str(OUTDIR / "ref" / f"{REF_BASENAME}.mmi"),
        str(OUTDIR / "ref" / f"{REF_BASENAME}.fai"),
        str(OUTDIR/"ref/GCF_000001405.40/promoters_1kb.pc.bed"),
        str(OUTDIR/"ref/GCF_000001405.40/promoters_1kb.pc.renamed.bed"),
        # FASTQC
        expand("{outdir}/qc/fastqc/{sample}_fastqc.html", outdir=str(OUTDIR), sample=SAMPLES),
        # Alignment
        expand("{outdir}/alignment/{sample}.merged.fastq", outdir=str(OUTDIR), sample=SAMPLES),
        expand("{outdir}/alignment/{sample}.bam",outdir=str(OUTDIR),sample=SAMPLES),
        expand("{outdir}/alignment/{sample}.sorted.bam", outdir=str(OUTDIR),sample=SAMPLES),
        expand("{outdir}/alignment/{sample}.sorted.bam.bai",outdir=str(OUTDIR),sample=SAMPLES),
        expand("{outdir}/alignment/{sample}.summary.txt", outdir=str(OUTDIR),sample=SAMPLES),
        # QC  AND  stats
        expand("{outdir}/coverage/{sample}.mosdepth.summary.txt",sample=SAMPLES,outdir=str(OUTDIR)),
        expand("{outdir}/coverage/{sample}.mosdepth.global.dist.txt",  sample=SAMPLES,outdir=str(OUTDIR)),
        expand("{outdir}/stats/{sample}.flagstat.txt",sample=SAMPLES,outdir=str(OUTDIR)),
        expand("{outdir}/stats/{sample}.stats.txt", sample=SAMPLES,outdir=str(OUTDIR)),
        expand("{outdir}/stats/{sample}.idxstats.txt", sample=SAMPLES,outdir=str(OUTDIR)),
        expand("{outdir}/coverage/{sample}.chrom_mean.tsv",sample=SAMPLES, outdir=str(OUTDIR)),
        expand("{outdir}/coverage/{sample}.run_mean.txt",sample=SAMPLES,outdir=str(OUTDIR)),
        # Methylation
        expand("{outdir}/mods/{sample}.5mC.CpG.bed",sample=SAMPLES,outdir=str(OUTDIR)),
        expand("{outdir}/mods/{sample}.promoter_methylation.csv",sample=SAMPLES,outdir=str(OUTDIR)),
        expand("{outdir}/mods/{sample}.promoter_methylation.clean.csv",sample=SAMPLES,outdir=str(OUTDIR)),
        expand("{outdir}/mods/{sample}.promoter_methylation.clean.annotated.csv", sample=SAMPLES, outdir=str(OUTDIR)),
        expand("{outdir}/mods/{sample}.promoter_methylation.clean.annotated.PDF", sample=SAMPLES, outdir=str(OUTDIR)),
        # # # variants
        # Clair3
        expand("{outdir}/variants/{sample}/clair3/merge_output.vcf.gz", sample=SAMPLES, outdir=str(OUTDIR)),
        expand("{outdir}/variants/{sample}/clair3/merge_output.norm_sort.vcf.gz", sample=SAMPLES, outdir=str(OUTDIR)),
        # # Longshot with clair3 as input
        # expand("{outdir}/variants/{sample}/longshot/longshot_clair3.raw.vcf.gz", sample=SAMPLES, outdir=OUTDIR),
        # expand("{outdir}/variants/{sample}/longshot/longshot_clair3.norm_sort.vcf.gz", sample=SAMPLES, outdir=OUTDIR),
        # sniffles2
        expand("{outdir}/variants/{sample}/sniffles/sniffles.snf", sample=SAMPLES, outdir=str(OUTDIR))

# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# MAPPING
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
'''
index_ref — Build reference indices once per project.
Uses minimap2 to create the .mmi and samtools to make the .fai, both written under OUTDIR/ref.
Note: faidx normally lands next to the FASTA, so we copy it into OUTDIR for tidy project-local artifacts.
'''
rule index_ref:
    input:
        fasta = str(REF)
    output:
        mmi = str(OUTDIR / "ref" / f"{REF_BASENAME}.mmi"),
        fai = str(OUTDIR / "ref" / f"{REF_BASENAME}.fai")
    params:
        outdir = str(OUTDIR)
    log:
        str(OUTDIR / "logs" / "index_ref.log")
    threads: 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}/ref {params.outdir}/logs
        minimap2 -d {output.mmi} {input.fasta} &> {log}
        samtools faidx {input.fasta} >> {log} 2>&1
        cp -f {input.fasta}.fai {output.fai}
        """

rule prepare_fastq_for_alignment:
    input:
        reads_dir = lambda wc: join(config["data_root"], wc.sample, "basecalling", "pass")
    output:
        merged_fastq="{outdir}/alignment/{sample}.merged.fastq"
    log:
        "{outdir}/logs/{sample}.prepare_merged_fastq.log"
    threads: 10
    shell:
        r"""
            set -euo pipefail

            mkdir -p "$(dirname "{output.merged_fastq}")" "$(dirname "{log}")"
        
            shopt -s nullglob
        
            FASTQ_FILES=( "{input.reads_dir}"/*.fastq "{input.reads_dir}"/*.fq )
        
            if [ "${{#FASTQ_FILES[@]}}" -gt 0 ]; then
              echo "Found ${{#FASTQ_FILES[@]}} plain FASTQ files. Merging..." >> "{log}"
              cat "${{FASTQ_FILES[@]}}" > "{output.merged_fastq}"
            else
              GZ_FILES=( "{input.reads_dir}"/*.fastq.gz "{input.reads_dir}"/*.fq.gz )
        
              if [ "${{#GZ_FILES[@]}}" -gt 0 ]; then
                echo "No plain FASTQ found. Found ${{#GZ_FILES[@]}} FASTQ.GZ files. Decompressing & merging..." >> "{log}"
                pigz -dc -p 8 "${{GZ_FILES[@]}}" > "{output.merged_fastq}"
              else
                echo "ERROR: No FASTQ or FASTQ.GZ files found in: {input.reads_dir}" >> "{log}"
                exit 1
              fi
            fi
        
            shopt -u nullglob
        
            LINES=$(wc -l < "{output.merged_fastq}")
            echo "Merged FASTQ: {output.merged_fastq}" >> "{log}"
            echo "Total lines: $LINES" >> "{log}"
        
            if [ $((LINES % 4)) -ne 0 ]; then
              echo "ERROR: Merged FASTQ line count ($LINES) is not divisible by 4." >> "{log}"
              exit 1
            fi


        """

rule fastqc_merged_fastq:
    input:
        fastq = rules.prepare_fastq_for_alignment.output.merged_fastq
    output:
        html = "{outdir}/qc/fastqc/{sample}_fastqc.html"
    threads: 4
    params:
        outdir=str(OUTDIR)
    log:
        "{outdir}/logs/{sample}.fastqc.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.outdir}/qc/fastqc" "$(dirname "{log}")"

        fastqc \
          -t {threads} \
          --outdir "{params.outdir}/qc/fastqc" \
          "{input.fastq}" \
          > "{log}" 2>&1
        """






'''
align_dorado — Map ONT reads and produce a sorted BAM in one pass.
Uses the minimap2 .mmi (map-ont with tuned params) and pipes to samtools sort → {outdir}/alignment/{sample}.sorted.bam.
Threaded; stderr goes to a per-sample log for debugging.
'''
rule align_dorado:
    input:
        ref = rules.index_ref.output.mmi,
        merged_fastq = rules.prepare_fastq_for_alignment.output.merged_fastq
    params:
        tmpdir    = "{outdir}/alignment/{sample}.dorado_tmp",
    output:
        bam        = "{outdir}/alignment/{sample}.bam",
        sorted_bam = "{outdir}/alignment/{sample}.sorted.bam",
        bai        = "{outdir}/alignment/{sample}.sorted.bam.bai",
        summary    = "{outdir}/alignment/{sample}.summary.txt"
    threads: 32
    log:
        "{outdir}/logs/{sample}.dorado_align.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.tmpdir} $(dirname {output.bam}) $(dirname {log})
     
        # 1) Align: when INPUT is a directory, dorado requires --output-dir
        minimap2 -t {threads} -a -Y -x map-ont {input.ref} "{input.merged_fastq}" 2>> {log} | samtools view -@ {threads} -b -o {output.bam}

        # 3) Sort & index
        samtools sort -@ {threads} -o {output.sorted_bam} {output.bam} >> {log} 2>&1
        samtools index -@ 8 {output.sorted_bam} >> {log} 2>&1

        # 4) Summary
        samtools stats {output.sorted_bam} > {output.summary} 2>> {log}

        # 5) Clean
        rm -rf "{params.tmpdir}"
        """

# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# STATS and QC
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================

'''
mosdepth_coverage — Compute per-sample coverage summaries with mosdepth.
Produces summary and distribution files (global/region) under {outdir}/coverage/, skipping huge per-base tracks by default (-n).
If you need windowed coverage, swap -n for --by <window> to emit regions.bed.gz.
'''
rule mosdepth_coverage:
    input:
        bam = "{outdir}/alignment/{sample}.sorted.bam",
        bai = "{outdir}/alignment/{sample}.sorted.bam.bai"
    output:
        summary     = "{outdir}/coverage/{sample}.mosdepth.summary.txt",
        globaldist  = "{outdir}/coverage/{sample}.mosdepth.global.dist.txt",
        chrom_mean  = "{outdir}/coverage/{sample}.chrom_mean.tsv",
        run_mean    = "{outdir}/coverage/{sample}.run_mean.txt"
    params:
        prefix = "{outdir}/coverage/{sample}"
    threads: 6
    log:
        "{outdir}/logs/{sample}.mosdepth.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.prefix}) $(dirname {log})

        # -n avoids huge per-base outputs
        mosdepth -t {threads} -n {params.prefix} {input.bam} > {log} 2>&1

        # Per-chromosome mean: skip header and any 'total' row
        awk 'NR>1 && tolower($1)!="total" {{print $1"\t"$4}}' \
            {output.summary} > {output.chrom_mean}

        # Length-weighted run-wide mean depth across chromosomes
        awk 'NR>1 && tolower($1)!="total" {{L+=$2; S+=$2*$4}} END{{if(L>0) printf("%.6f\n", S/L); else print "NA"}}' \
            {output.summary} > {output.run_mean}
        """

'''
mapping_stats — Generate mapping QC with samtools.
Writes flagstat, stats, and idxstats text files under {outdir}/stats/ for each sample (MultiQC-friendly).
Requires the BAM index and runs with modest threading.
'''
rule mapping_stats:
    input:
        bam = "{outdir}/alignment/{sample}.sorted.bam",
        bai = "{outdir}/alignment/{sample}.sorted.bam.bai"
    output:
        flagstat = "{outdir}/stats/{sample}.flagstat.txt",
        stats    = "{outdir}/stats/{sample}.stats.txt",
        idx      = "{outdir}/stats/{sample}.idxstats.txt"
    threads: 4
    log:
        "{outdir}/logs/{sample}.samtools_stats.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.flagstat}) $(dirname {log})
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat}
        samtools stats    -@ {threads} {input.bam} > {output.stats}
        samtools idxstats               {input.bam} > {output.idx}
        """
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# Methylation
# 5mC (CpG) in human WBCs is abundant and informative for regulation (especially around promoters/TSS).
# 6mA in mammals is rare/controversial at genomic scale—treat as exploratory.
# 4mC in human is generally not expected; likely negligible with standard models (watch for artifacts).
# Focus first on 5mC in CpG context (--cpg you already ran).
# Coverage filter: keep sites with N_valid_cov ≥ 5 (or ≥10 if you have depth to spare).
# Get a GRCh38.p14 gene annotation (GTF). (GENCODE is clean and current for hg38/GRCh38.p14.)
# GENCODE
# Build a promoter BED (e.g., TSS ±1 kb, strand-aware).
# Intersect CpG sites (filtered) with promoters.
# Compute coverage-weighted promoter methylation per gene:
# weighted_mean = sum(fraction * N_valid_cov) / sum(N_valid_cov)
# Flag genes with weighted_mean ≥ 0.6 and sum coverage ≥ 25 (example cutoffs; tune to your data).
# If you ran modkit pileup for --motif A 0 --mod-codes a (6mA) or --motif C 0 --mod-codes c (4mC), you’ll get analogous bedMethyl.
# Analyze the same way (filter by coverage, aggregate to promoters/genes).
# Interpret with caution in human WBC—bulk 6mA evidence is limited; many signals are near repeats or aligner artifacts.
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
rule make_promoters_1kb_and_rename:
    input:
        gtf = config["genome_gtf_file"],
        rename_map = config["rename_chr_file"]
    output:
        prom_bed = "{outdir}/ref/GCF_000001405.40/promoters_1kb.pc.bed",
        prom_bed_renamed = "{outdir}/ref/GCF_000001405.40/promoters_1kb.pc.renamed.bed"
    params:
        flank = 1000
    log:
        "{outdir}/logs/make_promoters_1kb.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.prom_bed}")" "$(dirname "{output.prom_bed_renamed}")" "$(dirname "{log}")"

        # 1) Build promoters (TSS +/- flank) from GTF -> BED6
        # Note: tss is 0-based for BED: for '+' TSS=$4-1; for '-' TSS=$5-1
        awk 'BEGIN{{FS=OFS="\t"}}
             $3=="gene"{{
               gene_id=""; gene_name="";
               n=split($9,a,";");
               for(i=1;i<=n;i++){{
                 gsub(/^ +/,"",a[i]);
                 if(a[i]~/(^| )gene_id /){{ split(a[i],b," "); gsub(/"/,"",b[2]); gene_id=b[2] }}
                 if(a[i]~/(^| )gene_name /){{ split(a[i],b," "); gsub(/"/,"",b[2]); gene_name=b[2] }}
               }}
               if($7=="+"){{ tss=$4-1 }} else {{ tss=$5-1 }}
               start=tss-{params.flank}; if(start<0)start=0; end=tss+{params.flank};
               print $1, start, end, gene_id";"gene_name, 0, $7
             }}' "{input.gtf}" \
        | sort -k1,1 -k2,2n \
        > "{output.prom_bed}"

        # 2) Rename contigs using mapping file (NC_... -> 1..24)
        awk 'BEGIN{{OFS="\t"}}
             NR==FNR {{map[$1]=$2; next}}
             {{
               if ($1 in map) $1=map[$1];
               print
             }}' "{input.rename_map}" "{output.prom_bed}" \
        | sort -k1,1 -k2,2n \
        > "{output.prom_bed_renamed}"

        echo "Wrote: {output.prom_bed}" >> "{log}"
        echo "Wrote: {output.prom_bed_renamed}" >> "{log}"
        """


rule modkit_pileup_5mc_cpg:
    input:
        bam = "{outdir}/alignment/{sample}.sorted.bam",
        bai = "{outdir}/alignment/{sample}.sorted.bam.bai",
        ref = config["reference"],
        promotors_bed_file= rules.make_promoters_1kb_and_rename.output.prom_bed
    output:
        bed = "{outdir}/mods/{sample}.5mC.CpG.bed"
    params:
        modkit_thrs = config.get("modkit_threshold_c", ""),
        workdir=os.getcwd(),
        isz= config.get("modkit_interval_bp", 100000),
        csize=config.get("modkit_chunk_intervals", 8),
        outdir="{outdir}/mods/",
        ref_dir=dirname(config["reference"]),
        image=config.get("modkit_docker","ontresearch/modkit:latest"),
    threads: 12
    log:
        "{outdir}/logs/{sample}.modkit_5mc_cpg.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bed}) $(dirname {log})
        docker run --rm -u $(id -u):$(id -g) \
          -v {params.workdir}:{params.workdir} \
          -v {params.ref_dir}:{params.ref_dir} \
          -v {params.outdir}:{params.outdir} \
          -w {params.workdir} \
          {params.image} \
          modkit pileup {input.bam} {output.bed} \
          --mod-thresholds {params.modkit_thrs} \
          --include-bed {input.promotors_bed_file} \
          --interval-size {params.isz} \
          --cpg --ref {input.ref} -t {threads} \
          > {log} 2>&1
        """

"""
    Compute coverage-weighted promoter CpG methylation from modkit pileup (--cpg).
    Output columns: gene, methylation_weighted_pct, num_cpgs, total_coverage
"""
shell.executable("bash")
rule gene_promoter_methylation:
    input:
        bed = "{outdir}/mods/{sample}.5mC.CpG.bed",
        fai = str(OUTDIR / "ref" / f"{REF_BASENAME}.fai"),
        promoters_bed_file = rules.make_promoters_1kb_and_rename.output.prom_bed_renamed
    output:
        csv = "{outdir}/mods/{sample}.promoter_methylation.csv"
    threads: 1
    log:
        "{outdir}/logs/{sample}.mod_promoter_methylation.log"
    shell:
        r"""
        mkdir -p $(dirname {output.csv})

        bedtools map -g {input.fai} -sorted \
          -a <(bedtools sort -g {input.fai} -i {input.promoters_bed_file}) \
          -b <(bedtools sort -g {input.fai} -i {input.bed}) \
          -c 1,12,10 -o count,sum,sum | awk 'BEGIN{{OFS=","; print "gene,methylation_weighted_pct,num_cpgs,total_coverage"}}
               {{count=$(NF-2); sum_m=$(NF-1); sum_c=$NF; pct=(sum_c>0?100*sum_m/sum_c:0);
                 print $4, sprintf("%.2f", pct), count, sum_c}}' > {output.csv}
        """

rule clean_mods:
    input:
        in_csv="{outdir}/mods/{sample}.promoter_methylation.csv"
    output:
        out_csv="{outdir}/mods/{sample}.promoter_methylation.clean.csv"
    threads: 1
    shell:
        r"""
        python scripts/annotate_mods.py {input.in_csv} -o {output.out_csv}
        """

rule annotate_mods:
    input:
        in_csv="{outdir}/mods/{sample}.promoter_methylation.clean.csv"
    output:
        out_csv="{outdir}/mods/{sample}.promoter_methylation.clean.annotated.csv"
    params:
        methylation_markers= config.get("methylation_markers.csv", "methylation_markers.csv")
    threads: 1
    shell:
        r"""
        python scripts/annotate_moded_genes.py {input.in_csv} {params.methylation_markers} -o {output.out_csv}
        """

rule visualize_annotated_mods:
    input:
        in_csv="{outdir}/mods/{sample}.promoter_methylation.clean.annotated.csv"
    output:
        out_plot="{outdir}/mods/{sample}.promoter_methylation.clean.annotated.PDF"
    params:
        show_mods_above=40
    threads: 1
    shell:
        r"""
         python scripts/visualize_moded_genes.py {input.in_csv} -o {output.out_plot} --only-biomarkers --label_threshold {params.show_mods_above}
        """


# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# Variant calling


#TODO
# sudo dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.2.0 /var/lib/minknow/data/SequencingData/SAM5/ --output-dir /var/lib/minknow/data/SequencingData/SAM5/basecalling/ -r -x cuda:all --min-qscore 7 --emit-fastq --batchsize 128
#
# bcftools annotate --rename-chrs /home/mateo/PycharmProjects/NanoporePipeline/rename_chr.txt merge_output.chr_only.vcf.gz -Oz -o merge_output.renamed_chr.vcf.gz
# bcftools index merge_output.renamed_chr.vcf.gz
# bcftools annotate -a /home/mateo/PycharmProjects/NanoporePipeline/clinvar.vcf.gz -c INFO merge_output.renamed_chr.vcf.gz -o annotated.vcf
# bcftools view -f 'PASS,.' annotated.vcf > filtered_annotated.vcf
# bcftools view -i 'QUAL>30 && DP>30' annotated.vcf -Oz -o FINAL.severe.vcf.gz
# https://myvariant.info/v1/variant/chr1:g.92840760G%3EA?assembly=hg38
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# Clair3 — small variant calling for ONT
rule clair3_call:
    input:
        bam = "{outdir}/alignment/{sample}.sorted.bam",
        bai = "{outdir}/alignment/{sample}.sorted.bam.bai",
        ref = config["reference"]
    output:
        vcf = "{outdir}/variants/{sample}/clair3/merge_output.vcf.gz",
        vcf_sorted = "{outdir}/variants/{sample}/clair3/merge_output.norm_sort.vcf.gz"
    params:
        image   = config.get("clair3_docker_image", "hkubal/clair3:latest"),
        gpus    = config.get("clair3_gpus", "0"),
        outdir  = "{outdir}/variants/{sample}/clair3/",
        workdir = os.getcwd(),
        ref_dir= dirname(config["reference"]),
        model   = config.get("clair3_model_path", "/opt/models/r1041_e82_400bps_sup_v500"),
        bed_file= config.get("target_bed", ""),
        min_mq=config.get("clair3_min_mq", 1),
        min_coverage=config.get("clair3_min_coverage", 1),
        snp_min_af=config.get("clair3_snp_min_af", 1),
        indel_min_af=config.get("clair3_indel_min_af", 1),
        qual=config.get("clair3_qual", 1),
        contigs=config.get("clair3_contigs", "")
    threads: 10
    log:
        "{outdir}/logs/{sample}.clair3.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {output.vcf}) $(dirname {log})

        docker run --rm --gpus "device=0" \
        -e CUDA_VISIBLE_DEVICES=0 \
        -u $(id -u):$(id -g) \
          -v {params.workdir}:{params.workdir} \
          -v {params.ref_dir}:{params.ref_dir} \
          -w {params.workdir} \
          {params.image} \
          bash -lc 'set -euo pipefail
                  source /opt/conda/etc/profile.d/conda.sh
                  conda activate clair3        

                    test -d "{params.model}" || {{ echo "Model not found: {params.model}"; exit 2; }}


                    /opt/bin/run_clair3.sh \
                      --bam_fn={input.bam} \
                      --ref_fn={input.ref} \
                      --remove_intermediate_dir \
                      --longphase_for_phasing \
                      --threads={threads} \
                      --min_mq={params.min_mq} \
                      --min_coverage={params.min_coverage} \
                      --snp_min_af={params.snp_min_af} --indel_min_af={params.indel_min_af} --qual={params.qual} \
                      --platform=ont \
                      --ctg_name='{params.contigs}' \
                      --model_path={params.model} \
                      --output {params.outdir}'
                      
        """

# Longshot — small variant calling for ONT
rule longshot_with_clair_input_call:
    input:
        bam = "{outdir}/alignment/{sample}.sorted.bam",
        bai = "{outdir}/alignment/{sample}.sorted.bam.bai",
        clair_vcf = "{outdir}/variants/{sample}/clair3/merge_output.vcf.gz",
        ref = config["reference"]
    output:
        vcf = "{outdir}/variants/{sample}/longshot/longshot_clair3.raw.vcf.gz",
        vcf_sorted = "{outdir}/variants/{sample}/longshot/longshot_clair3.norm_sort.vcf.gz"
    params:
        image    = config.get("longshot_docker_image", "quay.io/biocontainers/longshot:1.0.0--h8dc4d9d_3"),
        outdir   = "{outdir}/variants/{sample}/longshot",
        workdir  = os.getcwd(),
        ref_dir  = dirname(config["reference"]),
        min_mq   = int(config.get("longshot_min_mq", 20)),
        min_af   = config.get("longshot_min_af", 0.25),
        max_cov  = config.get("longshot_max_cov", 0),   # set 0 to disable
    threads: 10
    log:
        "{outdir}/logs/{sample}.longshot.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {output.vcf}) $(dirname {log})

        docker run --rm \
          -u $(id -u):$(id -g) \
          -v {params.workdir}:{params.workdir} \
          -v {params.ref_dir}:{params.ref_dir} \
          -w {params.workdir} \
          {params.image} \
          bash -lc 'set -euo pipefail
            # raw call
            longshot \
              --bam {input.bam} \
              --ref {input.ref} \
              -A -S \
              -v {input.clair_vcf} \
              --out {params.outdir}/longshot_clair3.raw.vcf.gz \
              --sample_id {wildcards.sample} \
              --min_mapq {params.min_mq} \
              --min_alt_frac {params.min_af} 
              > {log} 2>&1
          '
          bcftools norm -f {input.ref} -m -both {output.vcf} | bcftools sort -Oz -o {output.vcf_sorted}
          tabix -f -p vcf {output.vcf_sorted} >> {log} 2>&1
        """

# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# CNV
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# Sniffles2 — CNV caller for ONT
rule sniffles_call:
    input:
        bam = "{outdir}/alignment/{sample}.sorted.bam",
        bai = "{outdir}/alignment/{sample}.sorted.bam.bai",
        ref = config["reference"]
    output:
        vcf = "{outdir}/variants/{sample}/sniffles/sniffles.raw.vcf.gz",
        snf = "{outdir}/variants/{sample}/sniffles/sniffles.snf"
    params:
        image   = config.get("sniffles_docker_image",
                             "quay.io/biocontainers/sniffles:2.6.3--pyhdfd78af_0"),
        outdir  = "{outdir}/variants/{sample}/sniffles",
        workdir = os.getcwd(),
        ref_dir = dirname(config["reference"]),
        minsupport = config.get("sniffles_min_support", ""),
        minsvlen = config.get("sniffles_minsvlen", ""),
        qc_coverage = config.get("sniffles_qc_coverage", ""),
        mapq=config.get("sniffles_mapq", ""),
        sniffles_tandem_repeats_bed=config.get("sniffles_tandem_repeats_bed", "")
    threads: 12
    log:
        "{outdir}/logs/{sample}.sniffles.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {output.vcf}) $(dirname {log})


        docker run --rm \
          -u $(id -u):$(id -g) \
          -v {params.workdir}:{params.workdir} \
          -v {params.ref_dir}:{params.ref_dir} \
          -w {params.workdir} \
          {params.image} \
          bash -lc 'set -euo pipefail
            sniffles \
              --input "{input.bam}" \
              --vcf "{params.outdir}/sniffles.raw.vcf.gz" \
              --snf "{output.snf}" \
              --reference "{input.ref}" \
              --minsupport {params.minsupport} \
              --minsvlen {params.minsvlen} \
              --mapq {params.mapq} \
              --qc-coverage {params.qc_coverage} \
              --threads {threads} 
          ' > {log} 2>&1
        """

rule cuteSV_call:
    input:
        bam = "{outdir}/alignment/{sample}.sorted.bam",
        bai = "{outdir}/alignment/{sample}.sorted.bam.bai",
        ref = config["reference"]
    output:
        vcf = "{outdir}/variants/{sample}/cuteSV/cutesv.raw.vcf"
    params:
        outdir = "{outdir}/variants/{sample}/cuteSV",
        workdir = os.getcwd(),

        # cuteSV parameters
        min_support = config.get("cutesv_min_support", 10),
        min_svlen   = config.get("cutesv_minsvlen", 50),
        min_mapq    = config.get("cutesv_mapq", 20),
        max_cluster_bias = config.get("cutesv_max_cluster_bias", 100),
        diff_ratio_merging_INS = config.get("cutesv_diff_ratio_merging_INS", 0.3),
        diff_ratio_merging_DEL = config.get("cutesv_diff_ratio_merging_DEL", 0.3)
    threads: 12
    shell:
        r"""
        mkdir -p {params.outdir}

        cuteSV \
            {input.bam} \
            {input.ref} \
            {output.vcf} \
            {params.outdir} \
            --threads {threads} \
            --min_support {params.min_support} \
            --min_svlen {params.min_svlen} \
            --min_mapq {params.min_mapq} \
            --max_cluster_bias {params.max_cluster_bias} \
            --diff_ratio_merging_INS {params.diff_ratio_merging_INS} \
            --diff_ratio_merging_DEL {params.diff_ratio_merging_DEL}
        """

rule svim_call:
    input:
        bam = "{outdir}/alignment/{sample}.sorted.bam",
        bai = "{outdir}/alignment/{sample}.sorted.bam.bai",
        ref = config["reference"]
    output:
        vcf = "{outdir}/variants/{sample}/svim/svim.raw.vcf"
    params:
        outdir = "{outdir}/variants/{sample}/svim",
        svim_bin = config.get("svim_bin", "svim"),
        # common filters / knobs
        min_mapq = config.get("svim_mapq", 20),
        min_sv_size = config.get("svim_min_sv_size", 50),
        # SVIM has different “modes”: typically `alignment` for BAM input
        mode = config.get("svim_mode", "alignment"),
    threads: 12
    log:
        "{outdir}/variants/{sample}/svim/svim.log"
    shell:
        r"""
        mkdir -p {params.outdir}

        # SVIM writes results into {params.outdir}. We then copy/standardize the VCF path.
        {params.svim_bin} {params.mode} \
            {params.outdir} \
            {input.bam} \
            {input.ref} \
            --min_mapq {params.min_mapq} \
            --min_sv_size {params.min_sv_size} \
            --cores {threads} \
            > {log} 2>&1

        # Standardize output path for the pipeline:
        # Depending on SVIM version, the VCF name can vary. Most commonly it's:
        #   {params.outdir}/variants.vcf
        # We normalize it to {output.vcf}.
        if [ -f "{params.outdir}/variants.vcf" ]; then
            cp "{params.outdir}/variants.vcf" "{output.vcf}"
        elif [ -f "{params.outdir}/variants.vcf.gz" ]; then
            # In case SVIM already gzipped (rare), decompress to raw vcf
            zcat "{params.outdir}/variants.vcf.gz" > "{output.vcf}"
        else
            echo "ERROR: SVIM output VCF not found in {params.outdir}" >&2
            ls -lah "{params.outdir}" >&2 || true
            exit 2
        fi
        """



