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
        f"{base}/{s}/basecalling/pass",
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
        # Alignment
        expand("{outdir}/alignment/{sample}.bam",outdir=str(OUTDIR),sample=SAMPLES),
        expand("{outdir}/alignment/{sample}.sorted.bam", outdir=str(OUTDIR),sample=SAMPLES),
        expand("{outdir}/alignment/{sample}.sorted.bam.bai",outdir=str(OUTDIR),sample=SAMPLES),
        expand("{outdir}/alignment/{sample}.summary.tsv", outdir=str(OUTDIR),sample=SAMPLES),
        # QC  AND  stats
        expand("{outdir}/coverage/{sample}.mosdepth.summary.txt",sample=SAMPLES,outdir=OUTDIR),
        expand("{outdir}/coverage/{sample}.mosdepth.global.dist.txt",  sample=SAMPLES,outdir=OUTDIR),
        expand("{outdir}/stats/{sample}.flagstat.txt",sample=SAMPLES,outdir=OUTDIR),
        expand("{outdir}/stats/{sample}.stats.txt", sample=SAMPLES,outdir=OUTDIR),
        expand("{outdir}/stats/{sample}.idxstats.txt", sample=SAMPLES,outdir=OUTDIR),
        expand("{outdir}/coverage/{sample}.chrom_mean.tsv",sample=SAMPLES, outdir=OUTDIR),
        expand("{outdir}/coverage/{sample}.run_mean.txt",sample=SAMPLES,outdir=OUTDIR),
        # # # variants
        # Clair3
        expand("{outdir}/variants/{sample}/clair3/merge_output.vcf.gz", sample=SAMPLES, outdir=OUTDIR),
        expand("{outdir}/variants/{sample}/clair3/merge_output.norm_sort.vcf.gz", sample=SAMPLES, outdir=OUTDIR),
        # Longshot
        expand("{outdir}/variants/{sample}/longshot/longshot.raw.vcf.gz", sample=SAMPLES, outdir=OUTDIR),
        expand("{outdir}/variants/{sample}/longshot/longshot.norm_sort.vcf.gz", sample=SAMPLES, outdir=OUTDIR),
        # Longshot with clair3 as input
        expand("{outdir}/variants/{sample}/longshot/longshot_clair3.raw.vcf.gz", sample=SAMPLES, outdir=OUTDIR),
        expand("{outdir}/variants/{sample}/longshot/longshot_clair3.norm_sort.vcf.gz", sample=SAMPLES, outdir=OUTDIR),
        # sniffles2
        expand("{outdir}/variants/{sample}/sniffles/sniffles.snf", sample=SAMPLES, outdir=OUTDIR)

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
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}/ref {params.outdir}/logs
        minimap2 -d {output.mmi} {input.fasta} &> {log}
        samtools faidx {input.fasta} >> {log} 2>&1
        cp -f {input.fasta}.fai {output.fai}
        """

'''
align_dorado — Map ONT reads and produce a sorted BAM in one pass.
Uses the minimap2 .mmi (map-ont with tuned params) and pipes to samtools sort → {outdir}/alignment/{sample}.sorted.bam.
Threaded; stderr goes to a per-sample log for debugging.
'''
rule align_dorado:
    input:
        # Dorado aligner can use an mmi or fasta; leave your index rule as-is
        ref = rules.index_ref.output.mmi
    params:
        # point this to your per-sample basecalled reads directory
        # e.g. config.yaml: dorado_dir: "/mnt/.../UltraLong"
        reads_dir = lambda wc: join(config["data_root"], wc.sample, "basecalling", "pass"),
        tmpdir    = "{outdir}/alignment/{sample}.dorado_tmp",
        mm2_extra = "-Y",
        rg= lambda wc: f"@RG\tID:{wc.sample}\tSM:{wc.sample}\tPL:ONT"
    output:
        bam        = "{outdir}/alignment/{sample}.bam",
        sorted_bam = "{outdir}/alignment/{sample}.sorted.bam",
        bai        = "{outdir}/alignment/{sample}.sorted.bam.bai",
        summary    = "{outdir}/alignment/{sample}.summary.tsv"
    threads: 32
    log:
        "{outdir}/logs/{sample}.dorado_align.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.tmpdir} $(dirname {output.bam}) $(dirname {log})

        # 1) Align: when INPUT is a directory, dorado requires --output-dir
        dorado aligner {input.ref} "{params.reads_dir}" \
          --mm2-opts "{params.mm2_extra}" \
          --output-dir "{params.tmpdir}" >> {log} 2>&1

        # 2) Merge Dorado's chunked BAMs into a single BAM
        samtools merge -@ {threads} -o {output.bam} {params.tmpdir}/*.bam >> {log} 2>&1
        
        # 3) Sort & index
        samtools sort -@ {threads} -o {output.sorted_bam} {output.bam} >> {log} 2>&1
        samtools index -@ {threads} {output.sorted_bam} >> {log} 2>&1

        # 4) Summary
        dorado summary {output.sorted_bam} > {output.summary} 2>> {log}

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
    threads: 8
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
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================







# ===========================================================================================================================================================
# ===========================================================================================================================================================
# ===========================================================================================================================================================
# Variant calling
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
        model   = config.get("clair3_model_path", ""),
        modeldir= config.get("clair3_model_dir","models"),
        bed_file= config.get("target_bed", ""),
        min_mq=config.get("clair3_min_mq", 1),
        min_coverage=config.get("clair3_min_coverage", 1),
        snp_min_af=config.get("clair3_snp_min_af", 1),
        indel_min_af=config.get("clair3_indel_min_af", 1),
        qual=config.get("clair3_qual", 1)
    threads: 32
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
          -v {params.modeldir}:/models \
          -w {params.workdir} \
          {params.image} \
          bash -lc 'set -euo pipefail
                  source /opt/conda/etc/profile.d/conda.sh
                  conda activate clair3        

                    test -d "/models/{params.model}" || {{ echo "Model not found: /models/{params.model}"; exit 2; }}


                    /opt/bin/run_clair3.sh \
                      --bam_fn={input.bam} \
                      --ref_fn={input.ref} \
                      --remove_intermediate_dir \
                      --include_all_ctgs \
                      --longphase_for_phasing \
                      --threads={threads} \
                      --min_mq={params.min_mq} \
                      --min_coverage={params.min_coverage} \
                      --snp_min_af={params.snp_min_af} --indel_min_af={params.indel_min_af} --qual={params.qual} \
                      --platform=ont \
                      --device='cuda:0' \
                      --model_path=/models/{params.model} \
                      --output {params.outdir}'
                      
          bcftools norm -f {input.ref} -m -both {output.vcf} | bcftools sort -Oz -o {output.vcf_sorted}
          tabix -f -p vcf {output.vcf_sorted} >> {log} 2>&1
        """

# Longshot — small variant calling for ONT
rule longshot_call:
    input:
        bam = "{outdir}/alignment/{sample}.sorted.bam",
        bai = "{outdir}/alignment/{sample}.sorted.bam.bai",
        ref = config["reference"]
    output:
        vcf = "{outdir}/variants/{sample}/longshot/longshot.raw.vcf.gz",
        vcf_sorted = "{outdir}/variants/{sample}/longshot/longshot.norm_sort.vcf.gz"
    params:
        image    = config.get("longshot_docker_image", "quay.io/biocontainers/longshot:1.0.0--h8dc4d9d_3"),
        outdir   = "{outdir}/variants/{sample}/longshot",
        workdir  = os.getcwd(),
        ref_dir  = dirname(config["reference"]),
        min_mq   = int(config.get("longshot_min_mq", 20)),
        min_af   = config.get("longshot_min_af", 0.25),
        max_cov  = config.get("longshot_max_cov", 0),   # set 0 to disable
    threads:6
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
              --out {params.outdir}/longshot.raw.vcf.gz \
              --sample_id {wildcards.sample} \
              --min_mapq {params.min_mq} \
              --min_alt_frac {params.min_af} 
              > {log} 2>&1
          '
          
        bcftools norm -f {input.ref} -m -both {output.vcf} | bcftools sort -Oz -o {output.vcf_sorted}
        tabix -f -p vcf {output.vcf_sorted} >> {log} 2>&1
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
    threads:6
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
    threads: 32
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





