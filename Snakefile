# Minimal UL-ONT prep: index REF + per-sample FASTQ concatenation
# --------------------------------------------------------------
# config.yaml must define:
#   reference: "/abs/path/to/GRCh38.fa"
#   fastq_root: "/abs/path/to/UltraLong"
#   fastq_subpath: "basecalling/pass"     # subpath under each {sample}
#   outdir: "/abs/path/to/ProjectOut"     # where outputs go
import os
from pathlib import Path

from gitdb.util import dirname

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

# ---------------------- targets ----------------------
rule all:
    input:
        # reference indices (written into OUTDIR/ref/)
        str(OUTDIR / "ref" / f"{REF_BASENAME}.mmi"),
        str(OUTDIR / "ref" / f"{REF_BASENAME}.fai"),
        expand("{outdir}/reads/{sample}.all.fastq.gz",
               outdir=str(OUTDIR), sample=SAMPLES),
        expand("{outdir}/qc/{sample}", outdir=str(OUTDIR),sample=SAMPLES),
        expand("{outdir}/alignment/{sample}.sorted.bam",
           outdir=str(OUTDIR),sample=SAMPLES),
        expand("{outdir}/alignment/{sample}.sorted.bam.bai",
            outdir=str(OUTDIR),sample=SAMPLES),
        expand("{outdir}/coverage/{sample}.mosdepth.summary.txt",
            sample=SAMPLES,outdir=OUTDIR),
        expand("{outdir}/coverage/{sample}.mosdepth.global.dist.txt",
            sample=SAMPLES,outdir=OUTDIR),
        expand("{outdir}/stats/{sample}.flagstat.txt",
            sample=SAMPLES,outdir=OUTDIR),
        expand("{outdir}/stats/{sample}.stats.txt",
            sample=SAMPLES,outdir=OUTDIR),
        expand("{outdir}/stats/{sample}.idxstats.txt",
            sample=SAMPLES,outdir=OUTDIR),
        expand("{outdir}/coverage/{sample}.chrom_mean.tsv",
           sample=SAMPLES, outdir=OUTDIR),
        expand("{outdir}/coverage/{sample}.run_mean.txt",
            sample=SAMPLES,outdir=OUTDIR),
        expand("{outdir}/variants/{sample}/clair3.vcf.gz", sample=SAMPLES, outdir=OUTDIR)
        # expand("{outdir}/variants/{sample}/medaka.vcf.gz", sample=SAMPLES, outdir=OUTDIR),
        # expand("{outdir}/variants/{sample}/longshot.vcf.gz", sample=SAMPLES, outdir=OUTDIR),
        # expand("{outdir}/variants/{sample}/consensus.vcf.gz", sample=SAMPLES, outdir=OUTDIR),


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
 Concat_fastq — Merge all per-sample *.fastq.gz into a single file.
 Discovers the file list for each {sample} and concatenates them to {outdir}/reads/{sample}.all.fastq.gz (gzip concatenation is safe).
 Depends on the reference indices so mapping can start immediately after.
'''
rule concat_fastq:
    input:
        # force ref indexing to finish first (if you want strict ordering)
        ref_mmi = rules.index_ref.output.mmi,
        ref_fai = rules.index_ref.output.fai,
        fastqs  = fastq_list
    output:
        out = "{outdir}/reads/{sample}.all.fastq.gz"
    params:
        outdir = str(OUTDIR)
    log:
        "{outdir}/logs/{sample}.concat.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.out}) $(dirname {log})
        # Safe concat even with many files; gzip concatenation is valid.
        printf "%s\n" {input.fastqs} | xargs -r cat > {output.out}
        """

'''
nanoplot — Read-level QC on the merged FASTQ.
Runs NanoPlot to produce N50, length/quality plots, and an HTML report under {outdir}/qc/{sample}.
Outputs a tracked directory so all artifacts are captured.
'''
rule nanoplot:
    input:
        fq = "{outdir}/reads/{sample}.all.fastq.gz"
    output:
        qcdir = directory("{outdir}/qc/{sample}")
    threads: 4
    log:
        "{outdir}/logs/{sample}.nanoplot.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.qcdir} $(dirname {log})
        NanoPlot --fastq {input.fq} \
                 --outdir {output.qcdir} \
                 --threads {threads} \
                 --N50 --verbose > {log} 2>&1
        """

'''
align_minimap2 — Map ONT reads and produce a sorted BAM in one pass.
Uses the minimap2 .mmi (map-ont with tuned params) and pipes to samtools sort → {outdir}/alignment/{sample}.sorted.bam.
Threaded; stderr goes to a per-sample log for debugging.
'''
rule align_minimap2:
    input:
        ref = rules.index_ref.output.mmi,
        fq  = "{outdir}/reads/{sample}.all.fastq.gz"
    output:
        bam = "{outdir}/alignment/{sample}.sorted.bam"
    threads: 32
    log:
        "{outdir}/logs/{sample}.minimap2.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        minimap2 -t {threads} -K 200M -k 20 -w 25 -I50G -ax map-ont {input.ref} {input.fq} \
          2> {log} \
        | samtools sort -@ {threads} -o {output.bam}
        """
'''
index_bam — Create the .bai index for each sorted BAM.
Required for random access and for downstream tools (coverage and stats).
Lightweight, typically runs quickly.
'''
rule index_bam:
    input:
        "{outdir}/alignment/{sample}.sorted.bam"
    output:
        "{outdir}/alignment/{sample}.sorted.bam.bai"
    threads: 2
    shell:
        "samtools index -@ {threads} {input}"

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

# ==============================================
# Variant calling (3 callers) + consensus (≥2/3)
# ==============================================

# Clair3 — small variant calling for ONT
rule clair3_call:
    input:
        bam = "{outdir}/alignment/{sample}.sorted.bam",
        bai = "{outdir}/alignment/{sample}.sorted.bam.bai",
        ref = config["reference"]
    output:
        vcf = "{outdir}/variants/{sample}/clair3.vcf.gz",
        tbi = "{outdir}/variants/{sample}/clair3.vcf.gz.tbi"
    params:
        image   = config.get("clair3_docker_image", "hkubal/clair3:latest"),
        gpus    = config.get("clair3_gpus", "0"),
        outdir  = "{outdir}/variants/{sample}/clair3",
        workdir = os.getcwd(),
        ref_dir= dirname(config["reference"]),
        model   = config.get("clair3_model_path", ""),
        modeldir= config.get("clair3_model_dir","models"),
        bed_file= config.get("clair3_bed_file", ""),
        min_mq=5,
        min_coverage=2
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
                      --bed_fn={params.bed_file} \
                      --remove_intermediate_dir \
                      --include_all_ctgs \
                      --threads={threads} \
                      --min_mq={params.min_mq} \
                      --min_coverage={params.min_coverage} \
                      --snp_min_af=0.08 --indel_min_af=0.15 --qual=0 \
                      --platform=ont \
                      --model_path=/models/{params.model} \
                      --output {params.outdir}'

        # Normalize/sort/index on host (same paths as inside container)
        bcftools norm -f {input.ref} -m -both -Oz -o {params.outdir}/clair3.norm.vcf.gz {params.outdir}/merge_output.vcf.gz
        bcftools sort -Oz -o {output.vcf} {params.outdir}/clair3.norm.vcf.gz
        tabix -f -p vcf {output.vcf} >> {log} 2>&1
        """
