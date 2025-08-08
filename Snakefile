# Minimal UL-ONT prep: index REF + per-sample FASTQ concatenation
# --------------------------------------------------------------
# config.yaml must define:
#   reference: "/abs/path/to/GRCh38.fa"
#   fastq_root: "/abs/path/to/UltraLong"
#   fastq_subpath: "basecalling/pass"     # subpath under each {sample}
#   outdir: "/abs/path/to/ProjectOut"     # where outputs go

from pathlib import Path

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
        expand("{outdir}/coverage/{sample}.mosdepth.summary.txt", outdir=str(OUTDIR),sample=SAMPLES),
        expand("{outdir}/stats/{sample}.flagstat.txt",outdir=str(OUTDIR),sample=SAMPLES),
        expand("{outdir}/stats/{sample}.stats.txt",outdir=str(OUTDIR),sample=SAMPLES),
        expand("{outdir}/stats/{sample}.idxstats.txt",outdir=str(OUTDIR),sample=SAMPLES)

# ---------------------- reference index ----------------------
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

# ---------------------- per-sample concat ----------------------
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

# Read-level QC with NanoPlot (per sample)
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

# ---------------------- per-sample concat ----------------------
rule align_minimap2:
    input:
        ref = rules.index_ref.output.mmi,
        fq  = "{outdir}/reads/{sample}.all.fastq.gz"
    output:
        bam = "{outdir}/alignment/{sample}.sorted.bam"
    threads: 10
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

rule index_bam:
    input:
        "{outdir}/alignment/{sample}.sorted.bam"
    output:
        "{outdir}/alignment/{sample}.sorted.bam.bai"
    threads: 2
    shell:
        "samtools index -@ {threads} {input}"

# ---------- Coverage (mosdepth) ----------
rule mosdepth_coverage:
    input:
        bam = "{outdir}/alignment/{sample}.sorted.bam",
        bai = "{outdir}/alignment/{sample}.sorted.bam.bai"
    output:
        summary    = "{outdir}/coverage/{sample}.mosdepth.summary.txt",
        globaldist = "{outdir}/coverage/{sample}.mosdepth.global.dist.txt",
        regiondist = "{outdir}/coverage/{sample}.mosdepth.region.dist.txt"
    params:
        prefix = "{outdir}/coverage/{sample}"
    threads: 8
    log:
        "{outdir}/logs/{sample}.mosdepth.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.prefix}) $(dirname {log})
        # -n avoids gigantic per-base files; add --by 10000 for 10kb windows if you want regions.bed.gz
        mosdepth -t {threads} -n {params.prefix} {input.bam} > {log} 2>&1
        """

# ---------- Mapping stats (samtools) ----------
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

