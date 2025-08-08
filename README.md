# Nanopore UL Read Mapping & QC — Snakemake Pipeline

A lightweight Snakemake workflow to **prepare a reference**, **merge per-sample FASTQs**, and generate **alignments + QC** for Oxford Nanopore ultra-long reads.

- Builds **`minimap2`** index (`.mmi`) and **`samtools`** index (`.fai`)
- Concatenates all `*.fastq.gz` **per sample**
- Maps with **minimap2** → sorted **BAM** + **BAI**
- Per-sample **NanoPlot** read QC, **mosdepth** coverage, and **samtools** mapping stats
- Runs samples **in parallel**; all outputs live under your chosen `outdir`

---

## 1) Prerequisites

- Python 3.8+ and **Snakemake ≥ 7**
- Tools on `$PATH` (or install via conda): `minimap2`, `samtools`, `NanoPlot`, `mosdepth`, `gzip`, `xargs`
- Bash available at `/bin/bash`

**Optional (recommended):** create a conda env  
```bash
mamba create -n ont-map snakemake minimap2 samtools nanoplot mosdepth -c conda-forge -c bioconda
mamba activate ont-map
```

---

## 2) Input Layout

The pipeline **auto-discovers samples** from:

```
{fastq_root}/{sample}/basecalling/pass/*.fastq.gz
```

Example:
```
/data/UltraLong/
├── SAM1/basecalling/pass/run1.fastq.gz
│                          └─ run2.fastq.gz
└── SAM2/basecalling/pass/laneA.fastq.gz
                           └─ laneB.fastq.gz
```

---

## 3) Configuration

Create `config.yaml` next to your `Snakefile`:

```yaml
# Path to reference FASTA
reference: "/abs/path/to/GRCh38.fa"

# Where raw reads live (contains SAM1, SAM2, ...)
fastq_root: "/abs/path/to/UltraLong"

# Subpath under each sample containing reads
fastq_subpath: "basecalling/pass"

# Where all results, logs, and derived files go
outdir: "/abs/path/to/ProjectOut"
```

---

## 4) What Gets Produced

Under `outdir`:

```
ref/
  GRCh38.fa.mmi
  GRCh38.fa.fai
reads/
  SAM1.all.fastq.gz
  SAM2.all.fastq.gz
qc/
  SAM1/ NanoPlot-report.html, plots...
  SAM2/ ...
alignment/
  SAM1.sorted.bam
  SAM1.sorted.bam.bai
  SAM2.sorted.bam
  SAM2.sorted.bam.bai
coverage/
  SAM1.mosdepth.summary.txt
  SAM1.mosdepth.global.dist.txt
  SAM1.mosdepth.region.dist.txt
  SAM2...
stats/
  SAM1.flagstat.txt
  SAM1.stats.txt
  SAM1.idxstats.txt
  SAM2...
logs/
  *.log per rule & sample
```

---

## 5) Running

**Dry-run** (see the plan):
```bash
snakemake -n
```

**Run with 32 parallel jobs**:
```bash
snakemake -j 32
```

**Build only specific targets**:
```bash
# Reference indices only
snakemake -j 8 ref/GRCh38.fa.mmi ref/GRCh38.fa.fai

# Align + index BAMs for two samples
snakemake -j 16 alignment/SAM1.sorted.bam.bai alignment/SAM2.sorted.bam.bai
```

> **Parallelism tip:** `-j` controls how many jobs (rules) run in parallel; each rule has its own `threads:` setting (e.g., minimap2 uses 10 threads). Match both to your cores.

---

## 6) Rules (What Each Step Does)

- **`index_ref`** – Build reference indices once per project (minimap2 `.mmi`, samtools `.fai`) into `OUTDIR/ref`. `.fai` is copied so all artifacts remain project-local.

- **`concat_fastq`** – For each sample, discover all `*.fastq.gz` and concatenate to `{outdir}/reads/{sample}.all.fastq.gz`. Gzip concatenation is standards-compliant; downstream tools read it as one file.

- **`nanoplot`** – Run NanoPlot on the merged FASTQ to produce per-sample read QC (HTML + plots) in `{outdir}/qc/{sample}`.

- **`align_minimap2`** – Map with minimap2 (`map-ont` tuning) and stream to `samtools sort` → `{outdir}/alignment/{sample}.sorted.bam`. Stderr goes to a per-sample log.

- **`index_bam`** – Make `{sample}.sorted.bam.bai` for random access and downstream tools.

- **`mosdepth_coverage`** – Per-sample coverage summary + distributions via mosdepth. Default uses `-n` to skip huge per-base files.

- **`mapping_stats`** – `samtools` QC: `flagstat`, `stats`, and `idxstats` per sample (MultiQC-friendly).

---

## 7) Tips & Customization

- **Include `fail/` reads too**  
  Extend the helper that collects FASTQs to also search `basecalling/fail/**/*.fastq.gz`.

- **Windowed coverage**  
  Replace `mosdepth -n` with `--by 10000` to emit 10 kb windows and `regions.bed.gz`.

- **Minimap2 tuning**  
  Edit flags in `align_minimap2` (e.g., `-I`, `-k`, `-w`, `-K`) to suit your dataset and memory.

- **Logs**  
  Every rule writes to `outdir/logs/`. If something fails, check the corresponding `.log`.

- **Incremental re-runs**  
  Snakemake only rebuilds missing/out-of-date targets. Delete a file to force its rebuild or use `--force`.

---

## 8) Common Pitfalls

- **“Target rules may not contain wildcards.”**  
  You tried to run a wildcard rule directly (e.g., `snakemake concat_fastq`).  
  → Run plain `snakemake` (uses `rule all`) or target concrete files like `reads/SAM1.all.fastq.gz`.

- **`'Rules' object has no attribute 'index_ref'`**  
  You referenced `rules.index_ref` before defining it.  
  → Move `rule all` below `index_ref` or list concrete output paths.

- **Indexes created next to the reference instead of `outdir`**  
  `samtools faidx` writes `.fai` next to the FASTA; this workflow **copies** it into `outdir/ref`.

---

## 9) Extending

Add on as needed:
- **Small variants:** Clair3, Medaka, Longshot → consensus VCF  
- **SVs:** Sniffles2, cuteSV  
- **CNV:** CNVpytor or windowed coverage + segmentation  
- **Methylation:** Megalodon (needs FAST5)  
- **Phasing:** WhatsHap  
- **Reporting:** MultiQC to aggregate logs and metrics

---

## 10) Reproducibility

- Keep `config.yaml` and the `Snakefile` under version control
- Pin tool versions via conda (`environment.yml`) for long-lived projects
- The workflow is idempotent and resume-safe

---

## 11) License & Citation

Use and adapt freely. If this workflow supports a publication, please cite the tools you use: **Snakemake**, **minimap2**, **samtools**, **NanoPlot**, and **mosdepth**.
