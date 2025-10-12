#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
import pandas as pd


def normalize_gene(s: str) -> str:
    if pd.isna(s):
        return ""
    s = str(s).strip()
    # Remove parenthetical aliases: e.g., "CDKN2A (p16)" -> "CDKN2A"
    s = re.sub(r"\s*\([^)]*\)\s*$", "", s)
    # Collapse slashes/semicolons/“and” into commas for later splitting (when needed)
    s = re.sub(r"\s*/\s*|\s*;\s*|\s+\band\b\s+", ",", s, flags=re.I)
    # Keep only main token if multiple (here we only normalize, splitting is elsewhere)
    s = s.strip()
    # Uppercase for matching
    return s.upper()


def read_methylation_csv(path: Path) -> pd.DataFrame:
    """
    Supports:
      - With header: gene,biotype,weight,n_cpg,total_cov,methylation_class
      - Without header: same 6 columns in that order
    """
    df = pd.read_csv(path, dtype=str)
    # If no header (i.e., first row looks like data and columns are unnamed or numeric), rename:
    expected = ["gene", "biotype", "weight", "n_cpg", "total_cov", "methylation_class"]
    if list(df.columns) != expected:
        # Heuristic: if there are exactly 6 columns, assume there is no header.
        if df.shape[1] == 6:
            df.columns = expected
        else:
            # Try to align best-effort (fallback)
            cols = list(df.columns)
            # If the first column name literally looks like a gene name (e.g., "OR4F5"), treat as headerless
            if len(cols) >= 1 and not re.match(r"(?i)^gene$", str(cols[0])):
                # keep original cols but ensure we have the needed ones
                df = df.copy()
                df.columns = (cols + expected[len(cols):])[:len(cols)]
            # If still not aligned, raise a clear error
            missing = [c for c in expected if c not in df.columns]
            if missing:
                raise ValueError(f"Input methylation CSV must have columns {expected} "
                                 f"(with or without header). Missing: {missing}")
    # Convert numerics where possible
    for c in ["weight", "n_cpg", "total_cov"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def read_biomarkers_csv(path: Path) -> pd.DataFrame:
    """
    Expects columns at least: Group, Gene, Context, MethylationChange, Evidence
    Allows extra columns; splits multi-gene cells into multiple rows.
    """
    bm = pd.read_csv(path, dtype=str)
    required = ["Group", "Gene", "Context", "MethylationChange", "Evidence"]
    missing = [c for c in required if c not in bm.columns]
    if missing:
        raise ValueError(f"Biomarker CSV is missing required columns: {missing}")

    # Split rows where Gene has multiple entries (comma, slash, semicolon, 'and')
    gene_norm = (bm["Gene"]
                 .fillna("")
                 .str.replace(r"\s*/\s*|\s*;\s*|\s+\band\b\s+", ",", regex=True, case=False))
    bm = bm.assign(GeneSplit=gene_norm.str.split(",")).explode("GeneSplit")
    bm["GeneSplit"] = bm["GeneSplit"].fillna("").str.strip()
    # Normalize gene tokens (remove parentheses aliases, uppercase)
    bm["gene_norm"] = bm["GeneSplit"].map(normalize_gene)
    bm = bm[bm["gene_norm"].str.len() > 0].copy()

    # If a gene appears in multiple rows/categories, aggregate evidence/context/group
    agg = (bm.groupby("gene_norm", as_index=False)
             .agg({
                 "Group": lambda x: "; ".join(sorted(set([v for v in x if isinstance(v, str)]))),
                 "Context": lambda x: "; ".join(sorted(set([v for v in x if isinstance(v, str)]))),
                 "MethylationChange": lambda x: "; ".join(sorted(set([v for v in x if isinstance(v, str)]))),
                 "Evidence": lambda x: "; ".join(sorted(set([v for v in x if isinstance(v, str)]))),
             }))
    agg.rename(columns={
        "Group": "biomarker_group",
        "Context": "biomarker_context",
        "MethylationChange": "biomarker_change",
        "Evidence": "biomarker_evidence",
    }, inplace=True)
    return agg


def annotate(meth: pd.DataFrame, bm: pd.DataFrame) -> pd.DataFrame:
    meth = meth.copy()
    meth["gene_norm"] = meth["gene"].astype(str).map(normalize_gene)

    out = meth.merge(bm, how="left", on="gene_norm")
    out["is_biomarker"] = out["biomarker_group"].notna()
    # Reorder columns
    preferred = [
        "gene", "biotype", "weight", "n_cpg", "total_cov", "methylation_class",
        "is_biomarker", "biomarker_group", "biomarker_change",
        "biomarker_context", "biomarker_evidence"
    ]
    # Keep any extra columns at the end
    ordered = [c for c in preferred if c in out.columns] + [c for c in out.columns if c not in preferred]
    return out[ordered]


def main():
    ap = argparse.ArgumentParser(description="Annotate methylation CSV with biomarker info.")
    ap.add_argument("methylation_csv", help="Input methylation CSV (may or may not have header).")
    ap.add_argument("biomarkers_csv", help="Biomarker CSV with columns: Group,Gene,Context,MethylationChange,Evidence.")
    ap.add_argument("-o", "--out", default=None, help="Output CSV (default: <methylation_csv>.annotated.csv)")
    args = ap.parse_args()

    methylation_csv = Path(args.methylation_csv)
    biomarkers_csv = Path(args.biomarkers_csv)
    out_csv = Path(args.out) if args.out else methylation_csv.with_suffix(".annotated.csv")

    meth = read_methylation_csv(methylation_csv)
    bm = read_biomarkers_csv(biomarkers_csv)
    annotated = annotate(meth, bm)

    annotated.to_csv(out_csv, index=False)
    print(f"Wrote: {out_csv}")


if __name__ == "__main__":
    main()
