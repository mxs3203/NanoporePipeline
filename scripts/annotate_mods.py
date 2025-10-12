import argparse
from pathlib import Path
import pandas as pd


def classify(w):
    try:
        w = float(w)
    except Exception:
        return "NA"
    if w < 30:
        return "Low"
    elif w <= 70:
        return "Medium"
    else:
        return "High"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_csv", help="Input promoter_methylation.csv")
    ap.add_argument("-o", "--out", default=None, help="Output CSV (default: <input>.clean.csv)")
    args = ap.parse_args()

    infile = Path(args.input_csv)
    outfile = Path(args.out) if args.out else infile.with_suffix(".clean.csv")

    # Read as 4 columns with no header
    df = pd.read_csv(infile)
    df.columns = ["gene_biotype", "weight", "n_cpg", "total_cov"]

    # Split gene_biotype at ';;' into gene and biotype
    split = df["gene_biotype"].astype(str).str.split(";;", n=1, expand=True)
    df["gene"] = split[0].str.replace(";", "", regex=False).str.strip()
    df["biotype"] = split[1].fillna("").str.strip()

    # Convert numeric columns
    df["weight"] = pd.to_numeric(df["weight"], errors="coerce")
    df["n_cpg"] = pd.to_numeric(df["n_cpg"], errors="coerce", downcast="integer")
    df["total_cov"] = pd.to_numeric(df["total_cov"], errors="coerce", downcast="integer")

    # Classify methylation by weight
    df["methylation_class"] = df["weight"].apply(classify)

    # Select and order columns
    out = df[["gene", "biotype", "weight", "n_cpg", "total_cov", "methylation_class"]]

    out.to_csv(outfile, index=False)
    print(f"Wrote: {outfile}")

if __name__ == "__main__":
    main()