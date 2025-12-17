#!/usr/bin/env python3

import sys
import pandas as pd
import pyranges as pr
import gzip

if len(sys.argv) != 4:
    sys.exit('python3 %s <pairs> <enzyme.bed> <outFile(.gz)>' % (sys.argv[0]))

pairFile = sys.argv[1]
enzBED   = sys.argv[2]
outFile  = sys.argv[3]

# --- helper for gz/no-gz ---
def open_maybe_gz(path, mode):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


# --- read BED as PyRanges ---
bed_pr = pr.read_bed(enzBED)
bed_pr = bed_pr.drop(["Score", "Strand"])  # keep Chrom, Start, End, Name


CHUNK = 1_000_000
first_write = True


# ---------------------------------------------------------------------
# Map single side (A or B) → fragment name
# ---------------------------------------------------------------------
def map_side(df, chr_col, pos_col, bed_pr):

    if df.shape[0] == 0:
        return pd.Series([], dtype=object)

    q = pd.DataFrame({
        "Chromosome": df[chr_col].astype(str),
        "Start": df[pos_col].astype(int) - 1,
        "End": df[pos_col].astype(int),
        "idx": df.index,
    })

    prq = pr.PyRanges(q)
    inter = prq.join(bed_pr, strandedness=False).df

    out = pd.Series("NA", index=df.index, dtype=object)
    if inter.empty:
        return out

    inter = inter[["idx", "Name"]].drop_duplicates("idx")
    out.loc[inter["idx"]] = inter["Name"].values
    return out


# ---------------------------------------------------------------------
# Streaming Pairs
# ---------------------------------------------------------------------
with open_maybe_gz(pairFile, "rt") as f, open_maybe_gz(outFile, "wt") as fout:

    while True:
        lines = [f.readline() for _ in range(CHUNK)]
        lines = [l for l in lines if l]
        if not lines:
            break

        header = [l for l in lines if l.startswith("#")]
        body   = [l for l in lines if not l.startswith("#")]

        if first_write:
            for l in header:
                fout.write(l)
            first_write = False

        if not body:
            continue

        df = pd.DataFrame(
            [l.rstrip().split('\t') for l in body],
            columns=["readID","chrA","posA","chrB","posB","strandA","strandB","type"]
        )
        df["posA"] = df["posA"].astype(int)
        df["posB"] = df["posB"].astype(int)
        df.index = range(len(df))  # ensures stable mapping index

        # compute A, B fragments separately
        fragA = map_side(df, "chrA", "posA", bed_pr)
        fragB = map_side(df, "chrB", "posB", bed_pr)

        df["fragA"] = fragA
        df["fragB"] = fragB

        # write chunk
        df.to_csv(fout, sep="\t", header=False, index=False)
