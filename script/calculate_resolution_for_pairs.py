#!/usr/bin/env python3
import sys
import gzip
import numpy as np
from collections import defaultdict
import math

def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def read_positions(pairs_file, chrom1_idx=1, pos1_idx=2, chrom2_idx=3, pos2_idx=4):
    chrom_pos = defaultdict(list)
    chrom_sizes = {}
    with open_maybe_gz(pairs_file) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                if line.startswith("#chromsize:"):
                    _, chrom, size = line.strip().split()
                    chrom_sizes[chrom] = int(size)
            else:
                cols = line.rstrip("\n").split("\t")
                chrom_pos[cols[chrom1_idx]].append(int(cols[pos1_idx]))
                chrom_pos[cols[chrom2_idx]].append(int(cols[pos2_idx]))
    return chrom_pos, chrom_sizes

def compute_stats(chrom_pos, bin_size, total_bins):
    contact_cutoff = 1000
    fulfil_bins = 0

    for chrom, pos_list in chrom_pos.items():
        if not pos_list:
            continue
        arr = np.array(pos_list)
        bin_idx = arr // bin_size
        unique_bins, counts = np.unique(bin_idx, return_counts=True)

        fulfil_bins += (counts > contact_cutoff).sum()

    #print(fulfil_bins, total_bins)
    return {
        "fraction": fulfil_bins / total_bins,
    }

def total_bins_from_header(chrom_sizes, resolution):
    total = 0
    for chrom, size in chrom_sizes.items():
        total += math.ceil(size / resolution)
    return total

def gradient_scan(pairs_file, resolutions,
                  chrom1_idx=1, pos1_idx=2, chrom2_idx=3, pos2_idx=4):
    print("#Reading positions from pairs file...\n")
    chrom_pos, chrom_size = read_positions(pairs_file, chrom1_idx, pos1_idx, chrom2_idx, pos2_idx)
    print("#Finished reading. Running gradient resolution scan...\n")

    for r in resolutions:
        total = total_bins_from_header(chrom_size, r)
        #print('resolution\t%i\t%i' % (r, total))
        stats = compute_stats(chrom_pos, r, total)
        if stats is None:
            continue

        frac = stats["fraction"]

        print(f"Resolution {r:7d} bp | over 1k contact = {frac:.3f}")

        # ---- Choose your stopping condition ----
        perct_cutoff = 0.8
        if frac >= perct_cutoff:
            print("\n*** Recommended minimum resolution:", r, "bp ***")
            return r

    print("\nNo suitable resolution found in the scanned range.")
    return None


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit('python %s <pairs>' % (sys.argv[0]))
    
    pairs = sys.argv[1]
    resolutions = []

    # 1 kb → 10 kb step 1 kb
    resolutions.extend(range(1_000, 10_000 + 1, 1_000))
    # 20 kb → 100 kb step 10 kb
    resolutions.extend(range(20_000, 100_000 + 1, 10_000))
    # 200 kb → 1 Mb step 100 kb
    resolutions.extend(range(200_000, 1_000_000 + 1, 100_000))
    # 2 Mb → 10 Mb step 1 Mb
    resolutions.extend(range(2_000_000, 10_000_000 + 1, 1_000_000))

    gradient_scan(pairs, resolutions)
