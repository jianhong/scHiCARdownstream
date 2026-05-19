"""
prepare_scDeepLUCIA_epi.py
---------------------------
Builds epigenomic (epi) input arrays for scDeepLUCIA from BigWig files
produced by ArchR::getGroupBW().

ArchR names each file after the group (e.g. dpa value), so all BigWigs for
every sample live in ONE flat folder:

    <bw_root>/
        X0-TileSize-25-normMethod-nFrags-ArchR.bw
        X5-TileSize-25-normMethod-nFrags-ArchR.bw
        X7-TileSize-25-normMethod-nFrags-ArchR.bw
        …

The script:
  1. Parses the ArchR filename to extract the sample/group label
     (the token before the first '-TileSize-' — e.g. "X0", "X5", "X7").
  2. Reads each BigWig at 25 bp resolution (the ArchR tileSize) per chrom.
  3. Bins the signal into windows of `bin_size` bp  →  shape (n_bins, sub, 1)
     where sub = bin_size // 25.
  4. Since ArchR ATAC BigWigs are single-mark, the epi array per sample is
     shape (n_bins, sub, 1).  
  5. Saves one .npy per (sample, chromosome).

Output layout:
    <outdir>/
        X0/
            chr1_epi.npy        # shape: (n_bins, sub, n_marks)  float32
            chr2_epi.npy
            …
        X5/
            …

Usage — single R2 ATAC mark (typical ArchR output):
    python prepare_scDeepLUCIA_epi.py \
        --bw_roots  ./bigwigs \
        --marker_name r2_030M \
        --resolution 5000 \
        --outdir    ./scDeepLUCIA_inputs \
        [--chromosomes chr1 chr2 … chrX chrY] \
        [--chrom_sizes hg38.chrom.sizes] \
        [--nan_to 0.0]

Requirements:
    pip install pyBigWig numpy
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np

try:
    import pyBigWig
except ImportError:
    sys.exit(
        "pyBigWig is required.  Install with:\n"
        "  pip install pyBigWig"
    )


# ---------------------------------------------------------------------------
# Core helpers
# ---------------------------------------------------------------------------
BW_STEP = 25  # BigWig signal is assumed to be pre-computed at 25 bp bins


def read_signal_25bp(
    bw_path: str,
    chrom: str,
    chrom_len: int,
    nan_fill: float = 0.0,
) -> np.ndarray:
    """
    Extract values at 25 bp resolution for `chrom` from a BigWig file.

    Returns
    -------
    np.ndarray, shape (n_25bp_bins,), dtype float32
    """
    bw = pyBigWig.open(bw_path)

    if chrom not in bw.chroms():
        bw.close()
        raise KeyError(
            f"Chromosome '{chrom}' not found in {bw_path}.\n"
            f"  Available: {list(bw.chroms().keys())[:10]} …"
        )

    # stats() with type="mean" is far faster than values() for binned data
    values = bw.stats(chrom, 0, chrom_len, type="mean", nBins=chrom_len // BW_STEP)
    bw.close()

    arr = np.array(values, dtype=np.float32)
    nan_mask = np.isnan(arr)
    if nan_mask.any():
        arr[nan_mask] = nan_fill

    return arr


def build_epi(
    chrom: str,
    signal_25bp: np.ndarray,
    bin_size: int,
    n_bins: int,
) -> np.ndarray:
    """
    Reshape a 25 bp-resolution signal vector into (n_bins, sub, 1).

    Parameters
    ----------
    chrom       : chromosome name (used only for error messages)
    signal_25bp : 1-D float32 array of length >= n_bins * (bin_size // 25)
    bin_size    : genomic window size in bp (must be divisible by 25)
    n_bins      : number of bins for this chromosome

    Returns
    -------
    np.ndarray, shape (n_bins, sub, 1), dtype float32
        sub = bin_size // 25
    """
    if bin_size % BW_STEP != 0:
        raise ValueError(
            f"bin_size ({bin_size}) must be divisible by {BW_STEP} (the BigWig step)."
        )

    sub = bin_size // BW_STEP
    epi = np.zeros((n_bins, sub, 1), dtype=np.float32)

    required = n_bins * sub
    if len(signal_25bp) < required:
        # Pad the signal if the chrom length is not a perfect multiple
        pad = np.zeros(required - len(signal_25bp), dtype=np.float32)
        signal_25bp = np.concatenate([signal_25bp, pad])

    for i in range(n_bins):
        start = i * sub
        end   = start + sub
        epi[i, :, 0] = signal_25bp[start:end]

    return epi


def stack_marks(mark_arrays: list[np.ndarray]) -> np.ndarray:
    """
    Concatenate per-mark epi arrays along the last axis.

    Input : list of (n_bins, sub, 1) arrays
    Output: (n_bins, sub, n_marks) float32
    """
    return np.concatenate(mark_arrays, axis=-1)


# ---------------------------------------------------------------------------
# Chrom utilities
# ---------------------------------------------------------------------------
def default_chromosomes() -> list[str]:
    return [f"chr{i}" for i in range(1, 25)] + ["chrX", "chrY"]


def load_chrom_sizes(path: str) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def infer_chrom_len_from_bw(bw_path: str, chrom: str) -> int | None:
    """Return chromosome length recorded in a BigWig header, or None."""
    bw = pyBigWig.open(bw_path)
    length = bw.chroms().get(chrom)
    bw.close()
    return length


# ---------------------------------------------------------------------------
# Per-sample processing
# ---------------------------------------------------------------------------
def process_sample(
    sample_name: str,
    bw_paths: list[str],       # 
    marker: str, # marker_type
    chromosomes: list[str],
    chrom_sizes: dict[str, int],
    bin_size: int,
    outdir: str,
    nan_fill: float,
) -> None:
    sample_out = os.path.join(outdir, sample_name)
    os.makedirs(sample_out, exist_ok=True)

    for chrom in chromosomes:
        if chrom not in chrom_sizes:
            print(f"  [SKIP] {chrom} not in chrom_sizes — skipping.", file=sys.stderr)
            continue

        chrom_len = chrom_sizes[chrom]
        n_bins = (chrom_len + bin_size - 1) // bin_size  # ceiling

        print(f"  {chrom}  len={chrom_len:,}  n_bins={n_bins}", end="")

        mark_epis: list[np.ndarray] = []

        for bw_path in bw_paths:
            try:
                signal = read_signal_25bp(bw_path, chrom, chrom_len, nan_fill)
            except KeyError as exc:
                print(f"\n  [WARNING] {exc}", file=sys.stderr)
                # Fill with zeros for this mark
                sub = bin_size // BW_STEP
                signal = np.zeros(n_bins * sub, dtype=np.float32)

            epi_mark = build_epi(chrom, signal, bin_size, n_bins)  # (n_bins, sub, 1)
            mark_epis.append(epi_mark)

        epi = stack_marks(mark_epis)  # (n_bins, sub, n_marks)

        out_path = os.path.join(sample_out, f"{marker}_{chrom}_epi.npy")
        np.save(out_path, epi)

        print(f"  →  {epi.shape}  {out_path}")


# ---------------------------------------------------------------------------
# ArchR BigWig discovery
# ---------------------------------------------------------------------------
# ArchR getGroupBW() produces filenames of the form:
#   <group>-TileSize-<tileSize>-normMethod-<normMethod>-ArchR.bw
# e.g.  X0-TileSize-25-normMethod-nFrags-ArchR.bw
#
# We parse the group label as everything before the first "-TileSize-" token.
# ---------------------------------------------------------------------------
ARCHR_SUFFIX_MARKER = "-TileSize-"


def parse_archr_sample(filename: str) -> str | None:
    """
    Extract the group/sample label from an ArchR BigWig filename.

    'X0-TileSize-25-normMethod-nFrags-ArchR.bw'  →  'X0'
    'MySample-TileSize-25-normMethod-nFrags-ArchR.bw' → 'MySample'

    Returns None if the filename doesn't match the ArchR pattern.
    """
    stem = Path(filename).stem  # strip .bw / .bigwig
    idx = stem.find(ARCHR_SUFFIX_MARKER)
    if idx == -1:
        return None
    return stem[:idx]


def discover_archr_samples(bw_root: str) -> dict[str, str]:
    """
    Scan a flat directory of ArchR BigWigs and return {sample_label: bw_path}.

    All BigWig files that match the ArchR naming convention are collected;
    non-matching files are silently ignored.
    """
    root = Path(bw_root)
    samples: dict[str, str] = {}

    for f in sorted(root.iterdir()):
        if f.suffix.lower() not in {".bw", ".bigwig"}:
            continue
        label = parse_archr_sample(f.name)
        if label is None:
            continue
        if label in samples:
            print(
                f"[WARNING] Duplicate label '{label}' — keeping first match, "
                f"ignoring {f}",
                file=sys.stderr,
            )
            continue
        samples[label] = str(f)

    return samples  # {sample_label: bw_path}


def discover_multi_mark(
    bw_roots: list[str],
    marker_name: str,
) -> dict[str, list[str]]:
    """
    For each mark root, discover ArchR BigWigs.  Then intersect sample labels
    across all marks so every returned sample has a BigWig for every mark.

    Returns {sample_label: [bw_path_mark0, bw_path_mark1, …]}
    """
    per_mark: list[dict[str, str]] = []
    for root in bw_roots:
        found = discover_archr_samples(root)
        print(f"[discover] ({root}): {len(found)} samples found.")
        per_mark.append(found)

    # Intersect sample labels
    common = set(per_mark[0].keys())
    for m in per_mark[1:]:
        common &= set(m.keys())

    if not common:
        sys.exit(
            "No sample labels are shared across all mark directories. "
            "Check that filenames use the same group names."
        )

    result: dict[str, list[str]] = {}
    for label in sorted(common):
        result[label] = [m[label] for m in per_mark]

    dropped = set(per_mark[0].keys()) - common
    if dropped:
        print(
            f"[WARNING] {len(dropped)} sample(s) missing from ≥1 mark folder "
            f"and will be skipped: {sorted(dropped)}",
            file=sys.stderr,
        )

    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Build scDeepLUCIA epi arrays from ArchR getGroupBW() BigWig files.\n\n"
            "All BigWigs for all samples live in one flat folder per mark.\n"
            "Filenames must follow ArchR convention:\n"
            "  <group>-TileSize-<n>-normMethod-<method>-ArchR.bw"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--bw_roots", nargs="+", required=True,
        help=(
            "One or more flat directory per epigenomic mark."
        ),
    )
    p.add_argument(
        "--marker_name", type=str, required=True,
        help=(
            "Human-readable mark labels, "
            "(e.g. r2_30M.  Used only for logging."
        ),
    )
    p.add_argument("--resolution", type=int, default=5000,
                   help="Bin size in bp, must be divisible by 25 (default: 5000).")
    p.add_argument("--outdir", default="./scDeepLUCIA_inputs",
                   help="Output root directory (default: ./scDeepLUCIA_inputs).")
    p.add_argument("--chromosomes", nargs="+", default=None,
                   help="Chromosomes to process (default: chr1–chr25 + chrX + chrY).")
    p.add_argument("--chrom_sizes", default=None,
                   help="UCSC-style chrom.sizes file.  Auto-inferred from BigWig header if omitted.")
    p.add_argument("--nan_to", type=float, default=0.0,
                   help="Value to substitute for NaN signal (default: 0.0).")
    p.add_argument("--samples", nargs="+", default=None,
                   help="Optional allow-list of sample labels to process (e.g. X0 X5).  "
                        "Defaults to all discovered samples.")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    chromosomes = args.chromosomes or default_chromosomes()

    # ---- Discover samples across all mark folders --------------------------
    samples = discover_multi_mark(args.bw_roots, args.marker_name)

    if args.samples:
        keep = set(args.samples)
        missing = keep - set(samples)
        if missing:
            print(f"[WARNING] Requested samples not found: {sorted(missing)}", file=sys.stderr)
        samples = {k: v for k, v in samples.items() if k in keep}

    if not samples:
        sys.exit("No samples to process after filtering.")

    # ---- Resolve chrom sizes -----------------------------------------------
    if args.chrom_sizes:
        chrom_sizes = load_chrom_sizes(args.chrom_sizes)
    else:
        first_bw = next(iter(samples.values()))[0]
        bw = pyBigWig.open(first_bw)
        chrom_sizes = dict(bw.chroms())
        bw.close()
        print(f"[info] Inferred chrom sizes from {first_bw}")

    chromosomes = [c for c in chromosomes if c in chrom_sizes]
    if not chromosomes:
        sys.exit("None of the requested chromosomes found in chrom_sizes.")

    print("=" * 60)
    print("scDeepLUCIA epi preparation  (ArchR BigWig mode)")
    print(f"  Marks       : {args.marker_name}")
    print(f"  BigWig roots: {args.bw_roots}")
    print(f"  Resolution  : {args.resolution} bp  →  sub={args.resolution // BW_STEP}")
    print(f"  Chromosomes : {chromosomes}")
    print(f"  Samples     : {sorted(samples.keys())}  (n={len(samples)})")
    print(f"  Output dir  : {args.outdir}")
    print("=" * 60)

    os.makedirs(args.outdir, exist_ok=True)

    for sample_name, bw_paths in samples.items():
        print(f"\n[sample] {sample_name}")
        for bw_path in bw_paths:
            print(f"         {bw_path}")
        process_sample(
            sample_name=sample_name,
            bw_paths=bw_paths,
            marker=args.marker_name,
            chromosomes=chromosomes,
            chrom_sizes=chrom_sizes,
            bin_size=args.resolution,
            outdir=args.outdir,
            nan_fill=args.nan_to,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
