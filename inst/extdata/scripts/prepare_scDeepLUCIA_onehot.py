"""
prepare_scDeepLUCIA_inputs.py
------------------------------
Prepares one-hot encoded genome sequence input files for scDeepLUCIA.

Outputs (per chromosome, per resolution window):
  - <outdir>/<chrom>_seq.npy   : float32 array of shape (n_bins, bin_size, 4)

Usage:
    python prepare_scDeepLUCIA_onehot.py \
        --fa genome.fa \
        --resolution 5000 \
        --outdir ./scDeepLUCIA_inputs \
        [--chromosomes chr1 chr2 ... chr22 chrX chrY] \
        [--chroms_file chrom_sizes.txt]
"""

# python prepare_scDeepLUCIA_onehot.py --fa /home/jou34/work/iGenomes/Mus_musculus/mm10/Sequence/WholeGenomeFasta/genome.fa \
#  --resolution 5000 --outdir mm10_onehot
import argparse
import os
import sys
import numpy as np


# ---------------------------------------------------------------------------
# One-hot encoding table
# ---------------------------------------------------------------------------
MAPPING = {
    "A": [1, 0, 0, 0],
    "C": [0, 1, 0, 0],
    "G": [0, 0, 1, 0],
    "T": [0, 0, 0, 1],
    "N": [1, 1, 1, 1],   # scDeepLUCIA convention: N → all-one
    # Handle lowercase (softmasked repeats)
    "a": [1, 0, 0, 0],
    "c": [0, 1, 0, 0],
    "g": [0, 0, 1, 0],
    "t": [0, 0, 0, 1],
    "n": [1, 1, 1, 1],
}
# Pre-build a lookup array indexed by ord(base) for speed
_ORD_TABLE = np.zeros((128, 4), dtype=np.float32)
for _base, _vec in MAPPING.items():
    _ORD_TABLE[ord(_base)] = _vec


# ---------------------------------------------------------------------------
# FASTA parser (no BioPython dependency)
# ---------------------------------------------------------------------------
def parse_fasta(fa_path: str) -> dict[str, str]:
    """
    Stream-parse a FASTA file and return {chrom: sequence_str}.
    Handles gzipped (.gz) files transparently.
    """
    import gzip

    open_fn = gzip.open if fa_path.endswith(".gz") else open

    sequences: dict[str, str] = {}
    current_chrom: str | None = None
    chunks: list[str] = []

    print(f"[parse_fasta] Reading {fa_path} …")
    with open_fn(fa_path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_chrom is not None:
                    sequences[current_chrom] = "".join(chunks)
                    chunks = []
                # Take the first token of the header as the chrom name
                current_chrom = line[1:].split()[0]
            else:
                chunks.append(line)

    if current_chrom is not None:
        sequences[current_chrom] = "".join(chunks)

    print(f"[parse_fasta] Loaded {len(sequences)} sequences.")
    return sequences


# ---------------------------------------------------------------------------
# One-hot encoding helpers
# ---------------------------------------------------------------------------
def onehot_encode(seq: str) -> np.ndarray:
    """
    Encode a nucleotide string to a (len, 4) float32 array.
    Unknown characters map to [1,1,1,1].
    """
    arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    # Clamp to table bounds (handles chars > 127 gracefully)
    arr = np.clip(arr, 0, 127)
    return _ORD_TABLE[arr]  # (L, 4)


def encode_chromosome(seq: str, resolution: int) -> np.ndarray:
    """
    One-hot encode a chromosome sequence and reshape into bins.

    Returns
    -------
    np.ndarray, shape (n_bins, resolution, 4), dtype float32
        Last bin is zero-padded if len(seq) % resolution != 0.
    """
    n = len(seq)
    n_bins = (n + resolution - 1) // resolution  # ceiling division
    pad_len = n_bins * resolution - n

    encoded = onehot_encode(seq)  # (n, 4)

    if pad_len > 0:
        pad = np.zeros((pad_len, 4), dtype=np.float32)
        encoded = np.concatenate([encoded, pad], axis=0)

    return encoded.reshape(n_bins, resolution, 4)


# ---------------------------------------------------------------------------
# Main preparation routine
# ---------------------------------------------------------------------------
def prepare_inputs(
    fa_path: str,
    resolution: int,
    chromosomes: list[str],
    outdir: str,
) -> None:
    os.makedirs(outdir, exist_ok=True)

    # Parse genome
    seq_dict = parse_fasta(fa_path)

    missing = [c for c in chromosomes if c not in seq_dict]
    if missing:
        print(
            f"[WARNING] The following chromosomes were not found in the FASTA "
            f"and will be skipped: {missing}",
            file=sys.stderr,
        )

    for chrom in chromosomes:
        if chrom not in seq_dict:
            continue

        seq = seq_dict[chrom]
        print(f"[encode] {chrom}  length={len(seq):,}  …", end=" ", flush=True)

        encoded = encode_chromosome(seq, resolution)  # (n_bins, resolution, 4)

        out_path = os.path.join(outdir, f"{chrom}_seq.npy")
        np.save(out_path, encoded)

        print(f"→ {encoded.shape}  saved to {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def default_chromosomes() -> list[str]:
    return [f"chr{i}" for i in range(1, 25)] + ["chrX", "chrY"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare one-hot encoded genome inputs for scDeepLUCIA."
    )
    parser.add_argument("--fa", required=True, help="Path to genome FASTA (.fa or .fa.gz)")
    parser.add_argument(
        "--resolution",
        type=int,
        default=5000,
        help="Bin size in bp (default: 5000)",
    )
    parser.add_argument(
        "--outdir",
        default="./scDeepLUCIA_inputs",
        help="Output directory (default: ./scDeepLUCIA_inputs)",
    )
    parser.add_argument(
        "--chromosomes",
        nargs="+",
        default=None,
        help=(
            "Chromosome names to process (e.g. chr1 chr2 chrX). "
            "Defaults to chr1–chr25 + chrX + chrY."
        ),
    )
    parser.add_argument(
        "--chroms_file",
        default=None,
        help=(
            "Optional tab-separated chrom_sizes file (e.g. from UCSC). "
            "If provided, chromosome list is read from its first column "
            "and filtered by --chromosomes if that flag is also set."
        ),
    )
    return parser.parse_args()


def resolve_chromosomes(args: argparse.Namespace) -> list[str]:
    if args.chroms_file:
        with open(args.chroms_file) as fh:
            file_chroms = [line.split()[0] for line in fh if line.strip()]
        if args.chromosomes:
            # Keep only the intersection, preserving file order
            keep = set(args.chromosomes)
            return [c for c in file_chroms if c in keep]
        return file_chroms

    return args.chromosomes if args.chromosomes else default_chromosomes()


def main() -> None:
    args = parse_args()
    chromosomes = resolve_chromosomes(args)

    print("=" * 60)
    print("scDeepLUCIA onehot seq preparation")
    print(f"  FASTA      : {args.fa}")
    print(f"  Resolution : {args.resolution} bp")
    print(f"  Chromosomes: {chromosomes}")
    print(f"  Output dir : {args.outdir}")
    print("=" * 60)

    prepare_inputs(
        fa_path=args.fa,
        resolution=args.resolution,
        chromosomes=chromosomes,
        outdir=args.outdir,
    )

    print("\nDone. Output files:")
    for f in sorted(os.listdir(args.outdir)):
        p = os.path.join(args.outdir, f)
        size_mb = os.path.getsize(p) / 1e6
        print(f"  {f:30s}  {size_mb:8.2f} MB")


if __name__ == "__main__":
    main()
