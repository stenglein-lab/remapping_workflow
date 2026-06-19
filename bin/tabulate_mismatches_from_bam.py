#!/usr/bin/env python3
"""
tabulate_mismatches_from_bam.py

Compute per-reference-position match/mismatch counts and frequencies from a
BAM file, using pysam (a maintained Python wrapper around htslib).

For every position in the reference (within the requested region, or the
entire BAM if no region is given) -- including positions with ZERO coverage
-- this script reports:
    - the reference base at that position
    - total depth (number of aligned, non-deleted bases at that position)
    - number of bases matching the reference
    - raw per-base counts (A, C, G, T, N)
    - per-base MISMATCH counts (A, C, G, T, N), i.e. how many reads carried
      that base at this position AND that base differs from the reference
      (so the count for the reference base itself is always 0)
    - per-base mismatch FREQUENCIES (mismatch count / depth)

Requirements:
    - pysam
    - A coordinate-SORTED and INDEXED BAM file (required; see below)
    - A reference FASTA file with an accompanying .fai index
      (create with: samtools faidx ref.fa)

Why sorted + indexed is required:
    This script needs random access into the BAM (via region-based pileup)
    to enumerate every reference position, including those with zero
    coverage. That requires a coordinate-sorted, indexed BAM
    (samtools sort + samtools index).

Why a ref/BAM mismatch is a hard failure:
    If a reference position falls outside the bounds of the cached reference
    sequence for a contig, that almost always means the BAM was aligned
    against a different reference than the one provided here. Silently
    filling in "N" would hide a serious, easy-to-miss error, so the script
    aborts instead.

Example usage:
    python bam_mismatch_profile.py \
        --bam aligned.sorted.bam \
        --ref ref.fa \
        --out mismatch_profile.tsv

    # Restrict to a specific region:
    python bam_mismatch_profile.py \
        --bam aligned.sorted.bam \
        --ref ref.fa \
        --region chr1:1000-2000 \
        --out mismatch_profile.tsv

    # Only report positions with at least 10x depth (zero-coverage positions
    # below this threshold are still excluded, since min-depth filtering is
    # applied after zero-coverage positions are generated):
    python bam_mismatch_profile.py \
        --bam aligned.sorted.bam \
        --ref ref.fa \
        --min-depth 10 \
        --out mismatch_profile.tsv
"""

import argparse
import os
import sys
from collections import Counter

import pysam


BASES = ("A", "C", "G", "T", "N")


# --------------------------------------------------------------------------- #
# CLI argument parsing
# --------------------------------------------------------------------------- #

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-position match/mismatch counts and frequencies "
            "from a BAM file using pysam."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--bam",
        required=True,
        help="Path to input BAM file. Must be coordinate-sorted and indexed "
             "(.bai) -- this is required, not optional, because the script "
             "needs random access to enumerate zero-coverage positions.",
    )
    parser.add_argument(
        "--ref",
        required=True,
        help="Path to reference FASTA file. Must have a .fai index "
             "(create with `samtools faidx ref.fa` if missing).",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Path to output TSV file. If omitted, results are written to "
             "stdout.",
    )
    parser.add_argument(
        "--region",
        default=None,
        help="Optional region string to restrict analysis, e.g. "
             "'chr1:1000-2000' or just 'chr1' for a whole contig. "
             "If omitted, every contig in the BAM header is processed in "
             "full.",
    )
    parser.add_argument(
        "--min-depth",
        type=int,
        default=0,
        help="Minimum depth required to report a position. Positions with "
             "lower depth (including zero-coverage positions, if > 0) are "
             "skipped. Default 0 means all positions are reported, "
             "including zero-coverage ones.",
    )
    parser.add_argument(
        "--min-base-quality",
        type=int,
        default=20,
        help="Minimum base quality (Phred score) for a base to be counted. "
             "Bases below this threshold are excluded from all counts.",
    )
    parser.add_argument(
        "--min-mapping-quality",
        type=int,
        default=20,
        help="Minimum mapping quality (MAPQ) for a read to be included in "
             "the pileup at all.",
    )
    parser.add_argument(
        "--stepper",
        choices=["samtools", "all", "nofilter"],
        default="samtools",
        help="pysam pileup 'stepper' option controlling which reads are "
             "skipped. 'samtools' mimics samtools default filtering "
             "(skips unmapped/secondary/QC-fail/duplicate reads by "
             "default within pysam). 'all' / 'nofilter' include more reads; "
             "see pysam docs for exact semantics.",
    )

    return parser.parse_args()


# --------------------------------------------------------------------------- #
# Validation helpers
# --------------------------------------------------------------------------- #

def require_sorted_bam(bam):
    """
    Fail if the BAM is not coordinate-sorted (per its header) or is
    not indexed. Both are required because this script does region-based
    pileup queries to enumerate every position, including zero-coverage
    ones.
    """
    header_dict = bam.header.to_dict()
    sort_order = header_dict.get("HD", {}).get("SO")

    if sort_order != "coordinate":
        sys.stderr.write(
            "ERROR: input BAM does not declare itself coordinate-sorted "
            f"(found HD/SO = {sort_order!r}). This script requires a "
            "coordinate-sorted, indexed BAM.\n"
            "Sort with: samtools sort -o sorted.bam input.bam\n"
        )
        sys.exit(1)

    if not bam.has_index():
        sys.stderr.write(
            "ERROR: input BAM has no index. This script requires an "
            "indexed BAM for region-based pileup queries.\n"
            "Index with: samtools index sorted.bam\n"
        )
        sys.exit(1)


def parse_region(region_str, bam):
    """
    Parse a region string of the form 'chrom:start-stop' or just 'chrom'
    into (contig, start_0based, stop_exclusive). If no stop is given (or no
    region at all), defaults to the full contig length from the BAM header.

    Returns a list of (contig, start, stop) tuples to process. If no region
    is given, this is every contig in the BAM header, in full.
    """
    if region_str is None:
        # Process every contig in the BAM header, in full.
        return [
            (contig, 0, bam.get_reference_length(contig))
            for contig in bam.references
        ]

    if ":" in region_str:
        contig, coords = region_str.split(":", 1)
        try:
            start_str, stop_str = coords.split("-")
            # Convert from 1-based inclusive (user-facing) to 0-based
            # half-open (pysam/htslib internal convention).
            start = int(start_str) - 1
            stop = int(stop_str)
        except ValueError:
            sys.stderr.write(
                f"ERROR: could not parse region '{region_str}'. "
                "Expected format 'chrom:start-stop' or just 'chrom'.\n"
            )
            sys.exit(1)
    else:
        contig = region_str
        start = 0
        stop = bam.get_reference_length(contig)

    if contig not in bam.references:
        sys.stderr.write(
            f"ERROR: contig '{contig}' not found in BAM header.\n"
        )
        sys.exit(1)

    return [(contig, start, stop)]


# --------------------------------------------------------------------------- #
# Core logic
# --------------------------------------------------------------------------- #

def compute_mismatch_profile(args):
    """
    Generator yielding per-position match/mismatch statistics dictionaries,
    covering every position in the requested region(s) -- including
    zero-coverage positions.
    """

    # Open the BAM file for reading.
    bam = pysam.AlignmentFile(args.bam, "rb")

    # Hard requirement: coordinate-sorted + indexed BAM.
    require_sorted_bam(bam)

    # Open the reference FASTA. Requires a .fai index alongside it.
    try:
        ref = pysam.FastaFile(args.ref)
    except Exception as e:
        sys.stderr.write(
            f"ERROR: could not open reference FASTA '{args.ref}'. "
            f"Make sure a .fai index exists (samtools faidx {args.ref}).\n"
            f"Original error: {e}\n"
        )
        sys.exit(1)

    regions_to_process = parse_region(args.region, bam)

    pileup_common_kwargs = dict(
        stepper=args.stepper,
        min_base_quality=args.min_base_quality,
        min_mapping_quality=args.min_mapping_quality,
        truncate=True,
        ignore_overlaps=False,  # raw per-base counts, not overlap-collapsed
    )

    for contig, region_start, region_stop in regions_to_process:

        # Fetch the full contig sequence once. This is far more efficient
        # than fetching base-by-base, and lets us validate every position
        # in this contig's range up front.
        contig_seq = ref.fetch(contig)
        contig_len_in_ref = len(contig_seq)

        # If the requested range extends beyond what the reference FASTA
        # actually contains for this contig, the BAM and reference almost
        # certainly do not correspond to each other (e.g. different genome
        # builds, truncated FASTA, wrong file). Fail rather than
        # silently substituting "N".
        if region_stop > contig_len_in_ref:
            sys.stderr.write(
                f"ERROR: reference/BAM mismatch detected for contig "
                f"'{contig}'. BAM header reports length "
                f"{bam.get_reference_length(contig)}, but reference FASTA "
                f"'{args.ref}' only has {contig_len_in_ref} bases for this "
                f"contig. Requested range up to position {region_stop} "
                "exceeds the reference sequence length.\n"
                "This usually means the BAM was aligned against a "
                "different reference than the one provided. Aborting.\n"
            )
            sys.exit(1)

        # Build a lookup of position -> per-base counts from the pileup.
        # Only positions with at least one covering read appear in the
        # pileup output; positions absent from this dict have zero coverage
        # and are filled in below.
        position_counts = {}

        pileup_iter = bam.pileup(
            contig=contig,
            start=region_start,
            stop=region_stop,
            **pileup_common_kwargs,
        )

        for pileup_column in pileup_iter:
            pos0 = pileup_column.reference_pos

            base_counts = Counter()
            for pileup_read in pileup_column.pileups:
                # Skip deletions and reference-skips (introns); there is no
                # base to compare against the reference at these positions
                # for this particular read.
                if pileup_read.is_del or pileup_read.is_refskip:
                    continue

                aln = pileup_read.alignment
                query_pos = pileup_read.query_position
                query_name = aln.query_name
                if query_pos is None:
                    continue

                base = aln.query_sequence[query_pos].upper()
                #MARKMARK
                # if pos0 > 842 and pos0 < 848:
                    # print(f"{contig}\t{pos0}\t{query_name}\t{query_pos}\t{base}\n", file=sys.stderr)
                base_counts[base] += 1

            position_counts[pos0] = base_counts

        # --- Walk every position in the requested range, in order ------ #
        for pos0 in range(region_start, region_stop):

            # This should be unreachable given the up-front length check
            # above, but kept as a defensive double-check per-position.
            if pos0 >= contig_len_in_ref:
                sys.stderr.write(
                    f"ERROR: reference/BAM mismatch detected at "
                    f"{contig}:{pos0 + 1} (position beyond reference "
                    "sequence length). Aborting.\n"
                )
                sys.exit(1)

            ref_base = contig_seq[pos0].upper()

            base_counts = position_counts.get(pos0, Counter())

            # MARKMARK
            # if pos0 > 842 and pos0 < 848:
               # print(f"{contig}\t{pos0}\t{dict(base_counts)}\n", file=sys.stderr)

            depth = sum(base_counts.values())

            if depth < args.min_depth:
                continue

            matches = base_counts.get(ref_base, 0)

            # Per-base mismatch counts: count of each base observed at this
            # position, EXCLUDING the count for the reference base itself
            # (which is by definition a match, not a mismatch).
            mismatch_counts = {
                b: (base_counts.get(b, 0) if b != ref_base else 0)
                for b in BASES
            }

            # Per-base mismatch frequencies, relative to total depth at
            # this position. 0.0 when depth is 0 (zero-coverage position).
            mismatch_freqs = {
                b: (mismatch_counts[b] / depth) if depth > 0 else 0.0
                for b in BASES
            }

            yield {
                "chrom": contig,
                "pos": pos0 + 1,  # report 1-based position
                "ref_base": ref_base,
                "depth": depth,
                "matches": matches,
                # Raw per-base counts (includes the reference base's own count)
                "A": base_counts.get("A", 0),
                "C": base_counts.get("C", 0),
                "G": base_counts.get("G", 0),
                "T": base_counts.get("T", 0),
                "N": base_counts.get("N", 0),
                # Per-base mismatch counts (reference base's own count is 0)
                "mismatch_A": mismatch_counts["A"],
                "mismatch_C": mismatch_counts["C"],
                "mismatch_G": mismatch_counts["G"],
                "mismatch_T": mismatch_counts["T"],
                "mismatch_N": mismatch_counts["N"],
                # Per-base mismatch frequencies
                "mismatch_freq_A": mismatch_freqs["A"],
                "mismatch_freq_C": mismatch_freqs["C"],
                "mismatch_freq_G": mismatch_freqs["G"],
                "mismatch_freq_T": mismatch_freqs["T"],
                "mismatch_freq_N": mismatch_freqs["N"],
            }

    bam.close()
    ref.close()


# --------------------------------------------------------------------------- #
# Output writing
# --------------------------------------------------------------------------- #

def write_output(records, out_path):
    """
    Write per-position records to a TSV file, or to stdout if out_path is None.
    """
    # header = [
        # "chrom", "pos", "ref_base", "depth", "matches",
        # "A", "C", "G", "T", "N",
        # "mismatch_A", "mismatch_C", "mismatch_G", "mismatch_T", "mismatch_N",
        # "mismatch_freq_A", "mismatch_freq_C", "mismatch_freq_G",
        # "mismatch_freq_T", "mismatch_freq_N",
    # ]
    header = [
        "chrom", "pos", "ref_base", "depth", 
        "A", "C", "G", "T", "N"
    ]

    out_fh = open(out_path, "w") if out_path else sys.stdout

    try:
        # don't output header 
        # out_fh.write("\t".join(header) + "\n")
        n_rows = 0
        for rec in records:
            row = [
                rec["chrom"],
                str(rec["pos"]),
                rec["ref_base"],
                str(rec["depth"]),
                # str(rec["matches"]),
                str(rec["A"]),
                str(rec["C"]),
                str(rec["G"]),
                str(rec["T"]),
                str(rec["N"]),
                # str(rec["mismatch_A"]),
                # str(rec["mismatch_C"]),
                # str(rec["mismatch_G"]),
                # str(rec["mismatch_T"]),
                # str(rec["mismatch_N"]),
                # f"{rec['mismatch_freq_A']:.6f}",
                # f"{rec['mismatch_freq_C']:.6f}",
                # f"{rec['mismatch_freq_G']:.6f}",
                # f"{rec['mismatch_freq_T']:.6f}",
                # f"{rec['mismatch_freq_N']:.6f}",
            ]
            out_fh.write("\t".join(row) + "\n")
            n_rows += 1
        return n_rows
    finally:
        if out_path:
            out_fh.close()


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #

def main():
    args = parse_args()

    records = compute_mismatch_profile(args)
    n_rows = write_output(records, args.out)

    dest = args.out if args.out else "stdout"
    sys.stderr.write(f"Done. Wrote {n_rows} positions to {dest}.\n")


if __name__ == "__main__":
    main()
