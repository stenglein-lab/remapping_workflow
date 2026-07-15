#!/usr/bin/env python3
"""
tabulate_mismatches_from_bam.py

Compute per-reference-position match/mismatch counts and frequencies from a
BAM file, using pysam (a python wrapper around htslib).

For every position in every reference sequence, including positions with no coverage,
this script reports:
    - the reference sequence name
    - the position
    - the reference base at that position
    - read depth at that position (number of aligned, non-deleted bases at that position)
    - per-base counts (A, C, G, T, N)

Additionally, this script can optionally output mismatch type summaries
(e.g. how many C->A substitutions were observed):
    - --mismatch-types-per-ref-out: counts per (contig, ref_base, alt_base)
    - --mismatch-types-total-out: counts per (ref_base, alt_base), totalled
      across every contig/region processed

These type summaries are accumulated over every position considered,
independent of --min-depth filtering applied to the main per-position
report, and only classify A/C/G/T reference bases (positions with an 'N'
reference base are excluded, since "N->X" is not a meaningful substitution
type).

Requirements:
    - pysam
    - A sorted and indexed bam file with accompanying .bai index
      (created with: samtools sort and samtools index)
    - A reference FASTA file with accompanying .fai index
      (created with: samtools faidx ref.fa)

Why sorted + indexed is required:
    This script needs random access into the BAM (via region-based pileup)
    to enumerate every reference position, including those with zero
    coverage. That requires a coordinate-sorted, indexed BAM
    (samtools sort + samtools index).

It is possible to restrict output to a particular refseq/region

Example usage:
    python tabulate_mismatches_from_bam.py \
        --bam aligned.sorted.bam \
        --ref ref.fa \
        --out mismatch_profile.tsv

    # Restrict to a specific region:
    python tabulate_mismatches_from_bam.py \
        --bam aligned.sorted.bam \
        --ref ref.fa \
        --region chr1:1000-2000 \
        --out mismatch_profile.tsv

    # Only report positions with at least 10x depth 
    python tabulate_mismatches_from_bam.py \
        --bam aligned.sorted.bam \
        --ref ref.fa \
        --min-depth 10 \
        --out mismatch_profile.tsv

    # Also emit mismatch-type summaries (e.g. C->A counts) per contig and
    # totalled overall:
    python tabulate_mismatches_from_bam.py \
        --bam aligned.sorted.bam \
        --ref ref.fa \
        --out mismatch_profile.tsv \
        --mismatch-types-per-ref-out mismatch_types_per_ref.tsv \
        --mismatch-types-total-out mismatch_types_total.tsv
"""

# standard python libs
import argparse
import os
import sys
from collections import Counter, defaultdict

# pysam python lib
import pysam


BASES = ("A", "C", "G", "T", "N")


# --------------------------------------------------------------------------- #
# CLI argument parsing
# --------------------------------------------------------------------------- #

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-position match/mismatch counts and frequencies "
            "from a bam file using pysam."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--bam",
        required=True,
        help="Path to input bam file. Must be sorted and indexed (with .bai)"
             "(create with `samtools sort` and `samtools index` if missing)",
    )
    parser.add_argument(
        "--ref",
        required=True,
        help="Path to reference fasta file. Must have a .fai index "
             "(create with `samtools faidx ref.fa` if missing).",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Path to output TSV file. If omitted, results are written to stdout.",
    )
    parser.add_argument(
        "--prefix",
        default=None,
        help="Optional text to output as a 1st column of all tsv output.",
    )
    parser.add_argument(
        "--region",
        default=None,
        help="Optional region string to restrict analysis, e.g. "
             "'chr1:1000-2000' or just 'chr1' for a whole reference sequence. "
             "coordinates are as in `samtools`, and are 1-based",
    )
    parser.add_argument(
        "--min_depth",
        type=int,
        default=0,
        help="Minimum coverage depth required to report a position. Positions with "
             "Default 0 means all positions are reported.",
    )
    parser.add_argument(
        "--max_depth",
        type=int,
        default=20000,
        help="Maximum coverage depth allowed at a position. "
             "Passed to the pysam pileup() function max_depth parameter.",
    )
    parser.add_argument(
        "--min_base_quality",
        type=int,
        default=20,
        help="Minimum basecall quality score for a base to be counted. "
             "Bases below this threshold are excluded from all counts.",
    )
    parser.add_argument(
        "--min_mapping_quality",
        type=int,
        default=20,
        help="Minimum mapping quality (MAPQ) for a read to be included in analysis",
    )
    parser.add_argument(
        "--stepper",
        choices=["samtools", "all", "nofilter"],
        default="samtools",
        help="pysam pileup 'stepper' option controlling which reads are "
             "skipped. 'samtools' mimics samtools default filtering "
             "(skips unmapped/secondary/QC-fail/duplicate reads by "
             "default within pysam). 'all' / 'nofilter' include more reads; "
             "see pysam docs for more information, at: " 
             "https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile.pileup",
    )
    parser.add_argument(
        "--mismatch_types_per_ref_out",
        default=None,
        help="Optional path to write a TSV summarizing mismatch type counts "
             "(e.g. C->A, G->T) per reference sequence. One row "
             "per (contig, ref_base, alt_base) combination observed. "
             "Only A/C/G/T reference bases are classified -- positions "
             "where the reference base is 'N' are excluded from this "
             "summary, since 'N->X' is not a meaningful substitution type.",
    )
    parser.add_argument(
        "--mismatch_types_by_pos_out",
        default=None,
        help="Optional path to write a TSV summarizing mismatch type counts "
             "(e.g. C->A, G->T) per reference sequence as a function of position in read"
             "Position in read uses positive integers to indicate position in read (1-based)"
             "and negative integers to indicate distance from read end." 
             "Position 1 = 1st base in read; -1 = last base in read, -2 = penultimate base in read.", 
    )
    parser.add_argument(
        "--mismatch_types_total_out",
        default=None,
        help="Optional path to write a TSV summarizing mismatch type counts "
             "(e.g. C->A, G->T) totalled across all reference sequences in "
             "the processed region(s). One row per (ref_base, alt_base) "
             "combination observed. Same A/C/G/T-only restriction as "
             "--mismatch-types-per-ref-out.",
    )

    return parser.parse_args()


# --------------------------------------------------------------------------- #
# Validation helpers
# --------------------------------------------------------------------------- #

def require_sorted_bam(bam):
    """
    Fail if the BAM is not sorted (per its header) or is not indexed. 
    Both are required because this script does region-based
    pileup queries to enumerate every position, including zero-coverage
    positions.
    """
    header_dict = bam.header.to_dict()
    sort_order = header_dict.get("HD", {}).get("SO")

    if sort_order != "coordinate":
        sys.stderr.write(
            "ERROR: input BAM is not sorted "
            f"(found HD/SO = {sort_order!r}). This script requires a "
            "sorted, indexed BAM.\n"
        )
        sys.exit(1)

    if not bam.has_index():
        sys.stderr.write(
            "ERROR: input BAM is not indexed. This script requires an "
            "indexed BAM (.bai file)."
            "Index with: samtools index sorted.bam\n"
        )
        sys.exit(1)

def parse_region(region_str, bam):
    """
    Optional: Parse a region string of the form 'chrom:start-stop' or just 'chrom'
    into (contig, start_0based, stop_exclusive). 

    Returns a list of (contig, start, stop) tuples to process. 
    """
    if region_str is None:
        # Process every refseq in the bam, in full.
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


def compute_mismatch_profile(args, mismatch_types_per_ref, mismatch_types_total, mismatch_types_by_pos):
    """
    Generator that yields per-refseq/per-position match/mismatch counts

    As a side effect, this also accumulates mismatch type counts (e.g.
    C->A, G->T) into the caller-supplied accumulators:
        - mismatch_types_per_ref: defaultdict(Counter), keyed by
          contig -> Counter({(ref_base, alt_base): count, ...})
        - mismatch_types_total: Counter({(ref_base, alt_base): count, ...}),
          totalled across every contig/region processed.

    Both accumulators are mutated in place; the caller creates them before
    calling this generator and inspects them after the generator is fully
    consumed (e.g. after write_output() has iterated over all records).

    Only A/C/G/T reference bases are classified into mismatch types --
    positions where the reference base is 'N' are skipped for this
    particular accounting, since "N->X" is not a meaningful substitution
    type. This does NOT affect the main per-position output, which still
    reports counts at reference-N positions as usual.
    """

    # Open the BAM file for reading.
    bam = pysam.AlignmentFile(args.bam, "rb")

    # Confirm bam is sorted and indexed 
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

    # optionally limit to certain refseqs/regions
    regions_to_process = parse_region(args.region, bam)

    # define dictionary with pysam args
    pileup_common_kwargs = dict(
        stepper             = args.stepper,
        max_depth           = args.max_depth,
        min_base_quality    = args.min_base_quality,
        min_mapping_quality = args.min_mapping_quality,
        truncate            = True,   # Pysam doc: By default, the samtools pileup engine outputs all reads overlapping a region. If truncate is True and a region is given, only columns in the exact region specified are returned
        ignore_overlaps     = False,  # Pysam doc: If set to True, detect if read pairs overlap and only take the higher quality base.
    )

    # for each refseq/region
    for contig, region_start, region_stop in regions_to_process:

        # Fetch the reference sequence once. 
        contig_seq = ref.fetch(contig)
        contig_len_in_ref = len(contig_seq)

        # If the requested range extends beyond the reference sequence,
        # the bam and refseq likely do not correspond to each other 
        # (e.g. different genome builds, truncated FASTA, wrong file). 
        # Throw error 
        if region_stop > contig_len_in_ref:
            sys.stderr.write(
                f"ERROR: reference/BAM mismatch detected for reference sequnce "
                f"'{contig}'. BAM header reports length "
                f"{bam.get_reference_length(contig)}, but reference FASTA "
                f"'{args.ref}' only has {contig_len_in_ref} bases for this "
                f"sequence. Requested range up to position {region_stop} "
                "exceeds the reference sequence length.\n"
                "This usually means the BAM was aligned against a "
                "different reference than the one provided. Aborting.\n"
            )
            sys.exit(1)

        # create a "pileup" of the region and use to 
        # build a lookup of position -> per-base counts 
        # Only positions with at least one covering read appear in the
        # pileup output; positions absent from this dict have zero coverage
        # and are filled in below.
        position_counts = {}

        # call pysam bam.pileup() 
        pileup_iter = bam.pileup(
            contig = contig,
            start  = region_start,
            stop   = region_stop,
            **pileup_common_kwargs,
        )

        # for each position in the refseq
        for pileup_column in pileup_iter:
            # the position in the reference sequence
            pos0 = pileup_column.reference_pos

            base_counts = Counter()
            # for each read in the pileup
            for pileup_read in pileup_column.pileups:
                # Skip deletions and reference-skips; there is no
                # base to compare against the reference at these positions
                # for this read.
                if pileup_read.is_del or pileup_read.is_refskip:
                    continue

                # the "aligned segment" object
                aln = pileup_read.alignment
                # the position in the read (query)
                query_pos = pileup_read.query_position
                # the read (query) name
                query_name = aln.query_name
                if query_pos is None:
                    continue

                # get the base in this read at this position
                base = aln.query_sequence[query_pos].upper()

                # keep track of matches/mismatches as a function of position in read
                # get the reference sequence base at this position
                ref_base = contig_seq[pos0].upper()
                # use 1-based positions
                distance_from_begin = query_pos + 1
                # use negative numbers for distances from end of read
                distance_from_end   = 0-(aln.query_length - query_pos)

                # use a tuple (refseq_name, distance from begin (+ integers) or end(- integers)) as key for dictionary
                mismatch_types_by_pos[(contig, distance_from_begin)][(ref_base, base)] += 1
                mismatch_types_by_pos[(contig, distance_from_end)  ][(ref_base, base)] += 1

                # tabulate
                base_counts[base] += 1

            # store base counts at this position of this refseq
            position_counts[pos0] = base_counts


        # --- Walk every position in the requested range, in order ------ #
        for pos0 in range(region_start, region_stop):

            # double check position doesn't exceed refseq length
            if pos0 >= contig_len_in_ref:
                sys.stderr.write(
                    f"ERROR: reference/BAM mismatch detected at "
                    f"{contig}:{pos0 + 1} (position beyond reference "
                    "sequence length). Aborting.\n"
                )
                sys.exit(1)

            # the reference sequence base at this position
            ref_base = contig_seq[pos0].upper()

            # the count of mapped bases at this position
            base_counts = position_counts.get(pos0, Counter())

            # depth at this position (of reads passing filters, sufficient mapQ, basecall Q, etc)
            depth = sum(base_counts.values())

            # --- Accumulate mismatch type counts (e.g. C->A) in total and per refseq 
            # This happens before the --min-depth filter below, so the
            # mismatch-type summaries reflect every observed mismatch in
            # the processed region(s), independent of whatever depth
            # cutoff is applied to the main per-position report.
            # Only classify positions with an A/C/G/T reference base;
            # "N->X" is not a meaningful substitution type and is skipped.
            if ref_base in ("A", "C", "G", "T"):
                for alt_base in ("A", "C", "G", "T"):
                    # if alt_base == ref_base:
                        # continue
                    count = base_counts.get(alt_base, 0)
                    # if count == 0:
                        # continue
                    mismatch_types_per_ref[contig][(ref_base, alt_base)] += count
                    mismatch_types_total[(ref_base, alt_base)] += count


            # if a minimum depth to report is specified
            if depth < args.min_depth:
                continue

            # number of reads matching the refseq
            matches = base_counts.get(ref_base, 0)

            # Per-position mismatch counts: counts of each base observed at this
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
# Output counts
# --------------------------------------------------------------------------- #

def write_output(records, out_path, prefix):
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
            if prefix is not None:
                out_fh.write(f"{prefix}\t")
            out_fh.write("\t".join(row) + "\n")
            n_rows += 1
        return n_rows
    finally:
        if out_path:
            out_fh.close()



def write_mismatch_type_summary_per_ref(mismatch_types_per_ref, out_path, prefix):
    """
    Write a TSV of mismatch type counts (e.g. C->A, G->T) per reference
    sequence (contig). One row per (refseq, ref_base, alt_base) combination
    that was observed at least once. 
    """
    out_fh = open(out_path, "w")
    try:
        # out_fh.write("chrom\tref_base\talt_base\tcount\n")
        n_rows = 0
        # Sort for deterministic, reproducible output
        for contig in sorted(mismatch_types_per_ref.keys()):
            type_counts = mismatch_types_per_ref[contig]
            for (ref_base, alt_base) in sorted(type_counts.keys()):
                count = type_counts[(ref_base, alt_base)]
                if prefix is not None:
                    out_fh.write(f"{prefix}\t")
                out_fh.write(f"{contig}\t{ref_base}\t{alt_base}\t{count}\n")
                n_rows += 1
        return n_rows
    finally:
        out_fh.close()


def write_mismatch_type_by_pos(mismatch_types_by_pos, out_path, prefix):
    """
    Write a TSV of mismatch type counts (e.g. C->A, G->T) per reference
    sequence per position relative to the beginning and ends of 
    Distances from the ends of reads use negative numbers.
    """
    out_fh = open(out_path, "w")
    try:
        # out_fh.write("chrom\tdistance_from_end\tref_base\talt_base\tcount\n")
        n_rows = 0
        # sort keys
        sorted_mismatches = sorted(mismatch_types_by_pos.items(), key=lambda item: item[0])
        for (contig, distance), type_counts in sorted_mismatches:
            for (ref_base, alt_base) in sorted(type_counts.keys()):
                count = type_counts[(ref_base, alt_base)]
                if prefix is not None:
                    out_fh.write(f"{prefix}\t")
                out_fh.write(f"{contig}\t{distance}\t{ref_base}\t{alt_base}\t{count}\n")
                n_rows += 1
        return n_rows
    finally:
        out_fh.close()


def write_mismatch_type_summary_total(mismatch_types_total, out_path, prefix):
    """
    Write a TSV of mismatch type counts (e.g. C->A, G->T) totalled across
    every contig/region processed. One row per (ref_base, alt_base)
    combination observed at least once.
    """
    out_fh = open(out_path, "w")
    try:
        # out_fh.write("ref_base\talt_base\tcount\n")
        n_rows = 0
        for (ref_base, alt_base) in sorted(mismatch_types_total.keys()):
            count = mismatch_types_total[(ref_base, alt_base)]
            if prefix is not None:
                out_fh.write(f"{prefix}\t")
            out_fh.write(f"{ref_base}\t{alt_base}\t{count}\n")
            n_rows += 1
        return n_rows
    finally:
        out_fh.close()


# --------------------------------------------------------------------------- #
# main entry point
# --------------------------------------------------------------------------- #

def main():
    args = parse_args()

    # these accumulators will be populated as a side effect of compute_mismatch_profile().
    # these are created here, then mutated during iteration, and inspected after the generator
    # is fully consumed by write_output() below.
    mismatch_types_per_ref = defaultdict(Counter)
    mismatch_types_total   = Counter()
    mismatch_types_by_pos  = defaultdict(Counter)

    records = compute_mismatch_profile(args, mismatch_types_per_ref, mismatch_types_total, mismatch_types_by_pos)
    n_rows = write_output(records, args.out, args.prefix)

    dest = args.out if args.out else "stdout"
    sys.stderr.write(f"Wrote {n_rows} positions to {dest}.\n")

    # At this point the generator has been fully consumed (write_output
    # iterated over every record), so the accumulators are complete.
    if args.mismatch_types_per_ref_out:
        n_type_rows = write_mismatch_type_summary_per_ref(
            mismatch_types_per_ref, args.mismatch_types_per_ref_out, args.prefix
        )
        sys.stderr.write(
            f"Wrote {n_type_rows} per-reference mismatch-type rows to "
            f"{args.mismatch_types_per_ref_out}.\n"
        )

    if args.mismatch_types_by_pos_out:
        n_pos_rows = write_mismatch_type_by_pos(
            mismatch_types_by_pos, args.mismatch_types_by_pos_out, args.prefix
        )
        sys.stderr.write(
            f"Wrote {n_pos_rows} per-position mismatch-type rows to "
            f"{args.mismatch_types_by_pos_out}.\n"
        )

    if args.mismatch_types_total_out:
        n_total_rows = write_mismatch_type_summary_total(
            mismatch_types_total, args.mismatch_types_total_out, args.prefix
        )
        sys.stderr.write(
            f"Wrote {n_total_rows} total mismatch-type rows to "
            f"{args.mismatch_types_total_out}.\n"
        )


if __name__ == "__main__":
    main()

