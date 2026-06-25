#!/usr/bin/env python3
"""
bam_to_strand_specific_coverage_depth.py:

Reads a BAM file an classifies each alignment into forward/reverse strand
based on SAM flags and accumulates per-base/per-refseq strand-specific 
coverage values

Handles both single-end and paired-end datasets

Strand classification is auto-detected per read using the PAIRED flag (0x1):
    - If a read is paired (part of a paired-end fragment), strand is
      assigned using the standard fr-firststrand paired convention:
          read2 + not-reverse  -> fwd
          read1 + reverse      -> fwd
          read2 + reverse      -> rev
          read1 + not-reverse  -> rev
    - If a read is NOT paired (single-end), strand is assigned directly
      from its own orientation:
          not-reverse -> fwd
          reverse     -> rev

This means a single BAM containing a mix of SE and PE reads (e.g. merged
runs) is handled correctly 

If your paired-end library is fr-secondstrand (opposite convention),
flip PE_FWD_IS_READ2_FORWARD below to False.

Notes:
    - Coverage counts only reference-consuming, non-deletion cigar blocks
      (matches/mismatches), same convention as `samtools depth` default.
    - Uses a difference-array + cumsum (numpy) per contig instead of
      incrementing every base in a Python loop -> O(reads) work instead
      of O(total aligned bases), much faster for high-depth data.
    - Overlapping mate pairs are NOT deduplicated (each read counted
      independently), matching samtools depth's default behavior.
    - Requires: pip install pysam numpy
"""

import argparse
import os
import sys
from array import array

# import numpy as np
import pysam

# SAM flags
FLAG_PAIRED    = 0x1
FLAG_REVERSE   = 0x10
FLAG_READ1     = 0x40
FLAG_READ2     = 0x80
FLAG_SECONDARY = 0x100
FLAG_QCFAIL    = 0x200
FLAG_DUP       = 0x400
FLAG_SUPPLEMENTARY = 0x800
FLAG_UNMAPPED  = 0x4

def classify_strand(flag):
    """
    Classify one read based on SAM flags as corresponding to the refseq
    "forward" strand or the refseq "reverse" strand.

    This assumes the typical strand-specific interpretation:
    that read1 is in the reverse complement orientation relative to the 
    strand of origin.

    Interpretation depends on orientation of reference sequence too.

    Auto-detect SE vs PE per-read via the PAIRED flag, then classify
    into 'fwd' or 'rev'. Returns None for reads that can't be classified
    (e.g. paired but missing both READ1/READ2 flags).
    """
    # is this read mapped in the reverse orientation relative to the ref seq?
    is_reverse = bool(flag & FLAG_REVERSE)

    # this is a paired end read
    if bool(flag & FLAG_PAIRED):
        is_read2 = bool(flag & FLAG_READ2)
        is_read1 = bool(flag & FLAG_READ1)
        # this is a paired read that is neither read1 nor read2 (TODO: should error?)
        if not (is_read1 or is_read2):
            return None
        
        # handle read1 / read2 differently (oppositely)
        if is_read2:
            return "fwd" if is_reverse else "rev"
        else:  # read1
            return "rev" if is_reverse else "fwd"
    else:
        # single end read
        return "rev" if is_reverse else "fwd"


def main():

    # handle CLI
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bam", help="input BAM (sorted; index not required)")

    p.add_argument("--min_mapq", type=int, default=20, help="minimum MAPQ to count a read")

    p.add_argument("--skip_secondary", action="store_true", default=True,
                    help="skip secondary alignments (default: on)")

    p.add_argument("--skip_supplementary", action="store_true", default=True,
                    help="skip supplementary alignments (default: on)")

    # p.add_argument("--fr-secondstrand", action="store_true", default=False,
                   # help="non-default strand mapping. See: https://ccb.jhu.edu/software/tophat/manual.shtml (default: false)")

    p.add_argument("--skip_dup", action="store_true", help="skip PCR/optical duplicates (default: off)")

    p.add_argument("--all_positions", action="store_true",
                    help="write zero-depth positions too (equivalent to samtools depth -a); "
                         "default only writes covered positions")

    p.add_argument("--out", default=None,
                    help="Path to output TSV file. If omitted, results are written to stdout.")

    p.add_argument("--prefix", default=None,
                    help="Optional text to output as a 1st column of tsv output.")

    args = p.parse_args()

    # read in bam file
    bam = pysam.AlignmentFile(args.bam, "rb")

    # dictionary of refseq lengths
    ref_lengths = dict(zip(bam.references, bam.lengths))

    # diff is a counter of increments (positive or negative) to running coverage along refseqs
    # diff[contig] = (fwd_diff_array, rev_diff_array); length+1 to allow safe end+1 decrement
    diff = {}

    n_reads = 0
    n_skipped = 0

    # for each read in the bam
    for read in bam.fetch(until_eof=True):
        n_reads += 1

        # skip unmapped reads
        if read.flag & FLAG_UNMAPPED:
            n_skipped += 1
            continue
        # skip secondary alignments 
        if args.skip_secondary and (read.flag & FLAG_SECONDARY):
            n_skipped += 1
            continue
        # skip supplementary alignments
        if args.skip_supplementary and (read.flag & FLAG_SUPPLEMENTARY):
            n_skipped += 1
            continue
        # skip reads marked as PCR/optical duplicates
        if args.skip_dup and (read.flag & FLAG_DUP):
            n_skipped += 1
            continue
        # skip QCFAIL flagged reads
        if read.flag & FLAG_QCFAIL:
            n_skipped += 1
            continue

        # optional filter based on mapping quality
        if read.mapping_quality < args.min_mapq:
            n_skipped += 1
            continue

        # classify as being mapped to forward or reverse strand
        strand = classify_strand(read.flag)
        if strand is None:
            n_skipped += 1
            continue

        # what refseq is this?
        contig = read.reference_name
        if contig not in diff:
            length = ref_lengths.get(contig)
            if length is None:
                raise ValueError(
                    f"Could not find reference length for contig '{contig}' "
                    f"(read: {read.query_name}). This likely means the BAM "
                    f"header is missing or inconsistent with its alignments. "
                    f"Stopping execution."
                )

            # if not already present, initialize a dictionary, keyed on contig name, 
            # where the values are tuples containing arrays of forward/reverse read mapping counts
            # +1 length buffer so end==length decrements don't go out of bounds
            #
            # this is a dictionary of positional increments to running coverage
            diff[contig] = (
                array('i', [0]) * (length + 1),
                array('i', [0]) * (length + 1),
            )

        # identify which refseq/strand to increment counts for this read
        fwd_diff, rev_diff = diff[contig]
        target = fwd_diff if strand == "fwd" else rev_diff

        # increment *incremental* coverage counts for this read
        # get_blocks() returns 0-based half-open ref-consuming blocks (cigar M/=/X),
        # matching samtools depth's definition of "covered" positions.
        # note that we are adding 1 to incremental coverage total at read beginning
        # and subtracting 1 incremental coverage total at read end; this avoids incrementing
        # coverage totals at all positions in the middle of reads
        for start, end in read.get_blocks():
            target[start] += 1
            target[end] -= 1

    bam.close()

    sys.stderr.write(f"Processed {n_reads} reads, skipped {n_skipped}.\n")

    # now convert running coverage increments to actual coverage totals
    # coverage[contig] = (fwd_depth_list, rev_depth_list), 0-based positions
    coverage = {}
    for contig, (fwd_diff, rev_diff) in diff.items():
        fwd_depth = [0] * (len(fwd_diff) - 1)
        rev_depth = [0] * (len(rev_diff) - 1)

        running = 0
        for i in range(len(fwd_diff) - 1):
            # this is the running total coverage
            # the running total will be changed by the incremental counter at each position
            running += fwd_diff[i]
            fwd_depth[i] = running

        running = 0
        for i in range(len(rev_diff) - 1):
            # this is the running total coverage
            # the running total will be changed by the incremental counter at each position
            running += rev_diff[i]
            rev_depth[i] = running

        coverage[contig] = (fwd_depth, rev_depth)

    for contig, (fwd_depth, rev_depth) in coverage.items():
        out_fh = open(args.out, "w") if args.out else sys.stdout
        for pos in range(len(fwd_depth)):
            f, r = fwd_depth[pos], rev_depth[pos]
            # if no coverage and not outputing 0-depth positions
            if f == 0 and r == 0 and not args.all_positions:
                continue

            # optional 1st column prefix
            if args.prefix is not None:
                out_fh.write(f"{args.prefix}\t")

            # output forward and reverse coverage
            out_fh.write(f"{contig}\t{pos + 1}\t{f}\t{r}\n")  # 1-based pos, like samtools depth

        # close file handle if open
        if args.out: 
            out_fh.close()


if __name__ == "__main__":
    main()
