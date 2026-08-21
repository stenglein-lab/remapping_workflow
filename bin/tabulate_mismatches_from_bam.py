#!/usr/bin/env python3

"""
Quantify reference -> observed base substitutions in mapped reads.

Outputs are plain-text TSV files. Four independent output files can be
requested:

    --global-out
        Global counts of all 16 reference -> observed combinations.

    --refseq-out
        Frequencies of all 16 combinations for each reference sequence.

    --refpos-out
        Frequencies of all 16 combinations at each reference sequence
        and reference position.

    --readpos-out
        Frequencies of all 16 combinations as a function of read
        position, separately for R1/R2 and separately from the beginning
        and end of the read.

Read positions in readpos output are represented independently as:

    position_from_start
        1 = first sequenced base
        2 = second sequenced base
        ...

    position_from_end
        1 = last sequenced base
        2 = second-to-last sequenced base
        ...

Thus every eligible base contributes to BOTH the start-position and
end-position distributions.

Optional overlap collapsing:

    --collapse-overlaps

When enabled, paired reads are grouped by query name using a temporary
name-sorted BAM. Bases covered by both R1 and R2 are treated as one
molecular observation.

For an overlapping base:

    - if R1 and R2 agree, retain one observation;
    - if they disagree and one has higher base quality, retain the
      higher-quality observation;
    - if they disagree and base qualities are equal, discard the base.

This behavior can be changed with --overlap-disagreement.

Example without overlap collapsing:

    python mismatch_spectrum.py \
        --bam sample.bam \
        --ref reference.fa \
        --global-out global.tsv \
        --refseq-out refseq.tsv \
        --refpos-out refpos.tsv \
        --readpos-out readpos.tsv

Example with overlap collapsing:

    python mismatch_spectrum.py \
        --bam sample.bam \
        --ref reference.fa \
        --global-out global.tsv \
        --refseq-out refseq.tsv \
        --refpos-out refpos.tsv \
        --readpos-out readpos.tsv \
        --collapse-overlaps

"""

import argparse
import collections
import os
import shutil
import subprocess
import sys
import tempfile

import pysam


BASES = ("A", "C", "G", "T")

SUBSTITUTIONS = tuple(
    f"{ref}>{obs}"
    for ref in BASES
    for obs in BASES
)


# ----------------------------------------------------------------------
# Argument parsing
# ----------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Quantify reference-to-observed base substitutions in "
            "mapped reads from a coordinate-sorted BAM."
        )
    )

    # Required input arguments

    parser.add_argument(
        "--bam",
        required=True,
        help="Input coordinate-sorted and indexed BAM.",
    )

    parser.add_argument(
        "--ref",
        required=True,
        help=(
            "Reference FASTA. Must have an accompanying .fai index."
        ),
    )

    # Output files

    parser.add_argument(
        "--global-out",
        default=None,
        help=(
            "Output TSV for global substitution counts. "
            "If omitted, this table is not produced."
        ),
    )

    parser.add_argument(
        "--refseq-out",
        default=None,
        help=(
            "Output TSV for per-reference-sequence substitution "
            "counts. If omitted, this table is not produced."
        ),
    )

    parser.add_argument(
        "--refpos-out",
        default=None,
        help=(
            "Output TSV for per-reference-sequence/per-position "
            "substitution counts. If omitted, this table is "
            "not produced."
        ),
    )

    parser.add_argument(
        "--readpos-out",
        default=None,
        help=(
            "Output TSV for substitution counts by read position. "
            "If omitted, this table is not produced."
        ),
    )

    parser.add_argument(
        "--prefix",
        default=None,
        help=(
            "Optional text written as the first tab-delimited column "
            "of every output file."
        ),
    )

    parser.add_argument(
        "--output-headers",
        action="store_true",
        help="Output header lines as first line of output files. Default: do not.",
    )

    parser.add_argument(
        "--output-all-positions",
        action="store_true",
        help=(
            "Output positions in per-refseq-per-positition output even if they have "
            "zero counts for the indicated substitution (no observations of that substitution)."
            "Default: do not."
        ),
    )

    # Filtering

    parser.add_argument(
        "--min-mapq",
        type=int,
        default=20,
        help="Minimum mapping quality. Default: 20.",
    )

    parser.add_argument(
        "--min-baseq",
        type=int,
        default=30,
        help="Minimum base quality. Default: 30.",
    )

    parser.add_argument(
        "--include-duplicates",
        action="store_true",
        help="Include reads marked as duplicates. Default: exclude.",
    )

    parser.add_argument(
        "--include-secondary",
        action="store_true",
        help="Include secondary alignments. Default: exclude.",
    )

    parser.add_argument(
        "--include-supplementary",
        action="store_true",
        help="Include supplementary alignments. Default: exclude.",
    )

    # Read-position restriction

    parser.add_argument(
        "--max-read-position",
        type=int,
        default=None,
        help=(
            "If specified, only observations within this many bases "
            "of either read end are included in readpos output. "
            "Other output tables are unaffected.  Default: no limit."
        ),
    )

    # Overlap handling

    parser.add_argument(
        "--collapse-overlaps",
        action="store_true",
        help=(
            "Collapse overlapping R1/R2 observations into a single "
            "molecular observation. Requires temporary name sorting."
        ),
    )

    parser.add_argument(
        "--overlap-disagreement",
        choices=("higher-quality", "discard"),
        default="higher-quality",
        help=(
            "How to handle an R1/R2 overlap where the observed bases "
            "disagree. 'higher-quality' retains the observation with "
            "higher base quality; equal-quality disagreements are "
            "discarded. 'discard' discards all disagreements. "
            "Default: higher-quality."
        ),
    )

    return parser.parse_args()


# ----------------------------------------------------------------------
# Input validation
# ----------------------------------------------------------------------

def check_inputs(args):

    if not os.path.isfile(args.bam):
        raise FileNotFoundError(
            f"BAM not found: {args.bam}"
        )

    bai1 = args.bam + ".bai"
    bai2 = os.path.splitext(args.bam)[0] + ".bai"

    if not os.path.isfile(bai1) and not os.path.isfile(bai2):
        raise FileNotFoundError(
            f"BAM index not found for {args.bam}. "
            f"Expected {bai1} or {bai2}."
        )

    if not os.path.isfile(args.ref):
        raise FileNotFoundError(
            f"Reference FASTA not found: {args.ref}"
        )

    if not os.path.isfile(args.ref + ".fai"):
        raise FileNotFoundError(
            f"Reference FASTA index not found: {args.ref}.fai"
        )

    if args.min_mapq < 0:
        raise ValueError("--min-mapq must be >= 0")

    if args.min_baseq < 0:
        raise ValueError("--min-baseq must be >= 0")

    if args.max_read_position is not None:
        if args.max_read_position < 1:
            raise ValueError("--max-read-position must be >= 1")


# ----------------------------------------------------------------------
# Counters
# ----------------------------------------------------------------------

def empty_counter():
    return collections.Counter(
        {sub: 0 for sub in SUBSTITUTIONS}
    )


# ----------------------------------------------------------------------
# Counters
# ----------------------------------------------------------------------

def empty_counter():

    return collections.Counter(
        {sub: 0 for sub in SUBSTITUTIONS}
    )


# ----------------------------------------------------------------------
# Substitution
# ----------------------------------------------------------------------

def make_substitution(ref_base, obs_base):

    ref_base = ref_base.upper()
    obs_base = obs_base.upper()

    if ref_base not in BASES:
        return None

    if obs_base not in BASES:
        return None

    return f"{ref_base}>{obs_base}"


# ----------------------------------------------------------------------
# Read-level filtering
# ----------------------------------------------------------------------

# does this read pass filters?
def passes_read_filters(read, args):

    if read.is_unmapped:
        return False

    if read.is_secondary and not args.include_secondary:
        return False

    if read.is_supplementary and not args.include_supplementary:
        return False

    if read.is_duplicate and not args.include_duplicates:
        return False

    if read.mapping_quality < args.min_mapq:
        return False

    if read.query_sequence is None:
        return False

    return True

# ----------------------------------------------------------------------
# Read label
# ----------------------------------------------------------------------

def get_read_label(read):

    if not read.is_paired:
        return "R1"

    if read.is_read1:
        return "R1"

    if read.is_read2:
        return "R2"

    return "R1"


# ----------------------------------------------------------------------
# Extract observations from one read
# ----------------------------------------------------------------------

def get_read_observations(
    read,
    reference_sequence,
    args,
):
    """
    Return:

        {
            reference_position_0based: observation
        }

    Each observation contains:

        refpos
        ref_base
        obs_base
        substitution
        baseq
        query_pos
        position_from_start
        position_from_end

    Insertions and deletions are excluded.
    """

    observations = {}

    query_sequence = read.query_sequence

    if query_sequence is None:
        return observations

    read_length = len(query_sequence)

    query_qualities = read.query_qualities

    for query_pos, refpos in read.get_aligned_pairs(
        matches_only=False
    ):

        # Insertion relative to reference.
        if refpos is None:
            continue

        # Deletion relative to reference.
        if query_pos is None:
            continue

        if query_pos < 0 or query_pos >= read_length:
            continue

        # Base quality.

        if query_qualities is None:
            baseq = 0
        else:
            baseq = query_qualities[query_pos]

        if baseq < args.min_baseq:
            continue

        # Reference base.

        if refpos < 0 or refpos >= len(reference_sequence):
            continue

        ref_base = reference_sequence[refpos].upper()

        obs_base = query_sequence[query_pos].upper()

        sub = make_substitution(
            ref_base,
            obs_base,
        )

        if sub is None:
            continue

        observations[refpos] = {
            "refpos": refpos,
            "ref_base": ref_base,
            "obs_base": obs_base,
            "substitution": sub,
            "baseq": baseq,
            "query_pos": query_pos,

            # 1-based distance from the beginning.
            "position_from_start": query_pos + 1,

            # 1-based distance from the end.
            "position_from_end": (
                read_length - query_pos
            ),
        }

    return observations


# ----------------------------------------------------------------------
# Add one observation to all relevant counters
# ----------------------------------------------------------------------

def add_observation(
    read,
    observation,
    args,
    global_counts,
    refseq_counts,
    refpos_counts,
    readpos_counts,
):

    sub = observation["substitution"]

    refseq = read.reference_name

    refpos = observation["refpos"]

    # Global.

    global_counts[sub] += 1

    # Reference sequence.

    refseq_counts[refseq][sub] += 1

    # Reference position.
    #
    # Output positions are 1-based.

    refpos_counts[
        (refseq, refpos + 1)
    ][sub] += 1

    # Read-position distributions.
    #
    # These are independent distributions:
    #
    #   start 1, 2, 3, ...
    #
    # and:
    #
    #   end 1, 2, 3, ...
    #
    # The same base therefore contributes to both distributions.
    #
    # encode distances from end as negative integers and distances from 
    # start as positive integers

    label = get_read_label(read)

    start_position = observation[
        "position_from_start"
    ]

    end_position = observation[
        "position_from_end"
    ]

    if (
        args.max_read_position is None
        or start_position <= args.max_read_position
    ):

        key = (
            label,
            "start",
            start_position,
        )

        readpos_counts[key][sub] += 1

    if (
        args.max_read_position is None
        or end_position <= args.max_read_position
    ):

        key = (
            label,
            "end",
            -end_position,
        )

        readpos_counts[key][sub] += 1


# ----------------------------------------------------------------------
# Collapse one R1/R2 pair
# ----------------------------------------------------------------------

def collapse_pair(
    read1,
    read2,
    reference_sequence,
    args,
):
    """
    Collapse observations from one R1/R2 pair.

    For positions covered by only one read:
        retain that observation.

    For positions covered by both:

        same observed base:
            retain one observation, using the higher-BQ read's
            positional information.

        different observed bases:
            higher-BQ observation wins.

        equal BQ:
            discard.

    Returns:

        list of (chosen_read, observation)
    """

    obs1 = get_read_observations(
        read1,
        reference_sequence,
        args,
    )

    obs2 = get_read_observations(
        read2,
        reference_sequence,
        args,
    )

    result = []

    positions = sorted(
        set(obs1) | set(obs2)
    )

    for refpos in positions:

        a = obs1.get(refpos)
        b = obs2.get(refpos)

        # Only R1 covers the position.

        if a is not None and b is None:

            result.append(
                (read1, a)
            )

            continue

        # Only R2 covers the position.

        if b is not None and a is None:

            result.append(
                (read2, b)
            )

            continue

        # Both cover it.

        if a["ref_base"] != b["ref_base"]:
            # Should not happen with the same reference.
            continue

        # Same observed base.

        if a["obs_base"] == b["obs_base"]:

            if b["baseq"] > a["baseq"]:
                result.append(
                    (read2, b)
                )
            else:
                result.append(
                    (read1, a)
                )

            continue

        # Different observed bases.

        if args.overlap_disagreement == "discard":
            continue

        if a["baseq"] > b["baseq"]:

            result.append(
                (read1, a)
            )

        elif b["baseq"] > a["baseq"]:

            result.append(
                (read2, b)
            )

        else:
            # Equal-quality disagreement.
            continue

    return result


# ----------------------------------------------------------------------
# Standard processing
# ----------------------------------------------------------------------

def process_standard(
    bam,
    fasta,
    args,
    global_counts,
    refseq_counts,
    refpos_counts,
    readpos_counts,
):

    reads_examined = 0
    reads_passed = 0
    bases_counted = 0

    # Process reference sequences one at a time so that we only need
    # one reference sequence in memory at once.

    for refseq in bam.references:

        sys.stderr.write(
            f"Processing {refseq}...\n"
        )
        sys.stderr.flush()

        reference_sequence = fasta.fetch(
            refseq
        ).upper()

        for read in bam.fetch(refseq):

            reads_examined += 1

            if not passes_read_filters(
                read,
                args,
            ):
                continue

            reads_passed += 1

            observations = get_read_observations(
                read,
                reference_sequence,
                args,
            )

            for observation in observations.values():

                add_observation(
                    read,
                    observation,
                    args,
                    global_counts,
                    refseq_counts,
                    refpos_counts,
                    readpos_counts,
                )

                bases_counted += 1

    return (
        reads_examined,
        reads_passed,
        bases_counted,
    )


# ----------------------------------------------------------------------
# Overlap-collapsed processing
# ----------------------------------------------------------------------

def process_collapsed(
    bam,
    fasta,
    args,
    global_counts,
    refseq_counts,
    refpos_counts,
    readpos_counts,
):

    reads_examined = 0
    reads_passed = 0
    bases_counted = 0

    current_qname = None
    current_reads = []

    # Reference sequence cache.
    #
    # Because the BAM is query-name sorted, reads are not grouped by
    # reference sequence. We therefore cache reference sequences on
    # demand. For a very large genome this can consume substantial RAM.
    #
    # In practice, for very large references, the non-collapsed mode
    # is substantially more memory efficient.

    reference_cache = {}

    def get_reference(refseq):

        if refseq not in reference_cache:

            reference_cache[refseq] = fasta.fetch(
                refseq
            ).upper()

        return reference_cache[refseq]

    def process_group(reads):

        nonlocal reads_passed
        nonlocal bases_counted

        if not reads:
            return

        passed = [
            read
            for read in reads
            if passes_read_filters(
                read,
                args,
            )
        ]

        reads_passed += len(passed)

        if not passed:
            return

        r1 = None
        r2 = None
        other_reads = []

        for read in passed:

            if (
                read.is_paired
                and read.is_read1
                and r1 is None
            ):
                r1 = read
                continue

            if (
                read.is_paired
                and read.is_read2
                and r2 is None
            ):
                r2 = read
                continue

            other_reads.append(read)

        # Complete R1/R2 pair.

        if r1 is not None and r2 is not None:

            # The two mates should normally have the same reference
            # sequence for an ordinary overlapping pair. If they do
            # not, process them independently.

            if (
                r1.reference_name
                == r2.reference_name
            ):

                refseq = r1.reference_name

                reference_sequence = get_reference(
                    refseq
                )

                observations = collapse_pair(
                    r1,
                    r2,
                    reference_sequence,
                    args,
                )

                for chosen_read, observation in observations:

                    add_observation(
                        chosen_read,
                        observation,
                        args,
                        global_counts,
                        refseq_counts,
                        refpos_counts,
                        readpos_counts,
                    )

                    bases_counted += 1

            else:

                for read in (r1, r2):

                    refseq = read.reference_name

                    reference_sequence = get_reference(
                        refseq
                    )

                    observations = get_read_observations(
                        read,
                        reference_sequence,
                        args,
                    )

                    for observation in observations.values():

                        add_observation(
                            read,
                            observation,
                            args,
                            global_counts,
                            refseq_counts,
                            refpos_counts,
                            readpos_counts,
                        )

                        bases_counted += 1

        else:

            # No complete pair. Process available reads independently.

            for read in passed:

                refseq = read.reference_name

                reference_sequence = get_reference(
                    refseq
                )

                observations = get_read_observations(
                    read,
                    reference_sequence,
                    args,
                )

                for observation in observations.values():

                    add_observation(
                        read,
                        observation,
                        args,
                        global_counts,
                        refseq_counts,
                        refpos_counts,
                        readpos_counts,
                    )

                    bases_counted += 1

    for read in bam:

        reads_examined += 1

        qname = read.query_name

        if current_qname is None:
            current_qname = qname

        if qname != current_qname:

            process_group(
                current_reads
            )

            current_reads = []
            current_qname = qname

        current_reads.append(read)

    # Last group.

    process_group(
        current_reads
    )

    return (
        reads_examined,
        reads_passed,
        bases_counted,
    )


# ----------------------------------------------------------------------
# Temporary name sorting
# ----------------------------------------------------------------------

def make_name_sorted_bam(args):

    tmpdir = tempfile.mkdtemp(
        prefix="mismatch_spectrum_",
        dir=args.tmpdir,
    )

    output_bam = os.path.join(
        tmpdir,
        "name_sorted.bam",
    )

    sys.stderr.write(
        "Creating temporary query-name-sorted BAM...\n"
    )
    sys.stderr.flush()

    pysam.sort(
        "-n",
        "-o",
        output_bam,
        "-T",
        os.path.join(
            tmpdir,
            "sort_tmp",
        ),
        args.bam,
    )

    return output_bam, tmpdir


# ----------------------------------------------------------------------
# Prefix
# ----------------------------------------------------------------------

def prefix_columns(args, columns):

    if args.prefix is None:
        return columns

    return [
        args.prefix,
        *columns,
    ]


# ----------------------------------------------------------------------
# Output: global
# ----------------------------------------------------------------------

def write_global(
    path,
    counts,
    args,
):

    if path is None:
        return

    total = sum(
        counts.values()
    )

    with open(path, "w") as out:

        if args.output_headers:
           out.write(
               "\t".join(
                   prefix_columns(
                       args,
                       [
                           "ref_base",
                           "obs_base",
                           "count",
                       ],
                   )
               )
               + "\n"
           )

        for sub in SUBSTITUTIONS:

            ref_base, obs_base = sub.split(">")

            count = counts[sub]

            frequency = (
                count / total
                if total > 0
                else 0.0
            )

            out.write(
                "\t".join(
                    prefix_columns(
                        args,
                        [
                            ref_base,
                            obs_base,
                            str(count),
                        ],
                    )
                )
                + "\n"
            )


# ----------------------------------------------------------------------
# Output: reference sequence
# ----------------------------------------------------------------------

def write_refseq(
    path,
    counts_by_refseq,
    args,
):

    if path is None:
        return

    with open(path, "w") as out:

        if args.output_headers:
           out.write(
               "\t".join(
                   prefix_columns(
                       args,
                       [
                           "refseq",
                           "ref_base",
                           "obs_base",
                           "count",
                       ],
                   )
               )
               + "\n"
           )

        for refseq in sorted(
            counts_by_refseq
        ):

            counts = counts_by_refseq[
                refseq
            ]

            # Denominator = all observations for this reference
            # sequence.

            total = sum(
                counts.values()
            )

            for sub in SUBSTITUTIONS:

                ref_base, obs_base = sub.split(">")

                count = counts[sub]

                frequency = (
                    count / total
                    if total > 0
                    else 0.0
                )

                out.write(
                    "\t".join(
                        prefix_columns(
                            args,
                            [
                                refseq,
                                ref_base,
                                obs_base,
                                str(count),
                            ],
                        )
                    )
                    + "\n"
                )


# ----------------------------------------------------------------------
# Output: reference position
# ----------------------------------------------------------------------

def write_refpos(
    path,
    counts_by_refpos,
    args,
):

    if path is None:
        return

    with open(path, "w") as out:

        if args.output_headers:
           out.write(
               "\t".join(
                   prefix_columns(
                       args,
                       [
                           "refseq",
                           "position",
                           "ref_base",
                           "obs_base",
                           "count",
                       ],
                   )
               )
               + "\n"
           )

        for (
            refseq,
            position,
        ) in sorted(
            counts_by_refpos
        ):

            counts = counts_by_refpos[
                (refseq, position)
            ]

            total = sum(
                counts.values()
            )

            for sub in SUBSTITUTIONS:

                ref_base, obs_base = sub.split(">")

                count = counts[sub]

                frequency = (
                    count / total
                    if total > 0
                    else 0.0
                )

                # don't output zero count positions unless specified
                if count > 0 or args.output_all_positions:

                   out.write(
                       "\t".join(
                           prefix_columns(
                               args,
                               [
                                   refseq,
                                   str(position),
                                   ref_base,
                                   obs_base,
                                   str(count),
                               ],
                           )
                       )
                       + "\n"
                   )


# ----------------------------------------------------------------------
# Output: read position
# ----------------------------------------------------------------------

def write_readpos(
    path,
    counts_by_readpos,
    args,
):

    if path is None:
        return

    with open(path, "w") as out:

        if args.output_headers:
           out.write(
               "\t".join(
                   prefix_columns(
                       args,
                       [
                           "read",
                           # "end",
                           "position",
                           "ref_base",
                           "obs_base",
                           "count",
                       ],
                   )
               )
               + "\n"
           )

        for (
            read_label,
            end,
            position,
        ) in sorted(
            counts_by_readpos
        ):

            counts = counts_by_readpos[
                (
                    read_label,
                    end,
                    position,
                )
            ]

            total = sum(
                counts.values()
            )

            for sub in SUBSTITUTIONS:

                ref_base, obs_base = sub.split(">")

                count = counts[sub]

                frequency = (
                    count / total
                    if total > 0
                    else 0.0
                )

                # don't output zero count positions unless specified
                if count > 0 or args.output_all_positions: 
                   out.write(
                       "\t".join(
                           prefix_columns(
                               args,
                               [
                                   read_label,
                                   # end,
                                   str(position),
                                   ref_base,
                                   obs_base,
                                   str(count),
                               ],
                           )
                       )
                       + "\n"
                   )


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():

    args = parse_args()

    check_inputs(args)

    fasta = pysam.FastaFile(
        args.ref
    )

    global_counts = empty_counter()

    refseq_counts = collections.defaultdict(
        empty_counter
    )

    refpos_counts = collections.defaultdict(
        empty_counter
    )

    readpos_counts = collections.defaultdict(
        empty_counter
    )

    temporary_dir = None
    input_bam = args.bam

    try:

        # --------------------------------------------------------------
        # Prepare BAM
        # --------------------------------------------------------------

        if args.collapse_overlaps:

            input_bam, temporary_dir = (
                make_name_sorted_bam(args)
            )

        # --------------------------------------------------------------
        # Process
        # --------------------------------------------------------------

        sys.stderr.write(
            "Processing alignments...\n"
        )
        sys.stderr.flush()

        with pysam.AlignmentFile(
            input_bam,
            "rb",
        ) as bam:

            if args.collapse_overlaps:

                (
                    reads_examined,
                    reads_passed,
                    bases_counted,
                ) = process_collapsed(
                    bam,
                    fasta,
                    args,
                    global_counts,
                    refseq_counts,
                    refpos_counts,
                    readpos_counts,
                )

            else:

                (
                    reads_examined,
                    reads_passed,
                    bases_counted,
                ) = process_standard(
                    bam,
                    fasta,
                    args,
                    global_counts,
                    refseq_counts,
                    refpos_counts,
                    readpos_counts,
                )

        # --------------------------------------------------------------
        # Write outputs
        # --------------------------------------------------------------

        if args.global_out is not None:

            sys.stderr.write(
                f"Writing {args.global_out}\n"
            )

            write_global(
                args.global_out,
                global_counts,
                args,
            )

        if args.refseq_out is not None:

            sys.stderr.write(
                f"Writing {args.refseq_out}\n"
            )

            write_refseq(
                args.refseq_out,
                refseq_counts,
                args,
            )

        if args.refpos_out is not None:

            sys.stderr.write(
                f"Writing {args.refpos_out}\n"
            )

            write_refpos(
                args.refpos_out,
                refpos_counts,
                args,
            )

        if args.readpos_out is not None:

            sys.stderr.write(
                f"Writing {args.readpos_out}\n"
            )

            write_readpos(
                args.readpos_out,
                readpos_counts,
                args,
            )

        # --------------------------------------------------------------
        # Summary
        # --------------------------------------------------------------

        sys.stderr.write(
            "\nFinished.\n"
            f"  Reads examined:        "
            f"{reads_examined:,}\n"
            f"  Reads passing filters: "
            f"{reads_passed:,}\n"
            f"  Bases counted:         "
            f"{bases_counted:,}\n"
        )

    finally:

        fasta.close()

        if temporary_dir is not None:

            sys.stderr.write(
                "Removing temporary files...\n"
            )

            shutil.rmtree(
                temporary_dir,
                ignore_errors=True,
            )


if __name__ == "__main__":
    main()
