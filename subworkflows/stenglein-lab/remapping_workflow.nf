include { PARSE_MAPPING_SAMPLESHEET   } from '../../subworkflows/stenglein-lab/parse_mapping_samplesheet'
include { MARSHAL_FASTQ               } from '../../subworkflows/stenglein-lab/marshal_fastq'
include { BOWTIE2_BUILD_ALIGN         } from '../../subworkflows/stenglein-lab/bowtie2_build_align'
include { SPLIT_BAM_BY_REFSEQ         } from '../../subworkflows/stenglein-lab/split_bam_by_refseq'
include { MAPPING_STATS               } from '../../subworkflows/stenglein-lab/mapping_stats'
include { EXTRACT_INSERT_SIZES        } from '../../modules/stenglein-lab/extract_insert_sizes'
include { QUANTIFY_MISMATCHES         } from '../../modules/stenglein-lab/quantify_mismatches'
include { QUANTIFY_STRAND_BIAS        } from '../../subworkflows/stenglein-lab/quantify_strand_bias'
include { PROCESS_WORKFLOW_OUTPUT     } from '../../subworkflows/stenglein-lab/process_workflow_output'
include { GENERATE_CONSENSUS_SEQUENCE } from '../../subworkflows/stenglein-lab/generate_consensus_sequence'

// these save consolidated tidy output files to results directory
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_COVERAGE     } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_STATS        } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_DEPTH        } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_INSERT_SIZES } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_MISMATCHES   } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_STRAND_BIAS  } from '../../modules/stenglein-lab/save_output_file'

workflow REMAPPING_WORKFLOW {

 main:

  MARSHAL_FASTQ(params.fastq_dir, params.fastq_pattern)

  PARSE_MAPPING_SAMPLESHEET(params.mapping_samplesheet)

  // pull out just sample ID (ignoring single-end vs not) 
  // so we can join just based on sample ID 
  reads_ch       = MARSHAL_FASTQ.out.reads.map{meta, reads -> [meta.id, meta, reads]}
  samplesheet_ch = PARSE_MAPPING_SAMPLESHEET.out.sample_sheet.map{meta, fasta -> [meta.id, fasta]}

  // drop just sample ID, bring back in original meta
  mapping_ch = reads_ch.join(samplesheet_ch).map{id, meta, reads, fasta -> [meta, reads, fasta] }

  // run bowtie2 and align
  def save_unaligned = false
  def sort_bam = true
  BOWTIE2_BUILD_ALIGN (mapping_ch, save_unaligned, sort_bam)

  // split up bams by mapped-to refseq if necessary
  ch_split_bam       = Channel.empty()
  ch_split_bam_fasta = Channel.empty()
  if (params.tabulate_insert_sizes || params.generate_consensus_sequences || params.quantify_mismatches){

    // split up bam files into per-ref-seq bam files
    SPLIT_BAM_BY_REFSEQ(BOWTIE2_BUILD_ALIGN.out.bam_fasta)

    // assign output channels
    ch_split_bam           = SPLIT_BAM_BY_REFSEQ.out.per_refseq_bam
    ch_split_bam_fasta     = SPLIT_BAM_BY_REFSEQ.out.per_refseq_bam_fasta
    ch_split_bam_fasta_fai = SPLIT_BAM_BY_REFSEQ.out.per_refseq_bam_fasta_fai
  }

  // optionally extract insert sizes from mapped reads
  ch_insert_sizes = Channel.empty() 
  if (params.tabulate_insert_sizes) {

    // extract insert sizes from bam using samtools stats
    //
    // use split up bam files into per-ref-seq bam files
    // samtools stats quantifies insert sizes but not per refseq

    EXTRACT_INSERT_SIZES(ch_split_bam)

    ch_insert_sizes = ch_insert_sizes.mix(EXTRACT_INSERT_SIZES.out.insert_sizes)
  }

  // optionally quantify mismatches to reference sequences in mapped reads 
  ch_misincorporation = Channel.empty() 
  if (params.quantify_mismatches) {

    // extract mismatched bases from bam 
    // use per-ref-seq bam files

    QUANTIFY_MISMATCHES(ch_split_bam_fasta_fai)

    ch_misincorporation = ch_misincorporation.mix(QUANTIFY_MISMATCHES.out.misincorporation)
  }

  // optionally generate new consensus sequences
  ch_consensus_seqs = Channel.empty() 
  if (params.generate_consensus_sequences) {

     // params related to support needed for consensus calling
     ch_min_qual  = Channel.value(params.illumina_min_qual)
     ch_min_depth = Channel.value(params.illumina_min_depth)
     ch_min_freq  = Channel.value(params.illumina_min_freq)

    // create new consensus sequences 
    GENERATE_CONSENSUS_SEQUENCE(ch_split_bam_fasta, ch_min_depth, ch_min_qual, ch_min_freq)
    ch_consensus_seqs = GENERATE_CONSENSUS_SEQUENCE.out.merged_fasta

    ch_insert_sizes = ch_insert_sizes.mix(EXTRACT_INSERT_SIZES.out.insert_sizes)
  }

  // run optional workflow to quantify strand bias
  ch_strand_bias = Channel.empty()
  if (params.quantify_strand_bias) {
    QUANTIFY_STRAND_BIAS (BOWTIE2_BUILD_ALIGN.out.bam, params.R1_antisense_orientation)
    ch_strand_bias = ch_strand_bias.mix(QUANTIFY_STRAND_BIAS.out.strand_bias)
  }

  // tabulate mapping stats: samtools stats, coverage, and optionally per-base depth
  def per_base_coverage = !params.skip_per_base_coverage
  MAPPING_STATS(BOWTIE2_BUILD_ALIGN.out.bam_fasta, per_base_coverage)

  // save consolidated output files
  SAVE_COLLECTED_COVERAGE    (MAPPING_STATS.out.prepended_coverage.collectFile(name: "collected_per_refseq_coverage.tsv"){it[1]})
  SAVE_COLLECTED_STATS       (MAPPING_STATS.out.prepended_stats.collectFile(name: "collected_stats.tsv"){it[1]})
  SAVE_COLLECTED_DEPTH       (MAPPING_STATS.out.prepended_depth.collectFile(name: "collected_per_base_depth.tsv"){it[1]})
  SAVE_COLLECTED_INSERT_SIZES(ch_insert_sizes.collectFile(name: "collected_insert_sizes.txt"){it[1]})
  SAVE_COLLECTED_MISMATCHES  (ch_misincorporation.collectFile(name: "collected_mismatches.txt"){it[1]})
  SAVE_COLLECTED_STRAND_BIAS (ch_strand_bias.collectFile(name: "collected_strand_bias.txt"){it[1]})

  ch_coverage     = SAVE_COLLECTED_COVERAGE.out.file
  ch_stats        = SAVE_COLLECTED_STATS.out.file
  ch_depth        = SAVE_COLLECTED_DEPTH.out.file
  ch_insert_sizes = SAVE_COLLECTED_INSERT_SIZES.out.file
  ch_mismatches   = SAVE_COLLECTED_MISMATCHES.out.file
  ch_strand_bias  = SAVE_COLLECTED_STRAND_BIAS.out.file

  // optional workflow to further process/analyze the main output files
  if (params.process_workflow_output) {
    PROCESS_WORKFLOW_OUTPUT(ch_coverage, ch_stats, ch_depth, ch_insert_sizes, ch_mismatches, ch_strand_bias)
  }

 emit:

  bam              = BOWTIE2_BUILD_ALIGN.out.bam
  bowtie2_log      = BOWTIE2_BUILD_ALIGN.out.log
  coverage         = ch_coverage
  coverage_plots   = PROCESS_WORKFLOW_OUTPUT.out.coverage_plots
  mismatch_plots   = PROCESS_WORKFLOW_OUTPUT.out.mismatch_plots 
  stats            = ch_stats
  depth            = ch_depth
  insert_sizes     = ch_insert_sizes
  mismatches       = ch_mismatches 
  strand_bias      = ch_strand_bias
  consensus        = ch_consensus_seqs
  samples          = mapping_ch

}

