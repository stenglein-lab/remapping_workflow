include { PARSE_MAPPING_SAMPLESHEET   } from '../../subworkflows/stenglein-lab/parse_mapping_samplesheet'
include { MARSHAL_FASTQ               } from '../../subworkflows/stenglein-lab/marshal_fastq'
include { BOWTIE2_BUILD               } from '../../subworkflows/stenglein-lab/bowtie2_build_align'
include { BOWTIE2_ALIGN               } from '../../subworkflows/stenglein-lab/bowtie2_build_align'
include { SPLIT_BAM_BY_REFSEQ         } from '../../subworkflows/stenglein-lab/split_bam_by_refseq'
include { MAPPING_STATS               } from '../../subworkflows/stenglein-lab/mapping_stats'
include { EXTRACT_INSERT_SIZES        } from '../../modules/stenglein-lab/extract_insert_sizes'
include { QUANTIFY_MISMATCHES         } from '../../modules/stenglein-lab/quantify_mismatches'
include { DAMAGE_PROFILER             } from '../../modules/stenglein-lab/damage_profiler'
include { QUANTIFY_STRAND_BIAS        } from '../../subworkflows/stenglein-lab/quantify_strand_bias'
include { PROCESS_WORKFLOW_OUTPUT     } from '../../subworkflows/stenglein-lab/process_workflow_output'
include { GENERATE_CONSENSUS_SEQUENCE } from '../../subworkflows/stenglein-lab/generate_consensus_sequence'

// these save consolidated tidy output files to results directory
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_COVERAGE     } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_STATS        } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_DEPTH        } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_INSERT_SIZES } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_MISMATCHES   } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_MISMATCHES_PR} from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_MISMATCHES_BP} from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_MISMATCHES_T } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_STRAND_BIAS  } from '../../modules/stenglein-lab/save_output_file'


workflow REMAPPING_WORKFLOW {

 main:

  // marshal fastq 
  MARSHAL_FASTQ(params.fastq_dir, params.fastq_pattern)

  // parse sample sheet (ids -> refseq_fastas)
  PARSE_MAPPING_SAMPLESHEET(params.mapping_samplesheet)

  // wait for samplesheet parsing to finish and collect 
  samplesheet_list_ch = PARSE_MAPPING_SAMPLESHEET.out.sample_sheet
    .map { meta, fasta -> [meta.id, fasta] }
    .toList()                                    // emits a single list of [pat, fasta] pairs
    .map { it -> [it] }                          // wrap so combine sees it as one argument

  // merge fastq with refseq fasta, accounting for possible regular expressions in mapping sample sheet IDs
  mapping_ch = MARSHAL_FASTQ.out.reads
    .combine(samplesheet_list_ch)
    .map { meta, reads, pattern_rows ->

        def matches = pattern_rows.findAll { pat, fasta ->
            meta.id ==~ GlobUtils.globToRegex(pat)
        }

        if (matches.size() == 0) {
            log.warn "No mapping samplesheet entry matches sample id: ${meta.id} fastq: $reads — will not map."
            return null // nulls will be filtered out below
        }

        if (matches.size() > 1) {
            def matchedPats = matches.collect { pat, fasta -> "'${pat}'" }.join(', ')
            exit 1, "Sample '${meta.id}' matches more than one mapping samplesheet pattern. matched patterns: ${matchedPats}."
        }

        return [meta, reads, matches[0][1]]
    }
    .filter { it != null }

  // create one index per refseq fasta to avoid duplicate index building
  // keep track of the original fasta path (toString()) because path to fasta
  // can change once brought into work directories (e.g. in BOWTIE2_BUILD work dir)
  refseq_fasta_ch = mapping_ch
    .map { meta, reads, fasta -> [fasta.toString(), file(fasta)] }
    .unique { it[0] }          

  // build indexes, one per fasta
  BOWTIE2_BUILD(refseq_fasta_ch)        

  // merge in indexes to original mapping ch
  full_mapping_ch = mapping_ch
    .map { meta, reads, fasta -> [fasta.toString(), fasta, meta, reads] }
    .combine(BOWTIE2_BUILD.out.index.map{fasta_name, fasta, index_dir, index_base -> [fasta_name, index_dir, index_base]}, by: 0)
    .map { fasta_name, fasta, meta, reads, index_dir, index_base -> [meta, reads, fasta, index_dir, index_base] }

  // run bowtie2 and align
  def save_unaligned = false
  def bowtie2_options = Channel.value(params.bowtie2_options)
  BOWTIE2_ALIGN (full_mapping_ch, save_unaligned, bowtie2_options)

  // split up bams by mapped-to refseq if necessary
  ch_split_bam       = Channel.empty()
  ch_split_bam_fasta = Channel.empty()
  if (params.tabulate_insert_sizes || params.generate_consensus_sequences || params.quantify_mismatches){

    // split up bam files into per-ref-seq bam files
    SPLIT_BAM_BY_REFSEQ(BOWTIE2_ALIGN.out.bam_fasta)

    // assign output channels
    ch_split_bam               = SPLIT_BAM_BY_REFSEQ.out.per_refseq_bam
    ch_split_bam_fasta         = SPLIT_BAM_BY_REFSEQ.out.per_refseq_bam_fasta
    ch_split_bam_fasta_fai     = SPLIT_BAM_BY_REFSEQ.out.per_refseq_bam_fasta_fai
    ch_split_bam_fasta_fai_bai = SPLIT_BAM_BY_REFSEQ.out.per_refseq_bam_fasta_fai_bai
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
  ch_mismatches         = Channel.empty() 
  ch_mismatches_per_ref = Channel.empty() 
  ch_mismatches_by_pos  = Channel.empty() 
  ch_mismatches_total   = Channel.empty() 
  if (params.quantify_mismatches) {

    // extract mismatched bases from bam 
    // use per-ref-seq bam files
    // params related to support needed for consensus calling

    QUANTIFY_MISMATCHES(ch_split_bam_fasta_fai_bai, 
                        params.min_mismatch_depth, 
                        params.min_mismatch_base_quality, 
                        params.min_mismatch_mapping_quality)

    ch_mismatches         = ch_mismatches        .mix(QUANTIFY_MISMATCHES.out.per_refseq_per_position_mismatches)
    ch_mismatches_per_ref = ch_mismatches_per_ref.mix(QUANTIFY_MISMATCHES.out.per_refseq_mismatches)
    ch_mismatches_by_pos  = ch_mismatches_by_pos .mix(QUANTIFY_MISMATCHES.out.by_position_mismatches)
    ch_mismatches_total   = ch_mismatches_total  .mix(QUANTIFY_MISMATCHES.out.total_mismatches)
  }

  // optionally generate new consensus sequences
  ch_consensus_seqs = Channel.empty() 
  if (params.generate_consensus_sequences) {

     // params related to support needed for consensus calling
     ch_min_qual  = Channel.value(params.consensus_min_qual)
     ch_min_depth = Channel.value(params.consensus_min_depth)
     ch_min_freq  = Channel.value(params.consensus_min_freq)

    // create new consensus sequences 
    GENERATE_CONSENSUS_SEQUENCE(ch_split_bam_fasta, ch_min_depth, ch_min_qual, ch_min_freq)
    ch_consensus_seqs = GENERATE_CONSENSUS_SEQUENCE.out.merged_fasta

    ch_insert_sizes = ch_insert_sizes.mix(EXTRACT_INSERT_SIZES.out.insert_sizes)
  }

  // run optional workflow to quantify strand bias
  ch_strand_bias = Channel.empty()
  if (params.quantify_strand_bias) {
    QUANTIFY_STRAND_BIAS (BOWTIE2_ALIGN.out.bam, params.R1_antisense_orientation)
    ch_strand_bias = ch_strand_bias.mix(QUANTIFY_STRAND_BIAS.out.strand_bias)
  }

  // tabulate mapping stats: samtools stats, coverage, and optionally per-base depth
  def per_base_coverage = !params.skip_per_base_coverage
  MAPPING_STATS(BOWTIE2_ALIGN.out.bam_fasta, per_base_coverage)

  // save consolidated output files
  SAVE_COLLECTED_COVERAGE    (MAPPING_STATS.out.prepended_coverage.collectFile(name: "collected_per_refseq_coverage.tsv"){it[1]})
  SAVE_COLLECTED_STATS       (MAPPING_STATS.out.prepended_stats.collectFile(name: "collected_stats.tsv"){it[1]})
  SAVE_COLLECTED_DEPTH       (MAPPING_STATS.out.prepended_depth.collectFile(name: "collected_per_base_depth.tsv"){it[1]})
  SAVE_COLLECTED_INSERT_SIZES(ch_insert_sizes.collectFile(name: "collected_insert_sizes.txt"){it[1]})
  SAVE_COLLECTED_MISMATCHES    (ch_mismatches.collectFile(name: "collected_mismatches.txt"){it[1]})
  SAVE_COLLECTED_MISMATCHES_PR (ch_mismatches_per_ref.collectFile(name: "collected_per_refseq_mismatches.txt"){it[1]})
  SAVE_COLLECTED_MISMATCHES_BP (ch_mismatches_by_pos.collectFile(name: "collected_by_read_position_mismatches.txt"){it[1]})
  SAVE_COLLECTED_MISMATCHES_T  (ch_mismatches_total.collectFile(name: "collected_total_mismatches.txt"){it[1]})
  SAVE_COLLECTED_STRAND_BIAS (ch_strand_bias.collectFile(name: "collected_strand_bias.txt"){it[1]})

  ch_coverage     = SAVE_COLLECTED_COVERAGE.out.file
  ch_stats        = SAVE_COLLECTED_STATS.out.file
  ch_depth        = SAVE_COLLECTED_DEPTH.out.file
  ch_insert_sizes = SAVE_COLLECTED_INSERT_SIZES.out.file
  ch_mismatches   = SAVE_COLLECTED_MISMATCHES.out.file
  ch_mismatches_per_ref = SAVE_COLLECTED_MISMATCHES_PR.out.file
  ch_mismatches_by_pos  = SAVE_COLLECTED_MISMATCHES_BP.out.file
  ch_mismatches_total   = SAVE_COLLECTED_MISMATCHES_T.out.file
  ch_strand_bias  = SAVE_COLLECTED_STRAND_BIAS.out.file

  // optional workflow to further process/analyze the main output files
  ch_coverage_plots = Channel.empty()
  if (params.process_workflow_output) {
    PROCESS_WORKFLOW_OUTPUT(ch_coverage, ch_stats, ch_depth, ch_insert_sizes, ch_mismatches, ch_strand_bias)
    ch_coverage_plots = PROCESS_WORKFLOW_OUTPUT.out.coverage_plots
  }

 emit:

  bam              = BOWTIE2_ALIGN.out.bam
  bowtie2_log      = BOWTIE2_ALIGN.out.log
  coverage         = ch_coverage
  coverage_plots   = ch_coverage_plots 
  stats            = ch_stats
  depth            = ch_depth
  insert_sizes     = ch_insert_sizes
  mismatches       = ch_mismatches 
  mismatches_per_ref = ch_mismatches_per_ref
  mismatches_by_pos  = ch_mismatches_by_pos  
  mismatches_total   = ch_mismatches_total 
  strand_bias      = ch_strand_bias
  consensus        = ch_consensus_seqs
  samples          = mapping_ch

}

