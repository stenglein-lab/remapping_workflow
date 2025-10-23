include { PARSE_MAPPING_SAMPLESHEET       } from '../../subworkflows/stenglein-lab/parse_mapping_samplesheet'
include { MARSHALL_FASTQ                  } from '../../subworkflows/stenglein-lab/marshall_fastq'
include { BUILD_BWA_INDEX                 } from '../../subworkflows/stenglein-lab/build_bwa_index'
include { REMAP_TO_GENOMES                } from '../../subworkflows/stenglein-lab/remap_to_genomes'
include { MAPPING_STATS                   } from '../../subworkflows/stenglein-lab/mapping_stats'

include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_COVERAGE    } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_STATS       } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_DEPTH       } from '../../modules/stenglein-lab/save_output_file'


workflow REMAPPING_WORKFLOW {                                                    

  main: 

  MARSHALL_FASTQ(params.fastq_dir, params.fastq_pattern)

  PARSE_MAPPING_SAMPLESHEET(params.mapping_samplesheet)

  // pull out necessary info for building BWA Index 
  PARSE_MAPPING_SAMPLESHEET.out.sample_sheet.map{ meta, fasta -> [meta, fasta] }.set{fasta_ch}

  BUILD_BWA_INDEX(fasta_ch)

  // pull out just sample ID (ignoring single-end vs not) 
  // so we can join just based on sample ID 
  MARSHALL_FASTQ.out.reads.map{meta, reads -> [meta.id, meta, reads]}.set{reads_ch}
  BUILD_BWA_INDEX.out.index.map{meta, index -> [meta.id, index]}.set{index_ch}

  // drop just sample ID, bring back in original meta
  reads_ch.join(index_ch).map{id, meta, reads, index -> [meta, reads, index] }.set{ch_mapping}
  
  REMAP_TO_GENOMES(ch_mapping)


  // tabulate mapping stats
  def per_base_coverage = !params.skip_per_base_coverage

  // join bam with fasta using sample ID 
  REMAP_TO_GENOMES.out.bam.map{ meta, bam -> [meta.id, meta, bam]}.set{bam_id_ch}
  fasta_ch.map{ meta, fasta -> [meta.id, fasta]}.set{fasta_id_ch}
  bam_id_ch.join(fasta_id_ch).map{id, meta, bam, fasta -> [meta, bam, fasta] }.set{bam_fasta_ch}

  MAPPING_STATS(bam_fasta_ch, per_base_coverage)

  MAPPING_STATS.out.prepended_depth.view { it -> "prepended depth: $it" }

  // save consolidated output files
  SAVE_COLLECTED_COVERAGE(MAPPING_STATS.out.prepended_coverage.collectFile(name: "collected_per_refseq_coverage.tsv"){it[1]})
  SAVE_COLLECTED_STATS   (MAPPING_STATS.out.prepended_stats.collectFile(name: "collected_stats.tsv"){it[1]})
  SAVE_COLLECTED_DEPTH   (MAPPING_STATS.out.prepended_depth.collectFile(name: "collected_per_base_depth.tsv"){it[1]})

  emit:

  bam = REMAP_TO_GENOMES.out.bam
  sam = REMAP_TO_GENOMES.out.sam

}

