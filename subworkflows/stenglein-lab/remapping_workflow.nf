include { PARSE_MAPPING_SAMPLESHEET       } from '../../subworkflows/stenglein-lab/parse_mapping_samplesheet'
include { MARSHALL_FASTQ                  } from '../../subworkflows/stenglein-lab/marshall_fastq'
include { BUILD_BWA_INDEX                 } from '../../subworkflows/stenglein-lab/build_bwa_index'
include { REMAP_TO_GENOMES                } from '../../subworkflows/stenglein-lab/remap_to_genomes'

workflow REMAPPING_WORKFLOW {                                                    

  MARSHALL_FASTQ(params.fastq_dir, params.fastq_pattern)

  PARSE_MAPPING_SAMPLESHEET(params.mapping_samplesheet)

  // pull out necessary info for building BWA Index 
  PARSE_MAPPING_SAMPLESHEET.out.sample_sheet.map{ meta, fasta -> [meta, fasta] }.set{fasta_ch}

  BUILD_BWA_INDEX(fasta_ch)


/*
  // what's going on with reads channel?
  MARSHALL_FASTQ.out.reads
    .tap { reads_view_ch }

  reads_view_ch.view { "Reads: $it" }

  // what's going on with BWA build channel?
  BUILD_BWA_INDEX.out.index
    .tap { bwa_view_ch }

  bwa_view_ch.view { "BWA Build: $it" }
 */ 

  // pull out just sample ID (ignoring single-end vs not) 
  // so we can join just based on sample ID 
  MARSHALL_FASTQ.out.reads.map{meta, reads -> [meta.id, meta, reads]}.set{reads_ch}
  BUILD_BWA_INDEX.out.index.map{meta, index -> [meta.id, index]}.set{index_ch}

  // drop just sample ID, bring back in original meta
  reads_ch.join(index_ch).map{id, meta, reads, index -> [meta, reads, index] }.set{ch_mapping}
  
/*
  // combine matching fastq and fasta from previous outputs into a single channel
  ch_mapping = MARSHALL_FASTQ.out.reads
    .join(BUILD_BWA_INDEX.out.index, failOnMismatch: false)
*/

  // debug join: anything getting to ch_mapping?
  ch_mapping
    .tap { mapping_view_ch }

  mapping_view_ch.view { "Mapping View ch: $it" }
  
  REMAP_TO_GENOMES(ch_mapping)
}

