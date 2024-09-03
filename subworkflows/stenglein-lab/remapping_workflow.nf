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
  
  // combine matching fastq and fasta from previous outputs into a single channel
  ch_mapping = MARSHALL_FASTQ.out.reads
    .join(BUILD_BWA_INDEX.out.index)
  
  REMAP_TO_GENOMES(ch_mapping)
}

