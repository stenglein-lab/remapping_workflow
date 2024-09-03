include { BUILD_BWA_INDEX                            } from '../../subworkflows/stenglein-lab/build_bwa_index'
include { MAP_TO_GENOME                              } from '../../subworkflows/stenglein-lab/map_to_genome'
include { BAM_TO_SAM                                 } from '../../modules/stenglein-lab/bam_to_sam'
                                                                                
workflow REMAP_TO_GENOMES {
                                                                                
  take:
  input     // [meta, fastq, bwa_index]
 
  main:                                                                         

  // make input channels in typical nf-core format of [meta, files]
  reads = input.map{ meta, fastq, bwa -> [meta, fastq] }
  fasta = input.map{ meta, fastq, bwa -> [meta, bwa]   }

  MAP_TO_GENOME(reads, fasta)

  BAM_TO_SAM(MAP_TO_GENOME.out.bam.filter{it[1].size() > 0})                                             
}         
