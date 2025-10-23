include { PREPEND_TSV_WITH_ID                        } from '../../modules/stenglein-lab/prepend_tsv_with_id'
                                                                                
include { EXTRACT_STRAND_BIAS                        } from '../../modules/stenglein-lab/extract_strand_bias'
include { PROCESS_STRAND_BIAS_OUTPUT                 } from '../../modules/stenglein-lab/process_strand_bias'

include { SAVE_OUTPUT_FILE as SAVE_STRAND_BIAS_OUTPUT    } from '../../modules/stenglein-lab/save_output_file'
                                                                                
workflow QUANTIFY_STRAND_BIAS {
                                                                                
  take:
  input     // [meta, sam]
 
  main:                                                                         

  // make input channels in typical nf-core format of [meta, files]
  // reads = input.map{ meta, fastq, bwa, ori -> [meta, fastq] }
  // fasta = input.map{ meta, fastq, bwa, ori -> [meta, bwa]   }

  def r1_orientation = params.r1_orientation
  ori_ch = input.map{ meta, sam -> [meta, sam, r1_orientation]}

  EXTRACT_STRAND_BIAS(ori_ch)
  PREPEND_TSV_WITH_ID(EXTRACT_STRAND_BIAS.out.txt)                              
                                                                                
  // a directory with additional R scripts
  // R_script_dir_ch = Channel.fromPath(params.R_shared_script_dir)

  SAVE_STRAND_BIAS_OUTPUT(PREPEND_TSV_WITH_ID.out.tsv.collectFile(name: "collected_strand_bias.tsv"){it[1]})

  // PROCESS_STRAND_BIAS_OUTPUT(PREPEND_TSV_WITH_ID.out.tsv.collectFile(name: "collected_strand_bias.tsv"){it[1]}, COLLECT_METADATA.out.rds, virus_refseq_ch, R_lib_dir, R_script_dir_ch)

  // this will save main output file
  // SAVE_OUTPUT_FILE(PROCESS_STRAND_BIAS_OUTPUT.out.txt)
}         
