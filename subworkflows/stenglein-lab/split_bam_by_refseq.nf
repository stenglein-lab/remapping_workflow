/*
 * This process splits a bam file created by mapping to multiple reference sequences 
 * into one bam file per refseq
 */

workflow SPLIT_BAM_BY_REFSEQ {

 take:
  bam_fasta          // [meta, bam, fasta]

 main:

  // define some empty channels for keeping track of stuff
  ch_versions     = Channel.empty()

  // split up input channel into separate bam and fasta channels

  // this will create a channel of [meta, bam, fasta_sequence_id]
  split_fasta_ch = bam_fasta.map{meta, bam, fasta -> [meta, fasta]}
    .splitFasta( record: [id: true, seqString: true] )

  // combine creates the cartesian cross product of input channels
  // creates a channel of: [meta, bam, seq_id]
  split_bam_fasta_ch = bam_fasta.map{meta, bam, fasta -> [meta, bam]}
    .combine(split_fasta_ch, by: 0)

  // this splits the bam into one bam per refseq
  SPLIT_BAM_BY_ONE_REFSEQ(split_bam_fasta_ch)

  // count number of reads mapped to each split bam
  COUNT_MAPPING_READS(SPLIT_BAM_BY_ONE_REFSEQ.out.per_refseq_bam)

  // only proceed when split bam contain enough mapped reads, 
  // with cutoff defined by params.mapped_read_cutoff:
  ch_enough_mapped_reads = COUNT_MAPPING_READS.out.per_refseq_bam
   .filter{Integer.parseInt(it[3]) > params.mapped_read_cutoff}
   .map {meta, bam, refseq, num_mapped_reads -> [meta, bam, refseq] }

  // a little process to write the refseq in fasta format
  OUTPUT_REFSEQ_FASTA(ch_enough_mapped_reads)

 emit:

  per_refseq_bam               = ch_enough_mapped_reads
  per_refseq_bam_fasta         = OUTPUT_REFSEQ_FASTA.out.per_refseq_bam_fasta
  per_refseq_bam_fasta_fai     = OUTPUT_REFSEQ_FASTA.out.per_refseq_bam_fasta_fai
  per_refseq_bam_fasta_fai_bai = OUTPUT_REFSEQ_FASTA.out.per_refseq_bam_fasta_fai_bai

}

// outputs a fasta file for the refseq
process OUTPUT_REFSEQ_FASTA {
   tag "$meta.id"
   label 'process_low'

   conda "bioconda::samtools=1.16.1"
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
       'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

   input:
   tuple val(meta), path(bam), val(refseq)                

   output:
   // TODO: seems like switching from tuples to records (new in nextflow 26.04) would help here:
   tuple val(meta), path(bam), val(refseq), path ("*.fasta"),                                emit: per_refseq_bam_fasta, optional: true
   tuple val(meta), path(bam), val(refseq), path ("*.fasta"), path ("*.fai"),                emit: per_refseq_bam_fasta_fai, optional: true
   tuple val(meta), path(bam), val(refseq), path ("*.fasta"), path ("*.fai"), path("*.bai"), emit: per_refseq_bam_fasta_fai_bai, optional: true

   script:
   """
   # output a fasta file with this refseq
   printf ">%s\n%s\n" ${refseq.id}  ${refseq.seqString} > ${meta.id}.${refseq.id}.fasta

   # create a fai index of fasta file in case needed downstream
   samtools faidx ${meta.id}.${refseq.id}.fasta

   # index bam file in case needed downstream
   samtools index ${bam}
   """
}

process COUNT_MAPPING_READS {
   tag "$meta.id"
   label 'process_low'

   conda "bioconda::samtools=1.16.1"
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
       'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

   input:
   tuple val(meta), path(bam), val(refseq)

   output:
   tuple val(meta), path(bam, includeInputs: true), val(refseq), env(num_mapping_reads), emit: per_refseq_bam,       optional: true

   shell:
   """
   num_mapping_reads=\$(samtools stats ${bam} | grep "^SN" | grep "reads mapped:" | cut -f 3)
   """
}

process SPLIT_BAM_BY_ONE_REFSEQ {
   tag "$meta.id"
   label 'process_low'

   conda "bioconda::samtools=1.16.1"
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
       'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

   input:
   tuple val(meta), path(bam), val(refseq)

   output:
   tuple val(meta), path("*.bam", includeInputs: false), val(refseq), emit: per_refseq_bam, optional: true
   path  "versions.yml",                                              emit: versions

   when:
   task.ext.when == null || task.ext.when

   script:
   def new_bam_name = bam.name.replaceAll(/.bam$/, ".${refseq.id}.bam")
   """
   # first have to index bam
   samtools index $bam

   # pull out refseq of interest (a region in samtools parlance)
   samtools \
       view \
       -h \
       --threads ${task.cpus-1} \
       $bam \
       ${refseq.id} | 

   # the piped awk and second samtools are to only keep a single @SQ line in the bam header
   # because viral_consensus complains otherwise
   awk '/^[^@]/ || (/^@/ && !/^@SQ/) || /^@SQ\tSN:${refseq.id}[[:space:]]/' |
   samtools view \
       -o ${new_bam_name} 

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
   END_VERSIONS
   """
}

