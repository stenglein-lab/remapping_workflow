process VIRAL_CONSENSUS {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'https://depot.galaxyproject.org/singularity/viral_consensus:1.0.1--hcf1f8c1_0':
       'biocontainers/viral_consensus:1.0.1--hcf1f8c1_0' }"

    input:
    tuple val(meta), path(bam), val (refseq), path(ref_fasta)
    val(min_qual)
    val(min_depth)
    val(min_freq)

    output:
    tuple val(meta), path(ref_fasta), path("*.consensus.fasta") , emit: refseq_and_new
    tuple val(meta), path("*.consensus.fasta")                  , emit: fasta
    tuple val(meta), path("*.position_counts.txt")           , emit: position_counts
    tuple val(meta), path(ref_fasta)                         , emit: refseq
    path "versions.yml"                                      , emit: versions
  
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_base = ref_fasta.baseName
    """
    viral_consensus \
        -i $bam \
        -r $ref_fasta \
        -o ${ref_base}.consensus.not_renamed.fasta \
        -op ${ref_base}.position_counts.txt \
        -q $min_qual \
        -d $min_depth \
        -f $min_freq \
        $args 

    # rename fasta sequences: {sample_id}_{refseq_id}
    sed "s/^>.*/>${meta.id}_${refseq.id}/"  ${ref_base}.consensus.not_renamed.fasta > ${ref_base}.consensus.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(echo \$(viral_consensus --version 2>&1) | sed 's/^viral_consensus /viral_consensus: /')
    END_VERSIONS
    """
}
