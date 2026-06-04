/*
 * This process uses bowtie2 to build an index from the provided refseq fasta
 * and then maps reads to that index.
 *
 * There are separated build and align processes below
 *
 * It is based on nf-core bowtie2 module code:
 * https://github.com/nf-core/modules/tree/master/modules/nf-core/bowtie2/align
 * https://github.com/nf-core/modules/tree/master/modules/nf-core/bowtie2/build
 * 
 */
process BOWTIE2_BUILD_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b41b403e81883126c3227fc45840015538e8e2212f13abc9ae84e4b98891d51c/data' :
        'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:edeb13799090a2a6' }"

    input:
    tuple val(meta), path(reads), path(fasta)
    val   save_unaligned
    val   sort_bam

    output:
    tuple val(meta), path("*.bam"), path(fasta)   , emit: bam_fasta     
    tuple val(meta), path("*.bam")                , emit: bam     
    tuple val(meta), path("*.log")                , emit: log
    tuple val(meta), path("*fastq.gz")            , emit: fastq   , optional:true
    tuple val(meta), path(fasta)                  , emit: fasta
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    def unaligned = ""
    def reads_args = ""
    if (meta.single_end) {
        unaligned = save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ""
        reads_args = "-U ${reads}"
    } else {
        unaligned = save_unaligned ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ""
        reads_args = "-1 ${reads[0]} -2 ${reads[1]}"
    }

    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension = "bam"

    """
    mkdir bowtie2
    bowtie2-build $fasta bowtie2/${fasta.baseName}

    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
    [ -z "\$INDEX" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

    bowtie2 \\
        -x \$INDEX \\
        --no-unal \\
        $reads_args \\
        --threads $task.cpus \\
        $unaligned \\
        $args \\
        2> >(tee ${prefix}.bowtie2.log >&2) \\
        | samtools $samtools_command $args2 --threads $task.cpus -o ${prefix}.bt.${fasta.baseName}.${extension} -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

/*
 * This process uses bowtie2 to build a mapping index from one fasta
 * 
 * It is based on nf-core bowtie2 module code:
 * https://github.com/nf-core/modules/tree/master/modules/nf-core/bowtie2/build
 * */
process BOWTIE2_BUILD {
    tag   "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b41b403e81883126c3227fc45840015538e8e2212f13abc9ae84e4b98891d51c/data' :
        'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:edeb13799090a2a6' }"

    input:
    tuple val(fasta_name), path(fasta)

    output:
    tuple val(fasta_name), path(fasta), path("index/"), val("index/${fasta.baseName}") , emit: index
    path  "versions.yml"                                                               , emit: versions

    script:
    """
    mkdir index
    bowtie2-build \\
      --threads ${task.cpus} \\
      $fasta \\
      index/${fasta.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2-build: \$(echo \$(bowtie2-build --version 2>&1) | sed 's/^.*bowtie2-build-s version //; s/ .*\$//')
    END_VERSIONS
    """
}

/*
 * This process uses bowtie2 to map reads to a set of reference sequences
 * 
 * It is based on nf-core bowtie2 module code:
 * https://github.com/nf-core/modules/tree/master/modules/nf-core/bowtie2/align
 * */
process BOWTIE2_ALIGN {
    tag   "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b41b403e81883126c3227fc45840015538e8e2212f13abc9ae84e4b98891d51c/data' :
        'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:edeb13799090a2a6' }"

    input:
    tuple val(meta), path(reads), path(fasta), path(index_dir), val(index_base)
    val   save_unaligned
    val   bowtie2_options  // any additional options to bowtie2, e.g. --local or --end_to_end

    output:
    tuple val(meta), path("*.bam"), path(fasta)   , emit: bam_fasta     
    tuple val(meta), path("*.bam")                , emit: bam     
    tuple val(meta), path("*.log")                , emit: log
    tuple val(meta), path("*fastq.gz")            , emit: fastq   , optional:true
    tuple val(meta), path(fasta)                  , emit: fasta
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    def unaligned = ""
    def reads_args = ""
    if (meta.single_end) {
        unaligned = save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ""
        reads_args = "-U ${reads}"
    } else {
        unaligned = save_unaligned ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ""
        reads_args = "-1 ${reads[0]} -2 ${reads[1]}"
    }

    """
    bowtie2 \\
        -x $index_base \\
        --no-unal \\
        $reads_args \\
        --threads $task.cpus \\
        $unaligned \\
        $bowtie2_options \\
        2> >(tee ${prefix}.bowtie2.log >&2) \\
    | samtools sort \\
        $args2 \\
        --threads $task.cpus \\
        -o ${prefix}.bt.${fasta.baseName}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
