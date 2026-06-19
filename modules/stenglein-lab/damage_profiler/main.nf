/*
  run damage profiler program to assess mismatches in mapped reads:
  see: https://pubmed.ncbi.nlm.nih.gov/33890614/
 */
process DAMAGE_PROFILER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/damageprofiler:1.1--hdfd78af_2' :
        'biocontainers/damageprofiler:1.1--hdfd78af_2' }"

    input:
    tuple val(meta), path(bam), val(refseq), path (refseq_fasta), path (fai)

    when:
    task.ext.when == null || task.ext.when

    output:
    tuple val(meta), path("*_misincorporation.txt"), emit: misincorporation

    shell:
    def new_name = bam.name.replaceAll(/.bam$/, ".mismatches.txt")

    """
    damageprofiler -i $bam -o . -r $refseq_fasta 

    # grep: ignore comment lines (start with #) and the header line (starts with Chr\tEnd)
    # awk: prepend key output files with sample ID so will be tidier format
    # sed: fix apparent bug in DamageProfiler output where tab is missing between adjacent numbers
    grep -v -e '^#' -e "^Chr	End" misincorporation.txt | awk '{print "${meta.id}" "\t" \$0}' | sed 's/\\.00\\./\\.0	0\\./' >  ${meta.id}_${refseq.id}_misincorporation.txt
    """
}


/*

DamageProfiler usage info

> damageprofiler -h
usage: DamageProfiler [-h] [-v] [-i <INPUT>] [-o <OUTPUT>] [-r <REFERENCE>] [-t <THRESHOLD>] [-s <SPECIES>] [-sf <SPECIES FILE>]
       [-l <LENGTH>] [-title <TITLE>] [-yaxis_dp_max <MAX_VALUE>] [-color_c_t <COLOR_C_T>] [-color_g_a <COLOR_G_A>]
       [-color_insertions <COLOR_C_T>] [-color_deletions <COLOR_DELETIONS>] [-color_other <COLOR_OTHER>] [-only_merged] [-sslib]

Detailed description:

 -h,--help                            Shows this help page.
 -v,--version                         Shows the version of DamageProfiler.
 -i <INPUT>                           REQUIRED. The input sam/bam/cram file.
 -o <OUTPUT>                          REQUIRED. The output folder.
 -r <REFERENCE>                       The reference file (fasta format).
 -t <THRESHOLD>                       DamagePlot: Number of bases which are considered for plotting nucleotide misincorporations.
                                      Default: 25
 -s <SPECIES>                         Reference sequence name (Reference NAME flag of SAM record). For more details see
                                      Documentation.
 -sf <SPECIES FILE>                   List of species reference names (Reference NAME flag of SAM record). For more details see
                                      Documentation.
 -l <LENGTH>                          Number of bases which are considered for frequency computations. Default: 100.
 -title <TITLE>                       Title used for all plots. Default: input filename.
 -yaxis_dp_max <MAX_VALUE>            DamagePlot: Maximal y-axis value.
 -color_c_t <COLOR_C_T>               DamagePlot: Color (HEX code) for C to T misincoporation frequency.
 -color_g_a <COLOR_G_A>               DamagePlot: Color (HEX code) for G to A misincoporation frequency.
 -color_insertions <COLOR_C_T>        DamagePlot: Color (HEX code) for base insertions.
 -color_deletions <COLOR_DELETIONS>   DamagePlot: Color (HEX code) for base deletions.
 -color_other <COLOR_OTHER>           DamagePlot: Color (HEX code) for other bases different to reference.
 -only_merged                         Use only mapped and merged (in case of paired-end sequencing) reads to calculate the damage
                                      plot instead of using all mapped reads. The SAM/BAM entry must start with 'M_', otherwise it
                                      will be skipped. Default: false
 -sslib                               Single-stranded library protocol was used. Default: false

*/
