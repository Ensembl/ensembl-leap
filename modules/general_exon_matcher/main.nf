process GENERAL_EXON_MATCHER{
    publishDir 'outputs/exonMatched', mode: 'copy', overwrite: true
    memory '40 GB'
    cpus 2
    input:
    tuple val(id), path(human), path(capOrTail), path(fantom), path (longRead)
    val single_exon
    val direction
    output:
    tuple val(id), path("filtered_matched_human_exons.gff"), path(capOrTail), path("matched_fantom_blocks.gff"), path("matched_longread_blocks.gff") 
    """
    globalExonMatcher.py ${human} ${fantom} ${longRead} ${single_exon} ${direction} .
    """
}