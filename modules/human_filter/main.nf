process HUMAN_FILTER{
    input:
    tuple val(id), path(human), path(capOrTail), path(fantom), path(longRead)
    path readthroughs
    val single_exon
    output:
    tuple val(id), path('noReadthroughProteinCoding.gff3'), path(capOrTail), path(fantom), path (longRead)
    """
    humanFilter.py ${human} 'noReadthroughProteinCoding.gff3' ${readthroughs} ${single_exon}
    """
}