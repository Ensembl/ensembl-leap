//Returns split transcript, polyA and chr for next process
process SPLIT_CHROMOSOMES {
    errorStrategy 'ignore'
    input:
    tuple val(id), path(human), path(capOrTail), path(fantom), path(longRead), val(chr)
    val three
    output:
    tuple path('split_human_*'), path('split_capOrTail_*'), path('split_fantom_*'), path('split_longRead_*'), val(id), val(chr)
    """
    splitChromosomes.py  ${chr} ${fantom} ${longRead} ${capOrTail} ${human} ${three}
    """
}