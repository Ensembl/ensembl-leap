process START_OR_END_GRAB {
    input:
    tuple val(id), path(noReadThrough), path(capOrTail), path(fantom), path(longRead)
    val three
    output:
    tuple val(id), path ("grabbedhg38.gff"), path(capOrTail), path(fantom), path(longRead)
    """
    startOrEndGrab.py ${noReadThrough} ${three} grabbedhg38.gff
    """
}