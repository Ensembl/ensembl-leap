process PREP_NEXT {
    publishDir 'outputs/prep_next', mode: 'copy', overwrite: true
    input:
    path (csv)
    val (direction)
    tuple val(id), path(human), path(polyA), path(fantom), path(encode)
    output:
    path "${direction}_nextRun.gff"

    """
    prepNext.py ${csv} ${direction} ${human}
    """
}