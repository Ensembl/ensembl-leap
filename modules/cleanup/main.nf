process CLEANUP {
    publishDir 'outputs/cleanup', mode: 'copy', overwrite: true
    input:
    tuple val(id), path(result)
    val direction
    output:
    path "${id}_extended_transcripts_${direction}_final.csv", emit:csv
    val(id), emit:id
    path "${id}_extended_transcripts_${direction}_ExtendStats.txt"
    path "${id}_extended_transcripts_${direction}_ExtendPlot.png"

    """
    finalFilterandStats.py ${result} ${direction} "${id}_extended_transcripts.csv"
    """
}