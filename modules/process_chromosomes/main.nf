process PROCESS_CHROMOSOMES {
    publishDir 'outputs/processedChrs', mode: 'copy', overwrite: true
    input:
    tuple path(human), path(capOrTail), path(fantom), path(longRead), val(id), val(chr)
    val (direction)
    
    output:
    tuple val(id), val(chr), path('output_*')
    

    """
    globalTranscriptChecker.py ${human} ${capOrTail} ${fantom} ${longRead} ${direction} ${chr} output  
    """
}