process CAT_ALL {
    publishDir 'outputs/catted', mode: 'copy', overwrite: true
    input:
    tuple val(id), path(results)
    

    output:
    tuple val(id), path("${id}_result.csv")

    """
    cat ${results.join(' ')} > "${id}_result.csv"
    """

}