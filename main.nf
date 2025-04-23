/*
main.nf

This script defines the main workflow for the LEAP integrated pipeline, which processes genomic data 
for both 5' and 3' ends simultaneously. The pipeline uses a sample sheet (CSV file) as input and 
supports configurable parameters to control the workflow. It integrates subworkflows for processing 
3' and 5' ends and allows for iterative refinement of results.

Usage:
    nextflow run main.nf --input <input.csv> --outputDir <output_directory>

Arguments:
    input.csv           A CSV file containing input data with columns for identifiers and file paths.
    outputDir           Directory where the results will be published (default: "results").

Steps:
1. Load the input CSV file and parse it into channels for 3' and 5' processing based on the "End" column.
2. Split the data into separate channels for 3' and 5' workflows.
3. Call subworkflows (`THREE_PRIME_PIPELINE` and `FIVE_PRIME_PIPELINE`) to process the data.
4. Combine results from the first run of each subworkflow and refine the data for a second run.
5. Publish the results from both runs into the specified output directory.

Subworkflows:
    - THREE_PRIME_PIPELINE: Processes data for 3' ends.
    - FIVE_PRIME_PIPELINE: Processes data for 5' ends.

Output:
    - Results are organized into directories for each run (`run1` and `run2`) and for each type of processing 
      (`three_prime` and `five_prime`).

Dependencies:
    - Nextflow: For workflow orchestration.
    - Input CSV file: Contains the data to be processed.

Example:
    nextflow run main.nf --input input.csv --outputDir results
*/

chromosomes = Channel.from('1', '2',
'3', '4', '5', '6', '7', '8', '9', 
'10','11', '12', '13', '14', '15', 
'16', '17', '18', '19', '20','21', 
'22', 'X', 'Y' ) // ADD X AND Y BACK IN LATER  'X', 'Y'


include { THREE_PRIME_PIPELINE } from './subworkflows/three_prime_pipeline'
include { FIVE_PRIME_PIPELINE } from './subworkflows/five_prime_pipeline'

params.outputDir = 'results_DFbrainAndMixture'

workflow {
    prep_next = true
    // three prime is always single exon true, maybe not five prime 
    single_exon = true
    // Define the path to your CSV file
    csv_file = channel.fromPath("input_csvs/inputDFbrainAndMixture.csv").splitCsv( header: true, strip:true) //If no strip:true then any white space will void this function
    .view()
    .map { row ->
            def end = row.End.toLowerCase()
            def is_three = end.contains('three') || end.toLowerCase().contains('3prime') || end == '3'
            def is_five = end.contains('five') || end.toLowerCase().contains('5prime') || end == '5'
            [
                is_three ? 'three' : 'five',
                [
                    id: row.End,
                    human: row.human,
                    capOrTail: row.capOrTail,
                    fantom: row.fantom ,
                    longRead: row.longRead 
                ]
            ]
        }.view()

        // Split the channel into 'three' and 'five' channels
        three_ch = csv_file.filter { it[0] == 'three' }.map { it[1] }
        .view()
        five_ch = csv_file.filter { it[0] == 'five' }.map { it[1] }
        .view()

            // Call subworkflows
        three_prime_results_one = THREE_PRIME_PIPELINE(three_ch, prep_next, single_exon, chromosomes)
        five_prime_results_one = FIVE_PRIME_PIPELINE(five_ch, prep_next, single_exon, chromosomes)
        new_three_ch = three_ch
            .combine(five_prime_results_one)
            .map { original, new_human ->
                original.human = new_human
                original.id = original.id + "_modified" // Append "_modified" to the ID
                return original
            }
            .view()
        new_five_ch = five_ch
            .combine(three_prime_results_one)
            .map { original, new_human ->
                original.human = new_human
                original.id = original.id + "_modified" // Append "_modified" to the ID
                return original
            }
            .view()
        prep_next = false
        three_prime_results_two = THREE_PRIME_PIPELINE(new_three_ch, prep_next, single_exon, chromosomes)
        five_prime_results_two = FIVE_PRIME_PIPELINE(new_five_ch, prep_next, single_exon, chromosomes)

    // Publish outputs
    PUBLISH_RESULTS(three_prime_results_one, five_prime_results_one, three_prime_results_two, five_prime_results_two)
}

// New process to handle publishing results
process PUBLISH_RESULTS {
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    path three_prime_run1
    path five_prime_run1
    path three_prime_run2
    path five_prime_run2

    output:
    path "run1/three_prime/*"
    path "run1/five_prime/*"
    path "run2/three_prime/*"
    path "run2/five_prime/*"

    script:
    """
    mkdir -p run1/three_prime run1/five_prime run2/three_prime run2/five_prime
    mv ${three_prime_run1} run1/three_prime/
    mv ${five_prime_run1} run1/five_prime/
    mv ${three_prime_run2} run2/three_prime/
    mv ${five_prime_run2} run2/five_prime/
    """
}
