include {HUMAN_FILTER} from '../../modules/human_filter'
include {START_OR_END_GRAB} from '../../modules/start_or_end_grab'
include {SPLIT_CHROMOSOMES} from '../../modules/split_chromosomes'
include {PROCESS_CHROMOSOMES} from '../../modules/process_chromosomes'
include {CAT_ALL} from '../../modules/cat_all'
include {CLEANUP} from '../../modules/cleanup'
include {PREP_NEXT} from '../../modules/prep_next'
include {GENERAL_EXON_MATCHER} from '../../modules/general_exon_matcher'


include {HUMAN_FILTER as HUMAN_FILTER_2} from '../../modules/human_filter'
include {START_OR_END_GRAB as START_OR_END_GRAB_2} from '../../modules/start_or_end_grab'
include {SPLIT_CHROMOSOMES as SPLIT_CHROMOSOMES_2} from '../../modules/split_chromosomes'
include {PROCESS_CHROMOSOMES as PROCESS_CHROMOSOMES_2} from '../../modules/process_chromosomes'
include {CAT_ALL as CAT_ALL_2} from '../../modules/cat_all'
include {CLEANUP as CLEANUP_2} from '../../modules/cleanup'
include {PREP_NEXT as PREP_NEXT_2} from '../../modules/prep_next'
include {GENERAL_EXON_MATCHER as GENERAL_EXON_MATCHER_2} from '../../modules/general_exon_matcher'

workflow THREE_PRIME_PIPELINE {
    
    take:
    ch_three
    prep_next
    single_exon
    chromosomes
    
    main:
    readThroughs = file("/nfs/production/flicek/ensembl/havana/lucascortes/polyA-DB/data/readthroughList/readthroughList.txt")
    three = "threePrime"
    if (prep_next){
        humanOut = HUMAN_FILTER(ch_three, readThroughs, single_exon)
        threePrimeOut = START_OR_END_GRAB(humanOut, three).view()
        generalOut = GENERAL_EXON_MATCHER(threePrimeOut, single_exon, three)
        generalChromosomes = generalOut.combine(chromosomes).view()
        splitChrs = SPLIT_CHROMOSOMES(generalChromosomes, three).view()
        processChrOut = PROCESS_CHROMOSOMES(splitChrs, three)
        .groupTuple(by: 0)  // Group by ID (index 0)
        .map { id, chrs, files -> tuple(id, files.flatten()) }  // Flatten the list of files
        .view()
        catted = CAT_ALL(processChrOut)
        cleaned = CLEANUP(catted, three) 
        prepnext_out = PREP_NEXT(cleaned.csv, three, ch_three)
    } else {
        generalOut = GENERAL_EXON_MATCHER_2(ch_three, single_exon, three)
        generalChromosomes = generalOut.combine(chromosomes).view()
        splitChrs = SPLIT_CHROMOSOMES_2(generalChromosomes, three).view()
        processChrOut = PROCESS_CHROMOSOMES_2(splitChrs, three)
        .groupTuple(by: 0)  // Group by ID (index 0)
        .map { id, chrs, files -> tuple(id, files.flatten()) }  // Flatten the list of files
        .view()
        catted = CAT_ALL_2(processChrOut)
        cleaned = CLEANUP_2(catted, three) 
        prepnext_out = cleaned.csv
    }
    emit:
    prepnext_out
}