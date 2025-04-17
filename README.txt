Author: Lucas Cortes
Date: 13.03.25

I am going to try to create a list of filters that LEAP has. This may not be every single filter but hopefully this is instructive. 

PolyA - Only including PolyA sites that have an actual polyA signal within 50bp of the annotated site, this is done before the pipeline begins 
Cage - None
FANTOM transcripts - None
LongRead Transcripts - None

Human Genome -  Protein Coding Only 
                No Read throughs 
                No single exons in the 3' direction 
                No transcripts that do not have a 3' or 5' UTR 
                Only CDS
                MANE should only be used if there is no other CDS in the gene, in that case MANE should be duplicated 

Considerations: 

The pipeline is supposed to create the maximal 5' and 3' ends for any given set of transcripts. It will extend one end of a transcript first, and then go back 
and look at the other end. This means that multiple transcripts of the same gene can be extended. There are certain limitations to my pipeline, it will only find
transcripts that have matching final exons with the longread and fantom data. Therefore, if there are not matching exons it will not fall back to another transcript
that is extendable, it will simply not extend. It will also ONLY extend if there is FANTOM and LongRead support for the polyA or CAGE site. Therefore, there could be an 
appropriate extension that is not selected because one of the two data types are not fully present. The pipeline only considers a CAGE or polyA site to be valid if it 
is FULLY encased by both FANTOM and LongRead data, this mean that certain sites that may be valid will be missed because they are partially encased. 
These can be manually reviewed, but normally the pipeline falls back to another valid site that is not maximal. 

There is a post nextflow script called 'makeGFF.py' that will actually make the annotation file. There are a few considerations for this. The script will tag any
transcript that was originally a MANE as "MANE_copy". This is so that we don't ever extend MANE accidentally. THIS IS VERY IMPORTANT, maybe one of the most 
important parts of the entire pipeline. Furthermore, the script will add EVERY transcript extended by LEAP to both gencode_primary and gencode_basic. It is simple 
to remove this behaviour if it becomes unwanted. Lastly, it will add a LEAP tag to every transcript that has been extended. 

Queries: lucascortes96@outlook.com

