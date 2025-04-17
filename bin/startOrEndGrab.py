#!/usr/bin/env python3
import pandas as pd
import sys

'''
Author: Lucas Cortes
Date: 2020-10-15
Usage: python 3primeGrab.py <input_file> <fiveOrThreePrime?> <output_file>

This script is intended to find the furthest threePrime or fivePrime transcript for each gene in the Human Genome
The result will be a GTF with an additional column that contains the gene name associated with the transcript 
'''

def main(input_file, capOrTail, output_file):
    # Define the column names for the GFF file
    gff_column_names = [
        "seqname", "source", "feature", "Start", "End", "score", "Strand", "frame", "Attributes"
    ]

    # Read the GFF file with the specified column names
    df = pd.read_csv(input_file, sep='\t', names=gff_column_names, comment='#', header=None)
    print(df.head())
    # Initialize variables
    ensembl_gene_ids = []
    current_ensg = "NA"

    # Loop through the GFF file to extract gene IDs and add them to the DataFrame
    for idx, row in df.iterrows():
        if row['feature'] == 'gene':
            current_ensg = row['Attributes'].split('ID=gene:')[1].split(';')[0]
            
        ensembl_gene_ids.append(current_ensg)
        if row['Attributes'].startswith('###'):
            current_ensg = "NA"  # Reset the gene ID at the end of a block

    # Add the ensembl_gene_id column to the DataFrame
    df['ensembl_gene_id'] = ensembl_gene_ids
    print(df)

    # Group by GeneID and apply a lambda function to select the most 5' transcript
    def select_most_3_transcript(group):
        # Extract transcript IDs associated with both 5' and 3' UTRs
        five_prime_utrs = group[group['feature'] == 'five_prime_UTR']['Attributes'].str.extract('Parent=transcript:([^;]+)')[0]
        print(five_prime_utrs)
        three_prime_utrs = group[group['feature'] == 'three_prime_UTR']['Attributes'].str.extract('Parent=transcript:([^;]+)')[0]
        valid_transcripts = set(five_prime_utrs.dropna()) & set(three_prime_utrs.dropna())

        # Filter the group for valid transcripts
        valid_transcript_group = group[group['Attributes'].apply(lambda x: any(f'Parent=transcript:{tid}' in x for tid in valid_transcripts))]

        if not valid_transcript_group.empty:
            # Select the most 5' CDS from the valid transcripts
            three_prime_utrs = valid_transcript_group[valid_transcript_group['feature'] == 'exon']
            if not three_prime_utrs.empty:
                if group['Strand'].iloc[0] == '+':
                    most_3_transcript = three_prime_utrs.loc[three_prime_utrs['End'].idxmax()]
                else:
                    most_3_transcript = three_prime_utrs.loc[three_prime_utrs['Start'].idxmin()]

                if 'Parent=transcript:' in most_3_transcript['Attributes']:
                    transcript_id = most_3_transcript['Attributes'].split('Parent=transcript:')[1].split(';')[0]
                else:
                    return None
                transcript_group = group[group['Attributes'].str.contains(f'{transcript_id}')]
                print(transcript_group, "/n", "TRANSCRIPT GROUP")

                has_mane = any((transcript_group['feature'] == 'mRNA') & (transcript_group['Attributes'].str.contains('MANE_Select')))
                if has_mane:
                    print(has_mane, "/n" ,"HAS MANE")
                
                if has_mane:
                    most_3_transcript['Attributes'] = most_3_transcript['Attributes'].replace(transcript_id, transcript_id + '_MANE_copy')
                
                return most_3_transcript
            else:
                return None
        else:
            return None
    # Group by GeneID and apply a lambda function to select the most 5' transcript
    def select_most_5_transcript(group):
        # Extract transcript IDs associated with both 5' and 3' UTRs
        five_prime_utrs = group[group['feature'] == 'five_prime_UTR']['Attributes'].str.extract('Parent=transcript:([^;]+)')[0]
        three_prime_utrs = group[group['feature'] == 'three_prime_UTR']['Attributes'].str.extract('Parent=transcript:([^;]+)')[0]
        valid_transcripts = set(five_prime_utrs.dropna()) & set(three_prime_utrs.dropna())

        # Filter the group for valid transcripts
        valid_transcript_group = group[group['Attributes'].apply(lambda x: any(f'Parent=transcript:{tid}' in x for tid in valid_transcripts))]

        if not valid_transcript_group.empty:
            # Select the most 5' CDS from the valid transcripts
            five_prime_utrs = valid_transcript_group[valid_transcript_group['feature'] == 'exon']
            if not five_prime_utrs.empty:
                if group['Strand'].iloc[0] == '+':
                    most_5_transcript = five_prime_utrs.loc[five_prime_utrs['Start'].idxmin()]
                else:
                    most_5_transcript = five_prime_utrs.loc[five_prime_utrs['End'].idxmax()]

                if 'Parent=transcript:' in most_5_transcript['Attributes']:
                    transcript_id = most_5_transcript['Attributes'].split('Parent=transcript:')[1].split(';')[0]
                else:
                    return None
                transcript_group = group[group['Attributes'].str.contains(f'{transcript_id}')]
                print(transcript_group, "/n", "TRANSCRIPT GROUP")

                has_mane = any((transcript_group['feature'] == 'mRNA') & (transcript_group['Attributes'].str.contains('MANE_Select')))
                if has_mane:
                    print(has_mane, "/n" ,"HAS MANE")
                
                if has_mane:
                    most_5_transcript['Attributes'] = most_5_transcript['Attributes'].replace(transcript_id, transcript_id + '_MANE_copy')
                
                return most_5_transcript
            else:
                return None
        else:
            return None

    if capOrTail == 'threePrime':
        result = df.groupby('ensembl_gene_id').apply(select_most_3_transcript)
    elif capOrTail == 'fivePrime':
        result = df.groupby('ensembl_gene_id').apply(select_most_5_transcript)
    else:
        print("Invalid option for fiveOrThreePrime. Please use 'fivePrime' or 'threePrime'.")
        sys.exit(1)

    # Filter out None values
    result = result.dropna().reset_index(drop=True)

    # Save the result to a new file
    result.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python 3primeGrab.py <input_file> <fiveOrThreePrime> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    capOrTail = sys.argv[2]
    output_file = sys.argv[3]
    main(input_file, capOrTail, output_file)
