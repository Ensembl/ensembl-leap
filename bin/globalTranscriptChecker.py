#!/usr/bin/env python3
'''
Author: Lucas Cortes
Date: 2020-10-15

Usage: python globalTranscriptChecker.py <human_transcripts> <fantom> <longread_transcripts>  
<capOrTail/capOrTail_transcripts> <fiveprimeOrThreeprime?> <chromosome> <output_directory> 

This script will check in order:
1. If there is a capOrTail peak or capOrTail site 5' or 3' of the selected Human Transcript 
2. If there is a FANTOM transcript that matches the Human Transcript and consumes a capOrTail or capOrTail site
3. If there is a long read transcript that matches the Human Transcript and consumes a capOrTail or capOrTail site

It does this human transcript by human transcript. The files coming in must be exon matched already, 
so there should only be human, FANTOM, and long read transcripts that have matching acceptor sites. 
This script does an extra check for that to make sure that the transcripts are in the correct order.

'''


import pandas as pd
import sys

def strip_chr_prefix(df):
    df.iloc[:, 0] = df.iloc[:, 0].astype(str)
    if df.iloc[:, 0].str.contains('chr').any():
        df.iloc[:, 0] = df.iloc[:, 0].str.replace('chr', '')
    return df
# Convert numeric chromosomes to int, leave 'X' and 'Y' as str
def convert_chromosome(chromosome):
    try:
        return int(chromosome)
    except ValueError:
        return chromosome
def process_dataframe(df, column_names):
    #df.iloc[:, 0] = df.iloc[:, 0].apply(convert_chromosome)
    df.iloc[:, 3] = pd.to_numeric(df.iloc[:, 3])
    df.iloc[:, 4] = pd.to_numeric(df.iloc[:, 4])
    df.iloc[:, 6] = df.iloc[:, 6].astype(str)
    df.columns = column_names + list(df.columns[len(column_names):])
    df = df[column_names]
    return df

def importGffs(human_file, capOrTail_file, fantom_file, longRead_file):
    human = pd.read_csv(human_file, sep='\t', skiprows=1)
    capOrTail = pd.read_csv(capOrTail_file, sep='\t', header=None).iloc[:, :9]
    fantom = pd.read_csv(fantom_file, sep='\t', header=None).iloc[:, :9]
    longRead = pd.read_csv(longRead_file, sep='\t', header=None).iloc[:, :9]
    

    

    # Strip unwatned 'chr' prefix
    human = strip_chr_prefix(human)
    fantom = strip_chr_prefix(fantom)
    longRead = strip_chr_prefix(longRead)
    capOrTail = strip_chr_prefix(capOrTail)


    human_column_names = ['Chromosome', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes', 'gene_id']
    capOrTail_column_names = ['Chromosome', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes']
    fantom_column_names = ['Chromosome', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes']
    longRead_column_names = ['Chromosome', 'Source', 'Type', 'LONGREAD_Start', 'LONGREAD_End', 'Score', 'Strand', 'Phase', 'Attributes']


    human = process_dataframe(human, human_column_names)
    capOrTail = process_dataframe(capOrTail, capOrTail_column_names)
    fantom = process_dataframe(fantom, fantom_column_names) 
    longRead = process_dataframe(longRead, longRead_column_names)

    fantom['Name'] = fantom['Attributes'].str.extract('Name="([^"]*)')
    longRead['Name'] = longRead['Attributes'].str.extract('gene_id\s+"([^"]+)"')
    human['Name'] = human['Attributes'].str.extract('Parent=transcript:(.*?);')
    print(human.head(), capOrTail.head(), fantom.head(), longRead.head())
    return human, capOrTail, fantom, longRead

def findMatchesFivePrime(human, fantom, longRead, capOrTail):
    results = pd.DataFrame()
    # Loop over human exons
    for i, exon in human.iterrows():
        # Filter capOrTail based on the strand, chromosome, position, within 10000bp of the exon
        capOrTail_filtered = capOrTail[
            (capOrTail['Strand'] == exon['Strand']) &
            (capOrTail['Chromosome'] == exon['Chromosome']) &
            ((capOrTail['End'] < exon['Start']) if exon['Strand'] == '+' else (capOrTail['Start'] > exon['End'])) &
            ((capOrTail['Start'] >= (exon['Start'] - 10000)) if exon['Strand'] == '+' else (capOrTail['End'] <= (exon['End'] + 10000)))
        ]

        for j, capOrTail_site in capOrTail_filtered.iterrows():
            print(exon['Start'])
            print(fantom['Start'])
            # Filter FANTOM based on strand, chromosome, position, within 10000bp of the capOrTail site
            fantom_filtered = fantom[
                (fantom['Strand'] == exon['Strand']) &
                (fantom['Chromosome'] == exon['Chromosome']) &
                (fantom['End'] == exon['End'] if exon['Strand'] == '+' else fantom['Start'] == exon['Start']) &
                (fantom['Start'] < (capOrTail_site['Start']) if exon['Strand'] == '+' else fantom['End'] > (capOrTail_site['End']))
            ]
            fantom_filtered['capOrTail_Start'] = capOrTail_site['Start']
            fantom_filtered['capOrTail_End'] = capOrTail_site['End']

            for k, fantom_site in fantom_filtered.iterrows():
                # Filter LONGREAD based on strand, chromosome, position, within 10000bp of the capOrTail site
                longRead_filtered = longRead[
                    (longRead['Strand'] == exon['Strand']) &
                    (longRead['Chromosome'] == exon['Chromosome']) &
                    (longRead['LONGREAD_End'] == exon['End'] if exon['Strand'] == '+' else longRead['LONGREAD_Start'] == exon['Start']) &
                    (longRead['LONGREAD_Start'] < (fantom_site['capOrTail_Start']) if exon['Strand'] == '+' else longRead['LONGREAD_End'] > (fantom_site['capOrTail_End']))

                ].copy()  # Create a copy to avoid SettingWithCopyWarning ***
                
            # Add the capOrTail/fantom/longRead start & ends to the results
                if not longRead_filtered.empty:
                    result = exon.copy()
                    result['gene_id'] = exon['gene_id']
                    result['capOrTail_Start'] = capOrTail_site['Start']
                    result['capOrTail_End'] = capOrTail_site['End']
                    result['Transcript_Start'] = exon['Start']
                    result['Transcript_End'] = exon['End']
                    result['Transcript_Name'] = exon['Name']  # Assuming 'Name' column exists, otherwise default to 'Unknown'
                    results = pd.concat([results, pd.DataFrame([result])], ignore_index=True)
    return results

def findMatchesThreePrime(human, fantom, longRead, capOrTail):
    # Initialize an empty dataframe to store the results
    results = pd.DataFrame()
    # Ensure no NaN values in critical columns

    # Loop over human exons
    for i, exon in human.iterrows():
        capOrTail_filtered = capOrTail[
            (capOrTail['Strand'] == exon['Strand']) &
            (capOrTail['Chromosome'] == exon['Chromosome']) &
            ((capOrTail['Start'] > exon['End']) if exon['Strand'] == '+' else (capOrTail['End'] < exon['Start'])) &
            ((capOrTail['Start'] <= (exon['Start'] + 10000)) if exon['Strand'] == '+' else (capOrTail['End'] >= (exon['End'] - 10000)))
            
        ]
        # Debug filtered results
        print(capOrTail['Strand'], capOrTail['Start'], exon['End'])


        # Loop over filtered capOrTail sites
        for j, capOrTail_site in capOrTail_filtered.iterrows():
            # Filter FANTOM transcripts based on the strand, chromosome, and position
            fantom_filtered = fantom[
                (fantom['Strand'] == exon['Strand']) &
                (fantom['Chromosome'] == exon['Chromosome']) &
                ((fantom['Start'] == exon['Start']) if exon['Strand'] == '+' else (fantom['End'] == exon['End'])) &
                ((fantom['End'] >= capOrTail_site['Start']) if exon['Strand'] == '+' else (fantom['Start'] <= capOrTail_site['Start']))
            ]
            # Add the matching capOrTail Start to the fantom_filtered DataFrame
            fantom_filtered['capOrTail_Start'] = capOrTail_site['Start']
            fantom_filtered['capOrTail_End'] = capOrTail_site['End']

            # Filter longRead transcripts based on the same criteria as FANTOM transcripts
            longRead_filtered = longRead[
                (longRead['Strand'] == exon['Strand']) &
                (longRead['Chromosome'] == exon['Chromosome']) &
                ((longRead['LONGREAD_Start'] == exon['Start']) if exon['Strand'] == '+' else (longRead['LONGREAD_End'] == exon['End'])) &
                ((longRead['LONGREAD_End'] >= capOrTail_site['Start']) if exon['Strand'] == '+' else (longRead['LONGREAD_Start'] <= capOrTail_site['Start']))
            ]
            # Add the matching capOrTail Start to the fantom_filtered DataFrame
            longRead_filtered['capOrTail_Start'] = capOrTail_site['Start']
            # Append the filtered rows and drop duplicates
            combined_filtered = fantom_filtered[fantom_filtered['capOrTail_Start'].isin(longRead_filtered['capOrTail_Start'])]


            

            # If there are matching transcripts, add them to the results
            if not combined_filtered.empty:
                for k, transcript in combined_filtered.iterrows():
                    result = exon.copy()
                    result['gene_id'] = exon['gene_id']
                    result['capOrTail_Start'] = capOrTail_site['Start']
                    result['capOrTail_End'] = capOrTail_site['End']
                    result['Transcript_Start'] = exon['Start']
                    result['Transcript_End'] = exon['End']
                    result['Transcript_Name'] = exon['Name']  # Assuming 'Name' column exists, otherwise default to 'Unknown'
                    results = pd.concat([results, pd.DataFrame([result])], ignore_index=True)
    return(results)

def main():
    output_file = sys.argv[7]
    chromosome_value = sys.argv[6]
    #chromosome_value = int(chromosome_value)
    if len(sys.argv) > 4:
        arg = sys.argv[5].lower()
        if arg in ['fiveprime', '5', "5'"]:
            direction = 'fiveprime'
        elif arg in ['threeprime', '3', "3'"]:
            direction = 'threeprime'
            print("threeprime")
        else:
            raise ValueError("Invalid direction argument. Use 'fiveprime', 'threeprime', '5', '3', '5\' or '3\'.")
    else:
        direction = 'fiveprime'  # Default value
    imported = importGffs(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    human = imported[0]
    capOrTail = imported[1]
    fantom = imported[2]
    longRead = imported[3]
    # Filter dataframes based on the chromosome value
    human = human[human['Chromosome'] == chromosome_value]
    capOrTail = pd.DataFrame(capOrTail.loc[capOrTail['Chromosome'] == chromosome_value, :])  # Ensure capOrTail is a DataFrame
    fantom = fantom[fantom['Chromosome'] == chromosome_value]
    longRead = longRead[longRead['Chromosome'] == chromosome_value]
        # Determine which function to use based on the direction
    if direction == 'fiveprime':
        findMatches = findMatchesFivePrime
    elif direction == 'threeprime':
        print("threeprime")
        findMatches = findMatchesThreePrime
    else:
        raise ValueError("Invalid direction argument. Use 'fiveprime' or 'threeprime'.")
    matches = findMatches(human,fantom, longRead, capOrTail)

    output_file = f"{output_file}_matched_chr{chromosome_value}.csv"
    matches.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    main()