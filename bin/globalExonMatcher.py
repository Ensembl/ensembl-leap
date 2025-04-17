#!/usr/bin/env python3
'''
Author: Lucas Cortes
Date: 2020-10-15
Usage: python globalExonMatcher.py <human_exons.gff> <FANTOM_exons.gff> <long_read_exons.gff> <single_exon?> <fiveprimeOrThreeprime?> <output_directory>

This script is used to match exons of incoming files in both the 3' and 5' direction 
so that when the outputs are passed to the next script, we have matching acceptor 
splice sites in the 5' or 3' direction. 

The script takes in three GFF files, one for human exons, one for FANTOM exons, and one 
for the long read data (e.g. Nanopore or PacBio). The script then matches the exons and returns
a filtered_matched_human_exons file which contains the human exons that are matched for both 
the FANTOM and long read data.
'''

import pandas as pd
import re
import sys
import os
import subprocess
import argparse
import csv

def validate_gff(file_path):
    result = subprocess.run(['gffread', file_path, '-E'], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Validation failed for {file_path}:\n{result.stderr}")
    else:
        print(f"{file_path} is valid.")

def preprocess_gff(file_path):
    processed_file = 'processed_' + os.path.basename(file_path)
    with open(file_path, 'r') as infile, open(processed_file, 'w') as outfile:
        for line in infile:
            if not line.startswith('#'):
                outfile.write(line)
    return processed_file

def filter_by_ensg(df1, df2):
    ensg_set = set(df2.iloc[:, -2])
    filtered_df = df1[df1.iloc[:, -2].isin(ensg_set)]
    return filtered_df

def extract_transcript_id_and_exon_number(df):
    if 'Name=' in df.iloc[0, 8]:
        df['transcript_id'] = df[8].str.extract('Name="(.*?)\..*?"')
        df['block_num'] = df[8].str.extract('Name=".*?_block(.*?)"')
    else:
        df['transcript_id'] = df[8].str.extract('transcript_id "(.*?)"')
        df['block_num'] = df[8].str.extract('exon_number "(.*?)"')
    return df

def match_exons_with_blocks_threeprime(human_df, block_df, single_exon):
    # Check if 'chr' is present in any of the entries in the column
    if block_df[0].str.contains('chr').any():
        block_df[0] = block_df[0].str.replace('^chr', '', regex=True)
    print(block_df.head())
    human_df['transcript_id'] = human_df['Attributes'].str.extract('Parent=transcript:(.*?);')
    block_df = extract_transcript_id_and_exon_number(block_df)

    block_df['block_num'] = pd.to_numeric(block_df['block_num'], errors='coerce')
    block_df = block_df.dropna(subset=['block_num'])
    # FILTER OUT SINGLE EXON GENES
    if single_exon == True:
        max_block_num = block_df.groupby('transcript_id')['block_num'].transform('max')
        block_df = block_df[max_block_num > 1]
    blocks_forward = block_df[block_df[6] == '+']
    blocks_reverse = block_df[block_df[6] == '-']
    last_blocks_forward = blocks_forward.loc[blocks_forward.groupby('transcript_id')['block_num'].idxmax()]
    last_blocks_reverse = blocks_reverse.loc[blocks_reverse.groupby('transcript_id')['block_num'].idxmin()]

    last_blocks_forward[0] = last_blocks_forward[0].astype(str)
    last_blocks_reverse[0] = last_blocks_reverse[0].astype(str)

    matched_human_exons_forward = human_df[human_df['Strand'] == '+']
    matched_human_exons_forward = matched_human_exons_forward[matched_human_exons_forward['Start'].isin(last_blocks_forward[3])]
    matched_blocks_forward = last_blocks_forward[last_blocks_forward[3].isin(matched_human_exons_forward['Start'])]

    matched_human_exons_reverse = human_df[human_df['Strand'] == '-']
    matched_human_exons_reverse = matched_human_exons_reverse[matched_human_exons_reverse['End'].isin(last_blocks_reverse[4])]
    matched_blocks_reverse = last_blocks_reverse[last_blocks_reverse[4].isin(matched_human_exons_reverse['End'])]

    matched_human_exons_forward = matched_human_exons_forward.sort_values(by=['Start'])
    matched_human_exons_reverse = matched_human_exons_reverse.sort_values(by=['End'])
    matched_blocks_forward = matched_blocks_forward.sort_values(by=[3])
    matched_blocks_reverse = matched_blocks_reverse.sort_values(by=[4])

    matched_human_exons = pd.concat([matched_human_exons_forward, matched_human_exons_reverse])
    matched_blocks = pd.concat([matched_blocks_forward, matched_blocks_reverse])

    return matched_human_exons, matched_blocks

def match_exons_with_blocks_fiveprime(human_df, block_df, single_exon):
    # Check if 'chr' is present in any of the entries in the column
    if block_df[0].str.contains('chr').any():
        block_df[0] = block_df[0].str.replace('^chr', '', regex=True)
    human_df['transcript_id'] = human_df['Attributes'].str.extract('Parent=transcript:(.*?);')
    block_df = extract_transcript_id_and_exon_number(block_df)

    block_df['block_num'] = pd.to_numeric(block_df['block_num'], errors='coerce')
    block_df = block_df.dropna(subset=['block_num'])
    # FILTER OUT SINGLE EXON GENES
    if single_exon == True:
        max_block_num = block_df.groupby('transcript_id')['block_num'].transform('max')
        block_df = block_df[max_block_num > 1]
    blocks_forward = block_df[block_df[6] == '+']
    blocks_reverse = block_df[block_df[6] == '-']
    first_blocks_forward = blocks_forward.loc[blocks_forward.groupby('transcript_id')['block_num'].idxmin()]
    first_blocks_reverse = blocks_reverse.loc[blocks_reverse.groupby('transcript_id')['block_num'].idxmax()]
    first_blocks_forward[0] = first_blocks_forward[0].astype(str)
    first_blocks_reverse[0] = first_blocks_reverse[0].astype(str)
    matched_human_exons_forward = human_df[human_df['Strand'] == '+']
    matched_human_exons_forward = matched_human_exons_forward[matched_human_exons_forward['End'].isin(first_blocks_forward[4])]
    matched_blocks_forward = first_blocks_forward[first_blocks_forward[4].isin(matched_human_exons_forward['End'])]
    matched_human_exons_reverse = human_df[human_df['Strand'] == '-']
    matched_human_exons_reverse = matched_human_exons_reverse[matched_human_exons_reverse['Start'].isin(first_blocks_reverse[3])]
    matched_blocks_reverse = first_blocks_reverse[first_blocks_reverse[3].isin(matched_human_exons_reverse['Start'])]
    matched_human_exons_forward = matched_human_exons_forward.sort_values(by=['End'])
    matched_human_exons_reverse = matched_human_exons_reverse.sort_values(by=['Start'])
    matched_blocks_forward = matched_blocks_forward.sort_values(by=[4])
    matched_blocks_reverse = matched_blocks_reverse.sort_values(by=[3])
    matched_human_exons = pd.concat([matched_human_exons_forward, matched_human_exons_reverse])
    matched_blocks = pd.concat([matched_blocks_forward, matched_blocks_reverse])

    return matched_human_exons, matched_blocks

def write_file(file_path, data):
    data.to_csv(file_path, sep='\t', index=False, header=False)
def main():
    if len(sys.argv) < 5:
        print("Usage: python exonMatcher.py <human_file.gff> <fantom_file.gff> <longread_file.gff> <single_exon?>  <fiveprimeOrThreeprime?> <output_directory>")
        sys.exit(1)

    human_file = sys.argv[1]
    fantom_file = sys.argv[2]
    longread_file = sys.argv[3]
    single_exon = sys.argv[4].lower() == 'true' if len(sys.argv) > 4 else True
    if len(sys.argv) > 4:
        arg = sys.argv[5].lower()
        if arg in ['fiveprime', '5', "5'"]:
            direction = 'fiveprime'
        elif arg in ['threeprime', '3', "3'"]:
            direction = 'threeprime'
        else:
            raise ValueError("Invalid direction argument. Use 'fiveprime', 'threeprime', '5', '3', '5\' or '3\'.")
    else:
        direction = 'fiveprime'  # Default value
    output_dir = sys.argv[6] if len(sys.argv) > 6 else os.getcwd()

    validate_gff(human_file)
    validate_gff(fantom_file)
    validate_gff(longread_file)

    processed_human_file = preprocess_gff(human_file)
    processed_fantom_file = preprocess_gff(fantom_file)
    processed_longread_file = preprocess_gff(longread_file)

    human_df = pd.read_csv(processed_human_file, sep='\t', on_bad_lines='skip')
    fantom_df = pd.read_csv(processed_fantom_file, sep='\t', header=None, comment='#', on_bad_lines='skip')
    longread_df = pd.read_csv(processed_longread_file, sep='\t', header=None, comment='#', on_bad_lines='skip')

    # Determine which function to use based on the direction
    if direction == 'fiveprime':
        match_exons_with_blocks = match_exons_with_blocks_fiveprime
    elif direction == 'threeprime':
        match_exons_with_blocks = match_exons_with_blocks_threeprime
    else:
        raise ValueError("Invalid direction argument. Use 'fiveprime' or 'threeprime'.")

    # Call the appropriate function
    matched_human_exons_fantom, matched_fantom_blocks = match_exons_with_blocks(human_df, fantom_df, single_exon)
    matched_human_exons_longread, matched_longread_blocks = match_exons_with_blocks(human_df, longread_df, single_exon)

    output_dir = os.path.join(output_dir, '')
    matched_human_exons_fantom.to_csv(output_dir + 'matched_human_exons_fantom.gff', sep='\t', index=False, header=False)
    matched_fantom_blocks.to_csv(output_dir + 'matched_fantom_blocks.gff', sep='\t', index=False, header=False)
    matched_human_exons_longread.to_csv(output_dir + 'matched_human_exons_longread.gff', sep='\t', index=False, header=False)
    matched_longread_blocks.to_csv(output_dir + 'matched_longread_blocks.gff', sep='\t', index=False, header=False)

    # ENSG filtering
    filtered_human_exons_fantom = filter_by_ensg(matched_human_exons_fantom, matched_human_exons_longread)
    write_file(output_dir + 'filtered_matched_human_exons.gff' , filtered_human_exons_fantom)
    os.remove(processed_human_file)
    os.remove(processed_fantom_file)
    os.remove(processed_longread_file)

if __name__ == "__main__":
    main()