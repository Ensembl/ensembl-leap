#!/usr/bin/env python3
import sys
import pandas as pd

"""
splitChromosomes.py

This script processes input files containing genomic data and splits them by a specified chromosome. 
It supports multiple file formats, including GFF, GTF, GFF3, BED, and tab-delimited text files. The 
filtered data for the specified chromosome is written to separate output files for each input file.

Usage:
    python splitChromosomes.py <chr> <fantom> <longRead> <capOrTail> <human> <capOrTail_type>

Arguments:
    chr                 The chromosome to filter (e.g., "1", "X", "MT").
    fantom              Path to the FANTOM input file.
    longRead            Path to the long-read input file.
    capOrTail           Path to the cap or tail input file.
    human               Path to the human input file.
    capOrTail_type      Type of cap or tail data (not used in the current implementation).

Steps:
1. Parse the input files based on their format (GFF, GTF, BED, or tab-delimited).
2. Filter the data for the specified chromosome.
3. Convert BED files to GFF format if necessary.
4. Write the filtered data to output files with a prefix indicating the input file type.

Output:
    - Separate output files for each input file, containing only the data for the specified chromosome.

Dependencies:
    - pandas: For reading and processing tabular data.
    - sys: For command-line argument handling.

Example:
    python splitChromosomes.py 1 fantom.gff longRead.bed capOrTail.txt human.gtf cap
"""

def convert_to_gff(df):
    # Convert BAM-like DataFrame to GFF format
    gff_df = pd.DataFrame()
    gff_df['seqname'] = df[0]  # Chromosome
    gff_df['source'] = 'BAM'  # Source
    gff_df['feature'] = 'region'  # Feature type
    gff_df['start'] = df[1]  # Start position
    gff_df['end'] = df[2]  # End position
    gff_df['score'] = df[4]  # Score
    gff_df['strand'] = df[5]  # Strand
    gff_df['frame'] = '.'  # Frame (not applicable)
    gff_df['attribute'] = df[9]  # Attributes (e.g., TE or other info)
    return gff_df

def split_file(input_file, output_prefix, chr):
    file_extension = input_file.split('.')[-1]
    print(file_extension)
    
    if file_extension in ['gff', 'gtf', 'gff3']:
        df = pd.read_csv(input_file, sep='\t', comment='#', header=None, dtype=str)
        chr_col = 0
    elif file_extension == 'bed':
        df = pd.read_csv(input_file, sep='\t', header=None, dtype=str)
        print(df.head())
        chr_col = 0
    else:
        df = pd.read_csv(input_file, sep='\t', header=0, dtype=str)
        chr_col = df.columns[0]
    
    # Check if the chromosome column contains "chr" prefix
    if df[chr_col].astype(str).str.startswith('chr').any():
        chr_value = f"chr{chr}"
    else:
        chr_value = chr
    
    filtered_df = df[df[chr_col] == chr_value]
    if not filtered_df.empty:
        if file_extension == 'bed':
            filtered_df = convert_to_gff(filtered_df)
        output_file = f"{output_prefix}_{chr}.txt"
        filtered_df.to_csv(output_file, sep='\t', index=False, header=(file_extension not in ['gff', 'gtf', 'gff3', 'bed']))
        print(f"Written to {output_file}")
    else:
        print(f"No data found for chromosome {chr} in {input_file}")

def main(chr, fantom, longRead, capOrTail, human, capOrTail_type):
    # Process FANTOM input file
    split_file(fantom, "split_fantom", chr)

    # Process longRead input file
    split_file(longRead, "split_longRead", chr)

    # Process capOrTail input file
    # Process capOrTail input file
    split_file(capOrTail, "split_capOrTail", chr)

    # Process HUMAN input file
    split_file(human, "split_human", chr)

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python split_chromosomes.py <chr> <fantom> <longRead> <capOrTail> <human> <capOrTail_type>")
        sys.exit(1)

    chr = sys.argv[1]
    fantom = sys.argv[2]
    longRead = sys.argv[3]
    capOrTail = sys.argv[4]
    human = sys.argv[5]
    capOrTail_type = sys.argv[6]

    main(chr, fantom, longRead, capOrTail, human, capOrTail_type)


