#!/usr/bin/env python3
import re
import sys
import pandas as pd

"""
humanFilter.py

This script filters protein-coding genes from a GFF file based on specific criteria. It excludes genes 
that are marked as readthrough transcripts and optionally filters out single-exon genes. The filtered 
output is written to a new GFF file.

Usage:
    python humanFilter.py <input_file> <output_file> <readthrough_file> <single_exon>

Arguments:
    input_file          Path to the input GFF file.
    output_file         Path to the output filtered GFF file.
    readthrough_file    Path to a file containing a list of readthrough transcript stable IDs.
    single_exon         Boolean flag ('true' or 'false') indicating whether to include single-exon genes.

Steps:
1. Load the list of readthrough transcript stable IDs from the specified file.
2. Parse the input GFF file block by block (genes are separated by "###").
3. For each gene block:
   - Include the block if it contains a protein-coding gene and does not match any readthrough IDs.
   - Optionally exclude single-exon genes based on the `single_exon` flag.
4. Write the filtered gene blocks to the output file.

Output:
    - A filtered GFF file containing only the desired protein-coding genes.

Dependencies:
    - pandas: For reading the readthrough transcript list.
    - re: For regular expression matching.
    - sys: For command-line argument handling.

Example:
    python humanFilter.py input.gff output.gff readthrough.txt true
"""

def load_readthrough_list(readthrough_file):
    readthrough_df = pd.read_csv(readthrough_file, sep='\t')
    return set(readthrough_df['stable_id'])

def filter_protein_coding_genes(input_file, output_file, readthrough_file, single_exon):
    readthrough_ids = load_readthrough_list(readthrough_file)

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        gene_block = []
        keep_block = False
        if single_exon == True:
            exon_count = 0
        else:
            exon_count = 1 

        for line in infile:
            if line.startswith("###"):
                if keep_block and exon_count > 1:
                    outfile.write("".join(gene_block))
                    outfile.write(line)
                gene_block = []
                keep_block = False
                exon_count = 0
                
            else:
                gene_block.append(line)
                if "biotype=protein_coding" in line:
                    keep_block = True
                if any(str(readthrough_id) in line for readthrough_id in readthrough_ids):
                    keep_block = False
                if "\texon\t" in line:
                    exon_count += 1
                    

        # Write the last block if it should be kept
        if keep_block and exon_count > 1:
            outfile.write("".join(gene_block))

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python proteinCodingGeneFilterGff.py <input_file> <output_file> <readthrough_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    readthrough_file = sys.argv[3]
    single_exon = sys.argv[4].lower() == 'true'

    if single_exon:
        print("Single exon is True")
    else:
        print("Single exon is False")
    filter_protein_coding_genes(input_file, output_file, readthrough_file, single_exon)