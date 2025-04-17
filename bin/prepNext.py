#!/usr/bin/env python3
import sys
import os
import pandas as pd
import requests

"""
prepNext.py

This script processes an input GFF file and a GTF file to extract specific exons based on the provided 
identity (five_prime or three_prime). It identifies the first exon for the forward strand or the last 
exon for the reverse strand, depending on the specified identity. The output is written to a new GFF file.

Usage:
    python prepNext.py <input_file> <identity> <gtf_file>

Arguments:
    input_file          Path to the input GFF file containing transcript data.
    identity            Specifies whether to extract five_prime or three_prime exons. Acceptable values:
                        'five', '5', 'three', or '3'.
    gtf_file            Path to the GTF file containing exon information.

Steps:
1. Normalize the `identity` argument to determine whether to process five_prime or three_prime exons.
2. Read the input GFF file and extract transcript IDs, handling `_MANE_COPY` transcripts separately.
3. Parse the GTF file to filter for relevant transcript IDs and extract exon information.
4. Identify the first exon for the forward strand or the last exon for the reverse strand based on the 
   specified identity.
5. Write the selected exons to an output GFF file.

Output:
    - A GFF file containing the selected exons for each transcript.

Dependencies:
    - pandas: For reading and processing tabular data.
    - requests: For potential external requests (not used in the current implementation).
    - sys, os: For command-line argument handling and file operations.

Example:
    python prepNext.py input.gff five gtf_file.gtf
"""
# Get the input file name from the command line arguments
input_file = sys.argv[1]
identity = sys.argv[2]
gtf_file = sys.argv[3]

# Normalize the identity input
base_name, ext = os.path.splitext(input_file)
output_file = f"{identity}_nextRun.gff"
identity = identity.lower()  # Convert to lowercase for case-insensitive matching
if 'five' in identity or '5' in identity:
    identity = 'five_prime'
elif 'three' in identity or '3' in identity:
    identity = 'three_prime'
else:
    raise ValueError("Invalid identity input. Please use 'five', '5', 'three', or '3'.")
# Construct the output file name


# Read the input file
df = pd.read_csv(input_file, sep='\t', usecols=["Chromosome","Source","Type","Start","End","Score","Strand","Phase","Attributes","gene_id","Name","capOrTail_Start","capOrTail_End","Transcript_Start","Transcript_End","Transcript_Name"])
print(df.tail())

# Function to parse attributes from the GTF file
def parse_attributes(attribute_string):
    attributes = {}
    for attribute in attribute_string.split(';'):
        if attribute.strip():
            key, value = attribute.strip().split('=')
            attributes[key] = value.strip('"')
    return attributes
# Extract transcript IDs and handle _MANE_copy
transcript_ids = []
mane_transcripts = {}

for _, row in df.iterrows():
    attributes = parse_attributes(row['Attributes'])
    transcript_id = attributes.get('Parent').split(':')[1] if 'Parent' in attributes else row['Transcript_Name']
    if pd.notna(transcript_id):
        if '_MANE_COPY' in transcript_id.upper():
            base_id = transcript_id.upper().split('_MANE_COPY')[0]
            transcript_ids.append(base_id)
            mane_transcripts[base_id] = transcript_id
            print(mane_transcripts, "MANE")
        else:
            transcript_ids.append(transcript_id)

gene_ids = df.set_index('Transcript_Name')['gene_id'].to_dict()
print(gene_ids)
# Extract transcript IDs from the 'Transcript_Name' column
gtf_columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
gtf_df = pd.read_csv(gtf_file, sep="\t", comment='#', names=gtf_columns)

# Filter the GTF DataFrame for the relevant transcript IDs and "exon" feature
filtered_gtf_df = gtf_df[(gtf_df['attribute'].str.contains('|'.join(transcript_ids))) & (gtf_df['feature'] == 'exon')]
#print(filtered_gtf_df.tail())



# Extract transcript data from the filtered GTF DataFrame
transcript_data = []
for _, row in filtered_gtf_df.iterrows():
    attributes = parse_attributes(row['attribute'])
    transcript_id = attributes.get('Parent').split(':')[1] if 'Parent' in attributes else None
    exon_number = attributes.get('rank')  # Extract exon number from rank
    if transcript_id in transcript_ids:
        if transcript_id in mane_transcripts:
            transcript_id = mane_transcripts[transcript_id]
        gene_id = gene_ids.get(transcript_id)
        transcript_data.append({
            'seq_region_name': row['seqname'],
            'source': row['source'],
            'feature': row['feature'],
            'start': row['start'],
            'end': row['end'],
            'score': row['score'],
            'strand': row['strand'],
            'frame': row['frame'],
            'transcript_id': transcript_id if transcript_id not in mane_transcripts else mane_transcripts[transcript_id],
            'gene_id': gene_id,
            'exon_number': int(exon_number)  # Convert exon_number to integer for sorting
        })

# Sort transcript_data by transcript_id and exon_number
transcript_data.sort(key=lambda x: (x['transcript_id'], x['exon_number']))

# Filter to keep only the first exon for each transcript on the forward strand
# and the last exon for each transcript on the reverse strand
select_exons = {}
for data in transcript_data:
    transcript_id = data['transcript_id']
    strand = data['strand']
    start = int(data['start'])
    end = int(data['end'])
    if transcript_id not in select_exons:
        select_exons[transcript_id] = data
        #print(f"Adding first exon for transcript {transcript_id} on strand {strand}: {data}")
    else:
        if identity == 'five_prime':
            if strand == '+':
                if end > int(select_exons[transcript_id]['end']):
                    select_exons[transcript_id] = data
                    #print(f"Updating first exon for transcript {transcript_id} on strand {strand}: {data}")
            elif strand == '-':
                if start < int(select_exons[transcript_id]['start']):
                    select_exons[transcript_id] = data
        elif identity == 'three_prime':
            if strand == '+':
                if start < int(select_exons[transcript_id]['start']):
                    select_exons[transcript_id] = data
            elif strand == '-':
                if end > int(select_exons[transcript_id]['end']):
                    select_exons[transcript_id] = data
                    #print(f"Updating first exon for transcript {transcript_id} on strand {strand}: {data}")


# Convert the dictionary back to a list
select_exon_data = list(select_exons.values())

# Function to write the results to a GTF file
def write_gtf(transcript_data, output_file):
    with open(output_file, 'w') as f:
        f.write("seqname\tsource\tfeature\tStart\tEnd\tscore\tStrand\tframe\tAttributes\tgene_id\n")
        for data in transcript_data:
            f.write(f"{data['seq_region_name']}\t{data['source']}\t{data['feature']}\t{data['start']}\t{data['end']}\t{data['score']}\t{data['strand']}\t{data['frame']}\texon_number {data['exon_number']};Parent=transcript:{data['transcript_id']}; gene_id={data['gene_id']}\t\"{data['gene_id']}\" \n")

# Write the results to the output GTF file
write_gtf(select_exon_data, output_file)