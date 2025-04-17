import pandas as pd
import argparse

"""
makeGFF.py

This script processes four input CSV files and a reference GFF file to generate an extended GFF file. 
It merges the input data, consolidates overlapping columns, and appends additional tags to the GFF 
attributes based on specific conditions. The resulting GFF file includes extended transcript information 
with updated start and end positions, as well as standardized tags.

Usage:
    python makeGFF.py <file1> <file2> <file3> <file4> <reference_gff> <output_gff> <final_merge_file>

Arguments:
    file1               Path to the first input CSV file.
    file2               Path to the second input CSV file.
    file3               Path to the third input CSV file.
    file4               Path to the fourth input CSV file.
    reference_gff       Path to the reference GFF file.
    output_gff          Path to the output GFF file.
    final_merge_file    Path to save the final merged dataframe.

Steps:
1. Load the four input CSV files and merge them based on transcript IDs.
2. Consolidate overlapping columns in the merged dataframe to create a unified dataset.
3. Load the reference GFF file and process its rows to:
   - Extend transcript start and end positions based on CAGE and PolyA data.
   - Add standardized tags (`gencode_primary`, `gencode_basic`) to all GFF attributes.
   - Add the `MANE_copy` tag to transcripts marked as `MANE_Select`.
4. Write the processed GFF data to the output file.

Output:
    - A GFF file with extended transcript information and updated attributes.
    - A merged CSV file for inspection of the consolidated input data.

Dependencies:
    - pandas: For data manipulation.
    - argparse: For parsing command-line arguments.

Example:
    python makeGFF.py input1.csv input2.csv input3.csv input4.csv reference.gff output.gff merged.csv
"""

# Set up argument parser
parser = argparse.ArgumentParser(description="Process four input files and a reference GFF to generate an extended GFF file.")
parser.add_argument("file1", type=str, help="Path to the first (fiveprime1) input CSV file")
parser.add_argument("file2", type=str, help="Path to the second (threeprime1) input CSV file")
parser.add_argument("file3", type=str, help="Path to the third (fiveprime2) input CSV file")
parser.add_argument("file4", type=str, help="Path to the fourth (threeprime2) input CSV file")
parser.add_argument("reference_gff", type=str, help="Path to the reference GFF file")
parser.add_argument("output_gff", type=str, help="Path to the output GFF file")
parser.add_argument("final_merge_file", type=str, help="Path to save the final merged dataframe")
args = parser.parse_args()

# Load the four CSV files
file1 = args.file1
file2 = args.file2
file3 = args.file3
file4 = args.file4
reference_gff = args.reference_gff
output_gff = args.output_gff
final_merge_file = args.final_merge_file

df1 = pd.read_csv(file1, sep="\t")
df2 = pd.read_csv(file2, sep="\t")
df3 = pd.read_csv(file3, sep="\t")
df4 = pd.read_csv(file4, sep="\t")

# Merge the first two files on transcript ID
merged_df1 = pd.merge(df1, df2, left_on="Transcript_Name", right_on="Transcript_Name", how="outer", suffixes=("_fiveprime1", "_threeprime1"))

# Merge the second two files on transcript ID
merged_df2 = pd.merge(df3, df4, left_on="Transcript_Name", right_on="Transcript_Name", how="outer", suffixes=("_fiveprime2", "_threeprime2"))

# Merge the two resulting dataframes on transcript ID
final_merged_df = pd.merge(merged_df1, merged_df2, on="Transcript_Name", how="outer", suffixes=("_df1", "_df2"))
# Consolidate columns with similar names
def consolidate_columns(df, base_columns, exclude_columns=None):
    if exclude_columns is None:
        exclude_columns = []

    for base_col in base_columns:
        # Skip columns that are in the exclude list
        if base_col in exclude_columns:
            continue

        # Find all columns that match the base column pattern
        matching_cols = [col for col in df.columns if col.startswith(base_col)]
        if matching_cols:
            # Consolidate values by taking the first non-null value across matching columns
            df[base_col] = df[matching_cols].bfill(axis=1).iloc[:, 0]
            # Drop the original matching columns
            df.drop(columns=matching_cols, inplace=True)
    return df

# List of base column names to consolidate

base_columns = [
    "Chromosome", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes",
    "gene_id", "Name", "capOrTail_Start", "capOrTail_End", "Transcript_Start", "Transcript_End",
    "Transcript_Name", "Difference"
]

# Consolidate the columns in the final merged dataframe
# Columns to exclude from consolidation
exclude_columns = ["capOrTail_Start", "capOrTail_End"]
final_merged_df = consolidate_columns(final_merged_df, base_columns, exclude_columns)

# Save the final merged dataframe to a file for inspection
final_merged_df.to_csv(final_merge_file, sep="\t", index=False)

# Load the reference GFF into a DataFrame
gff_columns = [
    "Chromosome", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes"
]
gff_df = pd.read_csv(reference_gff, sep="\t", comment="#", names=gff_columns, header=None)

# Create a list to store the output GFF lines
output_lines = []
numberOfTranscripts = 0

# Open the output GFF file in write mode initially to add the header
with open(output_gff, "w") as gff_file:
    gff_file.write("##gff-version 3\n")

def process_transcripts(merged_df, gff_df):
    output_lines = []
    

    # Pre-index the GFF DataFrame by transcript name for faster lookups
    gff_df["Transcript_Name"] = gff_df["Attributes"].str.extract(r"transcript:([^;]+)")
    gff_index = gff_df.set_index("Transcript_Name")

    def process_row(row):
        numberOfTranscripts = 0
        output = []
        transcript_name = row["Name"] if not pd.isna(row["Name"]) else row["Name"]
        is_mane_copy = False

        # Check if the transcript was originally a MANE_copy
        if "MANE_copy" in transcript_name:
            is_mane_copy = True
            transcript_name = transcript_name.rstrip("_MANE_copy")


        chromosome = row["Chromosome"]

        # Safely convert chromosome to an integer if possible
        try:
            chromosome = int(float(chromosome))  # Handles cases where it's a float
        except (ValueError, TypeError):
            pass  # Leave it as is if it can't be converted

        strand = row["Strand"]

        # Retrieve GFF rows for the transcript
        if transcript_name in gff_index.index:
            transcript_gff = gff_index.loc[[transcript_name]]
        else:
            print(transcript_name, "not found in GFF")
            return output  # Skip if no matching GFF rows

        # Get the minimum and maximum Start and End from the GFF file
        gff_start = transcript_gff["Start"].min()
        gff_end = transcript_gff["End"].max()

        # Determine cage and polya values
        # Determine capOrTail (CAGE and PolyA) values for fiveprime and threeprime
        if strand == "+":
            # For fiveprime (CAGE)
            cage_start = (
                int(row["capOrTail_Start_fiveprime1"])
                if not pd.isna(row["capOrTail_Start_fiveprime1"])
                else (
                    int(row["capOrTail_Start_fiveprime2"])
                    if not pd.isna(row["capOrTail_Start_fiveprime2"])
                    else gff_start
                )
            )
            cage_end = (
                int(row["capOrTail_End_fiveprime1"])
                if not pd.isna(row["capOrTail_End_fiveprime1"])
                else (
                    int(row["capOrTail_End_fiveprime2"])
                    if not pd.isna(row["capOrTail_End_fiveprime2"])
                    else gff_start
                )
            )

            # For threeprime (PolyA)
            polya_start = (
                int(row["capOrTail_Start_threeprime1"])
                if not pd.isna(row["capOrTail_Start_threeprime1"])
                else (
                    int(row["capOrTail_Start_threeprime2"])
                    if not pd.isna(row["capOrTail_Start_threeprime2"])
                    else gff_end
                )
            )
            polya_end = (
                int(row["capOrTail_End_threeprime1"])
                if not pd.isna(row["capOrTail_End_threeprime1"])
                else (
                    int(row["capOrTail_End_threeprime2"])
                    if not pd.isna(row["capOrTail_End_threeprime2"])
                    else gff_end
                )
            )
        elif strand == "-":
            # For fiveprime (CAGE)
            cage_start = (
                int(row["capOrTail_Start_fiveprime1"])
                if not pd.isna(row["capOrTail_Start_fiveprime1"])
                else (
                    int(row["capOrTail_Start_fiveprime2"])
                    if not pd.isna(row["capOrTail_Start_fiveprime2"])
                    else gff_end
                )
            )
            cage_end = (
                int(row["capOrTail_End_fiveprime1"])
                if not pd.isna(row["capOrTail_End_fiveprime1"])
                else (
                    int(row["capOrTail_End_fiveprime2"])
                    if not pd.isna(row["capOrTail_End_fiveprime2"])
                    else gff_end
                )
            )

            # For threeprime (PolyA)
            polya_start = (
                int(row["capOrTail_Start_threeprime1"])
                if not pd.isna(row["capOrTail_Start_threeprime1"])
                else (
                    int(row["capOrTail_Start_threeprime2"])
                    if not pd.isna(row["capOrTail_Start_threeprime2"])
                    else gff_start
                )
            )
            polya_end = (
                int(row["capOrTail_End_threeprime1"])
                if not pd.isna(row["capOrTail_End_threeprime1"])
                else (
                    int(row["capOrTail_End_threeprime2"])
                    if not pd.isna(row["capOrTail_End_threeprime2"])
                    else gff_start
                )
            )

        # Determine extended start and end
        if strand == "+":
            extended_start = cage_start
            extended_end = polya_end
        elif strand == "-":
            extended_start = polya_start
            extended_end = cage_end

        # Process GFF rows
        for _, gff_row in transcript_gff.iterrows():
            gff_type = gff_row["Type"]
            gff_start = gff_row["Start"]
            gff_end = gff_row["End"]
            gff_attributes = gff_row["Attributes"]

            # Add MANE_copy tag if applicable
            if "tag=" in gff_attributes:
                tags = gff_attributes.split("tag=")[1].split(";")[0].split(",")
                if gff_type == "mRNA" and "MANE_Select" in tags:
                    if "MANE_copy" not in tags:
                        tags.append("MANE_copy")
                    if "LEAP" not in tags:
                        tags.append("LEAP")
                if "gencode_primary" not in tags:
                    tags.extend(["gencode_primary", "gencode_basic"])
                if "LEAP" not in tags:
                    tags.append("LEAP")
                gff_attributes = gff_attributes.split("tag=")[0] + "tag=" + ",".join(tags) + ";" + ";".join(gff_attributes.split("tag=")[1].split(";")[1:])

            if gff_type == "mRNA":
                output.append(
                    f"{chromosome}\t{gff_row['Source']}\tmRNA\t{extended_start}\t{extended_end}\t.\t{strand}\t.\t{gff_attributes}\n"
                )
            elif gff_type == "five_prime_UTR":
                if strand == "+":
                    if gff_start == transcript_gff["Start"].min():
                        output.append(
                            f"{chromosome}\t{gff_row['Source']}\tfive_prime_UTR\t{extended_start}\t{gff_end}\t.\t{strand}\t.\t{gff_attributes}\n"
                        )
                    else:
                        output.append(
                            f"{chromosome}\t{gff_row['Source']}\tfive_prime_UTR\t{gff_start}\t{gff_end}\t.\t{strand}\t.\t{gff_attributes}\n"
                        )
                else:
                    if gff_end == transcript_gff["End"].max():
                        output.append(
                            f"{chromosome}\t{gff_row['Source']}\tfive_prime_UTR\t{gff_start}\t{extended_end}\t.\t{strand}\t.\t{gff_attributes}\n"
                        )
                    else:
                        output.append(
                            f"{chromosome}\t{gff_row['Source']}\tfive_prime_UTR\t{gff_start}\t{gff_end}\t.\t{strand}\t.\t{gff_attributes}\n"
                        )
            elif gff_type == "three_prime_UTR":
                if strand == "+":
                    if gff_end == transcript_gff["End"].max():
                        output.append(
                            f"{chromosome}\t{gff_row['Source']}\tthree_prime_UTR\t{gff_start}\t{extended_end}\t.\t{strand}\t.\t{gff_attributes}\n"
                        )
                    else:
                        output.append(
                            f"{chromosome}\t{gff_row['Source']}\tthree_prime_UTR\t{gff_start}\t{gff_end}\t.\t{strand}\t.\t{gff_attributes}\n"
                        )
                else:
                    if gff_start == transcript_gff["Start"].min():
                        output.append(
                            f"{chromosome}\t{gff_row['Source']}\tthree_prime_UTR\t{extended_start}\t{gff_end}\t.\t{strand}\t.\t{gff_attributes}\n"
                        )
                    else:
                        output.append(
                            f"{chromosome}\t{gff_row['Source']}\tthree_prime_UTR\t{gff_start}\t{gff_end}\t.\t{strand}\t.\t{gff_attributes}\n"
                        )
            elif gff_type == "exon":
                if strand == "+":
                    if gff_start == transcript_gff["Start"].min():
                        gff_start = extended_start
                    if gff_end == transcript_gff["End"].max():
                        gff_end = extended_end
                elif strand == "-":
                    if gff_end == transcript_gff["End"].max():
                        gff_end = extended_end
                    if gff_start == transcript_gff["Start"].min():
                        gff_start = extended_start

                output.append(
                    f"{chromosome}\t{gff_row['Source']}\texon\t{gff_start}\t{gff_end}\t.\t{strand}\t.\t{gff_attributes}\n"
                )
            else:
                output.append(
                    f"{chromosome}\t{gff_row['Source']}\t{gff_type}\t{gff_start}\t{gff_end}\t.\t{strand}\t.\t{gff_attributes}\n"
                )

        numberOfTranscripts += 1
        print("processed transcript: ", numberOfTranscripts)
        return output

    # Use apply to process rows in a vectorized manner
    merged_df.apply(lambda row: output_lines.extend(process_row(row)), axis=1)

    return output_lines

if __name__ == "__main__":
    # Process transcripts and generate output lines
    output_lines = process_transcripts(final_merged_df, gff_df)

    # Write all lines to the output file at once
    with open(output_gff, "a") as gff_file:
        gff_file.writelines(output_lines)

    print(f"Extended GFF file created: {output_gff}")