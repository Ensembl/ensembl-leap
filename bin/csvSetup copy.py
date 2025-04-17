#!/usr/bin/env python3
import sys
import pandas as pd

# Get the input from the command line argument
input_data = sys.argv[1]

# Split the input string into a list of values
input_values = input_data.split(',')

# Define the number of columns
num_columns = 5

# Reshape the input values into a list of lists, each containing num_columns elements
input_data = [input_values[i:i + num_columns] for i in range(0, len(input_values), num_columns)]

# Define the column names
column_names = ['ID', 'human', 'polyA', 'fantom', 'encode']

# Create a DataFrame from the list with the specified column names
df = pd.DataFrame(input_data, columns=column_names)

# Replace 'na' values with the data from the first row in the respective column
for column in df.columns[1:-1]:  # Skip the first and last columns
    first_row_value = df[column].iloc[0]
    df[column] = df[column].replace('na', first_row_value)

# Define the output CSV file path
output_csv = 'updated_output.csv'

# Write the DataFrame to the CSV file without column names
df.to_csv(output_csv, index=False, header=False)

# Print the output CSV file path
print(output_csv)
