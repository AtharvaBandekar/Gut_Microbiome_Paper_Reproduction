import pandas as pd
import re

input_csv_file = 'meta_full.csv'
output_tsv_file = 'meta_full.tsv'

try:
    # Read the CSV file, explicitly telling pandas NOT to use the first column as an index.
    # Set header=0 to indicate the first row is the header.
    df = pd.read_csv(input_csv_file, sep=',', header=0, index_col=False) 

    # Write to TSV, ensuring tab separator, no index column, and no quoting.
    df.to_csv(output_tsv_file, sep='\t', index=False, quoting=3) # quoting=3 (csv.QUOTE_NONE)
    print(f"Successfully converted '{input_csv_file}' to '{output_tsv_file}'.")

except Exception as e:
    print(f"An error occurred: {e}")
    print("Please check your input CSV file for malformed lines or hidden characters.")
