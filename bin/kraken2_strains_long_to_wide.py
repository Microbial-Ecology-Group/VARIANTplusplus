import pandas as pd
import numpy as np
import argparse

def combine_kraken_reports(files):
    # Initialize an empty dictionary to store the taxon counts for each file
    taxa_dict = {}
    
    # Initialize an empty dictionary to store the taxonomic level for each taxon
    taxonomic_levels = {}
    
    # Loop through each file to populate the taxa_dict
    for file_path in files:
        sample_name = file_path.split('/')[-1].split('.')[0]  # Extract the sample name from the file path
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) != 6:
                        print(f"Skipping line in file {file_path} due to incorrect number of fields: {line.strip()}")
                        continue
                    percentage, count, count_of_taxon, taxon, taxaid, name = parts
                    if taxon not in ['U', 'D', 'K', 'P', 'C', 'O', 'F', 'G', 'S', 'S1']: #Skip all the potential super families and other weird taxa levels, include strain "S1"
                        continue
                    taxon_key = name.strip()
                    taxa_dict.setdefault(taxon_key, {})[sample_name] = count_of_taxon  # Use count_of_taxon instead of count
                    taxonomic_levels[taxon_key] = taxon
        except Exception as e:
            print(f"An error occurred while processing file {file_path}: {e}")

    # Convert the taxa_dict to a DataFrame
    combined_df = pd.DataFrame.from_dict(taxa_dict, orient='index')
    
    # Sort the DataFrame by index (Taxon names)
    combined_df.sort_index(inplace=True)
    
    # Fill NaN values with zeros
    combined_df.fillna(0, inplace=True)

    # Convert the taxonomic_levels dictionary to a DataFrame
    if taxonomic_levels:
        taxonomy_df = pd.DataFrame.from_dict(taxonomic_levels, orient='index', columns=['Level'])
        taxonomy_df.index.name = 'Taxa'
    else:
        taxonomy_df = pd.DataFrame(columns=['Level'])  # Empty DataFrame if no data
        print("Warning: No taxonomic levels found.")  # Debug print
    
    return combined_df, taxonomy_df

if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Combine multiple Kraken reports into a single matrix.')
    parser.add_argument('-i', '--input_files', nargs='+', required=True, help='List of Kraken report files to combine')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output file name for the combined matrix')
    args = parser.parse_args()
    
    # Combine Kraken reports into a single DataFrame
    combined_df, taxonomy_df = combine_kraken_reports(args.input_files)
    
    # Write the combined DataFrame to a CSV file
    combined_df.to_csv(args.output_file, index=True)
    
    # Write the taxonomic levels DataFrame to a CSV file
    taxonomy_df.to_csv(f'{args.output_file}_taxonomy.csv', index=True)
