import pandas as pd
import numpy as np
import argparse

def combine_kraken_reports(files):
    taxa_dict = {}
    taxonomic_levels = {}
    
    for file_path in files:
        sample_name = file_path.split('/')[-1].split('.')[0]
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) != 6:
                        print(f"Skipping line in file {file_path} due to incorrect number of fields: {line.strip()}")
                        continue
                    percentage, count, count_of_taxon, taxon, taxaid, name = parts
                    count_of_taxon = int(count_of_taxon)
                    taxon_key = name.strip()
                    taxa_dict.setdefault(taxon_key, {})[sample_name] = count_of_taxon
                    taxonomic_levels[taxon_key] = taxon
        except Exception as e:
            print(f"An error occurred while processing file {file_path}: {e}")

    combined_df = pd.DataFrame.from_dict(taxa_dict, orient='index')
    combined_df.sort_index(inplace=True)
    combined_df.fillna(0, inplace=True)
    combined_df = combined_df.loc[(combined_df != 0).any(axis=1)]
    
    taxonomy_df = pd.DataFrame.from_dict(taxonomic_levels, orient='index', columns=['Level']) if taxonomic_levels else pd.DataFrame(columns=['Level'])
    
    return combined_df, taxonomy_df

def calculate_summary_statistics(combined_df, taxonomy_df):
    summary_stats = {}
    for sample in combined_df.columns:
        total_reads = combined_df[sample].sum()
        unclassified_reads = combined_df.loc['unclassified', sample] if 'unclassified' in combined_df.index else 0
        mh_species_df = combined_df[combined_df.index.str.contains('Mannheimia haemolytica')]
        mh_species_reads = mh_species_df[sample].sum()
        mh_strain_indexes = taxonomy_df[(taxonomy_df.index.str.contains(r'Mannheimia[_\s]haemolytica')) & (taxonomy_df['Level'].isin(['S1', 'S2']))].index 
        mh_strain_indexes = mh_strain_indexes.intersection(combined_df.index)  # Only keep indexes that exist in combined_df
        mh_strain_reads = combined_df.loc[mh_strain_indexes, sample].sum() if not mh_strain_indexes.empty else 0
        number_of_unique_mh_strains = len(mh_strain_indexes)
        pasteurellaceae_reads = combined_df[combined_df.index.str.contains('Pasteurellaceae')][sample].sum()
        
        summary_stats[sample] = {
            'Number of Unclassified Reads': unclassified_reads,
            'Percent Unclassified Reads': (unclassified_reads / total_reads) * 100,
            'Percent Reads Classified to Mannheimia haemolytica': (mh_species_reads / total_reads) * 100,
            'Number of Strains Classified under Mannheimia haemolytica': number_of_unique_mh_strains,
            'Sum of Hits to Mannheimia haemolytica Strains': mh_strain_reads,
            'Percent Reads Classified to Pasteurellaceae': (pasteurellaceae_reads / total_reads) * 100,
            'Number of Reads Classified to Pasteurellaceae': pasteurellaceae_reads
        }
    
    summary_df = pd.DataFrame.from_dict(summary_stats, orient='index')
    return summary_df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combine multiple Kraken reports into a single matrix and compute summary statistics.')
    parser.add_argument('-i', '--input_files', nargs='+', required=True, help='List of Kraken report files to combine')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output file name for the combined matrix')
    args = parser.parse_args()
    
    combined_df, taxonomy_df = combine_kraken_reports(args.input_files)
    combined_df.to_csv(args.output_file, index=True)
    taxonomy_df.to_csv(f'taxonomy_{args.output_file}', index=True)
    
    summary_df = calculate_summary_statistics(combined_df, taxonomy_df)
    summary_df.to_csv(f'summary_{args.output_file}', index=True)
    
    dataset_summary_df = summary_df.describe().transpose()[['mean', 'min', '25%', '50%', '75%', 'max', 'std']]
    dataset_summary_df.to_csv(f'dataset_summary_{args.output_file}', index=True)
