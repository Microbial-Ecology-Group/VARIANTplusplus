import os
import glob
import pandas as pd

def parse_files(file_pattern, output_name):
    # Initialize dictionaries to store data
    reads_info = {}
    aligned_info = {}
    psv_data = {}

    # Match files using the general expression
    files = glob.glob(file_pattern)
    print(f"Found {len(files)} files matching the pattern")

    if not files:
        print("No files found. Please check the file pattern and ensure files are available.")
        return

    for file in files:
        # Extract sample name from file name
        sample_name = os.path.basename(file).replace('.txt', '')
        print(f"Processing file: {file}, Sample name: {sample_name}")

        # Initialize dictionaries to store PSV_group and relative abundance for this sample
        psv_groups = {}
        
        try:
            with open(file, 'r') as f:
                for line in f:
                    # Remove ^M characters and any trailing carriage return
                    line = line.replace('\r', '').replace('^M', '').strip()

                    if line.startswith("#num_reads"):
                        num_reads = int(line.split('\t')[1].strip())
                        reads_info[sample_name] = num_reads
                    elif line.startswith("#num_aligned"):
                        num_aligned = int(line.split('\t')[1].strip())
                        aligned_info[sample_name] = num_aligned
                    elif not line.startswith("#") and line.strip():
                        parts = line.split('\t')
                        if len(parts) == 2:
                            psv_group = parts[0].strip()
                            rel_abund = float(parts[1].strip())
                            if sample_name not in psv_data:
                                psv_data[sample_name] = {}
                            psv_data[sample_name][psv_group] = rel_abund
        except Exception as e:
            print(f"Error processing file {file}: {e}")
            continue

        print(f"Finished processing file: {file}")

    # Create a DataFrame for relative abundance
    print("Creating DataFrame for relative abundance")
    df = pd.DataFrame(psv_data).fillna(0)
    df.index.name = 'Sample'

    # Save the count matrix to a file
    count_matrix_file = f"{output_name}.csv"
    df.to_csv(count_matrix_file)
    print(f"Count matrix saved to {count_matrix_file}")

    # Create a DataFrame for reads and aligned info
    print("Creating DataFrame for reads info")
    reads_df = pd.DataFrame.from_dict({
        'num_reads': reads_info,
        'num_aligned': aligned_info
    }, orient='index')

    # Save the reads info to a file
    reads_info_file = f"reads_{output_name}.csv"
    reads_df.to_csv(reads_info_file)
    print(f"Reads info saved to {reads_info_file}")

if __name__ == "__main__":
    import argparse

    # Setup argument parser
    parser = argparse.ArgumentParser(description='Parse text files and create output matrices.')
    parser.add_argument('file_pattern', type=str, help='Pattern to match text files (e.g., "/path/to/files/*.txt")')
    parser.add_argument('output_name', type=str, help='Base name for the output files')

    # Parse arguments
    args = parser.parse_args()

    # Run the file parsing function
    parse_files(args.file_pattern, args.output_name)
