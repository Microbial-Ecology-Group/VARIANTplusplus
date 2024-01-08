import pandas as pd
import sys
import os
import glob

# Function: Summarize results from multiple files created by TreeCluster.py. 
# Example command: 
# python parse_ClusterTree_results.py '*.txt' output_results.csv

# This bash function can be run seperately to loop through multiple thresholds and all clustering methods

#!/bin/bash
# Define the list of threshold values, List of methods, and loop through each of them
#thresholds=("0.0001" "0.0002" "0.0003" "0.0004" "0.0005" "0.0006" "0.0007" "0.0008" "0.0009" "0.001" "0.002" "0.003" "0.004" "0.005" "0.006" "0.007" "0.008" "0.009" "0.01")
#methods=("avg_clade" "leaf_dist_avg" "leaf_dist_max" "leaf_dist_min" "length" "length_clade" "max" "max_clade" "med_clade" "root_dist" "single_linkage" "single_linkage_cut" "single_linkage_union" "sum_branch" "sum_branch_clade")
#for method in "${methods[@]}"; do
#    for threshold in "${thresholds[@]}"; do
#        output_file="output_TreeCluster_labels_${method}_t_${threshold}.txt"
#        command="TreeCluster.py -i Final_tree_min0.dnd -o Benchmarking/$output_file -t $threshold -m $method"
#        echo "Running: $command"
#        eval $command
#    done
#done

def analyze_clusters(file_path):
    try:
        # Read the data from the file
        df = pd.read_csv(file_path, sep="\t")
        
        # Calculate the total number of unique clusters
        total_clusters = df['ClusterNumber'].nunique()
        
        # Calculate the number of genomes labeled as "-1"
        singletons = (df['ClusterNumber'] == -1).sum()
        
        # Calculate the value for the cluster with the highest number of genomes in it
        max_cluster_size = df['ClusterNumber'].value_counts().max()
        
        # Calculate the average genomes per cluster
        avg_genomes_per_cluster = df['ClusterNumber'].value_counts().mean()
        
        return {
            'FileName': os.path.basename(file_path),
            'TotalUniqueClusters': total_clusters,
            'Singletons': singletons,
            'MaxClusterSize': max_cluster_size,
            'AvgGenomesPerCluster': avg_genomes_per_cluster
        }
    except Exception as e:
        print(f"Error processing file {file_path}: {str(e)}")
        return None

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py '<file_pattern>' <output_file>")
        sys.exit(1)
    
    file_pattern = sys.argv[1]
    output_file = sys.argv[2]

    results = []
    for file_path in glob.glob(file_pattern):
        result = analyze_clusters(file_path)
        if result is not None:
            results.append(result)
    
    if results:
        results_df = pd.DataFrame(results)
        results_df.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")
    else:
        print("No valid results to save.")