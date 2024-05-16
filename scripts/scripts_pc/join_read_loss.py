# AUTHOR: Darío González
# DESCRIPTION: Combine all read loss files into one table and save it

# Imports
import argparse
import os
import numpy as np
import pandas as pd

# Parse arguments
parser = argparse.ArgumentParser(description="Combine all read loss files into one table and save it")

parser.add_argument('-m', '--mode', type=int, choices=[1, 2], default=2, help="Save different tables depending on the mode")
parser.add_argument('--wd', help="Set script working directory")
parser.add_argument('runs_fofn', help="File with the absolute path to runs folders to be processed")
parser.add_argument('path_to_read_loss', help="Path to the reads loss tsv within the run folder")
parser.add_argument('write_file', help="Path of the file to be named")

args = parser.parse_args()

# Set working directory
os.chdir(args.wd)

# Get runs directories
with open(args.runs_fofn, 'r') as f:
    run_folders = f.readlines()
    run_folders = [run.rstrip("\n") for run in run_folders]

# Get file mapping bundle folder names to FASTQ name
bundle2filename = pd.read_csv("bundle_filename_10x.tsv", sep="\t")

# Combine all tables
counter = 0

for folder in run_folders:
    counter += 1

    # Get run name
    run = folder.split("/")[-1]

    # Input reads stats file
    file_path = folder + args.path_to_read_loss
    stats = pd.read_table(file_path, sep="\t")

    # Get the file name
    filename = bundle2filename.loc[bundle2filename['bundle_uuid'] == run, 'file_name'].values[0]

    if args.mode == 1:
        ## OPTION 1: Runs as columns names, append columns
        
        if counter == 1:
            new_cols = [ filename + "_" + colname for colname in stats.columns[1:] ]
            stats.rename(columns=new_cols)
            combined_stats = stats
        else:
            stats = stats.add_prefix(filename + "_")
            combined_stats = pd.DataFrame.join(combined_stats, stats)

    if args.mode == 2:
        ## OPTION 2: Runs added in separate column, append rows
        stats.insert(1, 'run', filename)
        stats.insert(1, 'bundle', run)

        # Join dataframes
        if counter == 1:
            combined_stats = stats
        else:
            combined_stats = pd.concat([combined_stats, stats], axis = 0)

# Save table
combined_stats.to_csv(args.write_file, sep="\t", index=False)