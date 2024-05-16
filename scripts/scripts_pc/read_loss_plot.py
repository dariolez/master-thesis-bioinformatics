# AUTHOR: Darío González
# DESCRIPTION: Generate plots of read loss per sample

# Imports
import argparse
import os
import pandas as pd
import plotly.graph_objects as go

# Parse arguments
parser = argparse.ArgumentParser(description="Generate plots of read loss per sample")

parser.add_argument('folders_fofn', help="File with the sample folders to be processed")
parser.add_argument('path_to_read_loss', help="Path to the reads loss tsv within the run folder")
parser.add_argument('write_name', help="Name of the plot file to be written")

args = parser.parse_args()

# Get run folders
with open(args.folders_fofn, 'r') as f:
    run_folders = f.readlines()
    run_folders = [run.rstrip("\n") for run in run_folders]

for folder in run_folders:
    # Input reads loss file
    file_path = folder + "/" + args.path_to_read_loss
    stats = pd.read_table(file_path, sep="\t")

    # Remove some columns from plotting
    if 'unique_aligned' in list(stats['step']):
        stats_plot = stats.loc[stats['step'].isin(['original', 'unique_aligned', 'corrected', 'spliced+unspliced+ambiguous'])]
    else:
        if 'ambiguous' in stats['step']:
            stats_plot = stats.loc[stats['step'].isin(['original', 'pseudoaligned', 'corrected', 'spliced+unspliced+ambiguous'])]
        else:
            stats_plot = stats.loc[stats['step'].isin(['original', 'pseudoaligned', 'corrected', 'spliced+unspliced'])]

    # Generate plot
    fig = go.Figure([go.Bar(x=stats_plot['step'], 
                            y=stats_plot['percentage'],
                            marker_color="cornflowerBlue",
                            text=[round(x, 2) for x in stats_plot['percentage']], 
                            textposition='inside', 
                            textfont=dict(family="Arial", size=15, color="black"),
                            width=0.3)]
                    )
    fig.update_yaxes(tick0=0, dtick=10)

    # Save plot
    write_path = folder + "/" + "/".join((args.path_to_read_loss).split("/")[:-1]) + "/" + args.write_name
    fig.write_html(write_path)
