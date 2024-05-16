# AUTHOR: Darío González
# DESCRIPTION: Generate plots for the joined read loss statistics

# Imports
import argparse
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# Parse arguments
parser = argparse.ArgumentParser(description="Generate plots for the joined read loss statistics")

parser.add_argument('table', help="Path to joined read loss table")
parser.add_argument('-m', '--mode', choices=["kb", "star"], help="Alignment program used. One of 'kb' (kallisto-bustools) or 'star'", required=True)
parser.add_argument('-o', '--write_path', help="Path of the folder to write files", required=True)
parser.add_argument('-n', '--filename', help="Name of the files that will be written", required=True)

args = parser.parse_args()

#######
# MAIN
#######

# Read table with read loss data
combined_stats = pd.read_table(args.table, sep="\t")

# Choose table processing according to the columns present
if args.mode == "kb":
    # Check if read stats come from kb-python or bustools
    if 'ambiguous' in list(combined_stats['step']):
        read_selection = 'count'
    else:
        read_selection = 'capture'

# GENERATE MEANS PER STEP TABLE

# Select columns that will be used to create the summary table, and remove the rest
cols_keep = ['step', 'percentage']
cols_drop = [col for col in combined_stats.columns if col not in cols_keep]

step_mean = combined_stats.drop(cols_drop, axis=1).groupby('step').mean()

# Reorder rows after means calculation (pandas returns an unordered dataframe)
if args.mode == "star":
    row_order = ['original', 'unique_aligned', 'corrected', 'spliced', 'unspliced', 'ambiguous', 'spliced+unspliced+ambiguous']
else:
    if read_selection == 'capture':
        row_order = ['original', 'pseudoaligned', 'corrected', 'spliced', 'unspliced', 'spliced+unspliced']
    elif read_selection  == 'count':
        row_order = ['original', 'pseudoaligned', 'corrected', 'spliced', 'unspliced', 'ambiguous', 'spliced+unspliced+ambiguous']

step_mean = step_mean.reindex(row_order).reset_index()

# Save means per step table
write_name = args.filename + "_step_means" + ".tsv"
step_mean.to_csv(args.write_path + write_name, sep = "\t", index = False)

########
# PLOTS
########

# Select columns to plot
if args.mode == "star":
    rows_to_plot = ['original', 'unique_aligned', 'corrected', 'spliced+unspliced+ambiguous']
else:
    if read_selection == 'capture':
        rows_to_plot = ['original', 'pseudoaligned', 'corrected', 'spliced+unspliced']
    elif read_selection == 'count':
        rows_to_plot = ['original', 'pseudoaligned', 'corrected', 'spliced+unspliced+ambiguous']

## Plot read loss per sample
# Create new dataframe with selected columns to plot
rows_drop = [ step not in rows_to_plot for step in combined_stats['step'] ]

combined_stats_plot = combined_stats.drop(combined_stats[rows_drop].index, axis=0).reset_index(drop=True)

# Generate plot with plotly
fig1 = px.line(data_frame=combined_stats_plot, 
               x="step", 
               y="percentage", 
               line_group="run", 
               color="run", 
               color_discrete_sequence=["cornflowerBlue"], 
               markers=True, 
               hover_data = ["bundle", "run", "percentage"], 
               width=1200)
fig1.update_traces(patch = {'marker': {'size': 6, 'color': 'blue'}})
fig1.update_yaxes(range=[0, None], dtick=10, minor={'showgrid': True, 'dtick': 2})

# Save plot
write_name = args.filename + ".html"
fig1.write_html(args.write_path + write_name)

## Plot read loss means per step
# Create new dataframe with selected columns to plot
rows_drop = [ step not in rows_to_plot for step in step_mean['step'] ]

step_mean_plot = step_mean.drop(step_mean[rows_drop].index, axis=0).reset_index(drop=True)

# Generate plot with plotly
fig2 = px.scatter(data_frame = combined_stats_plot, 
                  x="step", 
                  y="percentage", 
                  hover_data=["bundle", "run", "percentage"], 
                  width=1200)
fig2.update_traces(patch = {'marker': {'size': 6, 'color': 'blue'}})
fig2.add_trace(go.Bar(x=step_mean_plot['step'], 
                      y=step_mean_plot['percentage'],
                      marker_color="cornflowerBlue",
                      text=[round(x, 2) for x in step_mean_plot['percentage']], 
                      textposition='inside', 
                      textfont=dict(family="Arial", size=15, color="black"),
                      width=0.3,
                      name="mean"
                      )
                )
fig2.update_yaxes(tick0=0, dtick=10, minor={'showgrid': True, 'dtick': 2})

# Save means per step plot
write_name = args.filename + "_step_means" + ".html"
fig2.write_html(args.write_path + write_name)

# Save means per step table used for plotting, plus decrease in each step
# Add percentage decrease
decrease = [step_mean_plot.loc[i, 'percentage'] - step_mean_plot.loc[i-1, 'percentage'] for i in range(1, len(step_mean_plot))]
decrease.insert(0, 0.0)
step_mean_plot['percent_decrease'] = decrease

# Save file
write_name = args.filename + "_step_means_decrease" + ".tsv"
step_mean_plot.to_csv(args.write_path + write_name, sep = "\t", index = False)
