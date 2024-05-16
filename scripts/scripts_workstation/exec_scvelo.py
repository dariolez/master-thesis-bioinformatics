"""
DESCRIPTION: perform RNA velocity and save plots
AUTHOR: Darío González
"""

# Imports
import argparse
import os
import scanpy as sp
import scvelo as scv

# Parse arguments
parser = argparse.ArgumentParser(description="Perform RNA velocity and save plots")

parser.add_argument('-m', '--mode', help="Scvelo modes to use and save plots for", action='store', nargs='+', choices=['stochastic', 'deterministic', 'dynamical'], default=['stochastic', 'deterministic', 'dynamical'])
parser.add_argument('-d', '--directory', help="Folder into which to store the results. The script generates a folder for each mode in here.", action='store', default=".")
parser.add_argument('-w', '--write_name', help="Base name of the plots files to be written", action='store', required=True)
parser.add_argument('-f', '--format', help="Formats in which to save the plots. One or more of [svg,png,pdf]", action='store', nargs='+', choices=['svg', 'png', 'pdf'], default=['svg', 'png'])
parser.add_argument('-a', '--velocity', help="Perform velocity", action='store_true')
parser.add_argument('-u', '--umap', help="Calculate UMAP", action='store_true')
parser.add_argument('-c', '--confidence', help="Save confidence plots", action='store_true')
parser.add_argument('-r', '--rank', help="Save table with velocity gene rank per cluster group", action='store_true')
parser.add_argument('-s', '--save', help="File name to save anndata object with velocity data", action='store', default=None)
parser.add_argument('adata', help="File with Anndata object. It should have spliced and unspliced layers.")

args = parser.parse_args()

print(args)

# Funtions to quickly save plots
def save_stream(adata, file, format=["svg", "png"], **kwargs):
    for extension in format:
        file_name = file + "." + extension
        scv.pl.velocity_embedding_stream(adata, save = file_name, **kwargs)

def save_grid(adata, file, format=["svg", "png"], **kwargs):
    for extension in format:
        file_name = file + "." + extension
        scv.pl.velocity_embedding_grid(adata, save = file_name, **kwargs)

def save_embedding(adata, file, format=["svg", "png"], **kwargs):
    for extension in format:
        file_name = file + "." + extension
        scv.pl.velocity_embedding(adata, save = file_name, **kwargs)

# Read the data
adata = scv.read(args.adata)

if args.velocity:
    # Filter and normalize
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    
    # Calculate cell moments
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

if args.umap:
    scv.tl.umap(adata)

for mode in args.mode:
    save_folder = args.directory + "/scvelo_" + mode + "/"
    os.makedirs(save_folder, exist_ok=True)

    if args.velocity:
        # Estimate RNA velocity with stochastic mode
        if mode == 'dynamical':
            scv.tl.recover_dynamics(adata, n_jobs=4)
        
        scv.tl.velocity(adata, mode=mode)
        scv.tl.velocity_graph(adata, n_jobs=6)
        scv.tl.velocity_embedding(adata, basis='umap')

        # Plot velocity results
        save_stream(adata, file=(save_folder + args.write_name + "_stream"), format=args.format, basis="umap", color="Celltypes_2020_names", palette='tab20', figsize=(9, 7), legend_loc='right margin')

        save_grid(adata, file=(save_folder + args.write_name + "_grid"), format=args.format, basis='umap', color="Celltypes_2020_names", palette='tab20', figsize=(9, 7), legend_loc='right margin', arrow_length=5, arrow_size=3)

        save_embedding(adata, file=(save_folder + args.write_name + "_velocity"), format=args.format, basis="umap", color="Celltypes_2020_names", palette='tab20', figsize=(9, 7), legend_loc='right margin', arrow_length=5, arrow_size=3, dpi=120)

    # Calculate confidence
    if args.confidence:
        scv.tl.velocity_confidence(adata)
        scv.pl.scatter(adata, basis='umap', color=['velocity_length', 'velocity_confidence'], color_map='coolwarm', perc=[0, 100], save=(save_folder + args.write_name + "_confidence.svg"))
        scv.pl.scatter(adata, basis='umap', color=['velocity_length', 'velocity_confidence'], color_map='coolwarm', perc=[0, 100], save=(save_folder + args.write_name + "_confidence.png"))

    # Rank genes
    if args.rank:
        scv.tl.rank_velocity_genes(adata, groupby='Celltypes_2020_names', min_corr=0.3)
        df = scv.get_df(adata.uns['rank_velocity_genes']['names'])
        df.to_csv(save_folder + args.write_name + "_genes_rank.csv", index=False)

    # Save object
    if args.save is not None:
        save_name = args.save + "_" + mode + '.h5ad'
        adata.write(save_name, compression='gzip')
