import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

def read_gene_groups(file_path):
    """
    Reads gene group information from a file.
    """
    with open(file_path, 'r') as file:
        gene_groups = {i: line.strip() for i, line in enumerate(file, 1)}
    return gene_groups

def find_highest_expression(df, gene_groups, gene_column, expression_column):
    """
    Initializes a dictionary to store the gene with the highest expression in each group.
    """
    max_gene_expressions = {}
    for group_id, genes in gene_groups.items():
        group_genes = genes.split(' ')
        max_gene = df[df[gene_column].isin(group_genes)].nlargest(1, expression_column)
        if not max_gene.empty:
            max_gene_name = max_gene[gene_column].iloc[0]
            max_expression = max_gene[expression_column].iloc[0]
            max_gene_expressions[group_id] = (max_gene_name, max_expression)
    return max_gene_expressions

def calculate_percentiles(df, max_gene_expressions, expression_column):
    """
    Sorts the expressions and calculates the ranking percentiles for each representative gene's expression.
    """
    sorted_df = df.sort_values(by=expression_column, ascending=False)
    expression_percentiles = {}
    for group_id, (gene, max_expression) in max_gene_expressions.items():
        percentile = np.mean(sorted_df[expression_column] > max_expression) * 100
        expression_percentiles[group_id] = (gene, percentile)
    return expression_percentiles

def generate_heatmap_data(excel_files, gene_groups):
    """
    Generates data for the heatmap, setting the index as gene names.
    """
    heatmap_data = pd.DataFrame(index=gene_groups)
    for file in excel_files:
        df = pd.read_csv(file, sep='\t')
        gene_column, expression_column = df.columns[0], df.columns[4]
        max_gene_expressions = find_highest_expression(df, gene_groups, gene_column, expression_column)
        percentiles = calculate_percentiles(df, max_gene_expressions, expression_column)
        heatmap_data[file] = heatmap_data.index.map(lambda gene: percentiles.get(gene, None))
    heatmap_data = heatmap_data.transpose()
    return heatmap_data

def plot_heatmap(heatmap_data):
    """
    Plots a heatmap using seaborn's magma colormap.
    """
    magma = sns.color_palette("magma", as_cmap=True)
    colors = magma(np.linspace(0, 1, 20))
    colors = np.vstack((colors, np.array([0.75, 0.75, 0.75, 1])))  # Adding gray color
    custom_cmap = ListedColormap(colors)
    boundaries = np.arange(0, 22, 1)  # From 0 to 21
    norm = BoundaryNorm(boundaries, custom_cmap.N, clip=True)
    plt.figure(figsize=(20, 12))
    sns.heatmap(heatmap_data, annot=True, cmap=custom_cmap, norm=norm, linewidths=0.05, linecolor='black')
    plt.xticks(rotation=45)
    plt.xlabel('TF')
    plt.ylabel('Sample')
    plt.title('Gene Expression Percentiles Across Samples')
    plt.show()

# Example usage
gene_groups = read_gene_groups('motif_alt_ids_enh.txt')  # Replace with your gene groups file path

# Directory containing the files
directory = "./"
xls_files = [filename for filename in os.listdir(directory) if filename.endswith(".xls") or filename.endswith(".xlsx")]
xls_files.sort()

heatmap_data = generate_heatmap_data(xls_files, gene_groups)

# Flatten the dictionary into a single string for each group
column_names = {k: ' '.join(v.split()) for k, v in gene_groups.items()}

# Rename the columns of the DataFrame
heatmap_data.rename(columns=column_names, inplace=True)

plot_heatmap(heatmap_data)

# Save the heatmap_data INCLUDING ROW NAMES
heatmap_data.to_csv('heatmap_data.csv')
