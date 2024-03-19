import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patches as mpatches
import seaborn as sns



sns.set_palette("pastel")  # Use Seaborn's "pastel" color scheme

# Assuming the existence of "co_occurence_mat.csv"
matrix = pd.read_csv("co_occurence_mat.csv", header=None).values
np.fill_diagonal(matrix, 0)  # Remove self-loops

G = nx.from_numpy_array(matrix)

# Define node colors, sizes, and other attributes
colors = ['red', 'green', 'blue', 'yellow', 'cyan', 'magenta']
node_sizes = [2.90619987, 6.99241391, 2.434912567, 3.66670514, 3.01191436, 3.415088525]  # Node size for each degree increase by 100

# Normalize node sizes to [100, 500] range
min_size, max_size = min(node_sizes), max(node_sizes)
adjusted_sizes = [500 + ((size - min_size) / (max_size - min_size) * (2000 - 500)) for size in node_sizes]

# Node positioning and plotting
pos = nx.spring_layout(G)  # Node position

# Define light colors and node types for the legend
light_colors = ['#FBB4AE', '#B3CDE3', '#CCEBC5', '#DECBE4', '#FED9A6', '#FFFFCC']
node_types = ['CTCCTGACCTCA(BACH1)', 'CCTGTAATCCCAGC(OTX2)', 'CACACACACAC(ACAC Repeat)', 'AGGCAGGAG', 'CCTCTGC', 'CCAGGGC']

# Set custom node labels
node_labels = {node: f"Motif{node + 1}" for node in G.nodes()}

plt.figure(figsize=(10, 8))
nx.draw_networkx_nodes(G, pos, node_color=light_colors[:len(G.nodes())], alpha=0.8, node_size=adjusted_sizes, node_shape='s')

# Draw thicker edges for weights > 0.6
for (u, v, d) in G.edges(data=True):
    if d['weight'] > 0.6:
        nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=3, alpha=0.7)

# Draw edge labels
edge_labels = nx.get_edge_attributes(G, "weight")
edge_labels = {k: f"{v:.2f}" for k, v in edge_labels.items()}
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)

# Draw node labels
nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=12)

# Create and display the legend
patch_list = [mpatches.Patch(color=color, label=label) for color, label in zip(light_colors, node_types)]
plt.legend(handles=patch_list, title="Motif List", loc='upper right', bbox_to_anchor=(1.19, 1.0))

plt.title("Motif Co-occurrence Network of Sample A673_DMSO_rep1 in Enhancer Region")
plt.axis('off')
plt.subplots_adjust(right=0.85)  # Adjust layout to prevent legend from being clipped
plt.show()

