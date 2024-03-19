import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from tqdm import tqdm

def extract_sequences_from_fastq_gz(file_name):
    """
    Extracts sequences from a gzipped FASTQ file.
    """
    sequences = []
    with gzip.open(file_name, 'rt') as file:
        while True:
            identifier = file.readline().strip()
            sequence = file.readline().strip()
            file.readline()  # plus line (ignored)
            file.readline()  # quality line (ignored)
            if not identifier:
                break
            sequences.append(sequence)
    return sequences

def read_sequences_from_fasta(fasta_file):
    """
    Reads sequences from a FASTA file.
    """
    sequences = []
    with open(fasta_file, 'r') as file:
        current_seq = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = ''
            else:
                current_seq += line
        if current_seq:
            sequences.append(current_seq)
    return sequences

def hamming_distance(x, y):
    """
    Calculates the Hamming distance between two strings.
    """
    if len(x) != len(y):
        raise ValueError("Strings must be of the same length")
    return sum(xi != yi for xi, yi in zip(x, y))

def fit_curve(finite_distances):
    """
    Fits a curve to the given distances and plots the histogram and fitted curve.
    """
    def gaussian(x, mu, sigma):
        return 0.055 * np.exp(-2 * ((x - mu) / sigma) ** 2)

    counts, bin_edges = np.histogram(finite_distances, bins=30, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    smoothed_counts = savgol_filter(counts, 11, 3)

    popt, _ = curve_fit(gaussian, bin_centers, smoothed_counts, bounds=([20, 0], [70, 15]), maxfev=10000)

    plt.figure(figsize=(10, 5))
    plt.bar(bin_centers, counts, width=np.diff(bin_centers)[0], edgecolor='k', alpha=0.5, label='Original Histogram')
    plt.plot(bin_centers, gaussian(bin_centers, *popt), color='blue', label='Fitted Curve')
    plt.title('Smoothed Distance Distribution with Fitted Curve')
    plt.xlabel('Distance')
    plt.ylabel('Frequency')
    plt.legend()
    plt.show()

    mu, sigma = popt
    print(f"Parameters mu={mu}, sigma={sigma}")

def calculate_pairwise_distances(sequences):
    """
    Calculates the pairwise Hamming distances between sequences.
    """
    n = len(sequences)
    dist_matrix = np.zeros((n, n), dtype=int)
    for i in tqdm(range(n)):
        for j in range(i + 1, n):
            dist_matrix[i, j] = dist_matrix[j, i] = hamming_distance(sequences[i], sequences[j])
    return dist_matrix

def analyze_distances(distances):
    """
    Analyzes distances to find a suitable threshold for clustering.
    """
    k = 600
    nearest_neighbors_indices = np.argsort(distances, axis=1)[:, :k]
    nearest_neighbors_distances = np.take_along_axis(distances, nearest_neighbors_indices, axis=1)
    avg_distances = np.mean(nearest_neighbors_distances, axis=1)
    
    plt.hist(avg_distances, bins=30, edgecolor='k')
    plt.title('Distance Distribution')
    plt.xlabel('Distance')
    plt.ylabel('Frequency')
    plt.show()

    finite_distances = avg_distances[np.isfinite(avg_distances)]
    fit_curve(finite_distances)

fasta_path = './data/seqs_alined_rm_ref.fasta'
sequences = read_sequences_from_fasta(fasta_path)
distances = calculate_pairwise_distances(sequences)
analyze_distances(distances)
