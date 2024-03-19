def reverse_complement(sequence):
    """
    Returns the reverse complement of a DNA sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in reversed(sequence))

def extract_meme_consensus_with_alignment(file_path):
    """
    Extracts the Multilevel consensus sequence and its alignment from a file.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    characters1 = []
    characters2 = []

    for line in lines:
        if "Multilevel" in line:
            line = line.strip()
            characters1.extend(list(line))
        elif "consensus" in line:
            line = line.strip()
            characters2.extend(list(line))

    if len(characters2) < len(characters1):
        characters2 += [" " for _ in range(len(characters1) - len(characters2))]

    i = 10
    while characters1[i] == " ":
        i += 1

    for j in range(i - 1, len(characters2)):
        if characters2[j] == " ":
            characters2[j] = "-"

    consensus_with_dashes = "".join(characters2[i:])
    multilevel_sequence = "".join(characters1[i:])

    return (multilevel_sequence, consensus_with_dashes)

def longest_common_substring_with_alternatives(s1, s2, alternatives):
    """
    Finds the longest common substring between two sequences, considering alternative characters.
    """
    max_length = 0
    longest_sub = ""
    alt_dict = {i: set(chars.split('/')) for i, chars in enumerate(alternatives) if chars != '-'}

    for i in range(len(s1)):
        for j in range(len(s2)):
            length = 0
            while i + length < len(s1) and j + length < len(s2) and (s1[i + length] == s2[j + length] or s1[i + length] in alt_dict.get(j + length, set())):
                length += 1
                if length > max_length:
                    max_length = length
                    longest_sub = s1[i:i + length]

    return longest_sub

def extract_kmap_consensus(filepath):
    """
    Extracts sequences from a file.
    """
    sequences = []
    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('>'):
                sequences.append(line)

    return sequences

def main(filename1, filename2):
    """
    Main function to process files and calculate consensus.
    """
    filepath1 = filename1 + "/meme.txt"
    filepath2 = filename2 + "/top_chain_new_lc.txt"
    consensus1 = extract_meme_consensus_with_alignment(filepath1)
    consensus2 = extract_kmap_consensus(filepath2)
    result = []

    for seq in consensus2:
        overlap_po = longest_common_substring_with_alternatives(seq, consensus1[0], consensus1[1])
        recall_po = len(overlap_po) / len(consensus1[0])
        precision_po = len(overlap_po) / len(seq)
        f1_po = 2 * precision_po * recall_po / (precision_po + recall_po) if precision_po + recall_po != 0 else 0

        overlap_ne = longest_common_substring_with_alternatives(reverse_complement(seq), consensus1[0], consensus1[1])
        recall_ne = len(overlap_ne) / len(consensus1[0])
        precision_ne = len(overlap_ne) / len(seq)
        f1_ne = 2 * precision_ne * recall_ne / (precision_ne + recall_ne) if precision_ne + recall_ne != 0 else 0

        if f1_po > f1_ne:
            result.append((precision_po, recall_po, f1_po, seq, consensus1[0], consensus1[1]))
        else:
            result.append((precision_ne, recall_ne, f1_ne, seq, consensus1[0], consensus1[1]))

    max_result = max(result, key=lambda x: x[2])
    csv_output = f"{max_result[0]},{max_result[1]},{max_result[2]},{max_result[3]},{max_result[4]},{max_result[5]}"
    return csv_output

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Process sequences to find consensus.')
    parser.add_argument('--kmap_out', type=str, help='Path to the directory containing KMAP result')
    parser.add_argument('--meme_out', type=str, help='Path to the directory containing MEME result')
    args = parser.parse_args()
    print(main(args.meme_out, args.kmap_out))
