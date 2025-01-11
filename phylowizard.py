import numpy as np
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo
import argparse
import os

def remove_inner_labels(tree):
    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.name = None

def read_phy(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    n_sequences = int(lines[0].strip())
    
    labels = []
    for i in range(1, n_sequences + 1):
        labels.append(lines[i][:10].strip())
      
    matrix = []
    for line in lines[1:]:
        matrix.append(list(map(float, line[10:].strip().split())))
    
    distance_matrix = np.array(matrix)

    lower_triangle = []
    for i in range(len(distance_matrix)):
        row = distance_matrix[i][:i + 1].tolist()
        lower_triangle.append(row)

    return labels, lower_triangle

def construct_upgma_tree(labels, distance_matrix):
    dm = DistanceMatrix(labels, distance_matrix)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    
    # remove inner labels
    remove_inner_labels(tree)
    return tree

def construct_nj_tree(labels, distance_matrix):
    dm = DistanceMatrix(labels, distance_matrix)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    
    # remove inner labels if needed
    remove_inner_labels(tree)
    return tree

def save_tree(tree, output_file):
    # save newick
    with open(output_file + ".newick", "w") as f:
        Phylo.write(tree, f, "newick")
        print(f"Tree saved to {output_file}.newick")

def read_newick(filename):
    with open(filename, 'r') as file:
        newick_str = file.read().strip()
    tree = Phylo.read(filename, 'newick')
    return tree

def main(input_file, output_file):
    if not input_file.endswith('.phy'):
        raise ValueError("Unsupported file format. Please provide a .phy or .newick file.")
    
    labels, matrix = read_phy(input_file)
    tree1 = construct_upgma_tree(labels, matrix)
    save_tree(tree1, output_file+".upgma")
    tree2 = construct_nj_tree(labels, matrix)
    save_tree(tree2, output_file+".nj")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Construct a phylogenetic tree using UPGMA and save it as an image and Newick file.")
    parser.add_argument("input_file", help="Path to the input file containing the distances (phy format)")
    args = parser.parse_args()

    output_file_base = os.path.splitext(args.input_file)[0]

    main(args.input_file, output_file_base)