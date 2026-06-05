import argparse
import os

import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix


def remove_inner_labels(tree):
    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.name = None

def sort_tree(tree):
    def sort_clades(clade):
        if clade.is_terminal():
            return
        clade.clades.sort(key=lambda x: min([t.name for t in x.get_terminals()]))
        for child in clade.clades:
            sort_clades(child)
    
    sort_clades(tree.root)
    return tree

def read_phy(filename, normalize=False):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    n_sequences = int(lines[0].strip())
    
    labels = []
    for i in range(1, n_sequences + 1):
        labels.append(lines[i][:10].strip())
      
    matrix = []
    for line in lines[1:]:
        matrix.append(list(map(np.float64, line[10:].strip().split())))
    
    distance_matrix = np.array(matrix, dtype=np.float64)

    if normalize:
        max_value = np.max(distance_matrix)
        if max_value > 0:
            distance_matrix /= max_value
        else:
            raise ValueError("Distance matrix contains all zeros or invalid values.")

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

def construct_nj_tree(labels, distance_matrix, root_by="midpoint"):
    dm = DistanceMatrix(labels, distance_matrix)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    
    if root_by == "midpoint":
        tree.root_at_midpoint()
    elif root_by == "outgroup":
        outgroup_label = labels[0]
        outgroup = next(clade for clade in tree.find_clades() if clade.name == outgroup_label)
        tree.root_with_outgroup(outgroup)
    
    # Remove inner labels
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

def draw_tree(tree_str, output_file):
    """
    Draw tree using ete3.

    ete3 is imported only here, so users who only want Newick output
    do not need ete3, PyQt, or Qt working in their environment.
    """
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

    try:
        from ete3 import Tree, TreeStyle
    except ImportError as e:
        raise ImportError(
            "Drawing requires ete3 with TreeStyle support. "
            "Install ete3 and its Qt/PyQt dependencies, or run without --draw."
        ) from e

    ete_tree = Tree(tree_str, format=1, quoted_node_names=True)

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.branch_vertical_margin = 5

    png_file = output_file + ".png"
    ete_tree.render(png_file, w=1000, units="px", tree_style=ts)

    print(f"Tree image saved to {png_file}")

def main(input_file, mode, normalize, output_file, draw=False):
    if not input_file.endswith(".phy") and not input_file.endswith(".newick"):
        raise ValueError("Unsupported file format. Please provide a .phy or .newick")

    if input_file.endswith(".newick"):
        with open(input_file, "r") as f:
            tree = f.read().strip()

        if not tree:
            raise ValueError(f"Empty Newick file: {input_file}")

        if draw:
            draw_tree(tree, output_file)

        return
        
    if mode == "nj":
        labels, matrix = read_phy(input_file, normalize)
        tree = construct_nj_tree(labels, matrix)
        tree = sort_tree(tree)
        save_tree(tree, output_file+".nj")
    else: 
        labels, matrix = read_phy(input_file, normalize)
        tree = construct_upgma_tree(labels, matrix)
        tree = sort_tree(tree)
        save_tree(tree, output_file+".upgma")

    if draw:
        tree_str = tree.format("newick")
        draw_tree(tree_str, output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Construct a phylogenetic tree using UPGMA and save it as an image and Newick file.")
    parser.add_argument("input_file", help="Path to the input file containing the distances (phy format)")
    parser.add_argument('--nj', action='store_true')
    parser.add_argument('--upgma', action='store_true')
    parser.add_argument('--normalize', action='store_true')
    parser.add_argument( "--draw", action="store_true", help="Render a PNG image using ete3. Requires ete3 and Qt/PyQt.")
    args = parser.parse_args()

    output_file_base = os.path.splitext(args.input_file)[0]
    mode = "upgma" if args.upgma else "nj"
    
    main(args.input_file, mode, args.normalize, output_file_base, args.draw)