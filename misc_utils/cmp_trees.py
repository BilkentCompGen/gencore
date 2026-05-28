import argparse
import itertools
import dendropy
from dendropy.calculate import treecompare


def load_tree(path, taxon_namespace):
    return dendropy.Tree.get(
        path=path,
        schema="newick",
        taxon_namespace=taxon_namespace,
        preserve_underscores=True,
    )


def quartet_signature(tree, taxa4):
    """
    Returns the resolved quartet split, e.g. (('A','B'), ('C','D')).
    Returns None for unresolved/star quartet.
    """
    labels = list(taxa4)
    pdm = tree.phylogenetic_distance_matrix()

    # Three possible quartet splits
    splits = [
        ((labels[0], labels[1]), (labels[2], labels[3])),
        ((labels[0], labels[2]), (labels[1], labels[3])),
        ((labels[0], labels[3]), (labels[1], labels[2])),
    ]

    scores = []
    for (a, b), (c, d) in splits:
        ta = tree.taxon_namespace.get_taxon(label=a)
        tb = tree.taxon_namespace.get_taxon(label=b)
        tc = tree.taxon_namespace.get_taxon(label=c)
        td = tree.taxon_namespace.get_taxon(label=d)

        score = (
            pdm.distance(ta, tb)
            + pdm.distance(tc, td)
        )
        scores.append(score)

    min_score = min(scores)

    # If tied, quartet is unresolved
    if scores.count(min_score) > 1:
        return None

    best = splits[scores.index(min_score)]

    return tuple(sorted([
        tuple(sorted(best[0])),
        tuple(sorted(best[1]))
    ]))


def quartet_distance(tree1, tree2):
    taxa1 = {t.label for t in tree1.taxon_namespace}
    taxa2 = {t.label for t in tree2.taxon_namespace}
    common_taxa = sorted(taxa1 & taxa2)

    if len(common_taxa) < 4:
        raise ValueError("Need at least 4 shared taxa for quartet distance.")

    total = 0
    different = 0

    for taxa4 in itertools.combinations(common_taxa, 4):
        q1 = quartet_signature(tree1, taxa4)
        q2 = quartet_signature(tree2, taxa4)

        total += 1
        if q1 != q2:
            different += 1

    return different, total, different / total


def compare_newick_trees(path1, path2):
    taxon_namespace = dendropy.TaxonNamespace()

    tree1 = load_tree(path1, taxon_namespace)
    tree2 = load_tree(path2, taxon_namespace)

    tree1.encode_bipartitions()
    tree2.encode_bipartitions()

    rf = treecompare.symmetric_difference(tree1, tree2)

    max_rf = 2 * (len(taxon_namespace) - 3)
    normalized_rf = rf / max_rf if max_rf > 0 else 0

    weighted_rf = treecompare.weighted_robinson_foulds_distance(tree1, tree2)

    q_diff, q_total, q_norm = quartet_distance(tree1, tree2)

    return {
        "RF distance": rf,
        "Normalized RF distance": normalized_rf,
        "Weighted RF distance": weighted_rf,
        "Quartet distance": q_diff,
        "Total quartets": q_total,
        "Normalized quartet distance": q_norm,
    }

def main():
    parser = argparse.ArgumentParser(
        description="Compare two Newick trees using RF, normalized RF, weighted RF, and quartet distance."
    )

    parser.add_argument("tree1", help="Path to first Newick tree file")
    parser.add_argument("tree2", help="Path to second Newick tree file")

    args = parser.parse_args()

    results = compare_newick_trees(args.tree1, args.tree2)

    for metric, value in results.items():
        print(f"{metric}: {value}")


if __name__ == "__main__":
    main()