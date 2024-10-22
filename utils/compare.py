import sys
import dendropy
from dendropy.calculate import treecompare

taxon_namespace = dendropy.TaxonNamespace()

path1 = sys.argv[1]
path2 = sys.argv[2]

tree1 = dendropy.Tree.get(path=path1, schema="newick", taxon_namespace=taxon_namespace)
tree2 = dendropy.Tree.get(path=path2, schema="newick", taxon_namespace=taxon_namespace)

rf_distance = treecompare.weighted_robinson_foulds_distance(tree1, tree2)

print(f"{path1} vs. {path2}: {rf_distance:.4f}")
