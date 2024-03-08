from cluster_affinity import cluster_affinity_cost,cluster_support_cost
import pdb
import argparse
import random
from spr import add_node_above,calculate_spr_sites, calculate_spr_neighborhood
from dendropy import TreeList,Tree,Node,TaxonNamespace,Taxon
from tree import CustomTree,CustomNode

def sym_cluster_affinity_cost(t1,t2):
    t1t2 = cluster_affinity_cost(t1,t2)
    t2t1 = cluster_affinity_cost(t2,t1)
    return (t1t2 + t2t1)/2

def build_starting_tree(gene_trees,cost_function):
    taxa = set([i.label for i in gene_trees.taxon_namespace._taxa[:]])
    taxa = [Taxon(i) for i in taxa]
    random.shuffle(taxa)
    tree = CustomTree()
    seed = CustomNode(0)
    seed.taxon = taxa.pop()
    tree.seed_node = seed
    tree.is_rooted = True

    second_taxon = taxa.pop()
    second_leaf = CustomNode(1)
    second_leaf.taxon = second_taxon
    add_node_above(tree, seed, second_leaf) ## undefined
    tree.update_taxon_namespace()
    added_taxa = [seed.taxon, second_leaf.taxon]
    tree.taxon_namespace = TaxonNamespace(added_taxa)

    while len(taxa) > 0:
        next_taxon = taxa.pop()
        next_leaf = CustomNode(len(taxa))
        next_leaf.taxon = next_taxon
        added_taxa.append(next_leaf.taxon)
        add_node_above(tree,tree.seed_node, next_leaf)
        tree.reindex()
        tree.update_taxon_namespace()
        tree.taxon_namespace = TaxonNamespace(added_taxa)
        best_tree,distance = calculate_spr_sites(gene_trees, tree, next_leaf, cost_function)
        tree = best_tree
        tree.update_taxon_namespace()
    return tree


def load_trees(path):
    TreeList.DEFAULT_TREE_TYPE=CustomTree
    vals = TreeList.get(path=path,schema="newick",rooting="force-rooted",suppress_internal_node_taxa=True,suppress_leaf_node_taxa=True)
    for g in vals:
        g.regenerate_taxon()
    return vals

def find_median_tree(path,cutoff,stat=False):
    gene_trees = load_trees(path)
    start_tree = build_starting_tree(gene_trees, sym_cluster_affinity_cost)
    prev_dist = float('inf')
    prev_tree = None
    for i in range(cutoff):
        best_tree,dist = calculate_spr_neighborhood(
            gene_trees, start_tree, sym_cluster_affinity_cost)
        if prev_dist <= dist:
            break
        else:
            prev_dist = dist
            prev_tree = best_tree
        if stat:
            print("Current best distance {}".format(prev_dist))
    if stat:
        print("Best Distance Found: {}".format(prev_dist))
        print("Distance from each gene tree")
        for g in gene_trees:
            print("{} \t {} \t {}".format(g.as_string("newick").strip(),prev_tree.as_string("newick").strip(),sym_cluster_affinity_cost(g, prev_tree)))
    return prev_dist,prev_tree

def main():
    program_description = ""
    parser = argparse.ArgumentParser(description=program_description,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('treefile', help="The treefile containing the trees. Each tree must be on a newline ending in a semi-colon")
    parser.add_argument('cutoff', metavar="cutoff", help="The maximum spr neighborhood that should be explored", type=int)
    parser.add_argument('--stat', help="Reports the best tree found in each generation and the final cost", action="store_true")
    parser.add_argument('--output_tree_1','-o1', metavar="o1",help="The output file for the resultant tree 1, defaults to tree1_result.tre",default="tree1_result.tre")
    args = parser.parse_args()
    prev_dist,prev_tree = find_median_tree(args.treefile,args.cutoff,args.stat)
    prev_tree.write(path=args.output_tree_1,schema="newick")

if __name__=="__main__":
    main()

    

