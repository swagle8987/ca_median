from dendropy import Node,Tree,TreeList,TaxonNamespace
from copy import deepcopy
import pdb
from tree import CustomNode,CustomTree

def spr(initial_tree, prune_node, regraft_above_node, tree_indexer=None):

    if prune_node.parent_node is initial_tre.seed_node and regraft_above_node is initial_tree.seed_node:
        return initial_tree, prune_node, prune_node.sibling_nodes()[0]

    p_n = initial_tree.find_node(lambda x:x.id==prune_node.id)
    initial_tree.prune_subtree(p_n,update_bipartitions=True)
    initial_tree.update_taxon_namespace()
    add_node_above(initial_tree, p_n, regraft_above_node)
    return initial_tree,p_n, regraft_above_node


def add_node_above(tree, node, add_above_node):
    parent_node = add_above_node.parent_node
    if parent_node is None:
        parent_node = CustomNode(0)
        add_above_node.parent_node = parent_node  # setter in dendropy automatically adds "self" as a child
        parent_node.add_child(node)  # same here
        tree.seed_node = parent_node
    else:
        new_parent = parent_node.new_child()
        parent_node.remove_child(add_above_node)
        new_parent.add_child(node)
        new_parent.add_child(add_above_node)

    tree.reindex()
    tree.update_taxon_namespace()
    return add_above_node.parent_node

def calculate_spr_sites(gene_tree, base_tree, prune_node, cost_function):
    copy_tree = deepcopy(base_tree)
    pruned_gene_trees = TreeList()
    pn = copy_tree.find_node(lambda x:x.id==prune_node.id)
    copy_tree.prune_subtree(pn)
    taxa_lim = list(base_tree.taxon_namespace)
    taxn_ns = TaxonNamespace(taxa_lim)
    label_lim  = [i.label for i in taxa_lim]
    for g in gene_tree:
        pruned_gene_tree = deepcopy(g.extract_tree_with_taxa_labels(labels = label_lim))
        pruned_gene_tree.taxon_namespace = taxn_ns
        pruned_gene_trees.append(pruned_gene_tree)
    pruned_gene_trees.taxon_namespace = taxn_ns
    distance = float('inf')
    best_node = None
    for node in copy_tree.nodes():
        temp_tree = deepcopy(copy_tree)
        temp_pn = deepcopy(pn)
        temp_tree.update_taxon_namespace()
        temp_node = temp_tree.find_node(lambda x:x.id==node.id)
        try:
            add_node_above(temp_tree, temp_pn, temp_node)
        except :
            pdb.set_trace()
        cur_distance = 0
        for g in pruned_gene_trees:
            cur_distance += cost_function(temp_tree, g)
        if cur_distance < distance:
            distance = cur_distance
            best_node = node
    bn = copy_tree.find_node(lambda x: x.id==best_node.id)
    add_node_above(copy_tree, pn, bn)
    copy_tree.update_taxon_namespace()
    return copy_tree,distance

def calculate_spr_neighborhood(gene_tree, base_tree, cost_function):
    distance = float('inf')
    best_tree = None
    x = 0
    for node in base_tree.nodes():
        if node.parent_node:
            copy_tree =  deepcopy(base_tree)
            best_val_neighbor,best_dist = calculate_spr_sites(gene_tree, copy_tree, node, cost_function)
            if best_dist < distance:
                distance = best_dist
                best_tree = best_val_neighbor
    return best_tree, distance
            
    
