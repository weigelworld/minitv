"""
Handles phylogenetic trees
"""

from collections import defaultdict
import dendropy
import itertools
import pandas
from io import StringIO


def _deep_scan(node):
    """Recursively generates a dictionary tree structure.

    Credit
    ======

    https://github.com/AliTVTeam/AliTV-perl-interface/blob/master/lib/AliTV/Tree.pm

    :param node: a dendropy.Node object
    :type node: dendropy.Node
    :rtype: dict
    """
    if node.is_leaf():
        return {
            'name': str(node.taxon).strip('"\'')
        }
    else:
        current_stage = defaultdict(list)
        for child in node.child_node_iter():
            current_stage['children'].append(_deep_scan(child))

    return current_stage


def _balance_node_depth(tree):
    """Converts a proper dendropy Tree into a Tree that requires all leaf nodes to be the same depth.
    Additionally, adds a final intermediate node before all leaves, required for proper AliTV drawing.

    :param tree: a dendropy tree
    :type tree: a dendropy tree
    """

    required_level = max([leaf.level() for leaf in tree.leaf_node_iter()]) + 1

    # traverse through leaf nodes
    for node in tree.leaf_node_iter():
        new_nodes_required = required_level - node.level()

        parent = node.parent_node
        parent.remove_child(node)

        for i in range(0, new_nodes_required):
            intermediate_node = dendropy.Node()
            parent.add_child(intermediate_node)
            parent = intermediate_node

        parent.add_child(node)


def tree_to_json_structure(tree):
    """converts a dendropy Tree into a tree recognized by AliTV.

    Credit
    ======

    https://github.com/AliTVTeam/AliTV-perl-interface/blob/master/lib/AliTV/Tree.pm

    :param tree: a dendropy tree
    :type tree: dendropy.Tree
    :return: AliTV-format tree
    :rtype: dict
    """

    # convert dendropy tree into AliTV compatible tree with extra nodes
    _balance_node_depth(tree)

    return _deep_scan(tree.seed_node)


def distance_matrix_csv_from_tuple(pairs, distances):
    """Transforms a pair of tuples describing FASTA file name pairs and their respective distances
    into a distance matrix (csv format, stored in a StringIO object).

    :param pairs: iterable of pairs of names
    :type pairs: collections.Iterable[tuple(str, str)]
    :param distances: list of distances
    :type distances: list[int]
    :return: csv-formatted distance matrix
    :rtype: io.StringIO
    """

    names = set(itertools.chain(*pairs))

    df = pandas.DataFrame(index=names, columns=names).fillna(0)

    for i, pair in enumerate(pairs):
        df[pair[0]].loc[pair[1]] = distances[i]
        df[pair[1]].loc[pair[0]] = distances[i]

    csv_string = df.to_csv(index=True, header=True)

    csv_buffer = StringIO(csv_string)
    return csv_buffer


def distance_matrix_from_alignment_groups(alignment_groups):
    """Calculates a distance matrix from alignment files

    :param alignment_groups: dict with tuple keys in the form (target, query) and values of lists of
    minitv.align.Alignment objects
    :type alignment_groups: dict[tuple(str, str), list[minisyn.align.Alignment]]
    :return: distance matrix
    :rtype: dendropy.PhylogeneticDistanceMatrix
    """

    distances = []

    for target, query in alignment_groups:

        alignment_group = alignment_groups[target, query]

        # calculate simple "distance matrix" using count of alignment blocks
        #  more alignment blocks --> more diverged sequences (even with overlapping blocks)
        alignment_blocks_count = len(alignment_group)

        distances.append(alignment_blocks_count)

    return dendropy.PhylogeneticDistanceMatrix.from_csv(
        distance_matrix_csv_from_tuple(
            tuple([(target, query) for target, query in alignment_groups]),
            distances)
    )


def tree_from_alignment_groups(alignment_groups):
    """Calculates a tree from alignment groups.

    :param alignment_groups: dict with tuple keys in the form (target, query) and values of lists of
    minitv.align.Alignment objects
    :type alignment_groups: dict[tuple(str, str), list[minisyn.align.Alignment]]
    :return: tree calculated from alignment groups
    :rtype: dendropy.Tree
    """

    # cluster with UPGMA
    pdm = distance_matrix_from_alignment_groups(alignment_groups)
    tree = pdm.upgma_tree()

    # sort tree
    tree.ladderize()

    return tree
