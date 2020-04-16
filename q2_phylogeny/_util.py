# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import skbio


def midpoint_root(tree: skbio.TreeNode) -> skbio.TreeNode:
    return tree.root_at_midpoint()

def lca_root(tree: skbio.TreeNode, tip_labels: str = None) -> skbio.TreeNode:
    if tip_labels is not None:
        stripped_tip_labels = [tl.strip() for tl in tip_labels.split(',')]
        if len(stripped_tip_labels) < 2:
            raise ValueError("Provide at least two tip labels to determine "
                             "LCA.")
        else:
            tip_nodes_for_lca = []
            tip_label_set = set(stripped_tip_labels)
            for tl in tip_label_set:
                print('tip:', tl)
                try:
                    tip_nodes_for_lca.append(tree.find(tl))
                except:
                    print("Unable to find tip: %s" %tl)
                    continue
                    # would be nice to write a text file for these.
                #except ValueError:
                #    print("Unable to find tip: %s" %tl)
                #    continue
            if len(tip_nodes_for_lca) < 2:
                raise ValueError("After searching for the provided tree tips, "
                                 "less than two were found! Two tips are "
                                 "required to find the LCA node. Check your "
                                 "labels!")
            else:
                try:
                    lcan = tree.lowest_common_ancestor(tip_nodes_for_lca)
                    return tree.root_at(lcan)
                except ValueError:
                    print("Not enough tip labels to determine LCA! "
                          "Are you providing correct labels?")
    else:
        raise ValueError("Provide comma separated list of at least two tip "
                         "labels to determine LCA.")

def robinson_foulds(trees: skbio.TreeNode, labels: str = None,
                    missing_tips: str = 'error') -> skbio.DistanceMatrix:
    if labels is None:
        labels = ['tree_%d' % d for d in range(1, len(trees) + 1)]
    elif len(trees) != len(labels):
        raise ValueError("The number of trees and labels must match.")

    tips = [{t.name for t in tree.tips()} for tree in trees]
    shared_tips = set.intersection(*tips)
    if not shared_tips:
        raise ValueError("No tip names are shared between these trees.")
    if missing_tips == 'intersect-all':
        trees = [t.shear(shared_tips) for t in trees]
    elif missing_tips == 'error':
        all_tips = set.union(*tips)
        if shared_tips != all_tips:
            SKIP = 10
            miss = list(all_tips - shared_tips)
            missing_repr = ", ".join(map(repr, miss[:SKIP]))
            if len(miss) > SKIP:
                missing_repr += ", ...<%d ommitted>" % (len(miss) - SKIP)

            raise ValueError("Not all tips are shared between trees: "
                             + missing_repr)
    else:
        raise ValueError("Unknown argument for missing_tips=%r"
                         % missing_tips)

    return skbio.DistanceMatrix.from_iterable(
        trees, metric=skbio.TreeNode.compare_rfd,
        keys=labels, validate=False)
