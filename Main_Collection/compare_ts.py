"""
Code originally from Wilder Wohns
(Wilder wrote functions kc_distance_tree and kc_distance_ts)
but adapted by Marianne Aspbury for particular purposes
(added other functions)
Functions added September 2019, commenting improved Oct 2019.
"""
import tskit

import numpy as np
import itertools

def kc_distance_tree(tree1, tree2, topo_v_age=0):
    """
    Returns the Kendall-Colijn distance between the specified pair of trees.
    topo_v_age determines weight of topology vs branch lengths in calculating
    the distance. Set topo_v_age at 0 to only consider topology, set at 1 to
    only consider branch lengths. See Kendall & Colijn (2016):
    https://academic.oup.com/mbe/article/33/10/2735/2925548
    """
    samples = tree1.tree_sequence.samples()
    if not np.array_equal(samples, tree2.tree_sequence.samples()):
        raise ValueError("Trees must have the same samples")
    if not len(tree1.roots) == len(tree2.roots) == 1:
        raise ValueError("Trees must have one root")
    k = samples.shape[0]
    n = (k * (k - 1)) // 2
    m = [np.ones(n + k), np.ones(n + k)]
    M = [np.zeros(n + k), np.zeros(n + k)]
    for tree_index, tree in enumerate([tree1, tree2]):
        stack = [(tree.root, 0, tree.time(tree.root))]
        while len(stack) > 0:
            u, depth, time = stack.pop()
            children = tree.children(u)
            for v in children:
                stack.append((v, depth + 1, tree.time(v)))
            for c1, c2 in itertools.combinations(children, 2):
                for u in tree.samples(c1):
                    for v in tree.samples(c2):
                        a = min(u, v)
                        b = max(u, v)
                        pair_index = a * (a - 2 * k + 1) // -2 + b - a - 1
                        assert m[tree_index][pair_index] == 1
                        m[tree_index][pair_index] = depth
                        M[tree_index][pair_index] = tree.time(tree.root) - time
            if len(tree.children(u)) == 0:
                M[tree_index][u + n] = tree.branch_length(u)
    return np.linalg.norm((1 - topo_v_age) * (m[0] - m[1]) + topo_v_age * (M[0] - M[1]))

def kc_distance_ts(ts_1, ts_2, lambda_param):
    """
    Returns the Kendall-Colijn distance between two specified tree sequences:
    one overall KC distance considering all comparable trees in the tree seq.
    lambda_param determines weight of topology vs branch lengths in calculating
    the distance. Set lambda_param at 0 to only consider topology, set at 1 to
    only consider branch lengths. See Kendall & Colijn (2016):
    https://academic.oup.com/mbe/article/33/10/2735/2925548
    """
    ts_1_breakpoints = np.array(list(ts_1.breakpoints()))
    ts_2_breakpoints = np.array(list(ts_2.breakpoints()))
    comb_breakpoints = np.sort(np.unique(np.concatenate([ts_1_breakpoints, ts_2_breakpoints])))
    comb_isin_1 = np.isin(comb_breakpoints, ts_1_breakpoints)
    comb_isin_2 = np.isin(comb_breakpoints, ts_2_breakpoints)
    ts_1_trees = ts_1.trees()
    ts_2_trees = ts_2.trees()
    kc = 0
    last_breakpoint = 0

    for index in range(len(comb_breakpoints)):
        try:
            if comb_isin_1[index]:
                 tree_1 = next(ts_1_trees)
            if comb_isin_2[index]:
                 tree_2 = next(ts_2_trees)

        except StopIteration:
            last_breakpoint = comb_breakpoints[index]
            break

        if not len(tree_1.roots) == len(tree_2.roots) == 1:
            print('tree skipped')
            pass
        else:
#            print('tree kc used')
            kc += kc_distance_tree(tree_1, tree_2, lambda_param) * (comb_breakpoints[index + 1] - comb_breakpoints[index])
    kc /= (last_breakpoint-comb_breakpoints[0])
    print(index)
    return(kc)

def kc_distance_each_tree(ts_1, ts_2, topo_v_age):
    """
    Returns the Kendall-Colijn distances between each pair of trees in the two
    specified tree sequences. Returns a list of KC tree distances covering the
    tree sequence.
    topo_v_age determines weight of topology vs branch lengths in calculating
    the distance. Set topo_v_age at 0 to only consider topology, set at 1 to
    only consider branch lengths. See Kendall & Colijn (2016):
    https://academic.oup.com/mbe/article/33/10/2735/2925548
    """
    kc = np.full((ts_1.num_trees + ts_2.num_trees), fill_value=0.0)
    ts_1_breakpoints = np.array(list(ts_1.breakpoints()))
    ts_2_breakpoints = np.array(list(ts_2.breakpoints()))
    comb_breakpoints = np.sort(np.unique(np.concatenate([ts_1_breakpoints, ts_2_breakpoints])))
    comb_isin_1 = np.isin(comb_breakpoints, ts_1_breakpoints)
    comb_isin_2 = np.isin(comb_breakpoints, ts_2_breakpoints)
    ts_1_trees = ts_1.trees()
    ts_2_trees = ts_2.trees()
    last_breakpoint = 0

    for index in range(len(comb_breakpoints)):
        try:
            if comb_isin_1[index]:
                 tree_1 = next(ts_1_trees)
            if comb_isin_2[index]:
                 tree_2 = next(ts_2_trees)

        except StopIteration:
            last_breakpoint = comb_breakpoints[index]
            break

        if not len(tree_1.roots) == len(tree_2.roots) == 1:
#            print('tree skipped')
            pass
        else:
#            print('tree kc used')
            kc[index] = kc_distance_tree(tree_1, tree_2, topo_v_age)

#    print('index is ', index)
    return(kc)



def check_ts(ts_1, ts_2, lambda_param):
    ''' Something I used just to check functionality'''
    ts_1_breakpoints = np.array(list(ts_1.breakpoints()))
    ts_2_breakpoints = np.array(list(ts_2.breakpoints()))
    comb_breakpoints = np.sort(np.unique(np.concatenate([ts_1_breakpoints, ts_2_breakpoints])))
    comb_isin_1 = np.isin(comb_breakpoints, ts_1_breakpoints)
    comb_isin_2 = np.isin(comb_breakpoints, ts_2_breakpoints)
    ts_1_trees = ts_1.trees()
    ts_2_trees = ts_2.trees()

    num=0

    for index in range(len(comb_breakpoints)):
        try:
            if comb_isin_1[index]:
                 tree_1 = next(ts_1_trees)
            if comb_isin_2[index]:
                 tree_2 = next(ts_2_trees)

        except StopIteration:
            break

        if not len(tree_1.roots) == len(tree_2.roots) == 1:
            num+=1
            print(tree_1.roots, tree_2.roots)
            if num > 15:
                break
    return num
