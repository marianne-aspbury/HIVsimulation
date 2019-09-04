# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 12:40:19 2019

Truncating simulated data and forming missing data genomes

Updated truncated samples def so genotypes matrix works for inference steps
23/8/2019 by Yan Wong

@author: mazmysta taken from yan wong

"""
import msprime
import tsinfer
import numpy as np

def truncate_ts_samples(ts, average_span, random_seed, min_span=5):
    """
    Create a tree sequence that has sample nodes which have been truncated
    so that they span only a small region of the genome. The length of the
    truncated spans is given by a poisson distribution whose mean is average_span
    but which cannot go below a fixed min_span, or above the sequence_length
    
    mutations above removed edges are removed so doesn't break genotypes matrix

    Samples are truncated by removing the edges that connect them to the rest
    of the tree.
    """
    def keep_with_offset(keep, data, offset):
       """Copied from keep_intervals"""
       # We need the astype here for 32 bit machines
       lens = np.diff(offset).astype(np.int32)
       return (data[np.repeat(keep, lens)],
               np.concatenate([
                   np.array([0], dtype=offset.dtype),
                   np.cumsum(lens[keep], dtype=offset.dtype)]))
    np.random.seed(random_seed)
    # Make a list of (left,right) tuples giving the new limits of each sample
    # Keyed by sample ID.
    to_slice = {}
    # for simplicity, we pick lengths from a poisson distribution of av 300 bp
    for sample_id, span in zip(
            ts.samples(), np.random.poisson(average_span, ts.num_samples)):
        span = max(span, min_span)
        span = min(span, ts.sequence_length)
        start = np.random.uniform(0, ts.sequence_length-span)
        to_slice[sample_id] = (start, start+span)

    tables = ts.dump_tables()
    tables.edges.clear()
    for e in ts.tables.edges:
        if e.child not in to_slice:
            left, right = e.left, e.right
        else:
            if e.right <= to_slice[e.child][0] or e.left >= to_slice[e.child][1]:
                continue  # this edge is outside the focal region
            else:
                left = max(e.left, to_slice[e.child][0])
                right = min(e.right, to_slice[e.child][1])
        tables.edges.add_row(left, right, e.parent, e.child)
    # Remove mutations above isolated nodes
    mutations = tables.mutations
    keep_mutations = np.ones((mutations.num_rows,), dtype = bool)
    positions = tables.sites.position[:]
    for i, m in enumerate(mutations):
        if m.node in to_slice:
            if not to_slice[m.node][0] <= positions[m.site] < to_slice[m.node][1]:
                keep_mutations[i] = False
    new_ds, new_ds_offset = keep_with_offset(
        keep_mutations, mutations.derived_state, mutations.derived_state_offset)
    new_md, new_md_offset = keep_with_offset(
        keep_mutations, mutations.metadata, mutations.metadata_offset)
    mutations_map = np.append(np.cumsum(keep_mutations) - 1, [-1])
    mutations_map = mutations_map.astype(mutations.parent.dtype)
    # parent -1 always maps to parent -1
    tables.mutations.set_columns(
        site=mutations.site[keep_mutations],
        node=mutations.node[keep_mutations],
        derived_state=new_ds,
        derived_state_offset=new_ds_offset,
        parent=mutations_map[mutations.parent[keep_mutations]],
        metadata=new_md,
        metadata_offset=new_md_offset)
    return tables.tree_sequence()

# original def
#def truncate_ts_samples(ts, average_span, random_seed, min_span=5):
#        """
#        Create a tree sequence that has sample nodes which have been truncated
#        so that they span only a small region of the genome. The length of the
#        truncated spans is given by a poisson distribution whose mean is average_span
#        but which cannot go below a fixed min_span, or above the sequence_length
#
#        Samples are truncated by removing the edges that connect them to the rest
#        of the tree.
#        """
#        np.random.seed(random_seed)
#        # Make a list of (left,right) tuples giving the new limits of each sample
#        # Keyed by sample ID.
#        keep = {}
#        # for simplicity, we pick lengths from a poisson distribution of av 300 bp
#        for sample_id, span in zip(
#                ts.samples(), np.random.poisson(average_span, ts.num_samples)):
#            span = max(span, min_span)
#            span = min(span, ts.sequence_length)
#            start = np.random.uniform(0, ts.sequence_length-span)
#            keep[sample_id] = (start, start+span)
#
#        tables = ts.dump_tables()
#        tables.edges.clear()
#        for e in ts.tables.edges:
#            if e.child not in keep:
#                left, right = e.left, e.right
#            else:
#                if e.right <= keep[e.child][0] or e.left >= keep[e.child][1]:
#                    continue  # this edge is outside the focal region
#                else:
#                    left = max(e.left, keep[e.child][0])
#                    right = min(e.right, keep[e.child][1])
#            tables.edges.add_row(left, right, e.parent, e.child)
#        return tables.tree_sequence()