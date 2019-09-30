"""
Simulate a ts with ancient samples, infer ts with modern samples and add ancient samples on
From Wilder Wohns, adapted for HIV simulations, and for .py not jup notebk
Marianne Aspbury 27/9/2019
"""

##### Imports #####

from IPython.display import display, SVG
import numpy as np
import pandas as pd
import subprocess
import sys
import logging
from collections import namedtuple
# import vcf

import tsdate
import msprime
import tsinfer
import tskit

from tsinfer.eval_util import *

#sys.path.insert(1, 'C:\\Users\\mazmysta\\OneDrive - Nexus365\\BDI_proj\\scripts\\github\\tsdate')

#%% My HIV part
import pickle
import random # to pick from list of poss roots of subtree
from IndividualClass import Individual # might say unused, but not true, needed to read pickled file
from truncate_ts_samples import truncate_ts_samples as truncate_ts_samples

from Bio import SeqIO

### Import HIV genome ###
HIV_file_loc = "fasta/"

with open(HIV_file_loc + "HIV_CompleteGenome.txt", "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        #print(record.id)
        #print(len(record.seq))
        HIV_genome = record.seq


### import the trees ###
file_loc = ""

with open(file_loc + "pickled_data_all.pickle", "rb") as f:
    total_tree = pickle.load(f)

#proof it works
print('Pickle loaded, 1st root num kids: {}'.format(total_tree['-1_-1'].total_children_num))
#Pickle loaded, 1st root num kids: 16758

# Saving list of people with desired num kids beneath
samples_poss = []
desired_overall_children = 5 ## Seemingly my max is 30 pops with 10 samples, 100 pop size and 1000 length genome.

for key in total_tree:
    if total_tree[key].total_children_num == desired_overall_children:
        samples_poss.append(key)

random.seed(4) # repeatable choice(s) - sequence same
start = random.choice(samples_poss)
start='58582_0' #for 5 kids
print('starting sample: ', start)

def preorder_traversal(u):
    all_nodes = []
    stack = [u]
    while len(stack) > 0:
        v = stack.pop()
        if total_tree[v].direct_children_num > 0: #Returns True if the specified node is not a leaf
            stack.extend(total_tree[v].infected_others)
        all_nodes.append(v)
    return all_nodes

list_of_root_and_kids = preorder_traversal(start)

### msprime assimilation
# Find total num of generations (days) to model over
times_loop = []

for item in list_of_root_and_kids:
    #float to take as number not string, and access list element with [0]
    times_loop.append(float(total_tree[item].time_infected_by[0]))

#since floats, can do maths
total_time = (max(times_loop)-min(times_loop))*365

#Furthest back in time is the smallest date (lowest year), do everything in diffs
gens_list = []
for item in list_of_root_and_kids:
    #float to take as number not string, and access list element with [0]
    diff_calc = 365*(max(times_loop)-float(total_tree[item].time_infected_by[0]))
    gens_list.append(diff_calc + 100) # add 100 so most recent infection is '100 gens (days) in past'

#    print(float(total_tree[item].time_infected_by[0]))

# info on dates so far
for i in range(len(gens_list)):
    print((list_of_root_and_kids)[i],
          #times_loop[i],
          gens_list[i])

######## msprime variable populations setup #########

## source population
PopSource = msprime.PopulationConfiguration(initial_size = 1e8, growth_rate = 0)

## number of populations in present time of model, not including source pop...
final_num_pops = len(list_of_root_and_kids)

## sample_sizes for each population and effective pop sizes
sample_size = 10
infection_size = 1
stable_pop_size = 100

# ## subpops based on infected people, all subpops that want to exist at end of sim (present time) need stated here
pop_list = [PopSource]
sample_list = []
#Setting up the end, so all pops exist and at stable pop size, no death in simulation time
for pop in range(final_num_pops):
#    print(pop)
    pop_list.append(
        msprime.PopulationConfiguration(initial_size = stable_pop_size, growth_rate = 0)
                  )
    #historical samples rather than contemporaneous ones, 1 week after infection
    for sample in range(sample_size):
        if sample < sample_size//2:
            sample_list.append(msprime.Sample(population=(pop+1),
                                              time = gens_list[pop] - 30))
        else:
            sample_list.append(msprime.Sample(population=(pop+1),
                                              time = 0))

# no migration between sources accross time, only infection events,
# so migration matrix is zeros
M = np.zeros((final_num_pops+1,final_num_pops+1))


transfers_list = []
for entry in range(len(pop_list)):
#    print(entry)
    if entry == 0: # ignore 0 since this is the source pop
        pass

    elif entry == 1: # 1 is root so needs own bit
        entry_ID = list_of_root_and_kids[entry-1]
#        print(entry, entry_ID)
        dest_index = 0 #infected from source
        transfer_time = gens_list[entry-1] # time infected still stored.
#        print(transfer_time)
        transfers_list.append(msprime.MassMigration(time = transfer_time, source = entry, dest = dest_index, proportion = 1))

    elif entry > 1: # 1 is root so needs own bit
        # get the index of the infected_by population
        # index of current population is its index in pop_list (index in list_of_root... + 1)
        entry_ID = list_of_root_and_kids[entry-1]
#        print(entry, entry_ID)
        dest_ID = total_tree[entry_ID].infected_by[0]
        dest_index = list_of_root_and_kids.index(dest_ID) + 1
#        print(dest_index, dest_ID)
        transfer_time = gens_list[entry-1]
#        print(transfer_time)
        transfers_list.append(msprime.MassMigration(time = transfer_time, source = entry, dest = dest_index, proportion = 1))


#### Bottleneck list initiation and creation
Pop_bottleneck_ie_growth_list = []
for entry in range(len(pop_list)):
    if entry > 0: # ignore 0 since this is the source pop. Only need infection time for this so root doesn't need own case

        transfer_time = gens_list[entry-1]

        # setting 0 in population before infection
        pop_entry_null = msprime.PopulationParametersChange(
                time = (transfer_time+0.001), initial_size=1e-8, growth_rate=0, population = entry) #i.e. for epochs inf-> 100*entry, this is growth rate

        #infection size setting
        pop_entry_bottleneck_start = msprime.PopulationParametersChange(
                time = transfer_time, initial_size=infection_size, growth_rate=0, population = entry) #i.e. for epochs inf-> 100*entry, this is growth rate

        #growth after infection setting - trial that 20 gens in future (less time back) gives appropriate growth for these params of pop_size = 100, rate = 0.25
        pop_entry_bottleneck_end = msprime.PopulationParametersChange(
                time = transfer_time-20, growth_rate = 0.24, initial_size=stable_pop_size, population = entry) #i.e. for epochs inf-> 100*entry, this is growth rate

        #save to list for manip outside loop
        Pop_bottleneck_ie_growth_list.extend((pop_entry_null, pop_entry_bottleneck_start, pop_entry_bottleneck_end))

# put all events together then sort them
events = Pop_bottleneck_ie_growth_list + transfers_list
events_sorted = sorted(events, key=lambda x: x.time, reverse=False)

ts2 = msprime.simulate(population_configurations=pop_list, migration_matrix = M,
                       demographic_events = events_sorted,
                       length = len(HIV_genome),
                       random_seed = 17, recombination_rate = 0.7e-4,
                       mutation_rate=2.5e-5,
                       samples=sample_list,
                       end_time=500000) # end time incase doesn't coalesce, occassionally there are issues

print(ts2.first().roots, ts2.first().time(ts2.first().roots[0]))


# 5 pops, need 6 colours including fake source pop
colour_map = {0:"grey", 1:"blue", 2:"red", 3:"orange", 4:"green", 5:"cyan", 6:"black"}
node_colours = {u.id: colour_map[u.population] for u in ts2.nodes()}

# draw only one tree

X = 1
i = 0
for tree in ts2.trees():
    if i < X:
        display(SVG(tree.draw(height=500, width = 1500,
                              format='SVG',
                              node_colours=node_colours,
                              tree_height_scale='rank',
                              #max_tree_height=5000
                              )))
        print("Tree {} covers [{:.2f}, {:.2f}); TMRCA = {:.4f}".format(
            tree.index, *tree.interval, tree.time(tree.roots[0])))
#        print(tree.branch_length(46))
    else:
        break
    i+=1


#%% Dating/inference
#import tsdate

##### Function Definitions #####

def prune_sites(ts, keep_sites):
    # Remove sites and mutations not referenced in keep_sites
    def keep_with_offset(keep, data, offset):
        """Copied from keep_intervals"""
        # We need the astype here for 32 bit machines
        lens = np.diff(offset).astype(np.int32)
        return (
            data[np.repeat(keep, lens)],
            np.concatenate([
                np.array([0], dtype=offset.dtype),
                np.cumsum(lens[keep], dtype=offset.dtype)]))

    tables = ts.dump_tables()
    sites = ts.tables.sites
    mutations = ts.tables.mutations

    new_as, new_as_offset = keep_with_offset(
        keep_sites, sites.ancestral_state, sites.ancestral_state_offset)
    new_md, new_md_offset = keep_with_offset(
        keep_sites, sites.metadata, sites.metadata_offset)
    keep_mutations = keep_sites[mutations.site]
    tables.sites.set_columns( #replaces existing
        position=sites.position[keep_sites],
        ancestral_state=new_as,
        ancestral_state_offset=new_as_offset,
        metadata=new_md,
        metadata_offset=new_md_offset)

    new_ds, new_ds_offset = keep_with_offset(
        keep_mutations, mutations.derived_state, mutations.derived_state_offset)
    new_md, new_md_offset = keep_with_offset(
        keep_mutations, mutations.metadata, mutations.metadata_offset)

    site_map = np.cumsum(keep_sites, dtype=mutations.site.dtype) - 1
    parent_map = np.cumsum(keep_mutations, dtype=mutations.parent.dtype) - 1
    # parent -1 always maps to parent -1
    parent_map = np.append(parent_map, np.array([-1], dtype=mutations.parent.dtype))
    # Don't bother dealing with connected mutations as we are removing the whole site
    tables.mutations.set_columns(
        site=site_map[mutations.site[keep_mutations]],
        node=mutations.node[keep_mutations],
        derived_state=new_ds,
        derived_state_offset=new_ds_offset,
        parent=parent_map[mutations.parent[keep_mutations]],
        metadata=new_md,
        metadata_offset=new_md_offset)

    return tables.tree_sequence()

def index_map(a, b, missing=-1):
    """
    Return a vector of the same length as `a`, giving the index into `b`
    of matching values in `a`.
    """
    idx = np.searchsorted(b, a)
    overflow = idx == len(b)
    nomatch = b[idx[~overflow]] != a[~overflow]
    nomatch = np.concatenate((nomatch, np.ones(np.count_nonzero(overflow), bool)))
    idx[nomatch] = -1
    return idx


###### Simulation and corrections for dating/inference #######
# Define ancient time, this needs to come first
#ANCIENT_TIME = 5000/25

# GET modern samples and ancient samples
ts = ts2
sample_nodes = ts.tables.nodes.flags == 1

modern_sample_indices = np.where((ts.tables.nodes.time[sample_nodes] ==0) == True)[0]
ancient_sample_indices = np.where((ts.tables.nodes.time[sample_nodes] ==0) == False)[0]
#modern_sample_indices = [i for i, e in enumerate(sample_list) if e.time == 0]
#ancient_sample_indices = [i for i, e in enumerate(sample_list) if e.time != 0]


#def simulate_ts():
#    modern_samples = [msprime.Sample(population=0, time=0) for x in range(90)]
#    ancient_samples = [msprime.Sample(population=0, time=ANCIENT_TIME) for x in range(10)]
#    samples = modern_samples + ancient_samples
#    return(msprime.simulate(samples=samples, mutation_rate=1e-8,
#                            recombination_rate=1e-8, length=1e4,
#                            Ne=10000))
#
#ts = simulate_ts()

sample_data = tsinfer.formats.SampleData.from_tree_sequence(ts)

# Remove non modern samples from tree sequence
modern_samples = tsinfer.SampleData(path= "modern_only_hiv.samples",
                                    sequence_length=ts.sequence_length)

for individual in range(len(modern_sample_indices)):
    modern_samples.add_individual(
        ploidy=1, metadata={})

for v in ts.variants():
    modern_samples.add_site(
        position=v.site.position, alleles=v.alleles,
        genotypes=v.genotypes[modern_sample_indices])
modern_samples.finalise()

# make ancient sample data files
ancient_samples = tsinfer.SampleData(path= "ancient_only_hiv.samples", sequence_length=ts.sequence_length)
for individual in range(len(ancient_sample_indices)):
    ancient_samples.add_individual(
        ploidy=1, metadata={})

for v in ts.variants():
#    if np.random.randint(4) != 1:
        ancient_samples.add_site(
            position=v.site.position, alleles=v.alleles,
            genotypes=v.genotypes[ancient_sample_indices])
ancient_samples.finalise()

#%%
# Infer and date tree from modern samples

primary_ts = ts.simplify(modern_sample_indices, filter_sites=False)
primary_samples = tsinfer.SampleData.from_tree_sequence(primary_ts)

ancestors = tsinfer.generate_ancestors(primary_samples)
ancestors_ts = tsinfer.match_ancestors(primary_samples, ancestors) # This only has inference sites

primary_inferred_ts = tsinfer.match_samples(primary_samples, ancestors_ts, simplify=False)
primary_inferred_ts_simplified = primary_inferred_ts.simplify(np.where(primary_inferred_ts.tables.nodes.flags == 1)[0], keep_unary=True)

tsdate.date(primary_inferred_ts_simplified, Ne=stable_pop_size, mutation_rate=2.5e-5)

#%%
# rest of inference- augmenting older samples in

augment_samples = ancient_samples
## re-inserting older samples
augment_samples = augment_samples.copy()

# First we must check that the (inference) sites in augment_samples are a subset of
# those in ancestors_ts. Any inference sites unique to augment_samples cannot be used,
# so should be automatically marked as "not for inference" (and a warning given).
# Sites marked as not-for-inference in the augment_samples may still be used if they
# correspond to an inference site in ancestors_ts: here we should check that the
# allele -> index mapping in the augment_samples is compatible with variant.alleles

sample_variant = namedtuple("sample_variant", "id position inference alleles")
augment_samples_sites = (sample_variant(id, pos, inf, alleles)
    for id, (pos, inf, alleles) in enumerate(zip(
        # TODO - check if there is an efficient iterator over these zarr arrays
        # or if we should be doing e.g. augment_samples.sites_position[:] instead
        augment_samples.sites_position,
        augment_samples.sites_inference,
        augment_samples.sites_alleles)))

main_samples_sites = (sample_variant(id, pos, 1, alleles)
    for id, (pos, alleles) in enumerate(zip(
        # TODO - check if there is an efficient iterator over these zarr arrays
        # or if we should be doing e.g. augment_samples.sites_position[:] instead
        primary_samples.sites_position[:][primary_samples.sites_inference[:].view(bool)],
        primary_samples.sites_alleles[:][primary_samples.sites_inference[:].view(bool)])))


keep_ancestors_sites = np.zeros((ancestors_ts.num_sites, ), dtype=bool)
augment_site = next(augment_samples_sites)
for main_site in main_samples_sites:
    while main_site.position >= augment_site.position:
        pos = augment_site.position
        if main_site.position == augment_site.position:
            keep_ancestors_sites[main_site.id] = 1
            a_m = [m for m in main_site.alleles if m is not None]
            a_x = [m for m in augment_site.alleles if m is not None]

            if augment_site.inference:
                if a_x != a_m:
                    raise ValueError(
                        "Inference site at pos {} incompatible".format(pos))
            else:
                # site in augment_samples is not marked for inference (e.g.
                # because it is monomorphic in the augment_samples set) but
                # we might want to use it anyway, as it could help place
                # the sample
                augment_samples.sites_inference[augment_site.id] = 1
                if a_x != a_m[:len(a_x)]:
                    raise ValueError(
                        "Inference site at pos {} incompatible".format(pos))
                augment_samples.sites_alleles[augment_site.id] = a_m
        else:
            if augment_site.inference:
                pass
                # augment_site is an inference site not in the main inference.
#                   I've commented the warning cos it's annoying
#                logging.warning(
#                    "Variant at pos {} is marked for inference in the augment"
#                    " file but not in the inferred tree sequence".format(pos))
#                augment_samples.sites_inference[augment_site.id] = 0
        try:
            augment_site = next(augment_samples_sites)
        except StopIteration:
            # set the position to force the for loop to finish
            augment_site = sample_variant(None, np.inf, None, None)

# final samples to add
augment_samples.finalise()

augment_samples = augment_samples.copy()

# Faster alternative
aug_to_primary = index_map(augment_samples.sites_position[:], ancestors_ts.tables.sites.position[:])
aug_sites_alleles = augment_samples.sites_alleles[:]
pri_inference_sites_alleles = primary_samples.sites_alleles[:][primary_samples.sites_inference[:].view(bool)]
aug_sites_inference = augment_samples.sites_inference[:]

# Sites not marked for inference in augment_samples, but which are used in ancestors_ts
both_used = np.logical_and(aug_sites_inference == 1, aug_to_primary >= 0).view(bool)
# Sanity check
if not np.array_equal(
        aug_sites_alleles[both_used],
        pri_inference_sites_alleles[aug_to_primary[both_used]]):
    raise ValueError("Inference alleles in original and new sample data files incompatible")

only_in_primary = np.logical_and(aug_sites_inference == 0, aug_to_primary >= 0).view(bool)
pri_al = pri_inference_sites_alleles[aug_to_primary[only_in_primary]]
aug_al = aug_sites_alleles[only_in_primary]
tmp_pos = augment_samples.sites_position[:][only_in_primary]
non_matching_alleles = aug_al != pri_al
if np.any(non_matching_alleles):
    for aug, pri, pos, idx in zip(
            aug_al[non_matching_alleles],
            pri_al[non_matching_alleles],
            tmp_pos[non_matching_alleles],
            np.where(only_in_primary)[0][non_matching_alleles]):
        aug_alleles = [a for a in aug if a is not None]
        pri_alleles = [a for a in pri if a is not None]
        if aug_alleles != pri_alleles[:len(aug_alleles)]:
            raise ValueError(
                "Inference sites at position {} have incompatible alleles"
                .format(pos))
        else:
            augment_samples.sites_alleles[idx] = pri
aug_sites_inference[only_in_primary] = 1

only_used_in_augment = np.logical_and(
    augment_samples.sites_inference[:] == 1, aug_to_primary < 0)
if np.any(only_used_in_augment):
    logging.warning(
        "Inference sites at positions {} in the augmenting samples"
        " were not marked for inference originally".format(pos))
aug_sites_inference[only_used_in_augment] = 0

augment_samples.sites_inference = aug_sites_inference
augment_samples.finalise()

primary_to_aug = index_map(ancestors_ts.tables.sites.position[:], augment_samples.sites_position[:])
keep_ancestors_sites = primary_to_aug >= 0

assert np.sum(keep_ancestors_sites) == np.sum(augment_samples.sites_inference[:])

# Make ancestors_ts compatible by keeping only sites that are sample_data inference sites
pruned_ancestors_ts = prune_sites(ancestors_ts, keep_ancestors_sites)
assert np.array_equal(
    augment_samples.sites_position[:][augment_samples.sites_inference[:].view(bool)],
    pruned_ancestors_ts.tables.sites.position)

pruned_ancestors_ts, node_map = pruned_ancestors_ts.simplify(map_nodes=True, reduce_to_site_topology=True)
augmented_matches = tsinfer.match_samples(augment_samples, pruned_ancestors_ts, simplify=False)

# Check the new samples take indexes after the existing nodes in the map
assert pruned_ancestors_ts.num_nodes == np.sum(node_map!=-1) == min(augmented_matches.samples()) == np.max(node_map) + 1

# Now, augment the primary tree sequence with the matched ancient samples
# Add the new nodes, edges, sites, etc into first_match
tables1 = primary_inferred_ts.dump_tables()
tables2 = augmented_matches.tables

# TO DO: add in any sites that were present in extra_samples and not in main_inference_samples
# These can only be placed by parsimony, using map_mutations

# Add any new individuals
existing_individuals = tables1.individuals.num_rows
new_individuals = tables1.individuals.append_columns(
    flags=tables2.individuals.flags,
    # TO DO - add location, location_offset, metadata, and metadata_offset
)
new_individuals = np.arange(  # Need to do this until https://github.com/tskit-dev/tskit/issues/365 is fixed
    existing_individuals, tables1.individuals.num_rows,
    dtype=tables1.nodes.individual.dtype)
# Append -1 to the end, so that -1 always maps to -1
new_individuals = np.append(new_individuals, -1).astype(
    tables1.nodes.individual.dtype)

# Add new nodes.
augmented_nodes = np.arange(pruned_ancestors_ts.num_nodes, tables2.nodes.num_rows)
existing_nodes = tables1.nodes.num_rows
tables1.nodes.append_columns(
    flags=tables2.nodes.flags[:][augmented_nodes],
    time=tables2.nodes.time[:][augmented_nodes],
    population=tables2.nodes.population[:][augmented_nodes],
    individual=new_individuals[tables2.nodes.individual[:][augmented_nodes]] if len(new_individuals) else None
    # TO DO - add metadata & metadata_offset
)
# Need to do this until https://github.com/tskit-dev/tskit/issues/365 is fixed
augmented_nodes_in_pri = np.arange(existing_nodes, tables1.nodes.num_rows)


# === Add new edges ===
# We have the mapping from primary_ts to the ones in the new augmented ts
# We need to reverse this so that we can convert the augmented_ts values
# back to the primary_ts nodes. We also need to add the mapping from
# augmented_matches.samples() to the new augmented_sample_nodes_in_pri
# (we assumes sample order is maintained)

# Extend with nulls
node_map = np.concatenate((
    node_map,
    np.full(tables1.nodes.num_rows - len(node_map), tskit.NULL, dtype=node_map.dtype)))
# Include the new samples
node_map[augmented_nodes_in_pri] = augmented_nodes
# Add an extra item to the reverse map, so that -1 maps correctly
reverse_map = np.full(augmented_matches.num_nodes + 1, -1, dtype=node_map.dtype)
reverse_map[node_map] = np.arange(0,len(node_map))

# Here we detect newly added edges by finding those that come from an
# augmented child, but this will miss any path compressed edges
# so this can probably be improved to include any edges with a new
# child
new_edges = np.isin(tables2.edges.child[:], augmented_nodes)
tables1.edges.append_columns(
    left=tables2.edges.left[new_edges],
    right=tables2.edges.right[new_edges],
    parent=reverse_map[tables2.edges.parent[new_edges]],
    child=reverse_map[tables2.edges.child[new_edges]]
)

### BREAKS HERE - time[parent] must be greater than time[child]
tables1.sort()

# final tree sequence with older samples too
augmented_ts = tables1.tree_sequence()

#%%
##### Draw some tree(s) #####
## Pick tree sim
upper = 5
lower = 4
X = upper
i = lower
for tree in ts.trees():
    if i < X:
        display(SVG(tree.draw(height=500, width = 1500,
                              format='SVG',
                              #max_tree_height=5000
                              )))
        print("Tree {} covers [{:.2f}, {:.2f}); TMRCA = {:.4f}".format(
            tree.index, *tree.interval, tree.time(tree.roots[0])))
#        print(tree.branch_length(46))
    else:
        break
    i+=1

X = upper
i = lower
for tree in augmented_ts.trees():
    if i < X:
        display(SVG(tree.draw(height=1000, width = 3000,
                              format='SVG',
                              #max_tree_height=5000
                              )))
        print("Tree {} covers [{:.2f}, {:.2f}); TMRCA = {:.4f}".format(
            tree.index, *tree.interval, tree.time(tree.roots[0])))
#        print(tree.branch_length(46))
    else:
        break
    i+=1

augmented_ts_simp = augmented_ts.simplify(augmented_ts.samples(), keep_unary=False)

X = upper
i = lower
for tree in augmented_ts_simp.trees():
    if i < X:
        display(SVG(tree.draw(height=1000, width = 1500,
                              format='SVG',
                              #max_tree_height=5000
                              )))
        print("Tree {} covers [{:.2f}, {:.2f}); TMRCA = {:.4f}".format(
            tree.index, *tree.interval, tree.time(tree.roots[0])))
#        print(tree.branch_length(46))
    else:
        break
    i+=1
### end ###
