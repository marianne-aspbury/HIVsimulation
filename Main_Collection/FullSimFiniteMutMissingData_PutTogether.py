# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 12:33:36 2019

Combined new simulation using missing data and finite mutation model
Putting everything together

@author: marianne aspbury
"""
import sys
import os
import msprime
import tsinfer
import numpy as np
from IPython.display import display
from IPython.display import SVG

import pickle
import random # to pick from list of poss roots of subtree
from IndividualClass import Individual # might say unused, but not true, needed to read pickled file
from truncate_ts_samples import truncate_ts_samples as truncate_ts_samples

from Bio import SeqIO

base_loc = os.path.dirname(os.path.abspath(__file__))
### Import HIV genome ###
HIV_file_loc = os.path.join(base_loc, "fasta")

with open(os.path.join(HIV_file_loc, "HIV_CompleteGenome.txt"), "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        #print(record.id)
        #print(len(record.seq))
        HIV_genome = record.seq

#file_loc = 'C:\\Users\\mazmysta\\OneDrive - Nexus365\\BDI_proj\\scripts\\PopulationsHIV\\simulation1\\'
#with open(file_loc + "HIV_fullsim_fasta_file.txt", "r") as handle:
#    for record in SeqIO.parse(handle, "fasta"):
#        print(len(record.seq))


### import the trees ###
with open(os.path.join(base_loc, 'pickled_data_all.pickle'), 'rb') as f:
    total_tree = pickle.load(f)

#proof it works
print('Pickle loaded, 1st root num kids: {}'.format(total_tree['-1_-1'].total_children_num))
#Pickle loaded, 1st root num kids: 16758

# Saving list of people with desired num kids beneath
samples_poss = []
desired_overall_children = 20 ## Seemingly my max is 30 pops with 10 samples, 100 pop size and 1000 length genome.

for key in total_tree:
    if total_tree[key].total_children_num == desired_overall_children:
        samples_poss.append(key)

print(len(samples_poss))
# 8 for kids = 100

#pick only one random sample
random.seed(4) # repeatable choice(s) - sequence same
start = random.choice(samples_poss)
#for 5 kids: start='58582_0'
print(start)

print(total_tree[start].infected_others)

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

#print(preorder_traversal(start), len(preorder_traversal(start)))

#print(list(set(list_of_root_and_kids) - set(start)))

# all info required for mini-transmission chain
#for item in list_of_root_and_kids:
#    print(item, total_tree[item].infected_others,
#          total_tree[item].time_infected_others,
#          total_tree[item].age_birth,
#          total_tree[item].age_death)


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
sample_size = 20
infection_size = 1
stable_pop_size = 100
# ## subpops based on infected people, all subpops that want to exist at end of sim (present time) need stated here
pop_list = [PopSource]

#Setting up the end, so all pops exist and at stable pop size, no death in simulation time
for pop in range(final_num_pops):
#    print(pop)
    pop_list.append(
        msprime.PopulationConfiguration(sample_size = sample_size, initial_size = stable_pop_size, growth_rate = 0)
                  )

# no migration between sources accross time, only infection events,
    # so migration matrix is zeros
M = np.zeros((final_num_pops+1,final_num_pops+1))

# Now get transmission events from the data. Use index as population number, but +1 since have fake source pop at index 0.
#for i in list_of_root_and_kids:
#    print(list_of_root_and_kids.index(i) + 1)

####--- new version with sub-pops ---####

## a simple model where independent sub-pop is infection derived from source pop
# if infected by true pop, need to state when diverged from past pop if that's the case
# Oddly source is the destination, i.e. direction of migration is dest -> source if forwards in time view.
# backwards in time means that destination is destination (but it's where the migration has come from)

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

#check as expected
print(transfers_list)

# compare
#for i in range(len(gens_list)):
#    print((list_of_root_and_kids)[i],
#          #times_loop[i],
#          gens_list[i])

## now have set of populations, the transfers for pops (infection events)
    # still need the bottlenecks (pop growth & stabilisation)
    # then can order(sort) the complete demography list by time and simulate

## Bottlenecks: add population growth in so infections from source pop are only ~ 1-5 virions, which then balloons to e.g. ~1000

#### Bottleneck list initiation and creation
Pop_bottleneck_ie_growth_list = []
for entry in range(len(pop_list)):
    if entry > 0: # ignore 0 since this is the source pop. Only need infection time for this so root doesn't need own case

        transfer_time = gens_list[entry-1]

        #infection size setting
        pop_entry_bottleneck_start = msprime.PopulationParametersChange(
                time = transfer_time, initial_size=infection_size, growth_rate=0, population = entry) #i.e. for epochs inf-> 100*entry, this is growth rate

        #growth after infection setting - trial that 20 gens in future (less time back) gives appropriate growth for these params of pop_size = 100, rate = 0.25
        pop_entry_bottleneck_end = msprime.PopulationParametersChange(
                time = transfer_time-20, growth_rate = 0.23, initial_size=stable_pop_size, population = entry) #i.e. for epochs inf-> 100*entry, this is growth rate

        #save to list for manip outside loop
        Pop_bottleneck_ie_growth_list.extend((pop_entry_bottleneck_start, pop_entry_bottleneck_end))

# put all events together then sort them
events = Pop_bottleneck_ie_growth_list + transfers_list
events_sorted = sorted(events, key=lambda x: x.time, reverse=False)

#check
#for event in events_sorted:
#    print(event.__dict__) # just for easier digestion of output


my_history = msprime.DemographyDebugger(
    population_configurations=pop_list, migration_matrix = M,
    demographic_events = events_sorted)

#my_history.print_history()



### plot how pop changes ##
#
#time_steps= range(1,int(np.max(gens_list))+100,2)
## print('pop0:', my_history.population_size_trajectory(time_steps)[:,0])
## print('pop1:', my_history.population_size_trajectory(time_steps)[:,1])
## print('time:', np.array(time_steps))
## plot the populations, matplotlib understands array of y's as multiple y's so don't need to call individually
#
#import matplotlib.pyplot as plt
#plt.figure(1)
##plt.plot(time_steps, my_history.population_size_trajectory(time_steps)[:,0])
##plt.plot(time_steps, my_history.population_size_trajectory(time_steps)[:,1])
##plt.rc('axes', prop_cycle=(cycler(color=['r', 'g', 'b', 'y'])))
#fig, ax = plt.subplots(figsize=(15, 6), dpi=80)
#ax.set_prop_cycle(color=["green", "blue", "red", "orange", "grey", "cyan", "black"][1:])
#plt.plot(time_steps, my_history.population_size_trajectory(time_steps)[:,1:], '--', alpha=0.5) # this will plot each y (pop size var.) separately
#plt.xlim(np.max(time_steps),0) # switch the order of time so present (0) is RHS and past is LHS (max time step)
##plt.xlim(np.max(time_steps),1) # switch the order of time so present (0) is RHS and past is LHS (max time step)
##plt.ylim(np.log(0.5), np.log(150))
##plt.axvline(x=100, color='k', linestyle='-', alpha=0.5) # add a vertical line for migration step
##plt.legend(('1','2','3','4','5','6'),loc='best')
##box = ax.get_position()
##ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
## Put a legend to the right of the current axis
##plt.yscale("log")
##plt.xscale("log")
#ax.legend(('1: ' + list_of_root_and_kids[0],
#           '2: ' + list_of_root_and_kids[1],
#           '3: ' + list_of_root_and_kids[2],
#           '4: ' + list_of_root_and_kids[3],
#           '5: ' + list_of_root_and_kids[4],
#           '6: ' + list_of_root_and_kids[5]),
# loc='center left', bbox_to_anchor=(1.02, 0.5))
#plt.show()
## time = 0 is present, larger time is past


############ Simulation time ############

## simulate this extended simple model
# NO mutations (rate 2.5e-5 for HIV) since add in finite way later
ts2 = msprime.simulate(population_configurations=pop_list, migration_matrix = M,
                       demographic_events = events_sorted,
                       length = len(HIV_genome),
                       random_seed = 17, recombination_rate = 0.7e-4,
                       mutation_rate=2.5e-5,
                       end_time=100000) # end time incase doesn't coalesce, occassionally there are issues


# 5 pops, need 6 colours including fake source pop
#colour_map = {0:"grey", 1:"blue", 2:"red", 3:"orange", 4:"green", 5:"cyan", 6:"black"}
#node_colours = {u.id: colour_map[u.population] for u in ts2.nodes()}

# draw only one tree

X = 1
i = 0
for tree in ts2.trees():
    if i < X:
        display(SVG(tree.draw(height=800, width = 10000,
                              format='SVG',
                              tree_height_scale='rank')))
        print("Tree {} covers [{:.2f}, {:.2f}); TMRCA = {:.4f}".format(
            tree.index, *tree.interval, tree.time(tree.roots[0])))
#        print(tree.branch_length(46))
    else:
        break
    i+=1



########### Add finite site mutations - code from Yan Wong #############
import tskit
import msprime
import itertools
import random
import math
import numpy as np

# Construct a tree sequence with integerized breakpoints
#length = 1000
#recomb_map = msprime.RecombinationMap.uniform_map(length, rate=1, num_loci=length)
#ts = msprime.simulate(30, mutation_rate=1, recombination_map=recomb_map)

ts = ts2

null_char_allele = b"\0"
states = np.array(['A','T','G','C'], dtype='|S1')

tables = ts.dump_tables()
tables.sites.clear()
tables.mutations.clear()

## Using HIV fasta genome instead

print(HIV_genome[0], type((HIV_genome[0])))
HIV_genome_array = np.array(HIV_genome, dtype='|S1')
print(HIV_genome_array[0], type(HIV_genome_array[0]))

# INCLUDING ALL SITES NOT JUST VARIANT SITES
variable_sites = itertools.groupby(ts.sites(), lambda x: math.floor(x.position))
variant_pos, sites_at_vpos = next(variable_sites)

for pos in range(len(HIV_genome)):
  ancestral_state = HIV_genome_array[pos]
  site_id = tables.sites.add_row(pos, ancestral_state)
  if variant_pos == pos:
      # Order mutations by time
      mutations = []
      for site in sites_at_vpos:
        mutations.extend(site.mutations)
      mutations.sort(key = lambda x: ts.node(x.node).time, reverse=True)
      for m in mutations:
        # Assign mutations with null parents & derived_states (will be filled in later)
        tables.mutations.add_row(site_id, m.node, derived_state=null_char_allele, parent=tskit.NULL)
      try:
        variant_pos, sites_at_vpos = next(variable_sites)
      except StopIteration:
        variant_pos, sites_at_vpos = -1, []

# THIS IS ONLY VARIABLE SITES
#for integer_pos, sites_at_pos in itertools.groupby(ts.sites(), lambda x: math.floor(x.position)):
##    ancestral_state = random.choice(states)
#    ancestral_state = HIV_genome_array[integer_pos]
##    print(ancestral_state)
#    site_id = tables.sites.add_row(integer_pos, ancestral_state)
##  print(site_id)
#  # Order mutations by time
#    mutations = []
#    for site in sites_at_pos:
#        mutations.extend(site.mutations)
#    mutations.sort(key = lambda x: ts.node(x.node).time, reverse=True)
#    for m in mutations:
#        # Assign mutations with null parents & derived_states (will be filled in later)
#        tables.mutations.add_row(site_id, m.node, derived_state=null_char_allele, parent=tskit.NULL)
#

#print(tables.mutations)
# Assign parents
tables.compute_mutation_parents()
#print(tables.mutations)

# Assign derived states (use an array of chars)
ancestral_states = tables.sites.ancestral_state.view(dtype='|S1')
mutation_derived_state = np.full(tables.mutations.num_rows, null_char_allele, dtype='|S1')
for i, m in enumerate(tables.mutations):
    if m.parent == tskit.NULL:
        prev_state = ancestral_states[m.site]
    else:
        prev_state = mutation_derived_state[m.parent]
    # Pick a new state that is different from the old one
    new_state = random.choice([s for s in states if s != prev_state])
    mutation_derived_state[i] = new_state

tables.mutations.derived_state = mutation_derived_state.view(tables.mutations.derived_state.dtype)
finite_sites_ts = tables.tree_sequence()

# Try printing them out
#for h in finite_sites_ts.haplotypes():
#  print(len(h))

truncated_ts = truncate_ts_samples(finite_sites_ts, average_span=200, random_seed=123)
#truncated_ts_2 = truncate_ts_samples(ts, average_span=200, random_seed=123)

#print(truncated_ts_2.genotype_matrix()) # doesn't break
#print(truncated_ts.genotype_matrix()) # doesn't break with update

#for h in truncated_ts.haplotypes():
#  print(h)

# broken on truncated seqs with mutations on
#print(truncated_ts.genotype_matrix())

######## fasta file output - now works for missing (truncated) data ########
haps = []
for i in truncated_ts.haplotypes():
    haps.append(i)
sys.exit(0) #testing missing data part

sequence_IDs = []
for i in range(len(haps)):
    sequence_IDs.append(f'sample_{ts.samples()[i]}_pop_{ts.node(i).population}')

#fasta printing
#for i in range(len(haps)):
#    print(f'>{sequence_IDs[i]}\n{haps[i]}')

############# Saving things ############
sys.exit(0) # to not overwrite files unless intentional
#file_loc = "C:\\Users\\mazmysta\\OneDrive - Nexus365\\BDI_proj\\scripts\\PopulationsHIV\\fasta\\"
write_loc = os.path.join(base_loc, 'simulation2')

## write and save fasta file
with open(os.path.join(write_loc, 'HIV_truncated_sim_fasta_file.txt'), 'w') as f:
    for i in range(len(haps)):
        f.write(f'>{sequence_IDs[i]}\n{haps[i]}\n')

haps_full = []
for i in finite_sites_ts.haplotypes():
    haps_full.append(i)

with open(os.path.join(write_loc, 'HIV_fullsim_fasta_file.txt'), 'w') as f:
    for i in range(len(haps_full)):
        f.write(f'>{sequence_IDs[i]}\n{haps_full[i]}\n')

#file_loc = "C:\\Users\\mazmysta\\OneDrive - Nexus365\\BDI_proj\\scripts\\PopulationsHIV\\fasta\\"
#another_test = []
#
#for i in ts.haplotypes():
#    another_test.append(i)
#
#sequence_IDs = []
#for i in range(len(another_test)):
#    sequence_IDs.append(f'sample_{ts.samples()[i]}_pop_{ts.node(i).population}')
#
#with open(file_loc + 'test_fasta_file.txt', 'w') as f:
#    for i in range(len(another_test)):
#        f.write(f'>{sequence_IDs[i]}\n{another_test[i]}\n')

### save both .trees files !! ###

write_loc = os.path.join(base_loc, 'simulation2')

truncated_ts.dump(write_loc + 'truncated_simulation_tree.trees')
finite_sites_ts.dump(write_loc + 'full_simulation_tree.trees')

######### Pickling relevant files for sim #############

#how to save
with open(os.path.join(write_loc, 'pickled_pop_list.pickle'), 'wb') as f:
    pickle.dump(pop_list, f)

with open(os.path.join(write_loc, 'pickled_events_sorted.pickle'), 'wb') as f:
    pickle.dump(events_sorted, f)

with open(os.path.join(write_loc, 'pickled_M.pickle'), 'wb') as f:
    pickle.dump(M, f)

with open(os.path.join(write_loc, 'pickles_root_kids_list'), 'wb') as f:
    pickle.dump(list_of_root_and_kids, f)

## loading example
#with open(file_loc + 'pickled_pop_list.pickle', 'rb') as f:
#    pop_list = pickle.load(f)

########### inference on truncation part ###########

# ts infer

sd = tsinfer.SampleData.from_tree_sequence(truncated_ts, use_times=False)

ts_inferred = tsinfer.infer(sd, simplify=False)

ts_inferred = ts_inferred.simplify(filter_sites=False, keep_unary=True)

ts_inferred
ts_inferred.dump(os.path.join(write_loc, 'inferred_tree.trees'))
#Out[43]: <tskit.trees.TreeSequence at 0x1ed55ca9710>

X = 18
Y = 20
i = 0
for tree in ts_inferred.trees():
    if i > X and i <= Y:
        display(SVG(tree.draw(height=600, width=4000,
                              tree_height_scale='rank')))
        print("Tree {} covers [{:.2f}, {:.2f}); TMRCA = {:.4f}".format(
            tree.index, *tree.interval, tree.time(tree.roots[0])))
#        print(tree.branch_length(46))
    elif i > Y:
        break
    i+=1

print(ts_inferred.genotype_matrix())

