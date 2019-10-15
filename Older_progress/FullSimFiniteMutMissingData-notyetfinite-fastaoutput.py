# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 12:33:36 2019

Combined new simulation using missing data and finite mutation model

@author: marianne aspbury
"""

import msprime
import tsinfer
import numpy as np
from IPython.display import display
from IPython.display import SVG

import pickle
import random # to pick from list of poss roots of subtree
from IndividualClass import Individual # might say unused, but not true, needed to read pickled file
from truncate_ts_samples import truncate_ts_samples as truncate_ts_samples

### import the trees ###
file_loc = 'C:\\Users\\mazmysta\OneDrive - Nexus365\\BDI_proj\\HIV_Transmission_Networks_WillProbert\\19-08-08-first_example_network\\'

with open(file_loc + 'pickled_data_all.pickle', 'rb') as f:
    total_tree = pickle.load(f)

#proof it works
print('Pickle loaded, 1st root num kids: {}'.format(total_tree['-1_-1'].total_children_num))

# Saving list of people with desired num kids beneath
samples_poss = []
desired_overall_children = 5 ## Seemingly my max is 30 with 10 samples, 100 pop size and 1000 length genome.

for key in total_tree:
    if total_tree[key].total_children_num == desired_overall_children:
        samples_poss.append(key)

print(len(samples_poss))
# 8 for kids = 100

#pick only one random sample
random.seed(4) # repeatable choice(s) - sequence same
start = random.choice(samples_poss)
#start='58582_0'
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

#file_loc = 'C:\\Users\\mazmysta\\OneDrive - Nexus365\\BDI_proj\\scripts\\PopulationsHIV\\'
#
##how to save
#with open(file_loc + 'pickled_pop_list.pickle', 'wb') as f:
#    pickle.dump(pop_list, f)
#
#with open(file_loc + 'pickled_events_sorted.pickle', 'wb') as f:
#    pickle.dump(events_sorted, f)
#
#with open(file_loc + 'pickled_M.pickle', 'wb') as f:
#    pickle.dump(M, f)
#
## loading example
#with open(file_loc + 'pickled_pop_list.pickle', 'rb') as f:
#    pop_list = pickle.load(f)
#
## simulate this extended simple model
ts2 = msprime.simulate(population_configurations=pop_list, migration_matrix = M,
                       demographic_events = events_sorted,
                       length = 1000,
                       random_seed = 17, recombination_rate = 0.7e-4,
                       mutation_rate=2e-5, end_time=60000)


## simulate this extended simple model
#ts2 = msprime.simulate(population_configurations=pop_list, migration_matrix = M,
#                       demographic_events = events_sorted,
#                       length = 1000,
#                       random_seed = 17, recombination_rate = 0.7e-4)

# 5 pops
colour_map = {0:"grey", 1:"blue", 2:"red", 3:"orange", 4:"green", 5:"cyan", 6:"black"}
node_colours = {u.id: colour_map[u.population] for u in ts2.nodes()}


X = 1
i = 0
for tree in ts2.trees():
    if i < X:
        display(SVG(tree.draw(node_colours = node_colours,
                              height=800, width = 1000,
                              format='SVG',
                              tree_height_scale='log_time')))
        print("Tree {} covers [{:.2f}, {:.2f}); TMRCA = {:.4f}".format(
            tree.index, *tree.interval, tree.time(tree.roots[0])))
        print(tree.branch_length(46))
    else:
        break
    i+=1

#### truncation part

# Having run
# ReadingTransmissionNetworkCSV_OverallChildrenCalcs_save.py
# have a ts2 = simulated data


# ts2 from subsetting to msprime, length 100, picking just one population below
# simplification and subsetting

#ts3 = ts.simplify(samples=[3,8,0,4,7,2,6,5,1,9]
#)

truncated_ts3 = truncate_ts_samples(ts2, average_span=200, random_seed=123)

#print(truncated_ts3)
#<tskit.trees.TreeSequence object at 0x000001ED52232710>

#truncated_ts3.tables
#Out[34]: <tskit.tables.TableCollection at 0x1ed4f6c1860>

#truncated_ts3.tables.nodes
#Out[35]: <tskit.tables.NodeTable at 0x1ed52232908>

#print(truncated_ts3.tables)
#prints the tables - long
X = 10
i = 0
for tree in truncated_ts3.trees():
    if i < X:
        display(SVG(tree.draw(height=800, width = 1000,
                              format='SVG',
                              tree_height_scale='log_time')))
        print("Tree {} covers [{:.2f}, {:.2f}); TMRCA = {:.4f}".format(
            tree.index, *tree.interval, tree.time(tree.roots[0])))
        print(tree.branch_length(46))
    else:
        break
    i+=1

truncated_ts3.dump(file_loc + 'truncated_simulation_tree.trees')
# ts infer

sd = tsinfer.SampleData.from_tree_sequence(truncated_ts3, use_times=False)

ts_inferred = tsinfer.infer(sd, simplify=False)

ts_inferred = ts_inferred.simplify(filter_sites=False, keep_unary=True)

ts_inferred
ts_inferred.dump(file_loc + 'inferred_tree.trees')
#Out[43]: <tskit.trees.TreeSequence at 0x1ed55ca9710>

X = 10
i = 0
for tree in ts_inferred.trees():
    if i < X:
        display(SVG(tree.draw(height=400, width=600,
                              tree_height_scale='rank')))
        print("Tree {} covers [{:.2f}, {:.2f}); TMRCA = {:.4f}".format(
            tree.index, *tree.interval, tree.time(tree.roots[0])))
#        print(tree.branch_length(46))
    else:
        break
    i+=1

print(ts_inferred.genotype_matrix())
print(truncated_ts3.genotype_matrix())

haps = []
for i in ts2.haplotypes():
    haps.append(i)

sequence_IDs = []
for i in range(len(haps)):
    sequence_IDs.append(f'sample_{ts2.samples()[i]}_pop_{ts2.node(i).population}')

#fasta printing
#for i in range(len(haps)):
#    print(f'>{sequence_IDs[i]}\n{haps[i]}')

## write and save fasta file
with open(file_loc + 'test_fasta_file.txt', 'w') as f:
    for i in range(len(haps)):
        f.write(f'>{sequence_IDs[i]}\n{haps[i]}\n')

