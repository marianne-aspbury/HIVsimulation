# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 09:03:53 2019

Reading in csv transmission network data and saving to Python Data Structure
And enable manipulation of this info in python

Contains successful implementation of overall children counting for each person

@author: mazmysta
"""

## some tools for later
from IndividualClass import Individual
import numpy as np 
import csv

'''
# Using this class to store information on individuals
class Individual:
    def __init__(self):
        self.ID = None
        self.infected_others = []
        self.infected_by = []
        self.time_infected_others =[]
        self.time_infected_by = []
        self.age_birth = -1
        self.age_death = -1
        self.num_partners = -1
        self.direct_children_num = 0
        self.total_children_num = 0
        '''
    

## Some tests to understand functionality 
#my_dict = {}
#print(my_dict)
#ID_1 = 2135
#
#if ID_1 in my_dict:
#    print(my_dict)
#else:
#    my_dict[str(ID_1)] = Individual()
#
#print(my_dict[str(ID_1)].ID)
#
#my_dict[str(ID_1)].ID = ID_1
#print(my_dict[str(ID_1)].ID)
#
#indiv1 = Individual()
#indiv1.ID = 2135
#indiv1.infected_by.append(2134)
#
#print(indiv1)
#print(indiv1.ID)
#print(indiv1.infected_by)


### File locations to read in files 
file_loc = 'C:\\Users\\mazmysta\OneDrive - Nexus365\\BDI_proj\\HIV_Transmission_Networks_WillProbert\\19-08-08-first_example_network\\'
individual_metadata_file = 'phylogenetic_individualdata_all_100314_CL01_Za_B_V1.2_patchall_Rand10_Run1_PCseed0_0.csv'
transmission_datafile = 'phylogenetic_transmission_all_100314_CL01_Za_B_V1.2_patchall_Rand10_Run1_PCseed0_0.csv'

### reading in IDs of all Individuals and creating class instances, saving to dictionary for access
my_dict={}
with open(file_loc + individual_metadata_file, 'r') as f:
     reader = csv.reader(f, delimiter=',')
#     for row in islice(reader,2000,25000):
     for row in reader:
#         if firstline:
#             firstline = False
#             print(row)
#             continue
         if row[6] == '1': ## We are only building a list of HIV infected individuals
#             print('yes')
             my_dict[row[1]] = Individual()
             my_dict[row[1]].ID = row[1]
             my_dict[row[1]].age_birth = row[4]
             my_dict[row[1]].age_death = row[5]
#             print(row[1])
#             print(my_dict[row[1]].__dict__)
         
#example class instance from dictionary
print(my_dict['2369_0'].__dict__)

## Find who 's infected who and store infected by & infected from info.
with open(file_loc + transmission_datafile, 'r') as f:
     reader = csv.reader(f, delimiter=',')
     for row in reader:
         if row[2] in my_dict: # Want to see if infector ID is in the dictionary
#             print('Found Infected individual:', row[2], '. Infector:', row[3])
#             print('Time of infection:', row[6])
#             print('Num of infector partners:', row[-1])
             my_dict[row[2]].infected_by.append(row[3])
             my_dict[row[2]].time_infected_by.append(row[6])
             if row[3] in my_dict:
#                 my_dict[row[3]].ID = row[3]
                 pass
             else:
                 print('extra is', row[3]) ## there are two seed cases noted with IDs starting -1
                 my_dict[row[3]] = Individual()
                 my_dict[row[3]].ID = row[3]
             my_dict[row[3]].infected_others.append(row[2])
             my_dict[row[3]].time_infected_others.append(row[6])
             my_dict[row[3]].num_partners = row[-1]

# proof of function
print(my_dict['2369_0'].__dict__)
print(my_dict[my_dict['2369_0'].infected_by[0]].__dict__)
print(my_dict['8876_1'].__dict__)

print(len(my_dict))
#49011 : 49009 'true' cases and two seeds '-1_1' and '-1_2' infectors.
for key in my_dict:
    if my_dict[key].infected_by == []:
        print(key)

## find out the child number - direct and overall - for each person
## direct child number is just by counting num IDs in infected_others ##
for key in my_dict:
#    print(my_dict[key].__dict__)
#    print(len(my_dict[key].infected_others))
    my_dict[key].direct_children_num = len(my_dict[key].infected_others)

#%% 
## find out overall child number by iteratively summing direct_children_num of infected_others and their children
## but need to be careful that start accumulating from the bottom of the trees upwards - visit all bottom children before any parents etc.
## post-order tree traversal. Not sure how to go about this, will think

# reset, can += 1 to add to total_children_num for each true child, reset if wrong / repeat
for key in my_dict:
    my_dict[key].total_children_num = 0
    
# just something I was doing to find tips of tree and their direct parents but prob don't want this?
tips_zero_children = []
parents_of_zero_children = []

for key in my_dict:
    if my_dict[key].direct_children_num == 0:
        #print(key)
        tips_zero_children.append(key)
        parents_of_zero_children.extend(my_dict[key].infected_by)
        
#%% TEST CASE - simple tree - calculating overall children
        '''
        ---------0---------
        |                 |
      --1--------         12
      |         |
    --2--     --5-----
    |   |     |      |
    3   4     6    --7--
                   |    |
                 --8--  11
                 |   |
                 9   10
Should have 3=0, 4=0, 2=2, 
            9=0, 10=0, 8=2, 11=0, 7=4, 6=0, 5=6
            1 = 10
            12 = 0
            0 = 11
'''

# setting up test case
test_dict = {}
for i in range(13):
    test_dict[str(i)] = Individual()
    test_dict[str(i)].ID = str(i)
    #print(str(i))

zeroth = test_dict['0'] 
zeroth.infected_others = ['1', '12']
test_dict['1'].infected_others = ['2', '5']
test_dict['2'].infected_others = ['3', '4']
test_dict['3'].infected_others = []
test_dict['4'].infected_others = []
test_dict['5'].infected_others = ['6', '7']
test_dict['6'].infected_others = []
test_dict['7'].infected_others = ['8', '11']
test_dict['8'].infected_others = ['9', '10']
test_dict['9'].infected_others = []
test_dict['10'].infected_others = []
test_dict['11'].infected_others = []
test_dict['12'].infected_others = []

test_dict['1'].infected_by = ['0']
test_dict['2'].infected_by = ['1']
test_dict['3'].infected_by = ['2']
test_dict['4'].infected_by = ['2']
test_dict['5'].infected_by = ['1']
test_dict['6'].infected_by = ['5']
test_dict['7'].infected_by = ['5']
test_dict['8'].infected_by = ['7']
test_dict['9'].infected_by = ['8']
test_dict['10'].infected_by = ['8']
test_dict['11'].infected_by = ['7']
test_dict['12'].infected_by = ['0']

print(test_dict['1'].__dict__)

test_tree_tips = []
for i in test_dict:
    if test_dict[i].infected_others == []:
        test_tree_tips.append(i)
        #print(i)
print(test_tree_tips)
#['3', '4', '6', '9', '10', '11', '12']

accounted_for = []
parents_working_list = ['0']
for i in test_tree_tips:
    start_tip = '0' # root of the tree
    for j in test_dict[start_tip].infected_others:
        child_visit = j
        print(j)
        if len(test_dict[j].infected_others) > 0:
            to_visit = test_dict[j].infected_others
            parents_working_list.append(j)
            print(parents_working_list)
            print(to_visit)
        else:
            infector = test_dict[j].infected_by[0]
            test_dict[infector].total_children_num += 1
            accounted_for.append(j)
            print(accounted_for)
            print(infector)
            print(test_dict[infector].total_children_num)

for x in parents_working_list:
    if len(test_dict[x].infected_others) > 0: #i.e. has children
        for y in test_dict[x].infected_others:
            if len(test_dict[y].infected_others) > 0: #i.e. has children
                parents_working_list.append(y)

  
for i in range(13):
    test_dict[str(i)].total_children_num = 0     
         
x ='0'
accounted_for = []
tot_parents = []
i = 0
while len(accounted_for) < len(test_dict):
    print('while')
    if len(test_dict[x].infected_others) > 0 and not set(test_dict[x].infected_others).issubset(accounted_for): #i.e. x has kids
        print('len > 0, 1')
        temp_parent = x
        tot_parents.append(x)
        print(temp_parent, tot_parents)
        for kid in test_dict[x].infected_others:
            if kid in accounted_for:
                print(kid, 'accounted')
#                pass
            else:
                print(kid, 'not accounted')
                a = test_dict[x].infected_others.index(kid)
                x = test_dict[x].infected_others[a]
                break
    elif len(test_dict[x].infected_others) == 0:
        print('len = 0')
        for parent in tot_parents:
            test_dict[parent].total_children_num += 1
            print('x', x, 'parent', parent, 'tot_kids', test_dict[parent].total_children_num)
        accounted_for.append(x)
        tot_parents.remove(test_dict[x].infected_by[0])
        print('tot parents', tot_parents)
        x = test_dict[x].infected_by[0]
    elif set(test_dict[x].infected_others).issubset(accounted_for):
        print('check tot parents. x is', x, 'tot parents is', tot_parents)
        if tot_parents == []: #check not at the root, if at the root need to stop.
            break
        for parent in tot_parents:
            test_dict[parent].total_children_num += 1
            print('x', x, 'parent', parent, 'tot_kids', test_dict[parent].total_children_num)
        accounted_for.append(x)
        tot_parents.remove(test_dict[x].infected_by[0])
        print('tot parents', tot_parents)
        x = test_dict[x].infected_by[0]
#    if accounted_for == ['3', '4']:
#        if i > 8:
#            print('i is', i)
#            print('acc for', accounted_for)
#            break
#    print(x, 'x')
#    if i > 8:
#        print('acc for', accounted_for)
#        break
#    i += 1
#            
for i in range(13):
    print(i, test_dict[str(i)].total_children_num)

## this works!! :) try for true dict
    
#%% TRUE CSV TRANSMISSION NETWORK - calculating overall children
    ## For actual dict
i = 0
for key in my_dict:
    print(my_dict[key].total_children_num)
    i+=1
    if i >10:
        break

# reset, can += 1 to add to total_children_num for each true child, reset if wrong / repeat
for key in my_dict:
    my_dict[key].total_children_num = 0

# two seeds '-1_1' and '-1_2' infectors.
for key in my_dict:
    if my_dict[key].infected_by == []:
        print(key)

#1st seed         
x ='-1_-1'
accounted_for = []
tot_parents = []
i = 0
while len(accounted_for) < len(my_dict):
    #print('while')
    if len(my_dict[x].infected_others) > 0 and not set(my_dict[x].infected_others).issubset(accounted_for): #i.e. x has kids
        #print('len > 0, 1')
        temp_parent = x
        tot_parents.append(x)
        #print(temp_parent, tot_parents)
        for kid in my_dict[x].infected_others:
            if kid in accounted_for:
                #print(kid, 'accounted')
                pass
            else:
                #print(kid, 'not accounted')
                a = my_dict[x].infected_others.index(kid)
                x = my_dict[x].infected_others[a]
                break
    elif len(my_dict[x].infected_others) == 0:
        #print('len = 0')
        for parent in tot_parents:
            my_dict[parent].total_children_num += 1
            #print('x', x, 'parent', parent, 'tot_kids', my_dict[parent].total_children_num)
        accounted_for.append(x)
        tot_parents.remove(my_dict[x].infected_by[0])
        #print('tot parents', tot_parents)
        x = my_dict[x].infected_by[0]
    elif set(my_dict[x].infected_others).issubset(accounted_for):
        #print('check tot parents. x is', x, 'tot parents is', tot_parents)
        if tot_parents == []: #check not at the root, if at the root need to stop.
            break
        for parent in tot_parents:
            my_dict[parent].total_children_num += 1
            #print('x', x, 'parent', parent, 'tot_kids', my_dict[parent].total_children_num)
        accounted_for.append(x)
        tot_parents.remove(my_dict[x].infected_by[0])
        #print('tot parents', tot_parents)
        x = my_dict[x].infected_by[0]

#2nd seed        
x ='-1_2'
accounted_for = []
tot_parents = []
i = 0
while len(accounted_for) < len(my_dict):
    #print('while')
    if len(my_dict[x].infected_others) > 0 and not set(my_dict[x].infected_others).issubset(accounted_for): #i.e. x has kids
        #print('len > 0, 1')
        temp_parent = x
        tot_parents.append(x)
        #print(temp_parent, tot_parents)
        for kid in my_dict[x].infected_others:
            if kid in accounted_for:
                #print(kid, 'accounted')
                pass
            else:
                #print(kid, 'not accounted')
                a = my_dict[x].infected_others.index(kid)
                x = my_dict[x].infected_others[a]
                break
    elif len(my_dict[x].infected_others) == 0:
        #print('len = 0')
        for parent in tot_parents:
            my_dict[parent].total_children_num += 1
            #print('x', x, 'parent', parent, 'tot_kids', my_dict[parent].total_children_num)
        accounted_for.append(x)
        tot_parents.remove(my_dict[x].infected_by[0])
        #print('tot parents', tot_parents)
        x = my_dict[x].infected_by[0]
    elif set(my_dict[x].infected_others).issubset(accounted_for):
        #print('check tot parents. x is', x, 'tot parents is', tot_parents)
        if tot_parents == []: #check not at the root, if at the root need to stop.
            break
        for parent in tot_parents:
            my_dict[parent].total_children_num += 1
            #print('x', x, 'parent', parent, 'tot_kids', my_dict[parent].total_children_num)
        accounted_for.append(x)
        tot_parents.remove(my_dict[x].infected_by[0])
        #print('tot parents', tot_parents)
        x = my_dict[x].infected_by[0]

## inspect frequencies of diff numbers of kids
a = []
for item in my_dict:
    a.append(my_dict[item].total_children_num)
d = np.unique(a, return_counts=True)
#print('all times, frequency of diff pops', d)
print('requency of diff overall kids below', '\n', 'freqs:', d[0], 'counts:', d[1])
print(len(d[1]))

## compare numbers of individual children:
b = []
for item in my_dict:
    b.append(my_dict[item].direct_children_num)
e = np.unique(b, return_counts=True)
#print('all times, frequency of diff pops', d)
print('frequency of diff direct kids below', '\n', 'number:', e[0], '\n', 'counts:', e[1])
print(len(e[1]))
    
for key in my_dict:
    if my_dict[key].infected_by == []:
        print(key)
        print(my_dict[key].total_children_num)

## that did seem to work.... freq of total_child_num == 0 is same as freq direct_children_num == 0
        
'''class Individual:
    def __init__(self):
        self.ID = None
        self.infected_others = []
        self.infected_by = []
        self.time_infected_others =[]
        self.time_infected_by = []
        self.age_birth = -1
        self.age_death = -1
        self.num_partners = -1
        self.direct_children_num = 0
        self.total_children_num = 0'''
        
#%%  SAVE CALCULATED CLASS STRUCTURE AND INFO
# Pickle this custom class data to just load in in future
import pickle

#how to save
with open(file_loc + 'pickled_data_all.pickle', 'wb') as f:
    pickle.dump(my_dict, f)

# loading example    
with open(file_loc + 'pickled_data_all.pickle', 'rb') as f:
    cc = pickle.load(f)

#proof it works
print('Pickle loaded, 1st root num kids: {}'.format(cc['-1_-1'].total_children_num))

#further proof
ct = 0
for key in cc:
    if ct < 5:
        print(cc[key])
        print(cc[key].__dict__)
    else:
        break
    ct += 1
    

