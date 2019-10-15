# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 11:09:42 2019

Example for slides of truncation

@author: marianne aspbury
"""

import tskit
import msprime
import tsinfer

import sys
import numpy as np
from IPython.display import display
from IPython.display import SVG

import pickle
import random # to pick from list of poss roots of subtree
import math
import itertools

#from IndividualClass import Individual # might say unused, but not true, needed to read pickled file
from truncate_ts_samples import truncate_ts_samples as truncate_ts_samples

from Bio import SeqIO

ts_eg = msprime.simulate(10,
                         mutation_rate=1e-8,
                         recombination_rate=1e-10,
                         length=100, Ne=100,
                         random_seed = 17)

truncated_ts_eg = truncate_ts_samples(ts_eg, average_span=50, random_seed=123)
truncated_ts_eg2 = truncated_ts_eg.simplify()

for tree in truncated_ts_eg2.trees():
    display(SVG(tree.draw(height=200, width = 300,
                      format='SVG',
                      tree_height_scale='rank'
                      #max_tree_height=5000
                      )))
    print("Tree {} covers [{:.2f}, {:.2f})".format(
    tree.index, *tree.interval, tree.time(tree.roots[0])))
    #        print(tree.branch_length(46))