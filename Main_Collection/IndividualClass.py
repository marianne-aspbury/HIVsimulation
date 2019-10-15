# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 13:08:31 2019

Just the class of individual

@author: marianne aspbury
"""

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