# -*- coding: utf-8 -*-
"""
Created on Sun Jul 15 13:46:39 2018

@author: Administrator
"""

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter


hartford = nx.read_edgelist('hartford_drug.txt', create_using=nx.DiGraph(),nodetype=int)
N,K = hartford.order(), hartford.size()
nx.draw_spring(hartford)


##### Degree distribution
# This method returns a dictionary view object
in_degrees = hartford.in_degree() # A DegreeView object capable of iterating (node, degree) pairs)
# We can iterate Dictionary views to retrieve their respective data
# You can regard views as the data 
for d in in_degrees:
    print(d)
# We can list node degree by using the Python list function    
list(in_degrees([1, 2]))

# We first sort the degrees
in_values_sorted = sorted(d for n, d in in_degrees)
print(in_values_sorted)
# Then count the number of degrees
counts = Counter(d for n, d in in_degrees())
 # unpack a list of pairs into two tuples
labels, values = zip(*counts.items())
# Prepare for bar plot by evenly spacing values 
# within a given interval using numpy arrange()
indexes = np.arange(len(labels))
width = 1
plt.bar(indexes, values, width)
plt.xticks(indexes + width * 0.5, labels)
plt.show()


##### Degree centrality
in_degree_c = nx.in_degree_centrality(hartford)
# Centrality methods return a dictionaries type!!!
# You cannot print out the centrality value using the following loop
for c in in_degree_c:
    print(c)
    
for k, c in in_degree_c.items():
    print(k, c)
    
# Sorting is a bit complex in Python 3:
# https://stackoverflow.com/questions/613183/how-do-i-sort-a-dictionary-by-value    
# We need to use lambda function:
# https://www.python-course.eu/lambda.php
sorted_by_value = sorted(in_degree_c.items(), key=lambda kv: kv[1])

# You can also calculate out degree centrality 
c_out_degree = nx.out_degree_centrality(hartford)


