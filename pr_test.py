import os
import networkx as nx
import matplotlib.pyplot as plt
os.chdir('/Users/xueyiheng/Desktop')
filename = 'cit-HepTh_3.txt'
G=nx.DiGraph()
with open(filename) as file:
    for line in file:
        head, tail = [float(x) for x in line.split()]
        G.add_edge(head,tail)
pr=nx.pagerank(G,alpha=0.85)
x = 0;
for node, value in pr.items():
    x = x + value
print(x)

G=nx.binomial_graph(1000, 0.3, directed=True)
layout = nx.spring_layout(G)
plt.figure(1)
nx.draw(G, pos=layout, node_color='y')

pr=nx.pagerank(G,alpha=0.85)
print(pr)
for node, pageRankValue in pr.items():
    print("%d,%.4f" %(node,pageRankValue))

plt.figure(2)
nx.draw(G, pos=layout, node_size=[x * 6000 for x in pr.values()],node_color='m',with_labels=True)
plt.show()