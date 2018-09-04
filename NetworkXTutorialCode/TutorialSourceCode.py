import matplotlib.pyplot as plt
import networkx as nx




#### Adding nodes
g = nx.Graph ()
g.add_node(1)
g.add_nodes_from ( ['b','c','d'])
g.add_nodes_from ('xyz')


#### Adding edges

g = nx. Graph( [('a','b') ,('b','c') ,('c','a')] )
g.add_edge('a', 'd')
g.add_edges_from([( 'd', 'c'), ('d', 'b')])

#### Add node attributes
g = nx.Graph()
g.add_node (1, name ='Obrian')
g.add_nodes_from([2], name ='Quintana')
g.node[1]['name']
g.node[1] # Python dictionary

#### Add edge attributes
g.add_edge(1, 2, weight=4.0 )
g[1][2]['weight'] = 5.0  # edge already added
g[1][2]
g.add_edges_from([(3 ,4) ,(4 ,5)], w =3.0)
g.add_weighted_edges_from([(6 ,7 ,3.0)]) # adds third value in tuple as ‘weight ’ attr
g[6][7]
g.add_edge(5 ,6)
g[5][6]

#### Number of nodes
g = nx.Graph()
g.add_edges_from([(1 ,2) ,(1 ,3)])
g.add_node('a')
len(g)
g. number_of_nodes()
g. order ()

#### Number of edges
g.number_of_edges()

g.nodes()
g.edges()


#  access the adjacency dictionary as G[n]:
g[1]  
g.degree(1)


## node and edge iterators
g = nx.Graph()
g.add_edges_from([(1 ,2) ,(1 ,3)])
g.add_node('a')

for node in g.nodes():
    print(node, g.degree(node))

## edge iterators
g.add_edge(1,3,weight=2.5)
g[1][2]['weight'] = 1.5
for n1,n2,attr in g.edges(data=True):
    print(n1,n2,attr['weight'])

## Deal with weighted graphs
dg = nx.DiGraph()
dg.add_weighted_edges_from([(1,4,0.5), (3,1,0.75)])
dg.out_degree(1,weight='weight')
dg.degree(1,weight='weight')
dg.successors(1)
dg.predecessors(1)

## Deal with multigraphs
MG = nx.MultiGraph()
MG.add_weighted_edges_from([(1, 2, 0.5), (1, 2, 0.75), (2, 3, 0.5)])
dict(MG.degree(weight='weight'))
## To calcuate the shortest path, we have to convert it to a standard graph 
GG = nx.Graph()
for n, nbrs in MG.adjacency():
    for nbr, edict in nbrs.items():
        minvalue = min([d['weight'] for d in edict.values()])
        GG.add_edge(n, nbr, weight = minvalue)
nx.shortest_path(GG, 1, 3)


## Generating graphs
Gc = nx.complete_graph(5)
nx.draw(Gc)

# Chain (path) graph
Gp = nx.path_graph(5)
nx.draw(Gp)

# Grid graph
Gg = nx.grid_graph ([4 ,4 ,4]) # 3D, 4^3 nodes
nx.draw(Gg)

# Erdos-Reny random graphs
n = 50
p = 0.1
Gr = nx.erdos_renyi_graph(n,p)
nx.draw(Gr)

# Preferential Attachment networks
n = 50
m = 4
Gpa = nx.barabasi_albert_graph(n, m)
nx.draw(Gpa)

# Small world networks
n = 50
k = 6
p = 0.3
Gsw = nx.watts_strogatz_graph(n, k, p)
nx.draw(Gsw)


### Graph operations
# Induced subgraph view of G on nodes in nbunch
G = nx.path_graph(4)
nx.draw(G)
H = G.subgraph([0, 1, 2])
nx.draw(H)


# Graph union 
G1 = nx.path_graph(4)
G2 = nx.complete_graph(5)
U = nx.disjoint_union(G1,G2)   
nx.draw(U)

# Return Cartesian product graph
CP = nx.cartesian_product(G1,G2)
nx.draw(CP)

# Combine graphs identifying nodes common to both
G1 = nx.Graph()
G1.add_edges_from([(1 ,2) ,(1 ,3)])
G2 = nx.Graph()
G2.add_edges_from([(3 ,4) ,(4,5)])
CO = nx.compose(G1,G2)
nx.draw(CO)    

# graph complement 
Gc = nx.complete_graph(5)
nx.draw(Gc)
ebunch=[(1, 2), (2, 3)]
Gc.remove_edges_from(ebunch)
nx.draw(Gc)
GP = nx.complement(Gc)
nx.draw(GP)

#### Python dictionaries
# Keys and values can be of any data type
fruit_dict={"apple":1,"orange":[0.23,0.11],42:True}
# Can retrieve the keys and values as Python lists (vector)
fruit_dict.keys()
fruit_dict.values()
# Or return tuples of (key,value)
fruit_dict.items()

# To iterate:
for k, v in fruit_dict.items():
    print(k, v)


#### Drawing networks
G_florentine = nx.florentine_families_graph()
nx.draw(G_florentine, with_labels=True)

# We can use graph layouts to set node positions for graph drawing.
pos = nx.spring_layout(G_florentine)
nx.draw(G_florentine, pos, with_labels=True)

pos = nx.circular_layout(G_florentine)
nx.draw(G_florentine, pos, with_labels=True)


## Advanced graph drawing samples
G = nx.random_geometric_graph(200, 0.125)
# position is stored as node attribute data for random_geometric_graph
pos = nx.get_node_attributes(G, 'pos')

# find node near center (0.5,0.5)
dmin = 1
ncenter = 0
for n in pos:
    x, y = pos[n]
    d = (x - 0.5)**2 + (y - 0.5)**2
    if d < dmin:
        ncenter = n
        dmin = d

# color by path length from node near center
p = dict(nx.single_source_shortest_path_length(G, ncenter))

plt.figure(figsize=(8, 8))
nx.draw_networkx_edges(G, pos, nodelist=[ncenter], alpha=0.4)
nx.draw_networkx_nodes(G, pos, nodelist=list(p.keys()),
                       node_size=80,
                       node_color=list(p.values()),
                       cmap=plt.cm.Reds_r)

plt.xlim(-0.05, 1.05)
plt.ylim(-0.05, 1.05)
plt.axis('off')
plt.show()
























