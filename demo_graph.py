#python -m venv myenv    /// In terminal,  create the virtual env  
#source myenv/bin/activate        ///Activate the virtual environment
# deactivate

#graphml_file = '/mnt/d/M2/NET/demo_file/reactome-77-reaction_R-HSA-5696021.graphml'

#g = ig.Graph.Read_GraphML(graphml_file)

#layout = g.layout("auto")
#ig.plot(g, layout=layout, target="output_graph2.png")

#print(f"nb de arret: {len(g.vs)}")
#print(f"nb de sommet: {len(g.es)}")

# pip install networkx, pip install matplotlib

import igraph as ig

graph = ig.Graph(directed=True)


graph.add_vertices(10)  
graph.vs["name"] = ["Glycolysis", "", "ADPGK:Mg2+", "ADPGK:Mg2+ phosphorylates Glc to G6P",
                    "Mg2+", "ADPGK", "G6P", "Glc", "ADP", "AMP"]
graph.vs["label"] = graph.vs["name"]

edges = [
    (0, 3),  # Glycolysis -> ADPGK:Mg2+ phosphorylates Glc to G6P
    (1, 2),  # vide -> ADPGK:Mg2+
    (1, 3),  # vide -> ADPGK:Mg2+ phosphorylates Glc to G6P
    (2, 4),  # ADPGK:Mg2+ -> Mg2+
    (2, 5),  # ADPGK:Mg2+ -> ADPGK
    (3, 6),  # ADPGK:Mg2+ phosphorylates Glc to G6P -> G6P
    (3, 7),  # ADPGK:Mg2+ phosphorylates Glc to G6P -> Glc
    (3, 8),  # ADPGK:Mg2+ phosphorylates Glc to G6P -> ADP
    (3, 9)   # ADPGK:Mg2+ phosphorylates Glc to G6P -> AMP
]
graph.add_edges(edges)

#ig.plot(graph, target="reconstructed_graph.png")

import networkx as nx
import matplotlib.pyplot as plt

graphml_file = '/mnt/d/M2/NET/demo_file/reactome-77-reaction_R-HSA-5696021.graphml'

G = nx.read_graphml(graphml_file)


# Mapping node IDs to tag names
new_node_names = {}  # Create an empty dictionary for storing new node names

for node, data in G.nodes(data=True):
    if 'label' in data:  
        new_node_names[node] = data['label']  

G = nx.relabel_nodes(G, new_node_names)

# Create a dictionary for displaying the labels of the edges during plotting
edge_labels = {}
for u, v, data in G.edges(data=True):
    edge_labels[(u, v)] = data.get('label', '')

# plot
plt.figure(figsize=(10, 8))
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True, node_size=500, node_color="lightblue", font_size=10, font_weight="bold")
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8, font_color="black")
plt.show()

# Selecting Start and End Label Names
source = "Glycolysis"  # Starting Point 
target = "G6P"         # end point 

# Check the path and run Dijkstra's algorithm
if source in G and target in G:
    if nx.has_path(G, source, target):
        shortest_path = nx.shortest_path(G, source=source, target=target, weight='weight', method='dijkstra')
        print(f"de  {source} à {target} le plus court chemin: {shortest_path}")
    else:
        print(f"de {source} à {target} Pas de chemin")
else:
    print("can't find the noeud")

# Calculate and output the mediated centrality of each node, by using the label name
betweenness = nx.betweenness_centrality(G)
print("\n Intermediary centrality of nodes :")
for node, centrality in betweenness.items():
    print(f"{node}: {centrality}")

plt.show()

