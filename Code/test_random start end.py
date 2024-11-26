import networkx as nx
import random


graphml_file = '/Users/hakimbouazzaoui/Downloads/PathwayCommons12.netpath.BIOPAX.graphml'
G = nx.read_graphml(graphml_file)

# Renomme les nœuds avec leur nom s'il existe
new_names = {}
for node in G.nodes():
    if 'name' in G.nodes[node]:
        new_names[node] = G.nodes[node]['name']
G = nx.relabel_nodes(G, new_names)

# Trouve toutes les protéines dans le graphe
proteins = []
for node in G.nodes():
    if G.nodes[node].get('biopaxType') == 'Protein':
        proteins.append(node)

# selectionaléatoirement une protéine de départ et une d'arrivée
start = random.choice(proteins)
proteins.remove(start)  # Retire la protéine de départ de la liste
end = random.choice(proteins)

print(f"\nProtéine de départ: {start}")
print(f"Protéine d'arrivée: {end}")

