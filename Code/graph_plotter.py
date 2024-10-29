import networkx as nx
import matplotlib.pyplot as plt

class GraphPlotter:

    def __init__(self, graph, figsize = (15, 10), node_size = 2000,
                 node_color = 'white', edge_color = 'black',
                 font_size = 8, margins = 0.3):
        self.graph = graph
        self.figsize = figsize
        self.node_size = node_size
        self.node_color = node_color
        self.edge_color = edge_color
        self.font_size = font_size
        self.margins = margins

    @property
    def get_graph(self):
        return self.graph

    @property
    def get_figsize(self):
        return self.figsize

    @property
    def get_node_size(self):
        return self.node_size

    @property
    def get_node_color(self):
        return self.node_color

    @property
    def get_edge_color(self):
        return self.edge_color

    @property
    def get_font_size(self):
        return self.font_size

    @property
    def get_margins(self):
        return self.margins

    @get_graph.setter
    def set_graph(self, graph):
        self.graph = graph

    @get_figsize.setter
    def set_figsize(self, figsize):
        self.figsize = figsize

    @get_node_size.setter
    def set_node_size(self, node_size):
        self.node_size = node_size

    @get_node_color.setter
    def set_node_color(self, node_color):
        self.node_color = node_color

    @get_edge_color.setter
    def set_edge_color(self, edge_color):
        self.edge_color = edge_color

    @get_font_size.setter
    def set_font_size(self, font_size):
        self.font_size = font_size

    @get_margins.setter
    def set_margins(self, margins):
        self.margins = margins

    def plot(self, output_file = None):

        plt.figure(figsize = self.figsize)

        pos = nx.spring_layout(self.graph, k = 10, iterations = 100, scale =3)

        nx.draw_networkx_nodes(self.graph, pos,
                               node_shape = 'o',
                               node_size = self.node_size,
                               node_color = self.node_color,
                               edgecolors = self.edge_color)
        
        nx.draw_networkx_edges(self.graph, pos,
                               arrows = True,
                               arrowsize = 20,
                               width = 1,
                               connectionstyle = 'arc,rad=0.15')
        
        labels = nx.get_node_attributes(self.graph, 'label')
        for node in self.graph.nodes():
            if not labels.get(node):
                labels[node] = self.graph.nodes[node]['type']

        pos_attrs = {}
        for node, coords in pos.items():
            pos_attrs[node] = (coords[0], coords[1] + 0.08)

        nx.draw_networkx_labels(self.graph, pos, labels, font_size = self.font_size)

        edge_labels = {}
        for u, v, data in self.graph.edges(data = True):
            edge_labels[(u,v)] = data['type']

        nx.draw_networkx_edge_labels(self.graph, pos,
                                     edge_labels,
                                     font_size = self.font_size,
                                     label_pos = 0.3)
        
        plt.axis('off')
        plt.margins(self.margins)
        plt.tight_layout()
        plt.savefig(output_file, format = 'png')