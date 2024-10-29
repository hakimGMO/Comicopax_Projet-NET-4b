"""
Contient la classe GraphPlotter pour créer des
visualisations de graphes et les enregistrer
sous forme d'images.
"""

import networkx as nx
import matplotlib.pyplot as plt

class GraphPlotter:
    """
    Une classe pour créer des visualisations de graphes et les enregistrer sous forme d'images PNG.

    Attributes:
        graph (networkx.Graph): Le graphe à visualiser.
        figsize (tuple): La taille de la figure.
        node_size (int): La taille des nœuds.
        node_color (str): La couleur des nœuds.
        edge_color (str): La couleur des arêtes.
        font_size (int): La taille de la police.
        margins (float): Les marges de la figure.
    """

    def __init__(self, graph, figsize = (15, 10), node_size = 2000,
                 node_color = 'white', edge_color = 'black',
                 font_size = 8, margins = 0.3):
        """
        Initialise le GraphPlotter avec un graphe et des paramètres de visualisation.

        Args:
            graph (networkx.Graph): Le graphe à visualiser.
            figsize (tuple, optional): La taille de la figure. Par défaut (15, 10).
            node_size (int, optional): La taille des nœuds. Par défaut 2000.
            node_color (str, optional): La couleur des nœuds. Par défaut 'white'.
            edge_color (str, optional): La couleur des arêtes. Par défaut 'black'.
            font_size (int, optional): La taille de la police. Par défaut 8.
            margins (float, optional): Les marges de la figure. Par défaut 0.3.
        """
        self.graph = graph
        self.figsize = figsize
        self.node_size = node_size
        self.node_color = node_color
        self.edge_color = edge_color
        self.font_size = font_size
        self.margins = margins

    @property
    def get_graph(self):
        """
        Obtient le graphe à visualiser.

        Returns:
            networkx.Graph: Le graphe à visualiser.
        """
        return self.graph

    @property
    def get_figsize(self):
        """
        Obtient la taille de la figure.

        Returns:
            tuple: La taille de la figure.
        """
        return self.figsize

    @property
    def get_node_size(self):
        """
        Obtient la taille des nœuds.

        Returns:
            int: La taille des nœuds.
        """
        return self.node_size

    @property
    def get_node_color(self):
        """
        Obtient la couleur des nœuds.

        Returns:
            str: La couleur des nœuds.
        """
        return self.node_color

    @property
    def get_edge_color(self):
        """
        Obtient la couleur des arêtes.

        Returns:
            str: La couleur des arêtes.
        """
        return self.edge_color

    @property
    def get_font_size(self):
        """
        Obtient la taille de la police.

        Returns:
            int: La taille de la police.
        """
        return self.font_size

    @property
    def get_margins(self):
        """
        Obtient les marges de la figure.

        Returns:
            float: Les marges de la figure.
        """
        return self.margins

    @get_graph.setter
    def set_graph(self, graph):
        """
        Définit le graphe à visualiser.

        Args:
            graph (networkx.Graph): Le nouveau graphe à visualiser.
        """
        self.graph = graph

    @get_figsize.setter
    def set_figsize(self, figsize):
        """
        Définit la taille de la figure.

        Args:
            figsize (tuple): La nouvelle taille de la figure.
        """
        self.figsize = figsize

    @get_node_size.setter
    def set_node_size(self, node_size):
        """
        Définit la taille des nœuds.

        Args:
            node_size (int): La nouvelle taille des nœuds.
        """
        self.node_size = node_size

    @get_node_color.setter
    def set_node_color(self, node_color):
        """
        Définit la couleur des nœuds.

        Args:
            node_color (str): La nouvelle couleur des nœuds.
        """
        self.node_color = node_color

    @get_edge_color.setter
    def set_edge_color(self, edge_color):
        """
        Définit la couleur des arêtes.

        Args:
            edge_color (str): La nouvelle couleur des arêtes.
        """
        self.edge_color = edge_color

    @get_font_size.setter
    def set_font_size(self, font_size):
        """
        Définit la taille de la police.

        Args:
            font_size (int): La nouvelle taille de la police.
        """
        self.font_size = font_size

    @get_margins.setter
    def set_margins(self, margins):
        """
        Définit les marges de la figure.

        Args:
            margins (float): Les nouvelles marges de la figure.
        """
        self.margins = margins

    def plot(self, output_file = None):
        """
        Crée une visualisation du graphe et l'enregistre sous forme d'image PNG.

        Args:
            output_file (str): Le chemin du fichier de sortie pour l'image PNG.
        """
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