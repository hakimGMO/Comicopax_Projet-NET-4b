"""
Contient la classe GraphVisualizer pour charger
et manipuler des graphes à partir de fichiers GraphML.
"""

import networkx as nx
import os

class GraphVisualizer:
    """
    Une classe pour visualiser des graphes à partir de fichiers GraphML.

    Attributes:
        data_dir (str): Le répertoire contenant les fichiers GraphML.
        graph (networkx.Graph): Le graphe chargé.
    """

    def __init__(self, data_dir):
        """
        Initialise le GraphVisualizer avec un répertoire de données.

        Args:
            data_dir (str): Le répertoire contenant les fichiers GraphML.
        """
        self.data_dir = data_dir
        self.graph = None

    @property
    def data_dir(self):
        """
        Obtient le répertoire de données.

        Returns:
            str: Le répertoire contenant les fichiers GraphML.
        """
        return self.data_dir
        
    @data_dir.setter
    def data_dir(self, value):
        """
        Définit le répertoire de données.

        Args:
            value (str): Le nouveau répertoire contenant les fichiers GraphML.
        """
        self.data_dir = value

    @property
    def graph(self):
        """
        Obtient le graphe chargé.

        Returns:
            networkx.Graph: Le graphe chargé.
        """
        return self.graph
        
    @graph.setter
    def graph(self, value):
        """
        Définit le graphe chargé.

        Args:
            value (networkx.Graph): Le nouveau graphe à charger.
        """
        self.graph = value

    def load_graph(self, file_name):
        """
        Charge un fichier GraphML et initialise le graphe.

        Args:
            file_name (str): Le nom du fichier GraphML à charger.
        """
        file_path = os.path.join(self.data_dir, file_name)
        self.graph = nx.read_graphml(file_path)