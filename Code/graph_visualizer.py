import networkx as nx
import os

class GraphVisualizer:

    def __init__(self, data_dir):

        self.data_dir = data_dir
        self.graph = None

    @property
    def data_dir(self):
        return self.data_dir
        
    @data_dir.setter
    def data_dir(self, value):
        self.data_dir = value

    @property
    def graph(self):
        return self.graph
        
    @graph.setter
    def graph(self, value):
        self.graph = value

    def load_graph(self, file_name):
        file_path = os.path.join(self.data_dir, file_name)
        self.graph = nx.read_graphml(file_path)