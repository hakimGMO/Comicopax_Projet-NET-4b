import os
from graph_visualizer import GraphVisualizer
from graph_plotter import GraphPlotter

absolute_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(absolute_dir, '..', 'Data')
figures_dir = os.path.join(absolute_dir, '..', 'Figures')

os.makedirs(figures_dir, exist_ok = True)

visualizer = GraphVisualizer(data_dir)

graphml_files = [f for f in os.listdir(data_dir) if f.endswith('.graphml')]

for file_name in graphml_files:
    visualizer.load_graph(file_name)
    plotter = GraphPlotter(visualizer.graph)
    base_name = os.path.splitext(file_name)[0]
    output_file = os.path.join(figures_dir, f"{base_name}.png")
    plotter.plot(output_file)