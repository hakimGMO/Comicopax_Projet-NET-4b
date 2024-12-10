import networkx as nx
import random
import argparse
import time
import os
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict

class ProteinNetworkAnalyzer:
    """
    Analyzes protein interaction networks from GraphML files.
    
    Provides comprehensive analysis of protein networks including path analysis,
    node frequency analysis, centrality calculations, and report generation.
    The analyzer supports parallel processing for improved performance and
    can generate both HTML and text reports.
    
    Attributes:
        graph (networkx.Graph): NetworkX graph representing the protein network
        threshold (float): Threshold value for identifying ubiquitous nodes (0-1)
        node_types (dict): Maps nodes to their biological types
        proteins (list): List of protein nodes in the network
        shortest_paths (list): List of shortest paths found in analysis
        frequent_nodes (defaultdict): Frequency counter for nodes by type
        ubiquitous_nodes (set): Set of frequently occurring nodes
        undirected_graph (networkx.Graph): Undirected version of the network
    """
    
    def __init__(self, graphml_file: str, threshold: float = 0.9):
        """
        Initialize the ProteinNetworkAnalyzer with a GraphML file.

        Args:
            graphml_file (str): Path to the GraphML file containing the network
            threshold (float, optional): Threshold for ubiquitous node detection. Defaults to 0.9.
        """
        self.graph = nx.read_graphml(graphml_file)
        self.threshold = threshold
        self._initialize_network()
        
    def _initialize_network(self):
        """
        Initialize internal network data structures and attributes.
        
        Performs initial setup including node relabeling, type mapping,
        and creation of derived network representations.
        """
        self._relabel_nodes()
        
        self.node_types = {}
        for node in self.graph.nodes:
            self.node_types[node] = self.graph.nodes[node].get('biopaxType', 'Unknown')
        
        self.proteins = [node for node, node_type in self.node_types.items() 
                        if node_type == 'Protein']
        
        self.shortest_paths = []
        self.frequent_nodes = defaultdict(lambda: defaultdict(int))
        self.ubiquitous_nodes = set()
        
        self.undirected_graph = self.graph.to_undirected()
        
    def _relabel_nodes(self):
        """
        Relabel network nodes using their name attributes if available.
        
        Updates the graph in place by replacing node IDs with their corresponding
        names from the 'name' attribute when present.
        """
        new_names = {}
        for node in self.graph.nodes:
            if 'name' in self.graph.nodes[node]:
                new_names[node] = self.graph.nodes[node]['name']
        self.graph = nx.relabel_nodes(self.graph, new_names)

    def find_shortest_path(self, start: str, end: str) -> list:
        """
        Find the shortest path between two nodes in the network.

        Args:
            start (str): Starting node identifier
            end (str): Target node identifier

        Returns:
            list: Ordered list of nodes representing the shortest path if found
            None: If no path exists between the nodes
        """
        try:
            return nx.shortest_path(self.undirected_graph, source=start, target=end)
        except nx.NetworkXNoPath:
            return None

    def analyze_paths(self, num_iterations: int = 100) -> int:
        """
        Analyze random paths between protein pairs using parallel processing.

        Args:
            num_iterations (int, optional): Number of random protein pairs to analyze.
                Defaults to 100.

        Returns:
            int: Number of valid paths found
        """
        print("\n=== Starting Path Analysis ===")
        print(f"Analyzing {num_iterations} random protein pairs...")
        
        pairs = []
        while len(pairs) < num_iterations:
            p1, p2 = random.sample(self.proteins, 2)
            pairs.append((p1, p2))
        
        paths_found = 0
        
        with ProcessPoolExecutor(max_workers=4) as executor:
            futures = {executor.submit(self.find_shortest_path, p1, p2): (p1, p2) 
                    for p1, p2 in pairs}
            
            for future in as_completed(futures):
                path = future.result()
                if path:
                    paths_found += 1
                    self.shortest_paths.append({
                        'start': path[0],
                        'end': path[-1],
                        'path': path,
                        'length': len(path) - 1
                    })
                    
                    for node in path:
                        node_type = self.node_types[node]
                        self.frequent_nodes[node_type][node] += 1
        
        print(f"Found {paths_found} valid paths out of {num_iterations} attempts")
        print("=== Path Analysis Complete ===\n")
        return paths_found


    def calculate_centrality(self, k: int = 100) -> dict:
        """
        Calculate betweenness centrality for nodes in the network.

        Args:
            k (int, optional): Number of sample nodes for approximation. Defaults to 100.

        Returns:
            dict: Mapping of nodes to their betweenness centrality scores
        """
        print("\n=== Starting Centrality Calculation ===")
        print(f"Using k={k} samples for approximation...")
        
        k = min(k, len(self.graph))
        
        try:
            centrality = nx.betweenness_centrality(
                self.undirected_graph,
                k=k,
                normalized=True
            )
            print("Centrality calculation successful")
            print("=== Centrality Calculation Complete ===\n")
            return centrality
            
        except nx.NetworkXError:
            print("Warning: Graph is disconnected. Calculating centrality for each component in parallel...")
            centrality = {node: 0.0 for node in self.undirected_graph.nodes()}
            components = list(nx.connected_components(self.undirected_graph))
            
            def process_component(comp):
                if len(comp) <= 1:
                    return {node: 0.0 for node in comp}
                    
                subgraph = self.undirected_graph.subgraph(comp)
                try:
                    return nx.betweenness_centrality(
                        subgraph,
                        k=min(k, len(subgraph)),
                        normalized=True
                    )
                except:
                    print(f"Failed to calculate centrality for component of size {len(comp)}")
                    return {node: 0.0 for node in comp}

            with ProcessPoolExecutor(max_workers=4) as executor:
                futures = {
                    executor.submit(process_component, comp): i 
                    for i, comp in enumerate(components)
                }
                
                for future in as_completed(futures):
                    comp_index = futures[future]
                    try:
                        comp_centrality = future.result()
                        centrality.update(comp_centrality)
                        print(f"Processed component {comp_index + 1}/{len(components)}")
                    except Exception as e:
                        print(f"Error processing component {comp_index + 1}: {str(e)}")
            
            print("=== Centrality Calculation Complete ===\n")
            return centrality

    def identify_ubiquitous_nodes(self):
        """
        Identify nodes that appear frequently across analyzed paths.
        
        A node is considered ubiquitous if it appears in a proportion of paths
        greater than the threshold specified during initialization.
        
        Updates the ubiquitous_nodes attribute with the results.
        """
        print("\n=== Starting Ubiquitous Nodes Identification ===")
        print(f"Using threshold: {self.threshold}")
        
        total_paths = sum(
            sum(freq.values()) for freq in self.frequent_nodes.values()
        )
        
        if total_paths > 0:
            self.ubiquitous_nodes = {
                node for type_nodes, freq in self.frequent_nodes.items()
                for node, count in freq.items()
                if count / total_paths >= self.threshold
            }
            print(f"Found {len(self.ubiquitous_nodes)} ubiquitous nodes")
        else:
            print("No paths available for analysis")
        
        print("=== Ubiquitous Nodes Identification Complete ===\n")

    def _generate_path_rows(self):
        """
        Generates HTML table rows for shortest paths information.
        
        Each row contains:
        - Start node
        - End node
        - Path visualization with arrows (→) between nodes
        - Path length
        
        The method processes the self.shortest_paths list which should contain
        dictionaries with 'start', 'end', 'path', and 'length' keys.
        
        Returns:
            str: HTML string containing the generated table rows
        """
        rows = ""
        for info in self.shortest_paths:
            path_str = " → ".join(info['path'])
            rows += f"""
                <tr>
                    <td>{info['start']}</td>
                    <td>{info['end']}</td>
                    <td>{path_str}</td>
                    <td>{info['length']}</td>
                </tr>
            """
        return rows

    def _generate_frequent_nodes_section(self):
        """
        Generates HTML content displaying the most frequent nodes for each node type.
        
        For each node type in self.frequent_nodes, creates a section containing:
        - A header with the node type
        - A table showing the top 10 most frequent nodes and their occurrence counts
        
        The nodes are sorted by frequency in descending order and limited to the top 10.
        
        Returns:
            str: HTML string containing the formatted frequent nodes information
        
        """
        html = ""
        for node_type, occurrences in self.frequent_nodes.items():
            sorted_nodes = sorted(occurrences.items(), key=lambda x: x[1], reverse=True)[:10]
            html += f"""
                <h3>Type: {node_type}</h3>
                <table>
                    <tr><th>Node</th><th>Frequency</th></tr>
                    {"".join(f"<tr><td>{node}</td><td>{freq}</td></tr>"
                            for node, freq in sorted_nodes)}
                </table>
            """
        return html
    
    def _generate_shortest_paths_html(self) -> str:
        """
        Generate HTML content for the shortest paths section, organized by path length.
        
        This method creates a complete HTML section that includes:
        - Path information grouped by length, with each length in its own section
        
        The method processes self.shortest_paths which should contain dictionaries with
        path information including 'length', 'start', 'end', and 'path' keys.
        
        Returns:
            str: Formatted HTML string containing the complete shortest paths analysis. 
            If no paths are found, returns a simple "No paths found" message
        
        """
        if not self.shortest_paths:
            return "<p>No paths found.</p>"
        
        total_paths = len(self.shortest_paths)
        html = f"""
        <style>
            .paths-table {{ margin-top: 20px; }}
            .paths-stats {{ margin-bottom: 20px; padding: 10px; background-color: #f8f9fa; border-radius: 5px; }}
            .path-length-header {{ margin-top: 30px; margin-bottom: 10px; padding: 5px; background-color: #e9ecef; border-radius: 3px; }}
        </style>
        <div class="paths-stats">
            <p><strong>Total paths found: </strong>{total_paths}</p>
        </div>
        """
        
        paths_by_length = defaultdict(list)
        for path_info in self.shortest_paths:
            paths_by_length[path_info['length']].append(path_info)
        
        for length in sorted(paths_by_length.keys()):
            html += self._generate_path_length_section(length, paths_by_length[length])
            
        return html
    
    def _generate_path_length_section(self, length: int, paths: list) -> str:
        """
        Generate HTML content for a specific path length section, displaying all paths 
        of the given length in a formatted table.

        The method creates a section containing:
        - A header showing the path length and number of paths found
        - A sortable table with columns for:
            - Start node
            - End node
            - Complete path visualization (nodes connected by arrows)
        
        The paths are sorted by start node for consistent display.

        Args:
            length (int): The length of the paths in this section (number of edges)
            paths (list): List of path dictionaries, each containing:
                - 'start': Starting node identifier
                - 'end': Ending node identifier
                - 'path': List of nodes in the path order

        Returns:
            str: HTML string containing:
                - A section header with path length info
                - A formatted table of all paths
        """
        html = f"""
        <div class="path-length-header">
            <h3>Paths of length {length} ({len(paths)} paths found)</h3>
        </div>
        <table class="paths-table">
            <tr>
                <th>Start</th>
                <th>End</th>
                <th>Path</th>
            </tr>
        """
        
        sorted_paths = sorted(paths, key=lambda x: x['start'])
        for path_info in sorted_paths:
            path_str = " → ".join(path_info['path'])
            html += f"""
            <tr>
                <td>{path_info['start']}</td>
                <td>{path_info['end']}</td>
                <td>{path_str}</td>
            </tr>
            """
        
        html += "</table>"
        return html
    
    def _generate_ubiquitous_nodes_section(self) -> str:
        """
        Generate HTML content for displaying ubiquitous nodes analysis results.

        This method creates a formatted HTML table showing all nodes that appear frequently
        in the analyzed network paths (ubiquitous nodes). A node is considered ubiquitous
        if it appears in a proportion of paths greater than the threshold specified during
        initialization.

        Returns:
            str: HTML-formatted string containing either:
                - A table with columns for Node, Type, and Occurrences if ubiquitous nodes exist
                - A "No ubiquitous nodes found" message if no nodes meet the threshold
        """
        if not self.ubiquitous_nodes:
            return "<p>No ubiquitous nodes found.</p>"

        html = """
        <table>
            <thead>
                <tr>
                    <th>Node</th>
                    <th>Type</th>
                    <th>Occurrences</th>
                </tr>
            </thead>
            <tbody>
        """

        for node in sorted(self.ubiquitous_nodes):
            node_type = self.node_types.get(node, 'Unknown')
            occurrences = sum(node in path['path'] for path in self.shortest_paths)
            html += f"""
                <tr>
                    <td>{node}</td>
                    <td>{node_type}</td>
                    <td>{occurrences}</td>
                </tr>
            """

        html += """
            </tbody>
        </table>
        """
        return html

    def _generate_centrality_section(self, centrality_scores: dict) -> str:
        """
        Generate HTML content displaying betweenness centrality analysis results grouped by node type.
        
        This method creates a formatted HTML section showing the top 10 nodes with highest
        centrality scores for each node type in the network. The centrality scores indicate
        how important each node is in terms of its position in the network's paths.

        Args:
            centrality_scores (dict): Dictionary mapping node identifiers to their betweenness
                centrality scores. Format: {node_id: float_score}

        Returns:
            str: HTML-formatted string containing either:
                - Multiple tables showing top 10 nodes by type, each with their centrality scores
                - A "No centrality scores available" message if the input dictionary is empty
        """
        if not centrality_scores:
            return "<p>No centrality scores available.</p>"
    
        nodes_by_type = defaultdict(list)
        for node, score in centrality_scores.items():
            node_type = self.node_types.get(node, 'Unknown')
            nodes_by_type[node_type].append((node, score))
        
        for node_type in nodes_by_type:
            nodes_by_type[node_type].sort(key=lambda x: x[1], reverse=True)
            nodes_by_type[node_type] = nodes_by_type[node_type][:10]
        
        tables_html = ""
        for node_type, nodes in sorted(nodes_by_type.items()):
            tables_html += f"""
            <div class="centrality-type-section">
                <h3>Top 10 {node_type} Nodes</h3>
                <table>
                    <thead>
                        <tr>
                            <th>Node</th>
                            <th>Centrality Score</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            
            for node, score in nodes:
                tables_html += f"""
                    <tr>
                        <td>{node}</td>
                        <td>{score:.6f}</td>
                    </tr>
                """
                
            tables_html += """
                    </tbody>
                </table>
            </div>
            """
        
        return tables_html
    
    def generate_html_report(self, execution_time: float) -> str:
        """
        Generate a comprehensive HTML report of the network analysis.

        Args:
            execution_time (float): Total execution time of the analysis in seconds

        Returns:
            str: Complete HTML document containing analysis results
        """
        print("\n=== Starting HTML Report Generation ===")
        print("Compiling analysis results...")
        centrality_scores = self.calculate_centrality()
        
        template = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <title>Protein Network Analysis Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 2rem; color: #333; line-height: 1.6; }}
                .container {{ max-width: 1200px; margin: 0 auto; }}
                table {{ width: 100%; border-collapse: collapse; margin: 1rem 0; }}
                th, td {{ padding: 0.5rem; border: 1px solid #ddd; text-align: left; }}
                th {{ background-color: #f5f5f5; cursor: pointer; }}
                th:hover {{ background-color: #e5e5e5; }}
                tr:nth-child(even) {{ background-color: #f9f9f9; }}
                tr:hover {{ background-color: #f5f5f5; }}
                .section {{ margin: 2rem 0; padding: 1rem; border-radius: 5px; background-color: #fff; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
                h1, h2 {{ color: #2c3e50; margin-top: 0; }}
                .meta-info {{ color: #666; font-size: 0.9rem; }}
                .centrality-stats {{ margin-bottom: 20px; padding: 15px; background-color: #f8f9fa; border-radius: 5px; }}
                .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1rem; }}
                .stat-item {{ padding: 10px; }}
                .stat-item label {{ font-weight: bold; display: block; }}
                .search-input {{ width: 100%; padding: 12px; margin: 10px 0; border: 1px solid #ddd; border-radius: 4px; }}
                .centrality-table-container {{ margin-top: 20px; overflow-x: auto; }}
                .sticky-header {{ position: sticky; top: 0; background-color: white; z-index: 2; padding: 1rem 0; }}
            </style>
        </head>
        <body>
            <div class="container">
                <div class="sticky-header">
                    <h1>Protein Network Analysis Report</h1>
                    <div class="meta-info">
                        <p>Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                        <p>Total execution time: {execution_time:.2f} seconds</p>
                    </div>
                </div>
                
                <div class="section">
                    <h2>Network Overview</h2>
                    <p>Total proteins found: {len(self.proteins)}</p>
                    <p>Total nodes in network: {len(self.graph.nodes)}</p>
                </div>

                <div class="section">
                    <h2>Shortest Paths Analysis</h2>
                    {self._generate_shortest_paths_html()}
                </div>
                
                <div class="section">
                    <h2>Node Recurrence Analysis</h2>
                    {self._generate_frequent_nodes_section()}
                </div>

                <div class="section">
                    <h2>Ubiquitous Nodes</h2>
                    {self._generate_ubiquitous_nodes_section()}
                </div>

                <div class="section">
                    <h2>Centrality Analysis</h2>
                    {self._generate_centrality_section(centrality_scores)}
                </div>
            </div>
        </body>
        </html>
        """
        print("Report generation successful")
        print("=== HTML Report Generation Complete ===\n")
        return template

    def generate_txt_report(self, execution_time: float) -> str:
        """
        Generate a plain text report of the network analysis.

        Args:
            execution_time (float): Total execution time of the analysis in seconds

        Returns:
            str: Formatted text report containing analysis results
        """
        print("\n=== Starting TXT Report Generation ===")
        print("Compiling analysis results...")
        lines = [
            "=== Protein Network Analysis ===",
            f"Executed on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"Execution time: {execution_time:.2f} seconds\n",
            "=== General Statistics ===",
            f"Number of proteins: {len(self.proteins)}",
            f"Analyzed paths: {len(self.shortest_paths)}\n",
            "=== Shortest Paths ===",
        ]

        for info in self.shortest_paths:
            path_str = " -> ".join(info['path'])
            lines.append(
                f"From {info['start']} to {info['end']} "
                f"(length: {info['length']}):\n{path_str}\n"
            )

        lines.append("\n=== Frequent Nodes by Type ===")
        for node_type, occurrences in self.frequent_nodes.items():
            lines.append(f"\nType: {node_type}")
            for node, freq in sorted(occurrences.items(), 
                                   key=lambda x: x[1], 
                                   reverse=True)[:10]:
                lines.append(f"{node}: {freq} occurrences")

        print("Report generation successful")
        print("=== TXT Report Generation Complete ===\n")
        return "\n".join(lines)

def main():
    parser = argparse.ArgumentParser(description="Analyze protein interaction networks")
    parser.add_argument('-f', '--file', required=True, help="Path to GraphML file")
    parser.add_argument('-i', '--iterations', type=int, default=100,
                       help="Number of iterations for path analysis")
    parser.add_argument('-k', '--centrality-samples', type=int, default=100,
                       help="Number of samples for centrality calculation")
    parser.add_argument('-o', '--output', default='network_analysis.html',
                       help="Output file path (HTML or TXT)")
    parser.add_argument('--format', choices=['html', 'txt'], default='html',
                       help="Output report format")
    parser.add_argument('--no-open', action='store_true',
                       help="Don't automatically open the report")
    parser.add_argument('-t', '--threshold', type=float, default=0.9,
                       help="Threshold for ubiquitous nodes identification")

    args = parser.parse_args()
    

    print("\n=== Starting Protein Network Analysis ===")
    print(f"Input file: {args.file}")
    print(f"Number of iterations: {args.iterations}")
    print(f"Centrality samples: {args.centrality_samples}")
    print(f"Output format: {args.format}")

    start_time = time.time()
    
    print("\n=== Initializing Network ===")
    analyzer = ProteinNetworkAnalyzer(args.file, args.threshold)
    print(f"Loaded network with {len(analyzer.proteins)} proteins")
    print("=== Network Initialization Complete ===\n")
    
    analyzer.analyze_paths(args.iterations)
    analyzer.calculate_centrality(args.centrality_samples)
    analyzer.identify_ubiquitous_nodes()
    
    execution_time = time.time() - start_time
    
    print("\n=== Generating Final Report ===")
    if args.format == 'html':
        content = analyzer.generate_html_report(execution_time)
        extension = '.html'
    else:
        content = analyzer.generate_txt_report(execution_time)
        extension = '.txt'

    filename = os.path.splitext(args.output)[0] + extension
    
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(content)
    
    print(f"\n=== Analysis Summary ===")
    print(f"Report saved to: {filename}")
    print(f"Total execution time: {execution_time:.2f} seconds")
    
    absolute_path = os.path.abspath(filename)
    if not args.no_open:
        try:
            if os.name == 'nt':
                os.startfile(absolute_path)
            elif os.name == 'darwin':
                os.system(f'open "{absolute_path}"')
            else:
                os.system(f'xdg-open "{absolute_path}"')
            print(f"Opening report in default application...")
        except Exception as e:
            print(f"Error opening file: {e}")
            print("Please open the file manually.")
    
    print("=== Analysis Complete ===\n")

if __name__ == "__main__":
    main()