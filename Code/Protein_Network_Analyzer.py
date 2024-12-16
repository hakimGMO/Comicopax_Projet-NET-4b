import networkx as nx
import argparse
import time
import random
import os
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter, defaultdict
from itertools import combinations
from typing import List, Tuple

class ProteinNetworkAnalyzer:
    """
    A class for analyzing protein interaction networks from GraphML files.

    This analyzer processes protein networks to identify and analyze shortest paths,
    calculate node frequencies, and determine network centrality metrics.

    Attributes:
        DEFAULT_BLACKLIST (set): Default set of metabolites and small molecules to exclude
        graph (NetworkX.Graph): The loaded network graph
        threshold (float): Threshold value for filtering interactions
        blacklist (set): Current set of nodes to exclude from analysis
        node_types (dict): Mapping of nodes to their biological types
        proteins (list): List of valid protein nodes after filtering
        shortest_paths (list): List of all analyzed shortest paths
        pair_frequencies (dict): Frequencies of nodes in paths between specific pairs
        normalized_frequencies (dict): Overall normalized node frequencies
        analyzed_nodes (set): Set of all nodes encountered in analysis
        undirected_graph (NetworkX.Graph): Undirected version of the network

    Args:
        graphml_file (str): Path to the GraphML file containing the network
        threshold (float, optional): Threshold for filtering. Defaults to 0.9
        custom_blacklist (set, optional): Additional nodes to blacklist. Defaults to None
    """

    DEFAULT_BLACKLIST = {
        'ATP', 'ADP','NADH', 'NAD+','NADPH', 'NADP+','FADH2', 'FAD','Pyruvate','Pi', 'Phosphate',
        'PPi', 'Pyrophosphate','H+', 'Proton','CO2','H2O','O2'}

    def __init__(self, graphml_file: str, custom_blacklist: set = None):
        """
        Initialize the ProteinNetworkAnalyzer with a GraphML file and analysis parameters.
        
        Args:
            graphml_file (str): Path to the GraphML file containing the protein network
            threshold (float, optional): Threshold value for filtering interactions. Defaults to 0.9
            custom_blacklist (set, optional): Set of node names to exclude from analysis. Defaults to None
        """
        self.graph = nx.read_graphml(graphml_file)
        self.blacklist = self.DEFAULT_BLACKLIST.union(custom_blacklist or set())
        self._initialize_network()
        
    def _initialize_network(self):
        """
        Initialize the network by setting up internal data structures and attributes.
        
        This method:
        - Relabels nodes using their names
        - Identifies node types
        - Filters protein nodes against blacklist
        - Initializes data structures for path analysis
        - Creates an undirected version of the graph
        """
        self._relabel_nodes()
    
        self.node_types = {}
        for node in self.graph.nodes:
            self.node_types[node] = self.graph.nodes[node].get('biopaxType', 'Unknown')
        
        self.proteins = [node for node in self.graph.nodes 
                        if self.node_types[node] == 'Protein' 
                        and not self._is_blacklisted(str(node))]
        
        self.shortest_paths = []
        self.pair_frequencies = {}
        self.normalized_frequencies = defaultdict(lambda: defaultdict(float))
        self.analyzed_nodes = set()
        
        self.undirected_graph = self.graph.to_undirected()
    
    def _is_blacklisted(self, node_name: str) -> bool:
        """
        Check if a node is blacklisted using exact matching.
        
        Args:
            node_name: Name of the node to check
                
        Returns:
            bool: True if node is blacklisted, False otherwise
        """
        node_name = str(node_name).upper().strip()
        
        if any(node_name == item.upper().strip() for item in self.blacklist):
            return True
            
        node_words = set(node_name.split())
        for blacklisted in self.blacklist:
            blacklisted_words = set(blacklisted.upper().strip().split())
            if blacklisted_words and blacklisted_words.issubset(node_words):
                return True
                
        return False

    def _relabel_nodes(self):
        """
        Relabel network nodes using their 'name' attribute if available.
        
        Modifies the graph in place by replacing node IDs with their corresponding
        names from the 'name' attribute of each node.
        """
        new_names = {}
        for node in self.graph.nodes:
            if 'name' in self.graph.nodes[node]:
                new_names[node] = self.graph.nodes[node]['name']
        self.graph = nx.relabel_nodes(self.graph, new_names)

    def find_all_shortest_paths(self, start: str, end: str) -> list:
        """
        Find all shortest paths between two nodes in the network.

        Args:
            start (str): Starting node identifier
            end (str): Target node identifier
            
        Returns:
            list: List of paths, where each path is a list of node identifiers.
                Returns None if no path exists between the nodes.
        """
        try:
            return list(nx.all_shortest_paths(self.undirected_graph, source=start, target=end))
        except (nx.NetworkXNoPath, nx.NetworkXError):
            return None

    def select_nodes(self, num_nodes=None, node_list=None):
        """
        Select nodes for analysis, raising error if blacklisted nodes are found.
        
        Args:
            num_nodes: Number of nodes to randomly select from proteins list
            node_list: List of specific nodes to analyze
            
        Returns:
            list: Selected nodes if all are valid
            
        Raises:
            ValueError: If blacklisted nodes are found or no valid nodes remain
        """
        if num_nodes is None and node_list is None:
            raise ValueError("Either num_nodes or node_list must be provided")
        
        if num_nodes is not None:
            if not self.proteins:
                raise ValueError(
                    "No valid protein nodes found after filtering blacklist. "
                    "Use --disable-blacklist to include blacklisted nodes."
                )
            return random.sample(self.proteins, min(num_nodes, len(self.proteins)))
        
        blacklisted_nodes = []
        missing_nodes = []
        valid_nodes = []
        
        for node in node_list:
            if node not in self.graph.nodes():
                missing_nodes.append(node)
            elif self._is_blacklisted(str(node)):
                blacklisted_nodes.append(node)
            else:
                valid_nodes.append(node)
        
        error_msgs = []
        if missing_nodes:
            error_msgs.append(f"Nodes not found in network: {', '.join(missing_nodes)}")
        if blacklisted_nodes:
            error_msgs.append(
                f"Blacklisted nodes detected: {', '.join(blacklisted_nodes)}\n"
                "To proceed, either:\n"
                "1. Remove the blacklisted nodes from your input\n"
                "2. Use --disable-blacklist to disable the blacklist"
            )
        
        if error_msgs:
            raise ValueError("\n".join(error_msgs))
            
        if not valid_nodes:
            raise ValueError("No valid nodes remain for analysis")
            
        return valid_nodes


    def _compute_pair_frequencies(self, paths: List[List[str]], pair: Tuple[str, str]):
        """
        Calculate normalized frequencies for a specific pair of nodes.
        
        Args:
            paths (List[List[str]]): List of shortest paths for this pair of nodes
            pair (Tuple[str, str]): Tuple containing the node pair (start, end)
        """
        if not paths:
            return
        
        node_counts = Counter(node for path in paths for node in path)
        num_paths = len(paths)
        normalized_freqs = defaultdict(dict)
        
        for node, count in node_counts.items():
            node_type = self.node_types.get(node, 'Unknown')
            norm_freq = count / num_paths
            normalized_freqs[node_type][node] = norm_freq
            self.normalized_frequencies[node_type][node] += norm_freq
        
        self.pair_frequencies[pair] = normalized_freqs

    def analyze_paths(self, num_nodes=None, node_list=None):
        """
        Analyze shortest paths between selected proteins in the network and compute path frequencies.

        Args:
            num_nodes (int, optional): Number of random proteins to analyze
            node_list (List[str], optional): List of specific proteins to analyze

        Returns:
            int: Total number of shortest paths found in the analysis
            
        Raises:
            ValueError: If neither num_nodes nor node_list is provided, or if no valid nodes are found
        """
        print("\n=== Starting Path Analysis ===")
        
        try:
            selected_proteins = self.select_nodes(num_nodes, node_list)
            pairs = list(combinations(selected_proteins, 2))
            print(f"Analyzing all paths between {len(selected_proteins)} nodes...")
            print(f"Analyzing {len(pairs)} pairs")
            
            total_paths_found = 0
            self.normalized_frequencies = defaultdict(lambda: defaultdict(float))
            
            for p1, p2 in pairs:
                pair = (p1, p2)
                paths = self.find_all_shortest_paths(p1, p2)
                
                if paths:  # Only process if paths were found
                    total_paths_found += len(paths)
                    self._compute_pair_frequencies(paths, pair)
                    
                    for path in paths:
                        self.shortest_paths.append({
                            'start': path[0],
                            'end': path[-1],
                            'path': path,
                            'length': len(path) - 1
                        })
                        
                        for node in path:
                            self.analyzed_nodes.add(node)
            
            for pair_freqs in self.pair_frequencies.values():
                for node_type, nodes in pair_freqs.items():
                    for node, norm_freq in nodes.items():
                        self.normalized_frequencies[node_type][node] += norm_freq
            
            print(f"Found {total_paths_found} shortest paths analyzing {len(pairs)} pairs")
            print("=== Path Analysis Complete ===\n")
            return total_paths_found
            
        except Exception as e:
            print(f"Error during path analysis: {str(e)}")
            raise

    def calculate_centrality(self):
        """
        Calculate normalized betweenness centrality for all encountered nodes.
        
        Returns:
            dict: Mapping of nodes to their normalized centrality scores 
                (normalized by NetworkX's built-in normalization)
        """
        if not self.analyzed_nodes:
            return {}
            
        print("\n=== Starting Centrality Calculation ===")
        analyzed_subgraph = self.undirected_graph.subgraph(self.analyzed_nodes)
        normalized_centrality = defaultdict(float)
        
        try:
            if len(self.analyzed_nodes) < 1000:
                return nx.betweenness_centrality(analyzed_subgraph, normalized=True)
        except:
            pass
            
        components = list(nx.connected_components(analyzed_subgraph))
        total_nodes = len(analyzed_subgraph)
        chunk_size = max(1, len(components) // (os.cpu_count() or 4))
        component_chunks = [components[i:i + chunk_size] for i in range(0, len(components), chunk_size)]
        
        with ProcessPoolExecutor() as executor:
            futures = []
            for chunk in component_chunks:
                for comp in chunk:
                    if len(comp) <= 1:
                        normalized_centrality.update({node: 0.0 for node in comp})
                        continue
                        
                    subgraph = analyzed_subgraph.subgraph(comp)
                    futures.append(executor.submit(nx.betweenness_centrality, subgraph, normalized=True))
                    
            for future, comp in zip(as_completed(futures), (c for chunk in component_chunks for c in chunk if len(c) > 1)):
                try:
                    comp_size_factor = len(comp) / total_nodes
                    normalized_centrality.update({
                        node: score * comp_size_factor 
                        for node, score in future.result().items()
                    })
                except Exception as e:
                    print(f"Error in component: {e}")
                    normalized_centrality.update({node: 0.0 for node in comp})
        
        print("=== Centrality Calculation Complete ===\n")
        return dict(normalized_centrality)

    def get_normalized_frequencies(self):
        """
        Retrieve the combined normalized frequencies of nodes in all shortest paths.
        
        Returns:
            dict: A nested dictionary structure
        """
        return dict(self.normalized_frequencies)

    def _generate_frequent_nodes_section(self):
        """
        Generate HTML content displaying the normalized node frequencies.
        
        Creates a formatted section showing how often each node appears in the 
        shortest paths, normalized by the number of paths and grouped by node type.
        
        Returns:
            str: HTML-formatted string containing tables of node frequencies
        """
        html = """
        <p>These frequencies represent how often each node appears in the shortest paths, 
        normalized by the number of shortest paths for each protein pair.</p>
        """
        
        normalized_freqs = self.get_normalized_frequencies()
        
        for node_type, normalized_occurrences in normalized_freqs.items():
            sorted_nodes = sorted(normalized_occurrences.items(), 
                                key=lambda x: x[1], 
                                reverse=True)
            
            if sorted_nodes:  # If we have nodes of this type
                html += f"""
                    <h3>Type: {node_type}</h3>
                    <table>
                        <tr>
                            <th>Node</th>
                            <th>Normalized Frequency</th>
                        </tr>
                """
                
                for node, norm_freq in sorted_nodes:
                    html += f"""
                        <tr>
                            <td>{node}</td>
                            <td>{norm_freq:.4f}</td>
                        </tr>
                    """
                html += "</table>"
        
        return html

    def _generate_pair_frequencies_section(self):
        """
        Generate HTML content displaying the normalized frequencies for each pair of nodes.
        
        Creates sections for each analyzed node pair showing the frequency of 
        intermediate nodes in their shortest paths, organized by node type.
        
        Returns:
            str: HTML-formatted string containing pair frequency analysis,
                or a "No pair frequencies available" message if none exist
        """
        if not self.pair_frequencies:
            return "<p>No pair frequencies available.</p>"

        html = """
        <div class="pair-frequencies">
            <p>Each frequency shown below represents the number of times a node appears in the shortest paths 
            between a specific pair, normalized by the total number of shortest paths for that pair.</p>
        """
        
        for pair, freq_data in sorted(self.pair_frequencies.items()):
            start, end = pair
            html += f"""
            <div class="pair-section">
                <h3>Pair: {start} ↔ {end}</h3>
                <div class="scrollable-wrapper">
            """
            
            for node_type, frequencies in freq_data.items():
                if frequencies:  # Only show node types that have frequencies
                    html += f"""
                    <h4>Node Type: {node_type}</h4>
                    <table>
                        <tr>
                            <th>Node</th>
                            <th>Normalized Frequency</th>
                        </tr>
                    """
                    
                    for node, freq in sorted(frequencies.items(), key=lambda x: x[1], reverse=True):
                        html += f"""
                        <tr>
                            <td>{node}</td>
                            <td>{freq:.4f}</td>
                        </tr>
                        """
                    html += "</table>"
            
            html += """
                </div>
            </div>
            """
        
        return html + "</div>"
    
    def _generate_shortest_paths_html(self) -> str:
        """
        Generate HTML content for the shortest paths section.
        
        Creates a comprehensive overview of all shortest paths found, including:
        - Total number of paths
        - Paths organized by length
        - Path visualizations with node connections
        
        Returns:
            str: HTML-formatted string containing the complete shortest paths analysis,
                or a "No paths found" message if none exist
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
        
        Creates a complete HTML document containing:
        - Network overview statistics
        - Shortest paths analysis
        - Node frequency analysis
        - Centrality analysis
        - Execution metadata
        
        Args:
            execution_time (float): Total execution time of the analysis in seconds
        
        Returns:
            str: Complete HTML document string containing all analysis results
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
                body {{ 
                    font-family: Arial, sans-serif; 
                    margin: 2rem; 
                    color: #333; 
                    line-height: 1.6; 
                }}
                .container {{ 
                    max-width: 1200px; 
                    margin: 0 auto; 
                }}
                table {{ 
                    width: 100%; 
                    border-collapse: collapse; 
                    margin: 1rem 0;
                    table-layout: fixed;
                }}
                th, td {{ 
                    padding: 0.5rem; 
                    border: 1px solid #ddd; 
                    text-align: left;
                    min-width: 150px;
                    max-width: 300px;
                }}
                td {{ 
                    white-space: nowrap;
                    overflow-x: auto;
                    position: relative;
                }}
                td::-webkit-scrollbar {{
                    height: 8px;
                }}
                td::-webkit-scrollbar-track {{
                    background: #f1f1f1;
                    border-radius: 4px;
                }}
                td::-webkit-scrollbar-thumb {{
                    background: #888;
                    border-radius: 4px;
                }}
                td::-webkit-scrollbar-thumb:hover {{
                    background: #555;
                }}
                td.scrollable::after {{
                    content: '⟷';
                    position: absolute;
                    right: 4px;
                    top: 50%;
                    transform: translateY(-50%);
                    color: #888;
                    font-size: 12px;
                    opacity: 0.7;
                }}
                th {{ 
                    background-color: #f5f5f5; 
                    position: sticky;
                    top: 0;
                    z-index: 10;
                }}
                tr:nth-child(even) {{ 
                    background-color: #f9f9f9; 
                }}
                tr:hover {{ 
                    background-color: #f5f5f5; 
                }}
                .section {{ 
                    margin: 2rem 0; 
                    padding: 1rem; 
                    border-radius: 5px; 
                    background-color: #fff; 
                    box-shadow: 0 2px 4px rgba(0,0,0,0.1); 
                }}
                .scrollable-wrapper {{
                    overflow-x: auto;
                    margin: 1rem 0;
                    padding-bottom: 1rem;
                }}
                h1, h2 {{ 
                    color: #2c3e50; 
                    margin-top: 0; 
                }}
                .meta-info {{ 
                    color: #666; 
                    font-size: 0.9rem; 
                }}
                .sticky-header {{ 
                    position: sticky; 
                    top: 0; 
                    background-color: white; 
                    z-index: 20; 
                    padding: 1rem 0; 
                }}
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
                    <div class="scrollable-wrapper">
                        <p>Total proteins found: {len(self.proteins)}</p>
                        <p>Total nodes in network: {len(self.graph.nodes)}</p>
                    </div>
                </div>

                <div class="section">
                    <h2>Shortest Paths Analysis</h2>
                    <div class="scrollable-wrapper">
                        {self._generate_shortest_paths_html()}
                    </div>
                </div>
                
                <div class="section">
                    <h2>Node Frequencies</h2>
                    <div class="scrollable-wrapper">
                        {self._generate_frequent_nodes_section()}
                    </div>
                </div>

                <div class="section">
                    <h2>Centrality Analysis</h2>
                    <div class="scrollable-wrapper">
                        {self._generate_centrality_section(centrality_scores)}
                    </div>
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
    """
    Entry point for the protein network analysis script.
    
    Parses command line arguments to:
    - Load and analyze protein interaction networks from GraphML files
    - Compute shortest paths and centrality metrics
    - Generate analysis reports in HTML or TXT format
    - Display the report in the default system viewer
    
    Command line arguments:
        -f, --file: Path to GraphML file (required)
        -o, --output: Output file path (default: network_analysis.html)
        --format: Output format, 'html' or 'txt' (default: html)
        --no-open: Don't automatically open the report
        -r, --random-nodes: Number of nodes to select randomly
        -l, --node-list: List of specific nodes to analyze
        --disable-blacklist: Disable the default metabolite blacklist
    """
    parser = argparse.ArgumentParser(description="Analyze protein interaction networks")
    parser.add_argument('-f', '--file', required=True, help="Path to GraphML file")
    parser.add_argument('-o', '--output', default='network_analysis.html',
                       help="Output file path")
    parser.add_argument('--format', choices=['html', 'txt'], default='html',
                       help="Output report format")
    parser.add_argument('--no-open', action='store_true',
                       help="Don't automatically open the report")
    node_group = parser.add_mutually_exclusive_group()
    node_group.add_argument('-r', '--random-nodes', type=int,
                           help="Number of nodes to select randomly")
    node_group.add_argument('-l', '--node-list', nargs='+',
                           help="List of specific nodes to analyze")
    parser.add_argument('--disable-blacklist', action='store_true',
                       help="Disable the default metabolite blacklist")
    args = parser.parse_args()
    start_time = time.time()
    
    try:
        analyzer = ProteinNetworkAnalyzer(
            args.file,
            custom_blacklist=set() if args.disable_blacklist else None
        )
    except Exception as e:
        print(f"Error initializing network analyzer: {e}")
        return
    
    try:
        paths_found = analyzer.analyze_paths(
            num_nodes=args.random_nodes,
            node_list=args.node_list
        )
        print(f"Successfully analyzed {paths_found} paths")
    except ValueError as e:
        print(f"Error: {str(e)}")
        return
    except Exception as e:
        print(f"Error during path analysis: {e}")
        return
    
    execution_time = time.time() - start_time
    
    if args.format == 'html':
        content = analyzer.generate_html_report(execution_time)
        extension = '.html'
    else:
        content = analyzer.generate_txt_report(execution_time)
        extension = '.txt'

    filename = os.path.splitext(args.output)[0] + extension
    
    try:
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"\nReport saved to: {filename}")
    except Exception as e:
        print(f"Error saving report: {e}")
        return
    
    print(f"Total execution time: {execution_time:.2f} seconds")
    
    if not args.no_open:
        absolute_path = os.path.abspath(filename)
        try:
            if os.name == 'nt':
                os.startfile(absolute_path)
            elif os.name == 'darwin':
                os.system(f'open "{absolute_path}"')
            else:
                os.system(f'xdg-open "{absolute_path}"')
            print("Opening report in default application...")
        except Exception as e:
            print(f"Error opening file: {e}")
            print("Please open the file manually.")

if __name__ == "__main__":
    main()
