## Network graphs and paths for 36_37 and 20_39

The networks were created by running “create_network.py” with the input files “20_39 EC50s.csv” and “36_37 EC50s.csv”. For each input file, this generated files to describe the edges and nodes of a graph. The edges and nodes were imported into Cytoscape 3.8.0 and arranged into a graph containing the same number of mutations in each column. 

The paths in which each step increases by 16% were also calculated using a recursive algorithm by “create_network.py”. Each path that met the conditions was saved in the respective file, i.e. "20_39_paths.csv" and "36_37_paths.csv".

The final networks are saved as "Networks_20-39_36-37.cys"

**Files:**
* Networks_20-39_36-37.cys - the cytoscape file containing the completed networks
* create_network.py - python script used to create nodes and edges for Cytoscape, and calculate paths in which each step increases by 16%
* 20_39 EC50s.csv and 36_37 EC50s.csv - the input files for creating the graphe
* 20_39_nodes.csv and 36_37_nodes.csv - files containing nodes for Cytoscape
* 20_39_edges.csv and 36_37_edges.csv - files containing edges for Cytoscape
* 20_39_paths.csv and 36_37_paths.csv - all paths within the network in which each step increases by 16%
