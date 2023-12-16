import os
import matplotlib.pyplot as plt
import networkx as nx

script_directory = os.path.dirname(os.path.abspath(__file__))

# Directories
vertices_directory = os.path.join(script_directory, '../res/verticesData')
edges_directory = os.path.join(script_directory, '../res/MSTEdges')
plots_out_directory = os.path.join(script_directory, "../res/plots")



# Traversing through the files in directories
for file in os.listdir(vertices_directory):
    if file.endswith(".tsp"):
        vertices_file = os.path.join(vertices_directory, file)
        edges_file = os.path.join(edges_directory, file)

        # Create an empty graph
        G = nx.Graph()

        # Read vertex coordinates
        with open(vertices_file, 'r') as file:
            vertex_data = file.readlines()[8:-1]

        # Create a dictionary to store coordinates for each vertex
        vertices = {}
        for line in vertex_data:
            vertex, x, y = map(int, line.split())
            vertices[vertex] = (x, y)
            G.add_node(vertex)  # Add nodes to the graph
            G.nodes[vertex]['pos'] = (x, y)  # Assign positions to nodes

        # Read edges
        with open(edges_file, 'r') as file:
            edge_data = file.readlines()

        # Add edges to the graph
        for edge in edge_data:
            vertex1, vertex2 = map(int, edge.split())
            if vertex1 in vertices and vertex2 in vertices:  # Ensure both vertices exist
                G.add_edge(vertex1, vertex2)  # Add edges to the graph
            else:
                print(f"Issue with edge: {vertex1} - {vertex2}. One or both vertices are missing.")

        # Draw the graph
        pos = nx.get_node_attributes(G, 'pos')
        nx.draw(G, pos, with_labels=False, node_size=5, node_color='red', font_weight='bold')

        plt.suptitle(f'Local Search Cycle for {os.path.splitext(os.path.basename(file.name))[0]}')  # Set title for the current file

        plt.grid(True)
        plt.savefig(f'{plots_out_directory}/{os.path.splitext(os.path.basename(file.name))[0]}.png')
        plt.close()