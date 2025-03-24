import numpy as np
import sys
import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations

class DAG():
    
    def __init__(self, N):
        self.N = N
        
        self.order1 = np.random.permutation(N)
        self.order2 = np.random.permutation(N)
        self.rank1 = {node: rank for rank, node in enumerate(self.order1)}
        self.rank2 = {node: rank for rank, node in enumerate(self.order2)}
        
        # define graph structure:
        G = nx.DiGraph()
        for i in range(N):
            G.add_node(i, pos=(self.rank1[i], self.rank2[i]))
            for j in range(N):
                if self.rank1[i] < self.rank1[j] and self.rank2[i] < self.rank2[j]:
                    G.add_edge(i, j)
        
        self.G_tc = nx.transitive_closure(G)
        
        # count substructures:
        self.num_3chains = 0
        self.num_diamonds = 0
        self.bh_links = 0
        
        maximal_nodes = set()
        maximal_but_one_nodes = set()
        
        for i in self.G_tc.nodes:
            if len(list(self.G_tc.successors(i))) == 0:
                maximal_nodes.add(i)
            if len(list(self.G_tc.successors(i))) == 1:
                maximal_but_one_nodes.add(i)
        
        for i in self.G_tc.nodes:
            successors = list(self.G_tc.successors(i))
            for j in successors:
                self.num_3chains += len(set(self.G_tc.successors(j)))
            for j, k in combinations(successors, 2):
                if not self.G_tc.has_edge(j, k) and not self.G_tc.has_edge(k, j):
                    common_successors = set(self.G_tc.successors(j)) & set(self.G_tc.successors(k))
                    self.num_diamonds += len(common_successors)
        
        # define a random hyperplane:
        self.h = []
        self.slope = np.random.uniform(-1, 1)
        self.intercept = np.random.uniform(0, N)
        for i in range(N):
            x, y = self.rank1[i], self.rank2[i]
            if y > self.slope * x + self.intercept:
                self.h.append(1)
            else:
                self.h.append(0)
                
        for vertex in maximal_but_one_nodes:
            if self.h[vertex] == 0:
                successor = list(self.G_tc.successors(vertex))[0]
                if successor in maximal_nodes and self.h[successor] == 1:
                    self.bh_links += 1
                    
    def save_dag(self, filename):
        nx.write_gpickle(self.G_tc, filename)
                    
    def plot_dag(self, colors=None, slope=None, intercept=None, ax=None):
        if colors is not None:
            color_labels = [colors[node] for node in self.G_tc.nodes()]
                
        pos = nx.get_node_attributes(self.G_tc, 'pos')
        plt.figure(figsize=(8, 6))
        if colors is None:
            nx.draw(self.G_tc, pos, with_labels=True, node_size=500, node_color='skyblue', edge_color='gray', arrows=True, ax=ax)
        else:
            nx.draw(self.G_tc, pos, with_labels=True, node_size=500, node_color=color_labels, edge_color='gray', arrows=True, ax=ax)
        
        if slope is not None and intercept is not None:
            # Plot the random line
            if ax is not None:
                x_vals = np.array(ax.get_xlim())
            else:
                x_vals = np.array(plt.gca().get_xlim())
            y_vals = slope * x_vals + intercept
            
            if ax is not None:
                ax.plot(x_vals, y_vals, '--', color='red')
            else:
                plt.plot(x_vals, y_vals, '--', color='red')
            
        plt.xlabel("Rank in Order 1")
        plt.ylabel("Rank in Order 2")
        plt.title("DAG from Order Intersection")
        plt.show()
        
    def get_counts(self):
        return {"3-chains": self.num_3chains, "diamonds": self.num_diamonds, "bh-links": self.bh_links}
    
    def get_adjacency_matrix(self):
        adj_matrix = nx.to_numpy_array(self.G_tc, nodelist=range(self.N))
        adj_matrix = np.array(adj_matrix)
        return adj_matrix
    
    def color_graph(self, bitstring1, bitstring2=None, slope=None, intercept=None, ax=None):
        #print(bitstring)
        colors = ["lightblue" for node in range(self.N)]
        for k in range(len(bitstring1)):
            if bitstring1[k] == 1:
                colors[k] = "r"
        if bitstring2 is not None:
            for k in range(len(bitstring2)):
                if bitstring2[k] == 1:
                    colors[k] = "lightgreen"
        print(colors)
        
        self.plot_dag(colors, slope, intercept, ax)
    
    def get_slope_intercept(self):
        return self.slope, self.intercept
        
        
def generate_dag_from_order_intersection(N):
    """
    Generates a DAG from the intersection of two random linear orders.

    Parameters:
        N (int): Number of nodes

    Returns:
        G_tc (networkx.DiGraph): The DAG with explicit transitive closure
        counts (dict): Counts of 3-chains and diamonds
    """
    # Step 1: Generate two independent linear orders
    order1 = np.random.permutation(N)  # Random ranking 1 (e.g., time)
    order2 = np.random.permutation(N)  # Random ranking 2 (e.g., hidden factor)
    print(order1)
    print(order2)

    rank1 = {node: rank for rank, node in enumerate(order1)}
    rank2 = {node: rank for rank, node in enumerate(order2)}
    print(rank1)
    print(rank2)

    G = nx.DiGraph()

    # Step 2: Construct edges based on intersection of both orders
    for i in range(N):
        G.add_node(i, pos=(rank1[i], rank2[i]))  # Store 2D rank as position

        for j in range(N):
            if rank1[i] < rank1[j] and rank2[i] < rank2[j]:  # Intersection rule
                G.add_edge(i, j)

    # Step 3: Compute transitive closure
    G_tc = G#nx.transitive_closure(G)

    # Step 4: Count 3-chains and diamonds
    num_3chains = 0
    num_diamonds = 0

    for i in G_tc.nodes:
        successors = list(G_tc.successors(i))
        for j in successors:
            num_3chains += len(set(G_tc.successors(j)))
        for j, k in combinations(successors, 2):  # Check pairs of successors
            if not G_tc.has_edge(j, k) and not G_tc.has_edge(k, j):  # Ensure j and k are independent
                common_successors = set(G_tc.successors(j)) & set(G_tc.successors(k))
                num_diamonds += len(common_successors)  # Each shared successor forms a diamond

    return G_tc, {"3-chains": num_3chains, "diamonds": num_diamonds}

def plot_dag(G):
    """
    Plots the DAG using its 2D ranking positions.
    """
    pos = nx.get_node_attributes(G, 'pos')
    plt.figure(figsize=(8, 6))
    nx.draw(G, pos, with_labels=True, node_size=500, node_color='skyblue', edge_color='gray', arrows=True)
    plt.xlabel("Rank in Order 1")
    plt.ylabel("Rank in Order 2")
    plt.title("DAG from Order Intersection")
    plt.show()
    
# create a main funciton to run the code
if __name__ == "__main__":
    # take the number of nodes from the command line as an argument:
    N = int(sys.argv[1])
    
    # Generate a DAG from the intersection of two random linear orders
    G, counts = generate_dag_from_order_intersection(N)
    print(counts)
    
    # Plot the DAG using its 2D ranking positions
    plot_dag(G)