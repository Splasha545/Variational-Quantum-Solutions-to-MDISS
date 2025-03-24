import networkx as nx
import os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from RNG_DAG import DAG
from VQALib import CubicProgram, IsingHamiltonian
from QUBO_FUNCTIONS import QUBO_Problem

def get_new_folder(base_dir="./vis_results", size=5, method="three-chain"):
    """Creates a new folder with an incrementing index."""
    os.makedirs(base_dir, exist_ok=True)
    existing_folders = [f for f in os.listdir(base_dir) if f.startswith(f"{method}_size_{size}_")]
    existing_indices = [int(f.split('_')[-1]) for f in existing_folders if f.split('_')[-1].isdigit()]
    next_index = max(existing_indices, default=0) + 1
    new_folder = os.path.join(base_dir, f"{method}_size_{size}_{next_index}")
    os.makedirs(new_folder)
    return new_folder

class DAGVisualization:
    def __init__(self, num_nodes, method):
        self.dag = DAG(num_nodes)
        self.num_nodes = num_nodes
        self.solutions = self.dag.get_counts()
        self.adj_matrix = self.dag.get_adjacency_matrix()
        self.method = method
        self.bitstring1 = None
        self.bitstring2 = None
        
    def save_dag(self, folder, bitstring1=None, bitstring2=None, title="DAG Visualization", slope=None, intercept=None):
        """Saves the DAG visualization in the given folder."""
        fig, ax = plt.subplots(figsize=(8, 6))
        self.dag.color_graph(bitstring1=bitstring1, bitstring2=bitstring2, slope=slope, intercept=intercept, ax=ax)

        dag_path = os.path.join(folder, "dag.png")
        fig.savefig(dag_path, dpi=300, bbox_inches="tight")
        plt.close()
        
        # save the pickle of the graph
        with open(os.path.join(folder, "dag.pkl"), "wb") as f:
            pickle.dump(self, f)
        
    def save_histogram(self, counts, correct_bitstring, folder, title="Histogram of Bitstring Counts",
                   xlabel="Bitstrings", ylabel="Frequency", color="blue"):
        """Saves the histogram in the given folder."""

        fig, ax = plt.subplots(figsize=(7, 6))
        sorted_counts = sorted(counts.items(), key=lambda item: item[1], reverse=True)[:40]
        bitstrings, frequencies = zip(*sorted_counts)

        sns.barplot(x=list(bitstrings), y=list(frequencies), color=color, ax=ax)
        ax.set_xticklabels(bitstrings, rotation=90)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        
        correct_bitstring_str = "".join([str(int(x)) for x in correct_bitstring[::-1]])
        axe = fig.gca()
        for label in axe.get_xticklabels():
            if label.get_text() == correct_bitstring_str:
                label.set_fontweight('bold')

        hist_path = os.path.join(folder, "histogram.png")
        fig.savefig(hist_path, dpi=300, bbox_inches="tight")
        plt.close()
        
    def visualize_results(self, colors, counts, correct_bitstring):

        self.dag.plot_dag(colors=colors)
        self.plot_histogram(counts, correct_bitstring)

    def plot_histogram(self, counts, correct_bitstring=None, title="Histogram of Top 40 Bitstring Counts", xlabel="Bitstrings", ylabel="Frequency", color="blue"):
        sorted_counts = sorted(counts.items(), key=lambda item: item[1], reverse=True)[:50]
        bitstrings, frequencies = zip(*sorted_counts)
    
        plt.figure(figsize=(10, 6))
        sns.barplot(x=list(bitstrings), y=list(frequencies), color=color)
        plt.xticks(rotation=90)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        
        if correct_bitstring is not None:
            correct_bitstring_str = "".join([str(int(x)) for x in correct_bitstring[::-1]])
            ax = plt.gca()
            for label in ax.get_xticklabels():
                if label.get_text() == correct_bitstring_str:
                    label.set_fontweight('bold')
        
        plt.show()
    
    # def save_histogram(self, counts, filename="histogram.png", title="Histogram of Bitstring Counts",
    #                xlabel="Bitstrings", ylabel="Frequency", color="blue"):
    #     """Saves the histogram as an image."""
    #     fig, ax = plt.subplots(figsize=(7, 6))
    #     bitstrings = list(counts.keys())
    #     frequencies = list(counts.values())

    #     sns.barplot(x=bitstrings, y=frequencies, color=color, ax=ax)
    #     ax.set_xticklabels(bitstrings, rotation=90)
    #     ax.set_xlabel(xlabel)
    #     ax.set_ylabel(ylabel)
    #     ax.set_title(title)
    #     fig.savefig(filename, dpi=300, bbox_inches="tight")
    #     plt.close(fig) 

    def solve_qubo(self, method, qubit_mult, h=None):
        if h is not None:
            qubo_problem = QUBO_Problem(self.adj_matrix, self.num_nodes, method, h)
        else:
            qubo_problem = QUBO_Problem(self.adj_matrix, self.num_nodes, method)
        self.cubic_program = CubicProgram(qubo_problem.linear, qubo_problem.quadratic, [], qubit_mult*self.num_nodes)
        self.ising_hamiltonian = IsingHamiltonian(self.cubic_program)
        print(f"Number of linear terms: {len(self.cubic_program.linear)}, number of quadratic terms: {len(self.cubic_program.quadratic)}")

    def solve_exact(self):
        correct_bitstring, probability = self.ising_hamiltonian.run_numpy_eigensolver()
        return correct_bitstring, probability

    def solve_qaoa(self, correct_bitstring, p=4):
        correct_bitstring_str = "".join([str(int(x)) for x in correct_bitstring[::-1]])
        found_bitstring, probability, overlap = self.ising_hamiltonian.run_QAOA([correct_bitstring_str], p)
        found_bitstring = found_bitstring[::-1]
        print(f"Correct bitstring: {correct_bitstring}")
        print(f"Found bitstring:   {found_bitstring}, probability: {probability}, overlap: {overlap}")
        
        if self.method == "three-chain":
            xs = found_bitstring[:len(found_bitstring)//2]
            ys = found_bitstring[len(found_bitstring)//2:]
            self.bitstring1 = [1 if xs[i] and ys[i] else 0 for i in range(len(xs))]
            self.bitstring2 = None
        elif self.method == "three-chain-naive":
            self.bitstring1 = [1 if found_bitstring[i] else 0 for i in range(len(found_bitstring))]
            self.bitstring2 = None
        elif self.method == "two-link":
            xs = found_bitstring[:len(found_bitstring)//3]
            ys = found_bitstring[len(found_bitstring)//3:2*len(found_bitstring)//3]
            zs = found_bitstring[2*len(found_bitstring)//3:]
            self.bitstring1 = [1 if xs[i] else 0 for i in range(len(xs))]
            self.bitstring2 = [1 if ys[i] else 0 for i in range(len(ys))]
        elif self.method == "bh-links":
            xs = found_bitstring[:len(found_bitstring)//2]
            ys = found_bitstring[len(found_bitstring)//2:]
            self.bitstring1 = [1 if xs[i] else 0 for i in range(len(xs))]
            self.bitstring2 = [1 if ys[i] else 0 for i in range(len(ys))]
        else:
            print("Method not supported!")
            return
        
        #self.dag.color_graph(colors)
        
        counts = self.ising_hamiltonian.counts
        for bitstring in sorted(counts, key=counts.get, reverse=True)[:5]:
            print(f"{bitstring[::-1]}: {counts[bitstring]}")
        #self.visualize_results(colors, counts, correct_bitstring)

# Example usage
method = "bh-links"
size = 6


# dag_viz = DAGVisualization(num_nodes=size, method=method)
# while dag_viz.dag.bh_links < 2:
#     dag_viz = DAGVisualization(num_nodes=size, method=method)
# load dag from file
with open("vis_results/bh-links_size_6_13/dag.pkl", "rb") as f:
    dag_viz = pickle.load(f)


if method == "three-chain-naive-2":
    qubit_mult = 2
elif method == "three-chain-naive":
    qubit_mult = 1
elif method == "three-chain":
    qubit_mult = 2
elif method == "two-link":
    qubit_mult = 3*dag_viz.num_nodes
elif method == "bh-links":
    qubit_mult = 2
    h = dag_viz.dag.h

dag_viz.dag.plot_dag()
if method == "bh-links":
    dag_viz.solve_qubo(method, qubit_mult, h)
else:
    dag_viz.solve_qubo(method, qubit_mult)
    
correct_bitstring, probability = dag_viz.solve_exact()
dag_viz.solve_qaoa(correct_bitstring)
found = False
while not found:
    correct_bs1 = correct_bitstring[:len(correct_bitstring)//2]
    correct_bs2 = correct_bitstring[len(correct_bitstring)//2:]
    
    for i in range(len(correct_bs1)):
        if correct_bs1[i] != dag_viz.bitstring1[i]:
            break
        if correct_bs2[i] != dag_viz.bitstring2[i]:
            break
        if i == len(correct_bs1)-1:
            found = True
    if not found:
        dag_viz.solve_qaoa(correct_bitstring)

output_folder = get_new_folder(size=size, method=method)
print(dag_viz.bitstring1, dag_viz.bitstring2)
if method == "bh-links":
    dag_viz.save_dag(output_folder, dag_viz.bitstring1, dag_viz.bitstring2, slope=dag_viz.dag.slope, intercept=dag_viz.dag.intercept)
else:
    dag_viz.save_dag(output_folder, dag_viz.bitstring1, dag_viz.bitstring2)
dag_viz.save_histogram(dag_viz.ising_hamiltonian.counts, correct_bitstring, output_folder)
