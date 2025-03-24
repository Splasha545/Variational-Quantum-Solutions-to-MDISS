# This script contains the main functions for the experiments conducted in the report. 
# The main function is called at the end of the script and can be used to run the experiments.
# The functions in this script are used to generate the problems, run the experiments, and visualise the results.

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from collections import Counter

from VQALib import CubicProgram, IsingHamiltonian
from VQATests import tests

from qiskit.visualization import plot_histogram

# ------------ GRAPH FUNCTIONS -----------
def draw_graph(G, colors, pos, ax=None):
    default_axes = plt.axes(frameon=True) if ax is None else ax
    nx.draw_networkx(G, node_color=colors, node_size=600, alpha=0.8, ax=default_axes, pos=pos, font_size=14)
    edge_labels = nx.get_edge_attributes(G, "weight")
    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels)

def compute_matrix(n, G):
    # Computing the weight matrix from the random graph
    w = np.zeros([n, n])
    for i in range(n):
        for j in range(n):
            temp = G.get_edge_data(i, j, default=0)
            if temp != 0:
                w[i, j] = 1
    return w

def color_graph(bitstring, n, m, H, posH, ax=None):
    print(bitstring)
    colors = ["r" for node in range(m)]
    for k in range(n*m):
        i = k // m  # Integer division to retrieve the i value
        j = k % m   # Modulus operation to retrieve the j value
        if bitstring[k] == 1:
            colors[j] = "b"
    draw_graph(H, colors, posH, ax)

def create_graph(n, elistG, m, elistH):
    # Target graph:
    H = nx.DiGraph()
    H.add_nodes_from(np.arange(0, m, 1))
    H.add_edges_from(elistH)
    colorsH = ["r" for node in H.nodes()]
    posH = nx.spring_layout(H)

    # Pattern graph:
    G = nx.DiGraph()
    G.add_nodes_from(np.arange(0, n, 1))
    G.add_edges_from(elistG)
    colorsG = ["r" for node in G.nodes()]
    posG = nx.spring_layout(H)

    return G, colorsG, posG, H, colorsH, posH

def display_graph(n, elistG, m, elistH, bitstring):
    # Target graph:
    H = nx.DiGraph()
    H.add_nodes_from(np.arange(0, m, 1))
    H.add_edges_from(elistH)
    colorsH = ["w" for node in H.nodes()]
    posH = nx.spring_layout(H)

    # Pattern graph:
    G = nx.DiGraph()
    G.add_nodes_from(np.arange(0, n, 1))
    G.add_edges_from(elistG)
    colorsG = ["c" for node in G.nodes()]
    posG = nx.spring_layout(H)

    colors = ["r" for node in range(m)]
    for k in range(n*m):
        i = k // m  # Integer division to retrieve the i value
        j = k % m   # Modulus operation to retrieve the j value
        if bitstring[k] == '1':
            print(f"coloring {j}!")
            colors[j] = "c"

    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    draw_graph(G, colorsG, posG, axes[0])
    draw_graph(H, colors, posH, axes[1])
    plt.tight_layout()
    plt.savefig("graph.pdf")
    plt.show()

# ------------ ARITHMETIC FUNCTIONS ------
def expand_arithmetic_function(mat1, mat2, n, m):
    """
    Expands and groups terms in the CUBO cost function based on the given matrices and dimensions.

    Args:
        mat1 (list): The first matrix.
        mat2 (list): The second matrix.
        n (int): The number of vertices in the pattern graph.
        m (int): The number of vertices in the target graph.

    Returns:
        tuple: A tuple containing three lists: linear, quadratic, and cubic.
            - linear: A list of linear terms and their coefficients.
            - quadratic: A list of quadratic terms and their coefficients.
            - cubic: A list of cubic terms and their coefficients.
    """

    linear = []
    quadratic = []
    cubic = []

    penalty = n

    def translate(i, j):
        return (i*m+j)

    # FIRST TERM:
    sum = 0
    for i in range(n):
        for u in range(m):
            linear.append((translate(i, u), -1))

    # SECOND TERM (EXPANDED):
    sum2 = 0
    for v in range(m):
        inner = 0
        for i in range(n):
            for k in range(n):
                quadratic.append( ( (translate(i, v), translate(k, v) ), 1) )   # ( (var1, var2),  coeff )
        
        for i in range(n):
            linear.append((translate(i, v), -1))
    

    sum2 = 0
    # THIRD TERM (CHECKS THAT NUMBER OF OUTGOING EDGES OF EACH VERTEX U IS EQUAL TO OUTGOING EDGES OF EACH VERTEX I):
    for u in range(m):
        for i in range(n):
            for j in range(n):
                # Add the linear terms:
                linear.append((translate(i, u), penalty*mat1[i][j]*mat1[i][j]))

                for v in range(m):
                    # Add quadratic terms:
                    quadratic.append( ( (translate(j,v), translate(i, u) ), -1*2*penalty*mat1[i][j]*mat2[u][v]) )
                
                for v in range(m):
                    for k in range(m):
                        # Add cubic terms:
                        cubic.append( ( ( translate(j,v), translate(j, k), translate(i,u) ), penalty*mat2[u][v]*mat2[u][k] ) )
            

    # FOURTH TERM (CHECKS THAT NUMBER OF INCOMING EDGES OF EACH VERTEX U IS EQUAL TO OUTGOING EDGES OF EACH VERTEX I):
    for u in range(m):
        for i in range(n):
            for j in range(n):
                # Add the linear terms:
                linear.append((translate(i, u), penalty*mat1[j][i]*mat1[j][i]))

                for v in range(m):
                    # Add quadratic terms:
                    quadratic.append( ( (translate(j,v), translate(i, u) ), -1*2*penalty*mat1[j][i]*mat2[v][u]) )
                
                for v in range(m):
                    for k in range(m):
                        # Add cubic terms:
                        cubic.append( ( ( translate(j,v), translate(j, k), translate(i,u) ), penalty*mat2[v][u]*mat2[k][u] ) )
    
    return linear, quadratic, cubic

def clean_coefficients(linear, quadratic, cubic):
    """Removes all the terms with coefficient 0 from the function.

    Args:
        linear (_type_): linear terms
        quadratic (_type_): quadratic_terms
        cubic (_type_): cubic terms
    """    
    def drop(list):
        to_remove = []
        for idx, (term, coeff) in enumerate(list):
            if coeff == 0:
                to_remove.append((term, coeff))
        
        for item in to_remove:
            list.remove(item)
        return list
    
    linear = drop(linear)
    quadratic = drop(quadratic)
    cubic = drop(cubic)
    
    return linear, quadratic, cubic

def collect_terms(linear, quadratic, cubic):
    """Sums all identical terms together to reduce function to strictly unique terms with real coefficients.
    """    
    def combine(tuple_list):
        combined_tuples = {}
        
        for tpl, coeff in tuple_list:
            # Convert the tuple to a hashable format using sorted tuple of its Counter elements
            if isinstance(tpl, int):
                tpl = (tpl,)
            hashable_tpl = tuple(sorted(Counter(tpl).items()))
            # Check if the hashable tuple already exists in the dictionary
            if hashable_tpl in combined_tuples:
                # Accumulate the coefficient for the existing tuple
                combined_tuples[hashable_tpl] += coeff
            else:
                # Create a new entry in the dictionary
                combined_tuples[hashable_tpl] = coeff
        
        # Convert the dictionary back to the original format of list of tuples
        reduced_list = [( tuple(x for x, count in dict_key for _ in range(count)) , coeff) for dict_key, coeff in combined_tuples.items()]
        
        return reduced_list

    
    collected_linear_terms = combine(linear)
    collected_quadratic_terms = combine(quadratic)
    collected_cubic_terms = combine(cubic)

    return collected_linear_terms, collected_quadratic_terms, collected_cubic_terms

def produce_terms_from_graphs(mat1, mat2, n, m):
    linear, quadratic, cubic = expand_arithmetic_function(mat1, mat2, n, m)
    linear, quadratic, cubic = clean_coefficients(linear, quadratic, cubic)
    linear, quadratic, cubic = collect_terms(linear, quadratic, cubic)
    return linear, quadratic, cubic

def produce_cubic_problem_from_graph(n, elistG, m, elistH):
    # Create the CUBO formula from the problem instance
    G, colorsG, posG, H, colorsH, posH = create_graph(n, elistG, m, elistH)
    linear, quadratic, cubic = produce_terms_from_graphs(compute_matrix(n, G), compute_matrix(m, H), n, m)
    cubic_instance = CubicProgram(linear, quadratic, cubic, n*m)
    return cubic_instance, n, m, H, posH

# ------------ EXPERIMENT FUNCTIONS ------
import timeit
def run_experiment(n, elistG, m, elistH, VQE=False, QAOA=False, correct_bitstrings=[''], quadratic=False, runs=20, depth=4):
    """
    Runs an experiment for quantum variational algorithms (VQE or QAOA) or numpy eigensolver.

    Parameters:
    - n (int): Number of nodes in the pattern graph.
    - elistG (list): List of edges in the pattern graph.
    - m (int): Number of nodes in the target graph.
    - elistH (list): List of edges in the target graph.
    - VQE (bool): Flag indicating whether to run VQE algorithm. Default is False.
    - QAOA (bool): Flag indicating whether to run QAOA algorithm. Default is False.
    - correct_bitstrings (list): List of correct bitstrings. Default is [''].
    - quadratic (bool): Flag indicating whether to reduce the cubic instance to quadratic. Default is False.
    - runs (int): Number of runs for the experiment. Default is 20.
    - depth (int): Depth parameter for VQE or QAOA algorithms. Default is 4.

    Returns:
    - avg_overlap (float): Average performance indicating $p_{succ}$ of the experiment.
    - all_overlaps (list): List of performances $p_{succ}$ for each run of the experiment.
    """
    
    overlaps = []
    timing = []
    bitstrings = []
    all_overlaps = []
    
    avg_overlap = 0
    avg_timing = 0
    cubic_instance, n, m, H, posH = produce_cubic_problem_from_graph(n, elistG, m, elistH)

    if quadratic:
        cubic_instance.reduce_cubic()

    for i in range(runs):
        print("Current iteration: ", i)
        hamiltonian = IsingHamiltonian(cubic_instance)

        # Time the function execution
        start_time = timeit.default_timer()
        if VQE:
            likeliest_bitstring, highest_probability, overlap = hamiltonian.run_VQE(correct_bitstrings=correct_bitstrings, depth=depth)
        elif QAOA:
            likeliest_bitstring, highest_probability, overlap = hamiltonian.run_QAOA_custom(correct_bitstrings=correct_bitstrings, qaoa_reps=depth)
        else:
            likeliest_bitstring, overlap = hamiltonian.run_numpy_eigensolver()
        end_time = timeit.default_timer()
        
        timing.append(end_time - start_time)
        overlaps.append(overlap)
        bitstrings.append(likeliest_bitstring)
        all_overlaps.append(overlaps)
        #print(likeliest_bitstring)
    
    avg_overlap = np.average(overlaps)
    avg_timing = np.average(timing)
    print(f"Iteration likeliest bitstrings: {bitstrings}")
    # print(f"Iteration overlaps: {overlaps}")
    # print(f"Iteration timings: {timing}")
    print(f"Average overlap: {avg_overlap}, Average timing: {avg_timing}")
    return avg_overlap, all_overlaps

def basic_runthrough(n, elistG, m, elistH):
    """
    Runs a basic runthrough of the quantum algorithms for a given problem.

    Args:
        n (int): The number of nodes in the graph.
        elistG (list): The edge list of the graph G.
        m (int): The number of nodes in the hypergraph H.
        elistH (list): The edge list of the hypergraph H.

    Returns:
        None

    Prints the results of running the numpy eigensolver, QAOA, and VQE algorithms on the given problem.
    """
    
    cubic_instance, n, m, H, posH = produce_cubic_problem_from_graph(n, elistG, m, elistH)
    hamiltonian = IsingHamiltonian(cubic_instance)
    print("NUMPY EIGENSOLVER:")
    correct_bitstring = hamiltonian.run_numpy_eigensolver()
    print(''.join(map(lambda x: str(int(x)), correct_bitstring)))

    print(f"Energy: {hamiltonian.energy}")
    print(f"Optimal value: {hamiltonian.optimal_value}")
    print(f"Probability: {hamiltonian.max_prob}")
    print(f"Bitstring: {hamiltonian.bitstring}")

    print("QAOA:")
    bitstring, probabilities = hamiltonian.run_QAOA(correct_bitstring=''.join(map(lambda x: str(int(x)), correct_bitstring)))
    print(f"Energy: {hamiltonian.energy}")
    print(f"Probability of correct result: {hamiltonian.correct_prob}")
    print(f"Highest probability bitstring: {bitstring} and its probability: {hamiltonian.max_prob}")

    print("VQE:")
    hamiltonian.run_VQE()
    print(f"Energy: {hamiltonian.energy}")
    print(f"Optimal value: {hamiltonian.optimal_value}")
    print(f"Probability: {hamiltonian.max_prob}")
    print(f"Bitstring: {hamiltonian.bitstring}")

def graphs_run_through():
    """
    Runs through a set of test cases and displays the resulting graphs.
    """
    test_pkg = tests()
    for num in test_pkg.qubit_nums[7:8]:
        for problem in test_pkg.problems[num]:
            G, colorsG, posG, H, colorsH, posH = create_graph(problem['n'], problem['elistG'], problem['m'], problem['elistH'])
            for bitstring in problem['correct_bitstrings']:
                bitarray = [1 if x == '1' else 0 for x in bitstring]
                color_graph(bitarray[::-1], problem['n'], problem['m'], H, posH)
                plt.tight_layout()
                plt.show()

def difficulty_scaling(reps=20, VQE=False, QAOA=False, quadratic=False, low=0, high=9):
    """
    Runs the same algorithm on increasingly difficult problems.

    Args:
        reps (int): The number of times to run the experiment for each problem.
        VQE (bool): Whether to use the Variational Quantum Eigensolver (VQE) algorithm.
        QAOA (bool): Whether to use the Quantum Approximate Optimization Algorithm (QAOA).
        quadratic (bool): Whether to use a quadratic optimization problem.
        low (int): The lower bound of the range of qubit numbers to consider.
        high (int): The upper bound of the range of qubit numbers to consider.

    Returns:
        list: A list of tuples, where each tuple contains the average probability and the probability spread for a specific problem.

    """
    test_pkg = tests()
    # Run the same algorithm on increasingly difficult problems
    probability_scaling = []
    for num in test_pkg.qubit_nums[low:high]:
        print(f"Current qubit number: {num}")
        problem_avg_probabilities = []
        problem_probability_spread = []
        for problem in test_pkg.problems[num]:
            n = problem['n']
            elistG = problem['elistG']
            m = problem['m']
            elistH = problem['elistH']
            correct_bitstrings = problem['correct_bitstrings']
            avg_prob, all_probabilities = run_experiment(n, elistG, m, elistH, QAOA=QAOA, VQE=VQE, correct_bitstrings=correct_bitstrings, quadratic=quadratic, runs=reps)
            problem_avg_probabilities.append(avg_prob)
            problem_probability_spread.append(all_probabilities)
        probability_scaling.append((np.average(problem_avg_probabilities), problem_probability_spread))
    return probability_scaling

from qiskit_aer import AerSimulator
from qiskit.circuit.library import QAOAAnsatz
def complexity_estimation(n, elistG, m, elistH, correct_bitstrings, quadratic=False, iterations=10, low=0, high=9, depth=4):
    """
    Estimates the time complexity of QAOA by running QAOA on increasingly complex problem instances, while training on a simple problem provided.
    Follows the methodology established by Sami Boulebnane and Ashley Montanaro "Solving boolean satisfiability problems with the quantum 
    approximate optimization algorithm" (2022).

    Args:
        n (int): Number of nodes in the graph.
        elistG (list): List of edges in the graph.
        m (int): Number of edges in the graph.
        elistH (list): List of edges in the hypergraph.
        correct_bitstrings (list): List of correct bitstrings for the problem instances.
        quadratic (bool, optional): Flag indicating whether to reduce the problem to a quadratic form. Defaults to False.
        iterations (int, optional): Number of iterations to run the algorithm. Defaults to 10.
        low (int, optional): Lower bound of the range of qubit numbers to consider. Defaults to 0.
        high (int, optional): Upper bound of the range of qubit numbers to consider. Defaults to 9.
        depth (int, optional): Depth of the QAOA circuit. Defaults to 4.

    Returns:
        list: A list of tuples containing the average probability and all probabilities for each problem size.
    """
    
    test_pkg = tests()

    problem_instance, n, m, H, posH = produce_cubic_problem_from_graph(n, elistG, m, elistH)
    if quadratic:
        problem_instance.reduce_cubic()
    OG_hamiltonian = IsingHamiltonian(problem_instance)

    probability_scaling = []
    average_probabilities_accross_n_iterations = {}
    all_probabilities = {}
    for num in test_pkg.qubit_nums[low:high]:
        average_probabilities_accross_n_iterations[num] = []
        all_probabilities[num] = []

    # Repeat finding the optimal parameters and average probabilities n times:
    for i in range(iterations):
        print("Current iteration running: ", i)
        OG_hamiltonian.run_QAOA(correct_bitstrings=correct_bitstrings, qaoa_reps=depth)

       	# Memorise optimal parameters
        optimised_parameters = OG_hamiltonian.optimised_parameters

        # create the simulator:
        aer_simulator = AerSimulator(method="statevector")

        # Use these optimal parameters to find result for increasingly complex problems
        for num in test_pkg.qubit_nums[low:high]:
            print("Current problem size: ", num)
            print("\n\n")
            problem_probabilities = []
            for problem in test_pkg.problems[num][:2]:
                n = problem['n']
                elistG = problem['elistG']
                m = problem['m']
                elistH = problem['elistH']
                correct_bitstrings = problem['correct_bitstrings']

                cubic_instance, n, m, H, posH = produce_cubic_problem_from_graph(n, elistG, m, elistH)
                if quadratic:
                    cubic_instance.reduce_cubic()
                hamiltonian = IsingHamiltonian(cubic_instance)
                qaoa_ansatz = QAOAAnsatz(cost_operator=hamiltonian.qubit_op, reps=depth)
                qaoa_ansatz.measure_active()

                parameter_bindings = dict(zip(qaoa_ansatz.parameters, optimised_parameters))
                qaoa_with_parameters = qaoa_ansatz.assign_parameters(parameter_bindings)
                qaoa_with_parameters_decomposed = qaoa_with_parameters.decompose(reps=3)

                number_of_shots = 10*pow(2, hamiltonian.qubit_op.num_qubits)
                counts = aer_simulator.run(qaoa_with_parameters_decomposed, shots=number_of_shots).result().get_counts()
                probabilities = {bitstring: (count/number_of_shots) \
                                    for bitstring, count in counts.items()}
                correct_prob = 0
                for correct_bitstring in correct_bitstrings:
                    if correct_bitstring in probabilities.keys():
                        correct_prob += probabilities[correct_bitstring]
                if correct_prob == 0:
                    correct_prob = 1/(2*pow(2, hamiltonian.qubit_op.num_qubits))
                # print(f"Problem {num}, iteration: {i}, Prob: {correct_prob}")
                # Stores results for 3 problems at current number of qubits
                problem_probabilities.append(correct_prob)
            average_accross_3_problems = np.average(problem_probabilities)
            average_probabilities_accross_n_iterations[num].append(average_accross_3_problems)
            all_probabilities[num].append(problem_probabilities)
    
    # For each qubit number, find the average probability:
    for num in test_pkg.qubit_nums[low:high]:
        probability_scaling.append((np.average(average_probabilities_accross_n_iterations[num]), np.array(all_probabilities[num]).flatten()))
    return probability_scaling
        

def main():
    """
    This is the main function that serves as the entry point of the program.
    It contains various commented out code snippets that can be uncommented and executed to perform different tasks.
    """
    test_pkg = tests()
    # Uncomment and modify the code snippets below to perform different tasks
    
    # Example 1: Run an experiment with specific parameters
    # problem = test_pkg.problems[12][1]
    # n = problem['n']
    # elistG = problem['elistG']
    # m = problem['m']
    # elistH = problem['elistH']
    # correct_bitstrings = problem['correct_bitstrings']
    # print(run_experiment(n, elistG, m, elistH, QAOA=True, VQE=False, correct_bitstrings=correct_bitstrings, quadratic=False, runs=10, depth=6))
    
    # Example 2: Visualize a problem instance
    problem = test_pkg.problems[8][2]
    n = problem['n']
    elistG = problem['elistG']
    m = problem['m']
    elistH = problem['elistH']
    correct_bitstrings = problem['correct_bitstrings']
    cubic_instance, n, m, H, posH = produce_cubic_problem_from_graph(n, elistG, m, elistH)
    display_graph(n, elistG, m, elistH, correct_bitstrings[0][::-1])
    color_graph(correct_bitstrings[0], n, m, H, posH)
    
    # Example 3: Iterate over test package problems
    # for num in test_pkg.qubit_nums[:10]:
    #     print("Problem size: ", num)
    #     added_list = []
    #     for problem in test_pkg.problems[num]:
    #         n = problem['n']
    #         elistG = problem['elistG']
    #         m = problem['m']
    #         elistH = problem['elistH']
    #         correct_bitstrings = problem['correct_bitstrings']
    #         cubic_instance, n, m, H, posH = produce_cubic_problem_from_graph(n, elistG, m, elistH)
    #         print(f"QUBIT NUMBER Before: {cubic_instance.num_vars}")
    #         cubic_instance.reduce_cubic()
    #         print(f"QUBIT NUMBER After: {cubic_instance.num_vars}")
    #         added_list.append(cubic_instance.num_vars)
    #     print(f"Average number of qubits added for size {num}: {np.average(added_list)}")
    
    # Example 4: Run an experiment with VQE
    # problem = test_pkg.problems[14][2]
    # n = problem['n']
    # elistG = problem['elistG']
    # m = problem['m']
    # elistH = problem['elistH']
    # correct_bitstrings = problem['correct_bitstrings']
    # print(run_experiment(n, elistG, m, elistH, VQE=True, correct_bitstrings=correct_bitstrings, runs=1))
    
    # Example 5: Compare performance of cubic VQE vs quadratic VQE
    # print("Cubic:")
    # print(difficulty_scaling(reps=10, QAOA=True, quadratic=False, low=6, high=7))
    # print("Quadratic:")
    # print(difficulty_scaling(reps=5, VQE=True, quadratic=True, low=6, high=7))
    
    # Example 6: Estimate complexity of QAOA
    # print(complexity_estimation(2, [(0, 1)], 3, [(0, 1), (1, 2), (2, 0)], ['010001', '100010', '001100'], iterations=40, quadratic=True, low=0, high=7))
    # print(complexity_estimation(3, [(0, 1), (1, 2), (0, 2)], 4, [(0, 1), (1, 2), (0, 2), (3, 1), (3, 2)], ['010000100001', '010000101000'], iterations=40, quadratic=True, low=0, high=7))
    
if __name__ == "__main__":
    # Call the main function if the script is executed directly
    main()


