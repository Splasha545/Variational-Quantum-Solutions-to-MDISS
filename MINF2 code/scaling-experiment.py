# QAOA complexity scaling function:
from VQALib import CubicProgram, IsingHamiltonian
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error, thermal_relaxation_error, phase_damping_error
from qiskit.providers.fake_provider import FakeAthens
from qiskit.circuit.library import QAOAAnsatz
import numpy as np
from RNG_DAG import DAG
from QUBO_FUNCTIONS import QUBO_Problem
import pickle

def run_aer_simulator(graph, function, aer_simulator, optimised_parameters, depth, qubit_mult, h=None):
    problem = QUBO_Problem(graph.get_adjacency_matrix(), graph.N, function, h)
    program = CubicProgram(problem.linear, problem.quadratic, [], qubit_mult*graph.N)
    hamiltonian = IsingHamiltonian(program)
    
    correct_bitstring, _ = hamiltonian.run_numpy_eigensolver()
    correct_bitstrings = ["".join([str(int(x)) for x in correct_bitstring[::-1]])]
    
    qaoa_ansatz = QAOAAnsatz(hamiltonian.qubit_op, reps=depth)
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
        
    return counts, correct_prob

def pre_train_qaoa_parameters(function, qubit_mult, depth, h=None):
    # find the parameters that perfrom the best on average:
    training_graphs = []
    with open("./generated_graphs/dags_training_9.pkl", "rb") as f:
        training_graphs = pickle.load(f)

    best_parameters = []
    best_average_overlap = 0

    for iteration in range(10):
        print(f"    Starting training iteration {iteration+1}/10")
        random_dag = DAG(9)
        
        random_problem = QUBO_Problem(random_dag.get_adjacency_matrix(), random_dag.N, function, h)
        random_program = CubicProgram(random_problem.linear, random_problem.quadratic, [], qubit_mult*random_dag.N)
        random_ising = IsingHamiltonian(random_program)
        
        correct_bitstring, _ = random_ising.run_numpy_eigensolver()
        correct_bitstring_str = "".join([str(int(x)) for x in correct_bitstring[::-1]])
        random_ising.run_QAOA([correct_bitstring_str], depth)
        
        print(f"    Random DAG generated, QAOA run")
        
        parameters = random_ising.optimised_parameters
        aer_simulator = AerSimulator(method="statevector")
        
        average_overlap = 0
        for graph in training_graphs:
            print(f"            Processing graph {training_graphs.index(graph)+1}/{len(training_graphs)}")
            _, overlap = run_aer_simulator(graph, function, aer_simulator, parameters, depth, qubit_mult)
            average_overlap += overlap
            
        average_overlap /= len(training_graphs)
        
        if average_overlap > best_average_overlap:
            best_parameters = parameters
            best_average_overlap = average_overlap
        
    return best_parameters
    
    
    

def QAOA_scaling(function, iterations, qubit_mult, depth, h=None, noisy=False, parameters=None, graph_sizes=None):
    if graph_sizes is None:
        graph_sizes = [4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 18]

    overlaps_dict = {}
    for size in graph_sizes:
        overlaps_dict[str(size)] = []
    
    for i in range(iterations):
        print(f"Starting iteration {i+1}/{iterations}")
        # denote the beginning of a new iteration in the txt file:
        with open("./experiment_results/" + function + "_scaling_with_pre_training.txt", "a") as f:
            f.write(f"-------------------------------------------------- \n")
            f.write(f"Iteration: {i} \n")
            f.write(f"-------------------------------------------------- \n")
        
        # generate a graph that the QAOA parameters will be trained on
        # training_graph = DAG(9)
        
        # # create a QUBO problem out of the training graph
        # training_problem = QUBO_Problem(training_graph.get_adjacency_matrix(), training_graph.N, function)
        # training_program = CubicProgram(training_problem.linear, training_problem.quadratic, [], qubit_mult*training_graph.N)
        # training_ising = IsingHamiltonian(training_program)
        
        # # solve the problem
        # correct_bitstring, _ = training_ising.run_numpy_eigensolver()
        # correct_bitstring_str = "".join([str(int(x)) for x in correct_bitstring[::-1]])
        # training_ising.run_QAOA([correct_bitstring_str], depth)
        
        # # memorise optimal parameters
        # optimised_parameters = training_ising.optimised_parameters
        # if i == 0:
        #     optimised_parameters = [-2.54630331,  1.36924767, -2.86133025,  1.84223595, -0.44556164,  6.26947132, -0.28154989,  2.35246692]
        # else:
        if parameters is not None:
            optimised_parameters = parameters[i]
        else:
            optimised_parameters = pre_train_qaoa_parameters(function, qubit_mult, depth, h)
        # save parameters
        with open("./experiment_results/" + function + "_scaling_with_pre_training.txt", "a") as f:
            f.write(f"Optimised parameters: {optimised_parameters} \n")
        
        # create the simulator:
        if noisy:
            aer_simulator = AerSimulator(method="statevector", noise_model=noise_model)
        else:
            aer_simulator = AerSimulator(method="statevector")
        
        # use the optmised parameters to find the average overlap for each set of graphs:
        average_overlaps = []
        
        for size in graph_sizes:
            print(f"Processing graph size {size}")
            # retrieve the graph set
            graphs = []
            with open("./generated_graphs/dags_" + str(size) + ".pkl", "rb") as f:
                graphs = pickle.load(f)
            
            overlaps = []
            counter = 0
            for graph in graphs:
                counter += 1
                print(f"Processing graph {counter}/{len(graphs)}")
                counts, overlap = run_aer_simulator(graph, function, aer_simulator, optimised_parameters, depth, qubit_mult, h)
                
                overlaps.append(overlap)
                
            average_overlap = np.average(overlaps)
            
            average_overlaps.append(average_overlap)
            overlaps_dict[str(size)].append(average_overlap)
            
            # save the average overlaps to a txt file:
            with open("./experiment_results/" + function + "_scaling_with_pre_training.txt", "a") as f:
                f.write(f"Size: {str(size)}, Average Overlap: {average_overlap} \n")
        
        # save the average overlaps for this iteration to a txt file:
        with open("./experiment_results/" + function + "_scaling_with_pre_training.txt", "a") as f:
            f.write(f"Iteration: {i}, Average Overlaps: {average_overlaps} \n")
            f.write(f"-------------------------------------------------- \n\n")
    
    # save the average overlaps for each size to a txt file:
    with open("./experiment_results/" + function + "_scaling_with_pre_training.txt", "a") as f:
        f.write(f"Final overlaps dictionary: {overlaps_dict} \n")
        f.write(f"-------------------------------------------------- \n\n")
    return overlaps_dict

# backend = FakeAthens()
noise_model = NoiseModel() #NoiseModel.from_backend(backend)

# depolarizing noise on CNOTs:
error = depolarizing_error(0.05, 2)
noise_model.add_all_qubit_quantum_error(error, ['cx'])

# phase damping on Z gates
phase_damping = phase_damping_error(0.005)
noise_model.add_all_qubit_quantum_error(phase_damping, ['rz', "u1", "u2", "u3"])

# naive params
parameters_naive = [
    [-5.73883929, -1.73403527, -3.39915904,  1.25425761, -1.23849634,  2.04132594, -4.32219188,  4.19536475],
    [ 5.5736363,   6.29121998,  0.27725664, -6.43941542,  0.86606176, -1.7245623,  -0.57780193,  2.55813336] ,
    [ 4.80039407, -2.54129874,  3.72857594, -0.44184796, -5.63032638,  5.1546264, -0.21119081,  0.86423407] ,
    [ 6.16421064, -1.01846299,  4.44226166, -4.89053611, -3.79922156,  5.1551642,  2.29861678, -2.26376097] ,
    [ 3.17857094,  5.41633313, -0.39385122, -3.38095973, -1.52156093,  1.81377439,  0.52740284,  5.44745955] ,
    [-5.19234035,  3.84243972,  2.0430767,   2.24625403, -5.88352017,  1.42723849, -0.62469217, -5.05276646] ,
    [ 2.07788376, -2.63538882,  6.62831901, -2.78695967,  5.92955412,  5.28568278, -4.2554975,   2.89899599] ,
    [-5.75951177,  2.25046307, -4.02374824, -2.58451344, -3.46427135, -1.43149814,  3.62264527, -3.74157956] ,
    [ 0.5619902,   1.14410083,  4.74143008,  3.39250622, -1.41038779,  3.48792533,  2.77529294, -2.70180528] ,
    [-1.40686994, -1.85655089, -0.67687207, -1.8131774,  -2.77872093,  1.73568055,  1.25279447,  1.92920335]
]

parameters_naive_depth_6 = [
    [ 0.19005505,  2.94729113, -3.90736899,  3.48036523,  4.41356948,  1.49705241, -0.59451874, -5.33268887,  6.14635224, -0.57933956,  4.96683671,  3.20728206],
    [ 4.03568383,  1.5543011,  -0.61379774, -3.72219086,  2.42568595,  1.45692186, -6.2354162,  -4.95862339, -1.33580466, -5.38906659, -6.32547957, -1.40794064],
    [-3.58933847,  0.13512447, -1.75227329, -6.26069414,  0.39561861,  4.38412945,  0.39600439, -3.15683719, -1.7618681,   3.19692044, -4.60671771, -4.18612353],
    [ 0.24940229,  2.05280503,  6.83579662, -5.2731325,  -4.70250323, -4.48272821,  3.16387094, -6.56070228, -1.29459834, -1.89059551,  1.14698835,  1.87755104],
    [ 6.75219663, -1.82521645, -4.6291548,  -4.7030089,  -0.44222437, -0.87753625, -0.60596097, -3.71153676,  1.17398555,  0.38473744, -5.48596953,  3.83183494] 
]

# overlaps_dict = QAOA_scaling("three-chain", iterations=5, qubit_mult=2, depth=4, noisy=True)
# print(overlaps_dict)

# overlaps_dict = QAOA_scaling("three-chain-naive", iterations=10, qubit_mult=1, depth=4, noisy=False, parameters=parameters, graph_sizes=None)
# print(overlaps_dict)

# quadratic-params
parameters_quadratic = [
    [-3.14963899,  2.37372039,  6.21602451,  4.64048058, -1.84439194,  1.49296986, -2.65030634, -0.58557241],
    [-0.41604211, -6.34751201, -2.31432064,  0.09771495, -5.95360312,  2.93785202,  2.79834634, -3.64775963],
    [ 3.98113616, -2.36702703,  5.85399156,  2.21353208, -0.5851378,   2.9295248,  1.05013407, -0.99567032],
    [-0.46088283,  5.03640903,  2.61379914, -0.23354031, -3.52008869, -3.18199488,  5.52510862, -0.74357528],
    [ 0.7463327,  -0.05041994,  6.23389176, -3.0982337,  -0.37256875,  3.20531352,  3.85623815,  1.74194731]
]

parameters = None

# overlaps_dict = QAOA_scaling("bh-links", iterations=5, qubit_mult=2, depth=4, noisy=False, parameters=parameters, graph_sizes=[4, 5, 6, 7, 8, 9, 10])
# print(overlaps_dict)
# overlaps_dict = QAOA_scaling("three-chain-naive-2", iterations=5, qubit_mult=2, depth=4)
# print(overlaps_dict)

# noisy experiments:
# 3-chain naive with low p
overlaps_dict = QAOA_scaling("three-chain-naive", iterations=5, qubit_mult=1, depth=6, noisy=True, parameters=parameters_naive_depth_6, graph_sizes=[4, 5, 6, 7, 8, 9, 10, 12, 14])
print(overlaps_dict)

# 3-chain naive with high p
# overlaps_dict = QAOA_scaling("three-chain-naive", iterations=5, qubit_mult=1, depth=6, noisy=True, parameters=None, graph_sizes=[4, 5, 6, 7, 8, 9, 10, 12, 14])
# print(overlaps_dict)

# 3-chain quadratic with low p
# overlaps_dict = QAOA_scaling("three-chain", iterations=5, qubit_mult=2, depth=2, noisy=True, parameters=parameters_quadratic, graph_sizes=[4, 5, 6, 7])
# print(overlaps_dict)