# This file contains the implementation of the CubicProgram and IsingHamiltonian classes, which are used
# to represent a cubic program and its corresponding Ising Hamiltonian, respectively. The IsingHamiltonian 
# class also contains methods to run the VQE and QAOA algorithms on the Hamiltonian.

from typing import Dict, List, Optional, Tuple, Union, cast
import numpy as np
from collections import Counter

# Necessary imports for the VQE and QAOA algorithms
from qiskit.quantum_info import Pauli, SparsePauliOp
from qiskit.quantum_info import Statevector
from qiskit.result import QuasiDistribution

# VQE libraries
from qiskit.algorithms.optimizers import SPSA, SLSQP
from qiskit.circuit.library import EfficientSU2
from qiskit.algorithms.minimum_eigensolvers import NumPyMinimumEigensolver, VQE
from qiskit.primitives import Estimator
from qiskit import Aer, QuantumCircuit, transpile, assemble

# QAOA libraries
from qiskit_aer.primitives import Sampler as AerSampler
from qiskit.algorithms.optimizers import COBYLA
from qiskit.algorithms.minimum_eigensolvers import QAOA
from qiskit.opflow import OperatorStateFn
from qiskit.visualization import plot_histogram
from qiskit_aer import AerSimulator
from qiskit.circuit.library import QAOAAnsatz
from qiskit.opflow import PauliSumOp


class CubicProgram:
    """
    Represents a cubic program with linear, quadratic, and cubic terms.

    Attributes:
        linear (list): A list of tuples representing the linear terms of the program.
        quadratic (list): A list of tuples representing the quadratic terms of the program.
        cubic (list): A list of tuples representing the cubic terms of the program.
        num_vars (int): The number of variables in the program.

    Methods:
        get_num_vars(): Returns the number of variables in the program.
        set_num_vars(n): Sets the number of variables in the program.
        set_linear(linear_vars_list): Sets the linear terms of the program.
        set_quadratic(quadratic_vars_list): Sets the quadratic terms of the program.
        set_cubic(cubic_vars_list): Sets the cubic terms of the program.
        reduce_cubic(): Reduces the cubic terms of the program to quadratic terms by replacing common pairs.
    """

    def __init__(self, linear, quadratic, cubic, num_vars):
        """
        Initializes a CubicProgram object.

        Args:
            linear (list): A list of tuples representing the linear terms of the program.
            quadratic (list): A list of tuples representing the quadratic terms of the program.
            cubic (list): A list of tuples representing the cubic terms of the program.
            num_vars (int): The number of variables in the program.
        """
        self.linear = linear
        self.quadratic = quadratic
        self.cubic = cubic
        self.num_vars = num_vars

    def get_num_vars(self):
        return self.num_vars
    
    def set_num_vars(self, n):
        self.num_vars = n

    def set_linear(self, linear_vars_list):
        self.linear = linear_vars_list
    
    def set_quadratic(self, quadratic_vars_list):
        self.quadratic = quadratic_vars_list

    def set_cubic(self, cubic_vars_list):
        self.cubic = cubic_vars_list

    def reduce_cubic(self):
        """
        Reduces the cubic terms of the program by replacing common pairs using the substitution method.
        """
        # get the cubic terms and identify which must be replaced:
        cubic = self.cubic
        num_vars = self.get_num_vars()
        true_cubics = []
        remaining = []
        for tup, coeff in cubic:
            if len(set(tup)) == len(tup):  # check if all values in tuple are unique
                true_cubics.append((tup, coeff))  
            else:
                remaining.append((tup, coeff))
        self.cubic = remaining  

        # iteratively identify the most common pairs in the cubic terms and reduce all of the terms
        # containing these pairs until no cubic terms remain:
        while true_cubics != []:
            def find_most_common_pair(lst):
                pairs = []
                for tup, coeff in lst:
                    pairs.extend([(tup[i], tup[j]) for i in range(len(tup)) for j in range(i+1, len(tup))])
                most_common_pair = Counter(pairs).most_common(1)[0][0]
                sublist = [item for item in lst if most_common_pair[0] in item[0] and most_common_pair[1] in item[0]]
                for item in sublist:
                    lst.remove(item)
                return most_common_pair, sublist
            
            most_common_pair, mcp_sublist = find_most_common_pair(true_cubics)

            for i in range(len(mcp_sublist)):
                ((a, b, c), coeff) = mcp_sublist[i]
                if (a, b) == most_common_pair:
                    x1 = a
                    x2 = b
                    x3 = c
                if (a, c) == most_common_pair:
                    x1 = a
                    x2 = c
                    x3 = b
                if (b, c) == most_common_pair:
                    x1 = b
                    x2 = c
                    x3 = a
                # create the new terms:
                z = num_vars
                self.quadratic.append(((z, x3), coeff))
                self.quadratic.append(((x1, x2), 2*coeff))
                self.quadratic.append(((x1, z), -4*coeff))
                self.quadratic.append(((x2, z), -4*coeff))
                self.linear.append(((z, ), 6*coeff))

            # increment the number of qubits:
            num_vars += 1

        self.num_vars = num_vars


class IsingHamiltonian():
    """
    Represents an Ising Hamiltonian.

    Attributes:
        qubit_op (SparsePauliOp): The qubit operator representing the Ising Hamiltonian.
        offset (int): The constant value offset of the Ising Hamiltonian.
        random_generator (np.random): The random number generator.
        aer_simulator (AerSimulator): The Aer simulator.
        number_of_shots (int): The number of shots for the simulator.
        counts (dict): The counts of measurement outcomes.

    Methods:
        to_cubic_ising(cubic_instance): Converts a CubicProgram instance to a cubic Ising Hamiltonian.
    """

    def __init__(self, cubic_instance: CubicProgram) -> None:
        """
        Initializes an IsingHamiltonian object.

        Args:
            cubic_instance (CubicProgram): A CubicProgram instance.
        """
        self.to_cubic_ising(cubic_instance)
        self.ansatz = None
        self.energy = None
        self.optimal_value = None
        self.bitstring = []
        self.random_generator = np.random.default_rng()
        self.number_of_shots = 10*pow(2, self.qubit_op.num_qubits)
        self.aer_simulator = AerSimulator(method="statevector", shots=self.number_of_shots)
        self.counts = {}
    
    def to_cubic_ising(self, cubic_prog: CubicProgram) -> SparsePauliOp:
        """
        Converts a CubicProgram instance to a cubic Ising Hamiltonian.

        Args:
            cubic_prog (CubicProgram): A CubicProgram instance.

        Returns:
            SparsePauliOp: The qubit operator representing the cubic Ising Hamiltonian.
        """
        num_vars = cubic_prog.get_num_vars()
        zero = np.zeros(num_vars, dtype=bool)
        cubic_list = []

        # just to keep track of Rz, Rzz and Rzzz:
        r_z = []
        r_zz = []
        r_zzz = []

        offset = 0
        sense = 1

        # convert linear parts of the objective function into Hamiltonian.
        for idx, coef in cubic_prog.linear:
            z_p = zero.copy()
            weight = coef*sense / 2
            z_p[idx] = True
            cubic_list.append(SparsePauliOp(Pauli((z_p, zero)), -weight))
            r_z.append(idx[0])
            offset += weight

        # convert quadratic parts of the objective function into Hamiltonian.
        for (i, j), coeff in cubic_prog.quadratic:
            weight = coeff*sense / 4

            if (i == j):
                offset += weight
            else:
                z_p = zero.copy()
                z_p[i] = True
                z_p[j] = True
                cubic_list.append(SparsePauliOp(Pauli((z_p, zero)), weight))
                r_zz.append((i,j))
            
            z_p = zero.copy()
            z_p[i] = True
            cubic_list.append(SparsePauliOp(Pauli((z_p, zero)), -weight))
            r_z.append(i)

            z_p = zero.copy()
            z_p[j] = True
            cubic_list.append(SparsePauliOp(Pauli((z_p, zero)), -weight))
            r_z.append(j)

            offset += weight

        # convert cubic parts of the objective function into Hamiltonian.
        for (i, j, k), coeff in cubic_prog.cubic:
            weight = coeff*sense / 8

            # consider cubic terms:
            if i == j and j == k:
                # a*a*a = Ia = a
                z_p[i] = True
                cubic_list.append(SparsePauliOp(Pauli((z_p, zero)), -weight))
                r_z.append(i)
            else:
                if i == j and j != k:
                    # Then we have something like a^2c = Ic = c (quadratic case where i != j):
                    z_p = zero.copy()
                    z_p[k] = True
                    cubic_list.append(SparsePauliOp(Pauli((z_p, zero)), -weight))
                    r_z.append(k)
                else:
                    if i != j and i == k:
                        # a^2b = Ib = b:
                        z_p = zero.copy()
                        z_p[j] = True
                        cubic_list.append(SparsePauliOp(Pauli((z_p, zero)), -weight))
                        r_z.append(j)
                    else:
                        if i != j and j == k:
                            # aI = a
                            z_p = zero.copy()
                            z_p[i] = True
                            cubic_list.append(SparsePauliOp(Pauli((z_p, zero)), -weight))
                            r_z.append(i)
                        else:
                            # i != j != k: abc
                            z_p = zero.copy()
                            z_p[i] = True
                            z_p[j] = True
                            z_p[k] = True
                            cubic_list.append(SparsePauliOp(Pauli((z_p, zero)), -weight))
                            r_zzz.append((i,j,k))


            # consider quadratic terms:
            if (i == j):
                offset += weight
            else:
                z_p = zero.copy()
                z_p[i] = True
                z_p[j] = True
                cubic_list.append(SparsePauliOp(Pauli((z_p, zero)), weight))
                r_zz.append((i,j))
            
            if (i == k):
                offset += weight
            else:
                z_p = zero.copy()
                z_p[i] = True
                z_p[k] = True
                cubic_list.append(SparsePauliOp(Pauli((z_p, zero)), weight))
                r_zz.append((i,k))
            
            if (j == k):
                offset += weight
            else:
                z_p = zero.copy()
                z_p[j] = True
                z_p[k] = True
                cubic_list.append(SparsePauliOp(Pauli((z_p, zero)), weight))
                r_zz.append((j,k))
            
            # consider linear terms
            z_p = zero.copy()
            z_p[i] = True
            cubic_list.append(SparsePauliOp(Pauli((z_p, zero)), -weight))
            r_z.append(i)

            z_p = zero.copy()
            z_p[j] = True
            cubic_list.append(SparsePauliOp(Pauli((z_p, zero)), -weight))
            r_z.append(j)

            z_p = zero.copy()
            z_p[k] = True
            cubic_list.append(SparsePauliOp(Pauli((z_p, zero)), -weight))
            r_z.append(k)

            # consider constant term:
            offset += weight

        # combine all Pauli terms
        if cubic_list:
            # Remove paulis whose coefficients are zeros.
            qubit_op = sum(cubic_list).simplify(atol=0)
        else:
            # If there is no variable, we set num_nodes=1 so that qubit_op should be an operator.
            # If num_nodes=0, I^0 = 1 (int).
            num_vars = max(1, num_vars)
            qubit_op = SparsePauliOp("I" * num_vars, 0)

        connections = {"single": r_z, "double": r_zz, "triple": r_zzz}
        self.qubit_op = qubit_op
        self.offset = offset
        self.connections = connections
        return qubit_op    

    def show_ansatz(self):
        if self.ansatz is not None:
            self.ansatz.decompose(reps=3).draw(output="mpl", style="iqp")
        else:
            raise LookupError(
                "Ansatz has not been defined yet."
            )

    def sample_most_likely(self, state_vector, result):
        if isinstance(state_vector, QuasiDistribution):
            # THIS IS QAOA?
            probabilities = state_vector.binary_probabilities()
            bitstring = result.best_measurement['bitstring']
            bitstring = [1 if x == '1' else 0 for x in bitstring]
            return bitstring, result.best_measurement['probability']
            # binary_string = max(probabilities.items(), key=lambda kv: kv[1])[0]
            # x = np.asarray([int(y) for y in reversed(list(binary_string))])
            # return x, max(probabilities.items(), key=lambda kv: kv[1])[1]

        elif isinstance(state_vector, Statevector):
            # THIS IS VQE?
            probabilities = state_vector.probabilities()
            n = state_vector.num_qubits
            k = np.argmax(np.abs(probabilities))
            x = np.zeros(n)
            for i in range(n):
                x[i] = k % 2
                k >>= 1
            return x, np.max(np.abs(probabilities))
        else:
            raise ValueError(
                "state vector should be QuasiDistribution or Statevector "
                f"But it is {type(state_vector)}."
            )
    
    def run_VQE(self, correct_bitstrings, depth=2):
        """Runs the VQE algorithm on the graph encoded in the Hamiltonian.

        Args:
            correct_bitstrings (Array[string]): The solution bitstrings.
            depth (int, optional): Depth of the ansatz. Defaults to 2.

        Returns:
            Bitstring: the bitstring which has the highest probability to be the solution.
            Highest probability: Probability of the above bitstring
            Correct probability: Probability of the correct bitstring
        """
        # Set up:
        estimator = Estimator() 
        optimizer = COBYLA()#SLSQP(maxiter=500) 
        
        ansatz = EfficientSU2(self.qubit_op.num_qubits, su2_gates="ry", reps=depth, entanglement="linear")
        vqe = VQE(estimator, ansatz, optimizer)
        # Compute: 
        result = vqe.compute_minimum_eigenvalue(operator=self.qubit_op) 
        # Find the bitstring:
        optimal_parameters = result.optimal_parameters
        optimal_circuit = ansatz.assign_parameters(optimal_parameters)
        optimal_circuit.measure_active()
        sim_result = self.aer_simulator.run(transpile(optimal_circuit, self.aer_simulator)).result()

        self.counts = sim_result.get_counts(transpile(optimal_circuit, self.aer_simulator))
        bitstring, max_prob, correct_prob = self.bitstring_analysis(correct_bitstrings)
        self.ansatz = ansatz
        self.energy = result.eigenvalue
        self.optimal_value = result.eigenvalue + self.offset
        self.bitstring = bitstring
        self.max_prob = max_prob
        self.correct_prob = correct_prob
        return bitstring, max_prob, correct_prob
        
    
    def run_QAOA(self, correct_bitstrings, qaoa_reps=2):
        """
        Runs the Quantum Approximate Optimization Algorithm (QAOA) to find the minimum eigenvalue of the Hamiltonian.

        Args:
            correct_bitstrings (list): A list of correct bitstrings that the algorithm should aim to maximize the probability of.
            qaoa_reps (int): The number of repetitions (ansatz depth) for the QAOA algorithm. Default is 2.

        Returns:
            tuple: A tuple containing the bitstring with the highest probability, the maximum probability, and the cumulative probability of the correct bitstrings.

        """

        sampler = AerSampler()
        optimizer = COBYLA(maxiter=2000)
        #optimizer = SLSQP()
        qaoa = QAOA(sampler, optimizer, reps=qaoa_reps) # reps denote the ansatz depth!!
        result = qaoa.compute_minimum_eigenvalue(self.qubit_op)
        parameter_values_optimized = result.optimal_point

        self.optimised_parameters = parameter_values_optimized
        self.ansatz = qaoa.ansatz
        self.energy = result.best_measurement['value'] + self.offset

        optimal_parameters = result.optimal_point
        optimal_circuit = qaoa.ansatz.assign_parameters(optimal_parameters)
        sim_res = self.aer_simulator.run(transpile(optimal_circuit, self.aer_simulator)).result()
        self.counts = sim_res.get_counts(transpile(optimal_circuit, self.aer_simulator))
        #print(self.counts)

        bitstring, max_prob, correct_prob = self.bitstring_analysis(correct_bitstrings)
        if correct_prob == 0:
            correct_prob = 1/(2*pow(2, self.qubit_op.num_qubits))
        self.bitstring = bitstring
        self.max_prob = max_prob
        self.correct_prob = correct_prob
        return bitstring, max_prob, correct_prob
        
    def bitstring_analysis(self, correct_bitstrings):
        """
        Analyzes the bitstrings obtained from a quantum circuit measurement and calculates the probabilities of each bitstring.
        
        Args:
            correct_bitstrings (list): A list of correct bitstrings that are expected to be observed.
        
        Returns:
            tuple: A tuple containing the following elements:
                - A list of integers representing the bitstring with the highest probability, where '1' is represented as 1 and '0' is represented as 0.
                - The highest probability observed.
                - The cumulative probability of observing any of the correct bitstrings.
        """
        probabilities = {bitstring: (count/self.number_of_shots) \
                            for bitstring, count in self.counts.items()}
        correct_prob = 0
        for correct_bitstring in correct_bitstrings:
            if correct_bitstring in probabilities.keys():
                correct_prob += probabilities[correct_bitstring]

        bitstring_highest_prob = max(
            probabilities, key=lambda key: probabilities[key])
        highest_probability = probabilities[bitstring_highest_prob]
        return [1 if x == '1' else 0 for x in bitstring_highest_prob], highest_probability, correct_prob
    
    def simulate_QAOA_circuit(self, qaoa_ansatz=None, optimised_parameters=None, optimal_circuit=None):
        """
        Simulates a QAOA circuit.

        Args:
            qaoa_ansatz (Optional[QuantumCircuit]): The QAOA ansatz circuit to be used. If not provided, the default ansatz will be used.
            optimised_parameters (Optional[List[float]]): The optimized parameters for the QAOA circuit. If not provided, the default optimized parameters will be used.
            optimal_circuit (Optional[QuantumCircuit]): The optimal circuit to be used. If not provided, the circuit will be constructed using the optimized parameters.

        Returns:
            Result: The result of running the QAOA circuit.

        """
        if optimised_parameters is None:
            optimised_parameters = self.optimised_parameters
        if qaoa_ansatz is None:
            qaoa_ansatz = self.ansatz
        if optimal_circuit is None:
            # Find correct string chance by running the circuit with optimized params:
            parameter_bindings = dict(zip(qaoa_ansatz.parameters, optimised_parameters))
            qaoa_with_parameters = self.ansatz.assign_parameters(parameter_bindings)
            qaoa_with_parameters_decomposed = qaoa_with_parameters.decompose(reps=3)

            result = self.aer_simulator.run(qaoa_with_parameters_decomposed).result()
            return result
        else:
            qaoa_with_parameters_decomposed = optimal_circuit.decompose(reps=3)
            result = self.aer_simulator.run(qaoa_with_parameters_decomposed).result()
            return result
    
    def run_QAOA_custom(self, correct_bitstrings, qaoa_reps=2):
        """Runs the QAOA algorithm on the Hamiltonin without using the QAOA.compute_minimum_eigenvalue().

        Args:
            correct_bitstrings (Array[string]): The solution bitstrings.
            qaoa_reps (int, optional): Depth of the ansatz. Defaults to 2.

        Returns:
            Bitstring: the bitstring which has the highest probability to be the solution.
            Highest probability: Probability of the above bitstring
            Correct probability: Probability of the correct bitstring
        """
        qaoa_ansatz = QAOAAnsatz(cost_operator=self.qubit_op, reps=qaoa_reps, name='qaoa')
        qaoa_ansatz.measure_active()

        observable = OperatorStateFn(PauliSumOp(self.qubit_op), is_measurement=True)

        def energy_evaluation(parameter_values: List[float]):
            """parameter_values is expected to be of
            the form [beta_0, beta_1, ..., gamma_0, gamma_1, ...]"""
            parameter_bindings = dict(zip(qaoa_ansatz.parameters, parameter_values))
            qaoa_with_parameters = qaoa_ansatz.assign_parameters(parameter_bindings)
            qaoa_with_parameters_decomposed = qaoa_with_parameters.decompose(reps=3)
            result = self.aer_simulator.run(qaoa_with_parameters_decomposed).result()
            counts = result.get_counts()
            sqrt_probabilities = {bitstring: np.sqrt(count/self.number_of_shots) \
            for bitstring, count in counts.items()}
            expectation = np.real(observable.eval(sqrt_probabilities)) + self.offset
            return expectation
        
        # Find optimized parameters:
        betas_initial_guess = np.pi*self.random_generator.random(qaoa_reps)
        gammas_initial_guess = 2*np.pi*self.random_generator.random(qaoa_reps)
        parameter_values_initial_guess = [*betas_initial_guess, *gammas_initial_guess]

        cobyla_optimizer = COBYLA(maxiter=5000)

        result_optimization = cobyla_optimizer.minimize(
            fun=energy_evaluation, x0=parameter_values_initial_guess)

        parameter_values_optimized = result_optimization.x
        energy_optimized = result_optimization.fun
        # Number of evalutions of energy_evaluation
        number_function_evaluations = result_optimization.nfev

        self.optimised_parameters = parameter_values_optimized
        self.ansatz = qaoa_ansatz
        self.energy = energy_optimized
        self.function_evals = number_function_evaluations

        result = self.simulate_QAOA_circuit()
        self.counts = result.get_counts()

        bitstring, max_prob, correct_prob = self.bitstring_analysis(correct_bitstrings)
        if correct_prob == 0:
            correct_prob = 1/(2*pow(2, self.qubit_op.num_qubits))
        self.bitstring = bitstring
        self.max_prob = max_prob
        self.correct_prob = correct_prob
        return bitstring, max_prob, correct_prob
    
    def run_numpy_eigensolver(self):
        """
        Runs the NumPy eigensolver to compute the minimum eigenvalue and corresponding eigenvector.

        Returns:
            tuple: A tuple containing the computed bitstring and the maximum probability.
        """
        ee = NumPyMinimumEigensolver()
        result = ee.compute_minimum_eigenvalue(self.qubit_op)
        bitstring = self.sample_most_likely(result.eigenstate, result)
        
        self.energy = result.eigenvalue.real
        self.optimal_value = result.eigenvalue.real + self.offset
        self.bitstring = bitstring[0]
        self.max_prob = 1.0
        return bitstring[0], 1.0



