# Variational Solutions for the Maximal Disjoint Isomorphic Subgraph Selection

This folder contains the necessary files for running Variational Quantum Algorithms, namely VQE and QAOA, on the Maximal Disjoint Isomorphic Subgraph Selection (MDISS) directed graph problem.

## Files

- `VQALib.py`: This file defines the necessary classes to solve a problem instance with Variational Quantum Eigensolver (VQE) and Quantum Approximate Optimization Algorithm (QAOA).

- `VQATests.py`: This file defines 27 problem instances to be solved using VQE and QAOA. There are 3 problems defined for the following sizes (qubit numbers): 4, 6, 8, 10, 12, 14, 16, 18, 21.

- `VQA-experiments.py`: This file contains all the functions necessary to run experiments using VQE and QAOA.

- `experiment-results.ipynb`: This file contains experiment results and their visualisations.

### Demo files
- `QAOA-demo-2.ipynb` and `VQE-demo-2.ipynb` present walkthroughs for using the functions defined in VQALib and VQA-experiments to solve some problems defined in `VQATests.py`

## Usage

To use the files in this folder, follow these steps:

1. Import the necessary libraries, which are given in `requirements.txt`

2. Import the required classes and functions from the respective files.

3. Use the imported classes and functions to solve specific quantum computing problems or run experiments. Some examples are depicted in `QAOA-demo-2.ipynb` and `VQE-demo-2.ipynb`. Additionally, the `VQA-experiments.py` has some examples that can be run. 
