# Description: This file contains the class tests which is used to store the test cases (graphs) of different
# sizes (from 4 - 21 qubits) for solving with VQAs.

class tests:
    """
    A class that represents a collection of tests for a quantum algorithm.
    
    Attributes:
    - problems: A dictionary that stores the problems for different qubit numbers.
    - qubit_nums: A list that stores the qubit numbers for the problems.
    """
    
    problems: {}
    qubit_nums: []
    def __init__(self):
        """
        Initialize the MyClass object.

        This method initializes the object by setting up the `qubit_nums` and `problems` attributes.

        `qubit_nums` is a list of integers representing the number of qubits for each problem.
        `problems` is a dictionary that stores the problem data for each problem instance.

        Each problem is represented as a dictionary with the following keys:
        - `n`: an integer representing the number of vertices in the pattern graph
        - `elistG`: a list of tuples representing the edges of the pattern graph G
        - `m`: an integer representing the number of vertices in the target graph
        - `elistH`: a list of tuples representing the edges of the target graph H
        - `correct_bitstrings`: a list of correct bitstrings for the problem

        The `problems` dictionary is structured as follows:
        {
            <number of qubits>: [
                {
                    "n": <number of vertices in G>,
                    "elistG": <list of edges in G>,
                    "m": <number of vertices in H>,
                    "elistH": <list of edges in H>,
                    "correct_bitstrings": <list of correct bitstrings which find correct assignments from G to H>
                },
                ...
            ],
            ...
        }
        """
        self.qubit_nums = [4, 6, 8, 10, 12, 14, 16, 18, 21]
        self.problems = {
            4: [
                {"n": 1, "elistG": [], "m": 4, "elistH": [(0, 1), (1, 2), (2, 3)], "correct_bitstrings": ['1010', '0101', '1001']},
                {"n": 2, "elistG": [(0, 1)], "m": 2, "elistH": [(0, 1)], "correct_bitstrings": ['1001']},
                {"n": 2, "elistG": [(0, 1)], "m": 2, "elistH": [(0, 1)], "correct_bitstrings": ['1001']}
            ],
            6: [
                {"n": 2, "elistG": [(0, 1)], "m": 3, "elistH": [(0, 1), (1, 2), (2, 0)], "correct_bitstrings": ['010001', '100010', '001100']},
                {"n": 2, "elistG": [(0, 1)], "m": 3, "elistH": [(0, 1), (1, 2)], "correct_bitstrings": ['010001', '100010']},
                {"n": 2, "elistG": [(0, 1)], "m": 3, "elistH": [(0, 1), (0, 2)], "correct_bitstrings": ['100001', '010001']}
            ],
            8: [
                {"n": 2, "elistG": [(0, 1)], "m": 4, "elistH": [(0, 1), (2, 3)], "correct_bitstrings": ['10100101']},
                {"n": 2, "elistG": [(0, 1)], "m": 4, "elistH": [(0, 1), (1, 2), (2, 3)], "correct_bitstrings": ['00100001', '01000010', '10000100']},
                {"n": 2, "elistG": [(0, 1)], "m": 4, "elistH": [(0, 1), (1, 2), (2, 3), (3, 0)], "correct_bitstrings": ['00100001', '01000010', '10000100', '00011000']}
            ],
            10: [
                {"n": 2, "elistG": [(0, 1)], "m": 5, "elistH": [(0, 1), (1, 2), (2, 3), (3, 4)], "correct_bitstrings": ['1001001001']},
                {"n": 2, "elistG": [(0, 1)], "m": 5, "elistH": [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)], "correct_bitstrings": ['0001000001', '0010000010', '0100000100', '1000001000', '0000110000']},
                {"n": 2, "elistG": [(0, 1)], "m": 5, "elistH": [(0, 1), (1, 2), (2, 0), (1, 3), (3, 4)], "correct_bitstrings": ['1000101100']}
            ],
            12: [
                {"n": 2, "elistG": [(0, 1)], "m": 6, "elistH": [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)], "correct_bitstrings": ['010010001001', '100100010010', '100010010001']},
                {"n": 2, "elistG": [(0, 1)], "m": 6, "elistH": [(0, 1), (1, 2), (0, 2), (2, 3), (3, 4), (3, 5), (4, 5)], "correct_bitstrings": ['100010010001', '010010001001', '100010001001', '100100010001', '100100010010']},
                {"n": 3, "elistG": [(0, 1), (1, 2), (0, 2)], "m": 4, "elistH": [(0, 1), (1, 2), (0, 2), (3, 1), (3, 2)], "correct_bitstrings": ['010000100001', '010000101000']}
            ],
            14: [
                {"n": 2, "elistG": [(0, 1)], "m": 7, "elistH": [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)], "correct_bitstrings": ['00100100001001', '01000100010001', '10000100100001', '01001000010010', '10001000100010', '10010000100100']},
                {"n": 2, "elistG": [(0, 1)], "m": 7, "elistH": [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 0)], "correct_bitstrings": ['00100100001001', '01000100010001', '01001000010010', '10001000100010', '10010000100100', '00010011000100', '00100011001000']},
                {"n": 2, "elistG": [(0, 1)], "m": 7, "elistH": [(0, 1), (1, 4), (2, 3), (3, 4), (5, 4), (6, 4), (5, 6)], "correct_bitstrings": ['10010100100101']}
            ],
            16: [
                {"n": 2, "elistG": [(0, 1)], "m": 8, "elistH": [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)], "correct_bitstrings": ['1001001001001001']},
                {"n": 2, "elistG": [(0, 1)], "m": 8, "elistH": [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 0)], "correct_bitstrings": ['0001001000001001', '0010001000010001', '0100001000100001', '0010010000010010', '0100010000100010', '1000010001000010', '0100100000100100', '1000100001000100', '0000100110000100', '1001000001001000', '0001000110001000', '0010000110010000']},
                {"n": 2, "elistG": [(0, 1)], "m": 8, "elistH": [(0, 1), (0, 2), (2, 1), (1, 3), (3, 6), (6, 7), (4, 3), (4, 5)], "correct_bitstrings": ['1010001001010001', '1010010001010001', '1010001001010100']}
            ],
            18: [
                {"n": 3, "elistG": [(0, 1), (0, 2), (1, 2)], "m": 6, "elistH": [(0, 1), (1, 2), (0, 2), (2, 5), (3, 5), (3, 4), (4, 5)], "correct_bitstrings": ['000100000010000001', '100000010000001000']},
                {"n": 3, "elistG": [(0, 1), (0, 2), (1, 2)], "m": 6, "elistH": [(0, 1), (0, 3), (1, 2), (1, 3), (1, 4), (2, 4), (4, 3), (5, 3), (5, 4)], "correct_bitstrings": ['001000000010000001', '010000000100000010', '001000010000000010', '001000010000100000']},
                {"n": 3, "elistG": [(0, 1), (0, 2), (1, 2)], "m": 6, "elistH": [(0, 1), (1, 2), (0, 3), (1, 3), (2, 4), (2, 5), (4, 3), (5, 4)], "correct_bitstrings": ['001000000010000001', '010000100000000100']}
            ],
            21: [
                {"n": 3, "elistG": [(0, 1), (0, 2), (1, 2)], "m": 7, "elistH": [(0, 1), (0, 2), (1, 2), (2, 3), (3, 6), (4, 5), (4, 6), (5, 6)], "correct_bitstrings": ['100010001000100010001']},
                {"n": 3, "elistG": [(0, 1), (0, 2), (1, 2)], "m": 7, "elistH": [(0, 1), (1, 2), (2, 3), (0, 4), (1, 4), (1, 5), (2, 5), (2, 6), (3, 6), (4, 5), (5, 6)], "correct_bitstrings": ['001000000000100000001', '010000000001000000010', '100000000010000000100', '010000000100000000010', '100000001000000000100']},
                {"n": 3, "elistG": [(0, 1), (0, 2), (1, 2)], "m": 7, "elistH": [(0, 1), (1, 2), (0, 2), (1, 3), (3, 6), (6, 5), (3, 5), (5, 4), (4, 2)], "correct_bitstrings": ['000010000000100000001', '010000010000000001000']}
            ]
        }
