import numpy as np
from itertools import product

def generate_binary_variables(n):
    # Generate all combinations of binary variables of size n*m
    bin_vars = list(product([0, 1], repeat = 2*n))
    return bin_vars

def compute_arithmetic_function(x, y, a, n):
    
    term1 = 0
    for i in range(n):
        term1 += 0.5 * x[i] * y[i]
    
    term2 = 0
    for i in range(n):
        inner = 0
        for j in range(n):
            inner += a[i][j]
        term2 += inner * (1 - x[i])
    
    term3 = 0
    for i in range(n):
        inner = 0
        for j in range(n):
            inner += a[j][i]
        term3 += inner * (1 - y[i])
    
    return term1 + term2 + term3


def main(adj_matrix):
    """
    Finds the best combination of binary variables that minimizes the result of the CUBO cost function.

    Args:
        adj_matrix (list): The adjacency matrix of G.

    """
    n = len(adj_matrix)
    binary_variables = generate_binary_variables(n)

    min = 1000
    best = None
    for variables in binary_variables:
        xs = variables[:n]
        ys = variables[n:2*n]
        result = compute_arithmetic_function(xs, ys, adj_matrix, n)
        if result < min:
            min = result
            best = variables

    print(f"Final result:\n Minimum: {min}\n Best combination:\n {best}")
    return min, best
    
adjacency_matrix = [
                    [0, 1, 1, 1, 1, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 1, 0, 1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0]
                    ]

adjacency_matrix = np.array([[0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 1., 0., 0.],
       [1., 0., 0., 0., 1., 0., 0.],
       [1., 1., 0., 0., 1., 1., 0.],
       [0., 0., 0., 0., 0., 0., 0.],
       [1., 1., 0., 0., 1., 0., 0.],
       [1., 1., 0., 0., 1., 0., 0.]])

min, best = main(adjacency_matrix)