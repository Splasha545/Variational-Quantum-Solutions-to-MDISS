import numpy as np
from itertools import product

def generate_binary_variables(n):
    # Generate all combinations of binary variables of size n*m
    bin_vars = list(product([0, 1], repeat = 3*n))
    return bin_vars

def compute_arithmetic_function(x, y, z, h, a, n):
    
    term1 = 0
    for i in range(n):
        term1 += 2 * x[i] * (1 - h[i]) - x[i]
        for j in range(n):
            term1 += x[i] * a[i][j]
            
    term2 = 0
    for i in range(n):
        term2 += 2 * y[i] * h[i] - y[i]
    
    term3 = 0
    for i in range(n):
        inner = 0
        for j in range(n):
            inner += a[i][j]
        term3 += n * y[i] * (1 - inner + z[i] * inner)
    
    P = 2 * n + 7
    term4 = 0
    for i in range(n):
        inner = 0
        for j in range(n):
            inner += a[i][j]
        term4 += P * z[i] * (2 - inner)
    
    term5 = 0
    for i in range(n):
        term5 += -1 * (n + 3) * z[i]
    
    return term1 + term2 + term3 + term4 + term5

def expanded_compute_arithmetic_function(x, y, z, h, a, n):
    
    term1 = 0
    for i in range(n):
        term1 += x[i] * (2 * (1 - h[i]) - 1)
        term1 += y[i] * (2 * h[i] + n - 1)
        term1 += z[i] * (3 * n + 11)
        for j in range(n):
            term1 += x[i] * a[i][j]
            term1 += -1 * y[i] * n * a[i][j]
            term1 += y[i] * z[i] * n * a[i][j]
            term1 += -1 * (2*n + 7) * z[i] * a[i][j]
    
    return term1

def compute(adj_matrix1, adj_matrix2):
    """
    Computes the arithmetic function using the given adjacency matrices.

    Parameters:
    adj_matrix1 (list): The first adjacency matrix.
    adj_matrix2 (list): The second adjacency matrix.
    """
    
    n = len(adj_matrix1)
    m = len(adj_matrix2)

    print(compute_arithmetic_function([0, 1, 0, 0, 0], [1, 0, 0, 1, 0], [0, 0, 1, 0, 0], h, adj_matrix1, n))
                     

def main(adj_matrix, h):
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
        zs = variables[2*n:3*n]
        result = compute_arithmetic_function(xs, ys, zs, h, adj_matrix, n)
        if result < min:
            min = result
            best = variables

    print(f"Final result:\n Minimum: {min}\n Best combination:\n {best}")
    return min, best
    
# adjacency_matrix = np.array([[0., 1., 0., 1., 0.],
#        [0., 0., 0., 1., 0.],
#        [0., 0., 0., 1., 0.],
#        [0., 0., 0., 0., 0.],
#        [0., 0., 0., 0., 0.]])
h = [0, 0, 0, 1, 0]

adjacency_matrix = np.array([[0., 1., 0., 0., 0.],
       [0., 0., 0., 0., 0.],
       [1., 1., 0., 0., 1.],
       [0., 1., 0., 0., 0.],
       [1., 1., 0., 0., 0.]])
h = [0, 1, 0, 0, 0]

# min, best = main(adjacency_matrix, h)
compute(adjacency_matrix, adjacency_matrix)