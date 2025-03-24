import numpy as np
from itertools import product

def generate_binary_variables(n):
    # Generate all combinations of binary variables of size n*m
    bin_vars = list(product([0, 1], repeat = 2*n))
    return bin_vars

def compute_arithmetic_function(x, y, mat, n, h):
    """
    Computes the arithmetic function (CUBO) based on the given parameters (pattern and target graphs).

    Parameters:
    x (list): A 2D list representing the mapping between vertices of two graphs.
    adj_matrix1 (list): A 2D list representing the adjacency matrix of the first graph.
    adj_matrix2 (list): A 2D list representing the adjacency matrix of the second graph.
    n (int): The number of vertices in the first graph.
    m (int): The number of vertices in the second graph.

    Returns:
    float: The computed arithmetic function value.

    """

    # OBJECTIVE FUNCTION (MINIMIZATION):
    # sum1 = 0
    # for i in range(n):
    #     for j in range(n):
    #         sum1 += -1*y[i] * x[j] * adj_matrix[i][j]
    
    # # CONSTRAINT 1 (CHECKS THAT EACH VERTEX IN G1 IS MAPPED TO AT MOST ONE VERTEX IN G2):
    # sum2 = 0
    # for i in range(n):
    #     for j in range(n):
    #         sum2 += x[i] * adj_matrix[i][j]

    # sum3 = 0
    # # CONSTRAINT 2 (CHECKS THAT NUMBER OF OUTGOING EDGES OF EACH VERTEX U IS EQUAL TO OUTGOING EDGES OF EACH VERTEX I):
    # for i in range(n):
    #     inner = 0
    #     for j in range(n):
    #         inner += adj_matrix[i][j]
        
    #     sum3 += n * y[i] * (1 - inner + z[i] * inner)

    # sum4 = 0
    # # CONSTRAINT 3 (CHECKS THAT NUMBER OF INCOMING EDGES OF EACH VERTEX U IS EQUAL TO OUTGOING EDGES OF EACH VERTEX I):
    # for i in range(n):
    #     inner = 0
    #     for j in range(n):
    #         inner += adj_matrix[i][j]
    #     sum4 += (2 * n + 7) * z[i] * (2 - inner)
    
    # sum5 = 0
    # for i in range(n):
    #     sum5 += -1 * (n + 3) * z[i]
    P = n
    sum = 0
    for i in range(n):
        inner = 0
        for j in range(i+1, n):
            inner += mat[i][j]
        sum += x[i] * (P * inner + 0.5)
        sum += y[i] * (P * (inner - 1)**2 + 0.5)

        for j in range(i+1, n):
            sum += y[i] * x[j] * -1 * mat[i][j] * (1 - h[i]) * h[j]

    return sum

def simulate_arithmetic_function(x, y, z, out, n=5):
    term1 = -1 * x + x * out
    
    term2 = -1 * y
    
    term3 = n * y * (1 - out + z * out)
    
    term4 = (2 * n + 7) * z * (2 - out)
    
    term5 = -1 * (n + 3) * z
    
    return term1 + term2 + term3 + term4 + term5

def main(adj_matrix, h):
    """
    Finds the best combination of binary variables that minimizes the result of the CUBO cost function.

    Args:
        adj_matrix1 (list): The adjacency matrix 1.
        adj_matrix2 (list): The adjacency matrix 2.

    """
    n = len(adj_matrix1)
    binary_variables = generate_binary_variables(n)

    min = 1000
    best = None
    for variables in binary_variables:
        xs = variables[:n]
        ys = variables[n:]
        result = compute_arithmetic_function(xs, ys, adj_matrix, n, h)
        if result < min:
            min = result
            best = variables
    print(f"Final result:\n Minimum: {min}\n Best combination:\n {best}")
    
def compute_table():
    # compute the penalties for each combination of binary variables
    xs = [0, 1]
    ys = [0, 1]
    zs = [0, 1]
    outgoing = [0, 1, 2, 3]
    for out in outgoing:
        for x in xs:
            for y in ys:
                for z in zs:
                    print(f"out: {out}, x: {x}, y: {y}, z: {z}")
                    print(simulate_arithmetic_function(out, x, y, z))
    

def compute(adj_matrix1, h, solution):
    """
    Computes the arithmetic function using the given adjacency matrices.

    Parameters:
    adj_matrix1 (list): The first adjacency matrix.
    adj_matrix2 (list): The second adjacency matrix.
    """
    
    n = len(adj_matrix1)
    xs = solution[:n]
    ys = solution[n:]
    result = compute_arithmetic_function(xs, ys, adj_matrix1, n, h)

    print(result)
    

# Define your adjacency matrices here
# adj_matrix1 = [[0, 1, 0, 0],
#                [0, 0, 1, 0],
#                [0, 0, 0, 1],
#                [1, 0, 0, 0]]

adj_matrix1 = np.array([[0., 0., 0., 0., 1.],
       [0., 0., 0., 0., 1.],
       [0., 0., 0., 1., 0.],
       [0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0.]])
h = [0, 0, 0, 0, 1]

#compute(adj_matrix1, adj_matrix2)
#main(adj_matrix1, adj_matrix2)
main(adj_matrix1, h)
#compute(adj_matrix1, h, [0, 0, 0, 0, 1, 1, 1, 0, 0, 0])
# compute_table()