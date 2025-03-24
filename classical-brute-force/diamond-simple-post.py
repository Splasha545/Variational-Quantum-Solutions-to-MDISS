import numpy as np
from itertools import product

def generate_binary_variables(n):
    # Generate all combinations of binary variables of size n*m
    bin_vars = list(product([0, 1], repeat = 2*n))
    return bin_vars

def compute_arithmetic_function(x, y, a, n):
    
    term1 = 0
    for i in range(n):
        for j in range(i+1, n):
            if i != j:
                # print(f"The 1.5 term: {1.5 * x[i] * x[j]}")
                # if 1.5 * x[i] * x[j] > 0:
                #     print(f"i: {i}, j: {j}, x[i]: {x[i]}, x[j]: {x[j]}")
                term1 += 1.5 * x[i] * x[j]
    
    term2 = 0
    for i in range(n):
        for j in range(n):
            if i != j:
                inner = 0
                for k in range(n):
                    inner += a[i][k] * a[k][j]
                # if inner >= 2:
                #     print(f"i: {i}, j: {j}, k: {k}, inner: {inner}")
                #     print(f"x[i]: {x[i]}, x[j]: {x[j]}")
                # print(f"The large term: {inner * (1 - x[i]*x[j])}")
                # print(f"The right term: {inner * x[i] * (1 - x[j])}")
                # print(f"The left term: {inner * (1 - x[i]) * x[j]}")
                term2 += inner * (1 - x[i]*x[j])
                term2 += inner * x[i] * (1 - x[j])
                term2 += inner * (1 - x[i]) * x[j]
                # term2 += (inner - 1) * (1 - x[i]) * (1 - x[j])
    
    # term3 = 0
    # for i in range(n):
    #     for j in range(n):
    #         term3 += 2.5 * y[i]* y[j]
    #         term3 += 5 * y[i] * (1 - y[j])
    #         term3 += 5 * (1 - y[i]) * y[j]
            
    # term4 = 0
    # for i in range(n):
    #     for j in range(n):
    #         inner = 0
    #         for k in range(n):
    #             inner += a[i][k] * a[k][j]
    #         term4 += inner * (1 - y[i]*y[j])
    
    return term1 + term2 #+ term3 + term4


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
        ys = variables[n:]
        # # make matrices:
        # xs = np.array(xs).reshape(n, n)
        # ys = np.array(ys).reshape(n, n)
        result = compute_arithmetic_function(xs, ys, adj_matrix, n)
        if result <= min:
            min = result
            best = variables

    print(f"Final result:\n Minimum: {min}\n Best combination:\n {best}")
    return min, best

def compute(adj_matrix, bitstring):
    n = len(adj_matrix)
    xs = bitstring[:n]
    ys = bitstring[n:]
    result = compute_arithmetic_function(xs, ys, adj_matrix, n)
    return result
    
# adjacency_matrix = [[0, 1, 0, 0, 0, 0],
#  [0, 0, 0, 0, 0, 0],
#  [1, 1, 0, 1, 1, 0],
#  [0, 1, 0, 0, 0, 0],
#  [0, 1, 0, 0, 0, 0],
#  [0, 1, 0, 0, 1, 0]]
adjacency_matrix = [[0, 1, 1, 1, 1],
 [0, 0, 0, 1, 0],
 [0, 0, 0, 1, 0],
 [0, 0, 0, 0, 0],
 [0, 0, 0, 1, 0]]

# adjacency_matrix = [
#     [0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0],
#     [1, 1, 1, 0, 0, 0, 0],
#     [1, 0, 0, 0, 0, 0, 0],
#     [1, 0, 0, 0, 0, 0, 0],
#     [1, 1, 1, 0, 1, 1, 0]
#  ]
adjacency_matrix = [
    [0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [1, 1, 0, 1, 1, 1, 0],
    [0, 1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0],
    [1, 1, 0, 1, 1, 1, 0]]


print(compute(adjacency_matrix, [0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
# min, best = main(adjacency_matrix)

# best_xs = best[:len(best)//2]
# best_ys = best[len(best)//2:]
# # best_xs = np.array(best_xs).reshape(len(adjacency_matrix), len(adjacency_matrix))
# # best_ys = np.array(best_ys).reshape(len(adjacency_matrix), len(adjacency_matrix))

# n = len(adjacency_matrix)

# xs_sum = 0
# ys_sum = 0
# for i in range(n):
#     for j in range(n):
#         xs_sum += best_xs[i][j]
#         ys_sum += best_ys[i][j]
# print(f"Final solution: {xs_sum-ys_sum}")