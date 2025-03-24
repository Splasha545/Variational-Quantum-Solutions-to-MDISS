import numpy as np
from itertools import product

def generate_binary_variables(n, m, k):
    # Generate all combinations of binary variables of size n*m
    bin_vars_list = list(product([0, 1], repeat = k*n*m))
    # Reshape the list of variables into k separate 2D arrays for easy addressing:
    return bin_vars_list
    

def compute_arithmetic_function(x, adj_matrix1, adj_matrix2, z, w, n, m):
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
    objective = 0
    energy = 0

    # OBJECTIVE FUNCTION (MINIMIZATION):
    objective = 0
    for i in range(n):
        for u in range(m):
            objective -= (x[i][u])
    
    # CONSTRAINT 1 (CHECKS THAT THERE IS AT LEAST ONE CORRESPONDING OUTGOING EDGE SATISFYING THE MAPPING):
    P_1 = 0
    for u in range(m):
        for i in range(n):
            if (x[i][u] == 1):
                for j in range(n):
                    if adj_matrix1[i][j] != 1:
                        continue
                    out = 0
                    out2 = 0
                    # finding number of outgoing edges:
                    for v in range(m):
                        if (adj_matrix2[u][v] == 1):
                            out += x[j][v]
                            
                    # finding number of outgoing edges if more than 2 edges:
                    if z[i][u] == 1:
                        for v in range(m):
                            if (adj_matrix2[u][v] == 1):
                                out2 += x[j][v]
                        out2 -= 1
                    
                    # penalty for "at least 1":
                    P_1 += 2 * (1 - out + out2)

    # CONSTRAINT 2 (CHECKS THAT THERE IS AT LEAST ONE CORRESPONDING INCOMING EDGE SATISFYING THE MAPPING):
    P_2 = 0
    for u in range(m):
        for i in range(n):
            if (x[i][u] == 1):
                for j in range(n):
                    if adj_matrix1[j][i] != 1:
                        continue
                    out = 0
                    out2 = 0
                    # finding number of outgoing edges:
                    for v in range(m):
                        if (adj_matrix2[v][u] == 1):
                            out += x[j][v]
                            
                    # finding number of outgoing edges if more than 2 edges:
                    if w[i][u] == 1:
                        for v in range(m):
                            if (adj_matrix2[v][u] == 1):
                                out2 += x[j][v]
                        out2 -= 1
                    
                    # penalty for "at least 1":
                    P_2 += 2 * (1 - out + out2)

    # CONSTRAINT 3 (ENFORCES z[i] TO BE 1 IF THERE ARE MORE THAN 2 OUTGOING EDGES):
    P_3 = 0
    for i in range(n):
        for u in range(m):
            if z[i][u] != 1:
                continue
            for j in range(n):
                if adj_matrix1[i][j] != 1:
                    continue
                out = 0
                for v in range(m):
                    if adj_matrix2[u][v] == 1:
                        out += x[j][v]
                P_3 += 3 * (2 - out)
    
    # CONSTRAINT 4 (ENFORCES w[i] TO BE 1 IF THERE ARE MORE THAN 2 INCOMING EDGES):
    P_4 = 0
    for i in range(n):
        for u in range(m):
            if w[i][u] != 1:
                continue
            for j in range(n):
                if adj_matrix1[j][i] != 1:
                    continue
                out = 0
                for v in range(m):
                    if adj_matrix2[v][u] == 1:
                        out += x[j][v]
                P_4 += 3 * (2 - out)
    
    # CONSTRAINT 5 (ENSURES WE FIND ALL Z[I][U])
    P_5 = 0
    for i in range(n):
        for u in range(m):
            for j in range(n):
                
                P_5 -= 3 * z[i][u] * adj_matrix1[i][j]
    
    # CONSTRAINT 6 (ENSURES WE FIND ALL W[I][U])
    P_6 = 0
    for i in range(n):
        for u in range(m):
            for j in range(n):
                P_6 -= 3 * w[i][u] * adj_matrix1[j][i]
            
    energy += objective + P_1 + P_2 + P_3 + P_4 + P_5 + P_6

    return energy

def compute_arithmetic_function_2(x, a, b, z, w, n, m):
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
    objective = 0
    energy = 0
    
    P_1 = 5*m + m
    P_2 = 5*m*m + m
    P_3 = 5*m
    P_4 = 1 # P_3 - 0.5
    P_5 = 6*m

    # outgoing term 1:
    out_1 = 0
    out_3 = 0
    out_4 = 0
    out_5 = 0
    # for i in range(n):
    #     for u in range(m):
    #         for j in range(n):
    #             for v in range(m):
    #                 out_1 -= x[i][u] * a[i][j] * b[u][v] * x[j][v]
    #             # # outgoing term 3:
    #             # for v in range(m):
    #             #     out_3 += P_2 * w[u][j] * b[u][v] * x[j][v] * a[i][j]
    #             # # outgoing term 4:
    #             # out_4 -= P_3 * w[u][j] * a[i][j]
    #             # # outgoing term 5:
    #             # inner = 0
    #             # for k in range(n):
    #             #     for v in range(m):
    #             #         inner += b[u][v] * x[k][v]
    #             # out_5 += P_4 * (1 - inner) * w[u][j] * a[i][j]
                    
    # # outgoing term 2:
    # out_2 = 0
    # for i in range(n):
    #     for u in range(m):
    #         for j in range(n):
    #             out_2 += P_1 * x[i][u] * a[i][j] * w[u][j]
    
    # # # outgoing term 5:
    # out_5 = 0
    # for i in range(n):
    #     for u in range(m):
    #         for j in range(n):
    #             inner = 0
    #             for k in range(n):
    #                 for v in range(m):
    #                     inner += b[u][v] * x[k][v]
    #             out_5 += P_4 * (1 - inner) * w[u][j] * a[i][j]
    
    out_1 = 0
    out_2 = 0
    out_3 = 0
    out_4 = 0
    out_5 = 0
    for i in range(n):
        for u in range(m):
            #print(f"x[{i}][{u}]: {x[i][u]}")
            for j in range(n):
                k = 0
                for v in range(m):
                    k += b[u][v] * x[j][v]
                    out_1 -= x[i][u] * a[i][j] * b[u][v] * x[j][v] * P_5
                #print(f"Sum b[{u}][v] x[{j}][v]: {k}")
                #print(f"w[{u}][{j}]: {w[u][j]}")
                # outgoing term 2:
                out_2 += P_1 * x[i][u] * a[i][j] * w[u][j]
                # outgoing term 3:
                # for v in range(m):
                #     out_3 += P_2 * w[u][j] * b[u][v] * x[j][v] * a[i][j]
                # # outgoing term 5:
                # inner = 0
                # for k in range(n):
                #     for v in range(m):
                #         inner += b[u][v] * x[k][v]
                # # outgoing term 4:
                # out_4 -= w[u][j] * a[i][j] * P_3# (inner + 1)
                #print(f"Sum_k Sum_v b[{u}][{v}] x[{k}][{v}]: {inner}")
                #out_5 += P_4 * (1 - inner) * w[u][j] * a[i][j]
                energy += x[i][u] * (1 - w[u][j])
                #print(f"Penalties: {out_1} + {out_2} + {out_3} + {out_4} + {out_5} = {out_1 + out_2 + out_3 + out_4 + out_5} ------------ \n")
    
    # # outgoing term 3:
    out_3 = 0
    for u in range(m):
        for j in range(n):
            for v in range(m):
                out_3 += P_2 * w[u][j] * b[u][v] * x[j][v]

    # # outgoing term 4:
    out_4 = 0
    for u in range(m):
        for j in range(n):
            out_4 -= P_3 * w[u][j]
            
    out_5 = 0
    for u in range(m):
        for j in range(n):
            inner = 0
            for v in range(m):
                inner += b[u][v] * x[j][v]
            out_5 += P_5 * (1 - w[u][j]) * (1 - inner)
            
    # # incoming term 1:
    # in_1 = 0
    # for i in range(n):
    #     for u in range(m):
    #         for j in range(n):
    #             for v in range(m):
    #                 in_1 -= x[i][u] * a[j][i] * b[v][u] * x[j][v]
    
    # # incoming term 2:
    # in_2 = 0
    # for i in range(n):
    #     for u in range(m):
    #         for j in range(n):
    #             in_2 += P_1 * x[i][u] * a[j][i] * z[u][j]
    
    # # incoming term 5:
    # in_5 = 0
    # for i in range(n):
    #     for u in range(m):
    #         for j in range(n):
    #             inner = 0
    #             for k in range(n):
    #                 for v in range(m):
    #                     inner += b[v][u] * x[k][v]
    #             in_5 += P_4 * (1 - inner) * z[u][j] * a[j][i]
    
    #print("\nI1NCOMING TERMS: -------------------\n")
    #print(f"Energy so far: {out_1 + out_2 + out_3 + out_4 + out_5}")
    
    in_1 = 0
    in_2 = 0
    in_3 = 0
    in_4 = 0
    in_5 = 0
    for i in range(n):
        for u in range(m):
            #print(f"x[{i}][{u}]: {x[i][u]}")
            for j in range(n):
                k = 0
                for v in range(m):
                    k += b[v][u] * x[j][v]
                    in_1 -= x[i][u] * a[j][i] * b[v][u] * x[j][v] * P_5
                #print(f"Sum b[{u}][v] x[{j}][v]: {k}")
                #print(f"w[{u}][{j}]: {z[u][j]}")
                # outgoing term 2:
                in_2 += P_1 * x[i][u] * a[j][i] * z[u][j]
                # # outgoing term 3:
                # for v in range(m):
                #     in_3 += P_2 * z[u][j] * b[v][u] * x[j][v] * a[j][i]
                # outgoing term 5:
                # inner = 0
                # for k in range(n):
                #     for v in range(m):
                #         inner += b[v][u] * x[k][v]
                # # outgoing term 4:
                # in_4 -= z[u][j] * a[j][i] * P_3# (inner + 1)
                #print(f"Sum_k Sum_v b[{u}][{v}] x[{k}][{v}]: {inner}")
                #in_5 += P_4 * (1 - inner) * z[u][j] * a[j][i]
                energy += x[i][u] * (1 - z[u][j])
                #print(f"Penalties: {in_1} + {in_2} + {in_3} + {in_4} + {in_5} = {in_1 + in_2 + in_3 + in_4 + in_5} ------------ \n")    
    
    # incoming term 3:
    in_3 = 0
    for u in range(m):
        for j in range(n):
            for v in range(m):
                in_3 += P_2 * z[u][j] * b[v][u] * x[j][v]
    
    # incoming term 4:
    in_4 = 0
    for u in range(m):
        for j in range(n):
            in_4 -= P_3 * z[u][j]
    
    in_5 = 0
    for u in range(m):
        for j in range(n):
            inner = 0
            for v in range(m):
                inner += b[v][u] * x[j][v]
            in_5 += P_5 * (1 - z[u][j]) * (1 - inner)
    #print(f"Energy so far: {in_1 + in_2 + in_3 + in_4 + in_5}")
    
    energy += out_1 + out_2 + out_3 + out_4 + out_5 + in_1 + in_2 + in_3 + in_4 + in_5
    
    # for i in range(n):
    #     for u in range(m):
    #         energy -= 0.1 * x[i][u]
    return energy

def main(adj_matrix1, adj_matrix2):
    """
    Finds the best combination of binary variables that minimizes the result of the CUBO cost function.

    Args:
        adj_matrix1 (list): The adjacency matrix 1.
        adj_matrix2 (list): The adjacency matrix 2.

    """
    n = len(adj_matrix1)
    m = len(adj_matrix2)
    binary_variables = generate_binary_variables(n, m, 3)

    min = 1000
    best = None
    for variables in binary_variables:
        xs = np.reshape(variables[:n*m], (n, m))
        zs = np.reshape(variables[n*m:2*n*m], (m, n))
        ws = np.reshape(variables[2*n*m:], (m, n))
        # if np.array_equal(xs, [[1, 0, 1],[1, 1, 1]]) and np.array_equal(zs, [[1, 1, 1],[1, 1, 1]]) and np.array_equal(ws, [[1, 1, 1],[1, 1, 1]]):
        #     print("Found")
        if np.array_equal(xs, [[0, 0, 0],[0, 0, 0]]):
            continue
        result = compute_arithmetic_function_2(xs, adj_matrix1, adj_matrix2, zs, ws, n, m)
        if result < min:
            min = result
            best = {'xs': xs, 'zs': zs, 'ws': ws}
    print(f"Final result:\n Minimum: {min}\n Best combination:\n {best['xs']}")
    print(f"Z: {best['zs']}")
    print(f"W: {best['ws']}")

def compute(adj_matrix1, adj_matrix2):
    """
    Computes the arithmetic function using the given adjacency matrices.

    Parameters:
    adj_matrix1 (list): The first adjacency matrix.
    adj_matrix2 (list): The second adjacency matrix.
    """
    
    n = len(adj_matrix1)
    m = len(adj_matrix2)
    # z = np.array([[1, 1, 1], 
    #               [1, 1, 1]])
    # w = np.array([[1, 1, 1], 
    #               [1, 1, 1]])
    # print(compute_arithmetic_function_2(
    #                             [[1, 0, 1],
    #                              [1, 1, 1]], adj_matrix1, adj_matrix2, z, w, n, m))
    z = np.array([[1, 1], 
                  [0, 1],
                  [0, 0]])
    w = np.array([[0, 0], 
                  [1, 0],
                  [1, 1]])
    x = [[1, 1, 0], 
         [0, 1, 1]]
    # z = np.array([[0, 0], 
    #               [1, 0],
    #               [0, 0]])
    # w = np.array([[0, 1], 
    #               [0, 0],
    #               [0, 0]])
    # x = [[0, 1, 0], [1, 0, 1]]
    print(compute_arithmetic_function_2(x, adj_matrix1, adj_matrix2, z, w, n, m))
    
    P_1 = m + 1
    P_2 = m*m + m
    P_3 = m
    P_4 = 1# P_3 - 0.5
    
    energy = 0
    
    a = adj_matrix1
    b = adj_matrix2
    for i in range(n):
        for u in range(m):
            print(f"x[{i}][{u}]: {x[i][u]}")
            for j in range(n):
                out_1 = 0
                out_2 = 0
                out_3 = 0
                out_4 = 0
                out_5 = 0 
                k = 0
                for v in range(m):
                    k += b[u][v] * x[j][v]
                    out_1 -= x[i][u] * a[i][j] * b[u][v] * x[j][v]
                print(f"Sum b[{u}][v] x[{j}][v]: {k}")
                print(f"w[{u}][{j}]: {w[u][j]}")
                # outgoing term 2:
                out_2 += P_1 * x[i][u] * a[i][j] * w[u][j]
                # outgoing term 3:
                for v in range(m):
                    out_3 += P_2 * w[u][j] * b[u][v] * x[j][v] * a[i][j]
                # outgoing term 5:
                inner = 0
                for k in range(n):
                    for v in range(m):
                        inner += b[u][v] * x[k][v]
                # outgoing term 4:
                out_4 -= w[u][j] * a[i][j] * P_3 #(inner + 1)
                print(f"Sum_k Sum_v b[{u}][{v}] x[{k}][{v}]: {inner}")
                #out_5 += P_4 * (1 - inner) * w[u][j] * a[i][j]
                print(f"Penalties: {out_1} + {out_2} + {out_3} + {out_4} + {out_5} = {out_1 + out_2 + out_3 + out_4 + out_5} ------------ \n")
                energy += out_1 + out_2 + out_3 + out_4 + out_5
            
    print("\nI1NCOMING TERMS: -------------------\n")
            
    in_1 = 0
    in_2 = 0
    in_3 = 0
    in_4 = 0
    in_5 = 0
    a = adj_matrix1
    b = adj_matrix2
    for i in range(n):
        for u in range(m):
            print(f"x[{i}][{u}]: {x[i][u]}")
            for j in range(n):
                in_1 = 0
                in_2 = 0
                in_3 = 0
                in_4 = 0
                in_5 = 0
                k = 0
                for v in range(m):
                    k += b[v][u] * x[j][v]
                    in_1 -= x[i][u] * a[j][i] * b[v][u] * x[j][v]
                print(f"Sum b[{u}][v] x[{j}][v]: {k}")
                print(f"w[{u}][{j}]: {z[u][j]}")
                # outgoing term 2:
                in_2 += P_1 * x[i][u] * a[j][i] * z[u][j]
                # outgoing term 3:
                for v in range(m):
                    in_3 += P_2 * z[u][j] * b[v][u] * x[j][v] * a[j][i]
                # outgoing term 5:
                inner = 0
                for k in range(n):
                    for v in range(m):
                        inner += b[v][u] * x[k][v]
                # outgoing term 4:
                in_4 -= z[u][j] * a[j][i] * P_3# (inner + 1)
                print(f"Sum_k Sum_v b[{u}][{v}] x[{k}][{v}]: {inner}")
                #in_5 += P_4 * (1 - inner) * z[u][j] * a[j][i]
                print(f"Penalties: {in_1} + {in_2} + {in_3} + {in_4} + {in_5} = {in_1 + in_2 + in_3 + in_4 + in_5} ------------ \n")   
                energy += in_1 + in_2 + in_3 + in_4 + in_5
                
    print(f"Energy: {energy}")
    

# Define your adjacency matrices here
adj_matrix1 = np.array([[0, 1], 
                        [0, 0]])
# adj_matrix2 = np.array([
#     [0, 1, 1],
#     [0, 0, 0],
#     [0, 0, 0]
# ])
adj_matrix2 = np.array([
    [0, 1, 0],
    [0, 0, 1],
    [0, 0, 0]
])
# adj_matrix2 = np.array([[0, 1], [1, 0]])

compute(adj_matrix1, adj_matrix2)

#main(adj_matrix1, adj_matrix2)