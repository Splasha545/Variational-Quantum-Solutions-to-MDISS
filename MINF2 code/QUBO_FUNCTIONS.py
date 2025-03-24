from collections import Counter

class QUBO_Problem:
    def __init__(self, matrix, n, function_name, h=None):
        self.function_name = function_name
        self.dag_matrix = matrix
        self.n = n
        self.linear, self.quadratic, self.cubic = produce_terms_from_graphs(matrix, n, function=function_name, h=h)
        self.num_terms = len(self.linear) + len(self.quadratic) + len(self.cubic)
        

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
    """Sums all identical terms together to reduce function to strictly unique terms with real coefficient.
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

def produce_terms_from_graphs(mat, n, function, h=None):
    
    linear = []
    quadratic = []
    cubic = []
    
    if function == "three-chain":
        linear, quadratic = three_chain_quadratic_to_terms(mat, n)
        linear, quadratic, cubic = clean_coefficients(linear, quadratic, cubic)
        linear, quadratic, cubic = collect_terms(linear, quadratic, cubic)
        
    elif function == "three-chain-naive":
        linear, quadratic = three_chain_naive_to_terms(mat, n)
        linear, quadratic, cubic = clean_coefficients(linear, quadratic, cubic)
        linear, quadratic, cubic = collect_terms(linear, quadratic, cubic)
        
    elif function == "three-chain-naive-2":
        linear, quadratic = three_chain_naive_2_to_terms(mat, n)
        linear, quadratic, cubic = clean_coefficients(linear, quadratic, cubic)
        linear, quadratic, cubic = collect_terms(linear, quadratic, cubic)
        
    elif function == "two-link":
        linear, quadratic = two_link_simple_post_processing_to_terms(mat, n)
        linear, quadratic, cubic = clean_coefficients(linear, quadratic, cubic)
        linear, quadratic, cubic = collect_terms(linear, quadratic, cubic)
    
    elif function == "bh-links":
        linear, quadratic = bh_links_naive_to_terms(mat, h, n)
        linear, quadratic, cubic = clean_coefficients(linear, quadratic, cubic)
        linear, quadratic, cubic = collect_terms(linear, quadratic, cubic)
        
    else:
        raise ValueError("Function not recognized.")
    
    #print(f"Size: {len(linear)}, linear terms: {linear}")
    #print(f"Size: {len(quadratic)}, quadratic terms: {quadratic}")
    #print(f"Size: {len(cubic)}, cubic terms: {cubic}")
    return linear, quadratic, cubic

## ----------------------------- ARITHMETIC EXPANSIONS --------------------------------- ##

def three_chain_quadratic_to_terms(mat, n):
    
    linear = []
    quadratic = []
    
    def x(i):
        return (i)
    def y(i):
        return (n + i)
    
    for i in range(n):
        quadratic.append(((x(i), y(i)), 0.5))
        
    for i in range(n):
        inner1 = 0
        inner2 = 0
        for j in range(n):
            inner1 += mat[i][j]
            inner2 += mat[j][i]
            
        linear.append((x(i), -1*inner1))
        linear.append((y(i), -1*inner2))
        
    return linear, quadratic

def three_chain_naive_to_terms(mat, n):
    
    linear = []
    quadratic = []
    
    def x(i):
        return (i)
    
    for i in range(n):
        linear.append((x(i), 0.5))
        for j in range(n):
            inner = 0
            for k in range(n):
                inner += mat[i][k] * mat[k][j]
            quadratic.append(((x(i), x(j)), -1 * inner))
    
    return linear, quadratic

def three_chain_naive_2_to_terms(mat, n):
    
    linear = []
    quadratic = []
    
    def x(i):
        return (i)
    def y(i):
        return (n + i)
    
    for i in range(n):
        linear.append((x(i), 0.5))
        linear.append((y(i), 0.5))
        for j in range(n):
            inner = 0
            for k in range(n):
                inner += mat[i][k] * mat[k][j]
            quadratic.append(((x(i), y(j)), -2 * inner))
    
    return linear, quadratic
                

def two_link_simple_post_processing_to_terms(mat, n):
    
    linear = []
    quadratic = []
    
    def x(i, j):
        return (i*n+j)
    def y(i, j):
        return (n*n + i*n+j)
    def z(i, j):
        return (2*n*n + i*n+j)
    
    for i in range(n):
        for j in range(n):
            inner = 0
            for k in range(n):
                inner += mat[i][k] * mat[k][j]
            linear.append((x(i,j), 15 - 10 * inner))
            linear.append((y(i,j), 25 - 10 * inner))
            linear.append((z(i,j), 10 * n - 5))
            
            quadratic.append(((z(i,j), x(i,j)), -1))
            quadratic.append(((z(i,j), y(i,j)), 10 * n))
    
    return linear, quadratic

def bh_links_to_terms(mat, h, n):
    
    linear = []
    quadratic = []
    cubic = []

    P = n

    def x(i):
        return (i)
    def y(i):
        return (n + i)
    def z(i):
        return (2*n + i)

    for i in range(n):
        linear.append((x(i), 2 * (1 - h[i]) - 1))
        linear.append((y(i), 2 * h[i] + n - 1))
        linear.append((z(i), 3 * n + 11))
        
        for j in range(n):
            linear.append((x(i), mat[i][j]))
            linear.append((y(i), -n * mat[i][j]))
            quadratic.append(((y(i), z(i)), n * mat[i][j]))
            linear.append((z(i), -1 * (2 * n + 7) * mat[i][j]))
    
    return linear, quadratic

def bh_links_naive_to_terms(mat, h, n):
    
    linear = []
    quadratic = []
    cubic = []

    P = n

    def x(i):
        return (i)
    def y(i):
        return (n + i)
    
    for i in range(n):
        inner = 0
        for j in range(n):
            inner += mat[i][j]
        linear.append((x(i), P * inner))
        linear.append((y(i), P * (inner - 1)**2))

        for j in range(n):
            quadratic.append(((y(i), x(j)), -1 * mat[i][j] * (1 - h[i]) * h[j]))
    
    return linear, quadratic



## ----------------------------- MAIN --------------------------------- ##

# def produce_cubic_problem_from_graph(n, matrix):
#     linear, quadratic, cubic = produce_terms_from_graphs(matrix, n, function="3-chain")
#     cubic_instance = CubicProgram(linear, quadratic, cubic, n*3)
#     return cubic_instance, n, G, posG
