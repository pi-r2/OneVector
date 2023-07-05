"""
This file contains an implementation of the attack described in the paper "One Vector to rule them all", Pebereau 2023.
It only uses native sagemath functions for improved clarity.
"""



def attack(G,v) :
    """
    TODO: FAIRE MARCHER

    This function takes as input a UOV public key G, and a vector v in the secret subspace O.
    It returns a basis of O in polynomial time.
    """

    k = len(G)
    n = G[0].dimensions()[0]
    Current = matrix.identity(n)
    while n > 3*k :
        J = matrix([v*g for g in G]) 
        B = matrix(J.right_kernel().basis())        
        Ghat = [B*g*B.transpose() for g in G] #restriction of G to the kernel of J
        vhat = B.solve_left(v)
        
        G = Ghat
        v = vhat
        n = n-k 
        Current = B*Current
    return one_vector_to_rule_them_all(G,v)*Current

def one_vector_to_rule_them_all(G,v) :
    """ 
    This function is the one that completes the attack.

    Given G a set of k quadratic forms admitting a common totally isotropic subspace O
    of dimension at least (n-k)/2, and v in O, find a basis of O as a whole.
    """
    
    J = matrix([v*g for g in G]) 
    B = matrix(J.right_kernel().basis())
    
    Ghat = [B*g*B.transpose() for g in G] #restriction of G to the kernel of J
    B2 = []
    for g in Ghat :
        for b in g.kernel().basis() :
            if len(B2) == 0 or b not in span(B2) :
                B2.append(b)

    B3 = matrix(span(B2).basis())
        
    C = B3*B
    return C 
    
 
q,k,v= 256, 44, 68
load("UOV.sage")
print('Parameters : q=',q,' k=',k,' v=',v)
print('Key Generation')
(A,F), G = KeyGen(q, k, v)
print('Query the oracle for a vector in O.')
v = oracle(A,F)
print("Run the attack")
B = one_vector_to_rule_them_all(G,v)
print('Dimension of the subspace we found ',B.dimensions()[0], ' vs dimension of O', k)
print('Maximum rank of B^TGB ?')
print(max([(B*g*B.transpose()).rank() for g in G]))




