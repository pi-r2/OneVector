"""
This file contains an implementation of the attack described in the paper "One Vector to rule them all", Pebereau 2023.
It only uses native sagemath functions for improved clarity.
"""



def attack(G,v) :
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
        if len(B2) == k :
            break
    B3 = matrix(B2)
        
    C = B3*B
    return C
    
def in_secret_subspace(G, v) :
    """
    This function takes as input a vector x and a UOV public key G, and returns True if the vector belongs to the secret subspace,
    and False otherwise. 
    """
    n = G[0].dimensions()[0]
    k = len(G)

    for g in G : #Sanity check
        if v*g*v != 0 :
            return False
    
    J = matrix([v*g for g in G]) 
    B = matrix(J.right_kernel().basis())
    
    for g in G :
        ghat = B*g*B.transpose() #Restriction to K(v)
        if ghat.rank() > 2*(n-2*k) :
            return False 
    return True 
 

    
 
q,k,v= 256, 44, 68

load("UOV.sage")
print('Parameters : q=',q,' k=',k,' v=',v)
print('Key Generation')
(A,F), G = KeyGen(q, k, v)
print('Query the oracle for a vector in O.')
v = oracle(A,F)
print("Run the attack")
B = attack(G,v)
print('Dimension of the subspace we found ',B.dimensions()[0], ' vs dimension of O', k)
print('Maximum rank of B^TGB ?')
print(max([(B*g*B.transpose()).rank() for g in G]))

tests = [GF(2).random_element() for _ in range(10)]
test_vectors = []
for i in range(10) :
    if tests[i] :
        test_vectors.append(oracle(A,F))
    else :
        test_vectors.append( Sign((A,F), vector([0 for _ in range(k)])))
print("Now we test the in_secret_subspace algorithm: ")
for i in range(10) :
    val = in_secret_subspace(G, test_vectors[i])
    print("The tests yields ",val, " and the answer is ", bool(tests[i]))
