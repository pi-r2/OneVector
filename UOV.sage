""" 
This file contains the code used to generate a UOV instance, and an oracle that returns a secret vector.
I also included the code that allows to sign a message and verify a message for completion.
"""

###UOV###
def complete_basis(B, E) :
    """
    This function takes as input a list of linearily independant vectors B in E and 
    completes them into a basis naively.
    If the list is empty, a random basis of E is returned
    """
    res = list(B) #Avoid modifying the input.
    if len(res) == 0 :
        b = E.random_element()
        while b.is_zero() :
            b = E.random_element()
        res.append(b)

    while len(res) != E.dimension() :
        b = E.random_element()
        if not(b in span(res)) :
            res.append(b)
    return res


def ROMKeyGen(q,k,v) :
    """
    This function generates a fake UOV public key by generating a collection of random invertible matrices, 
    symmetrizing them in odd characteristic.
    It returns a pair (A,F), G where (A,F) are identity matrices to be consistent with the non-ROM KeyGen. 
    """
    FF = GF(q)
    n = k+v
    G = [ matrix(FF,complete_basis([], FF**(n))) for _ in range(k)]
    if (q%2) != 0 : #In even characteristic, we represent the public key cannonically with a triangular matrix. We can also represent it 
                    #with an arbitrary sum of this triangular matrix and any symmetric matrix. This property holds in any field.
        for i in range(k):
            G[i] = FF(2)**(-1)*(G[i] + G[i].transpose())
    return (matrix.identity(n), [matrix.identity(n) for _ in range(k)]), G 

def KeyGen(q,k,v, verbose = False):
    """
    This function generates a key pair for UOV parameters k,v where k is the dimension of the oil subspace and v the dimension of the vinegar subspace.
    To do this, we sample random matrices with a block of zero of size to produce the private key, and impose that they must be invertible.
    We sample a random invertible change of variables A, compute G = F\circ A, and return the pair (A,F), G.
    (A,F) is the UOV private key, G is the UOV public key. The corresponding systems are obtained by evaluating the quadratic forms over FF[x1,...,xn].
    Important notice: This code is a demonstration tool and should not used for applications where security matters. 
    """

    n = k+v 
    FF = GF(q)
    RR = None 
    if verbose :
        RR = PolynomialRing(FF, 'x', n)
    #Define A
    A = complete_basis([], FF**(n))
    A = matrix(FF, A)
    
    #Define F, G
    F = []
    G = []
    for e in range(k):
        Fe = matrix(FF, n, n)
        while Fe.determinant() == 0 :
            B0 = matrix(FF,k,k)
            B1 = random_matrix(FF,k,v)
            B2 = random_matrix(FF,v,k)
            B3 = random_matrix(FF,v,v)
            Fe = block_matrix([[B0,B1],[B2,B3]])
        if FF.characteristic() != 2 : 
            Fe = FF(2)**(-1)*(Fe+Fe.transpose())
        Ge = A.transpose()*Fe*A
        F.append(Fe)
        G.append(Ge)
        
    if verbose :
        print("The public key system is:")
        x = vector(RR.gens())
        for Ge in G :
            print(x*Ge*x)
    return ((A,F), G) 
    
def Sign(PrK, M, depth = 0, verbose = False): 
    
    max_depth = 10 #Fail after 10 retries (arbitrary value). Probability : (p_failure)^10 
    if depth > max_depth:
        raise Exception("The signature procedure has failed.")
    
    A,F = PrK
    q = A.base_ring().cardinality()
    k = len(F)
    n = F[0].dimensions()[0]
    v = n- k
    A_1 = A.inverse()
    FF = GF(q)
    RR = PolynomialRing(FF, 'x', k)
    X = RR.gens()

    #Generate random values for vinegar variables.
    Y_2 = [FF.random_element() for _ in range(v)]

    Y = vector(RR, list(X) + Y_2)

    #Use linear algebra to solve for oil variables.
    system = []
    for e in range(k):
        system.append(Y*F[e]*Y)
    if verbose:
        print("The signer solves the system: ", system)
    #This system is linear in x_1, ..., x_k 
    mat = [ [ 0 for _ in range(k)] for _ in range(k)]
    for i in range(k) :
        for j in range(k) :
            mat[i][j] = system[i].coefficient(X[j])

    S = matrix(FF, mat)
    B = vector([system[i]([0 for _ in range(k)]) for i in range(k)])
    target = M-B 

    try: 
        Y_1 = S.solve_right(target)
    except:
        return Sign(PrK, M, depth+1, verbose) 
    Y_hat = vector(FF, list(Y_1) + Y_2) #Recombination
    
    if verbose :
        print("Number of retries : ", depth)
    
    #Use Secret Key to return an 'X' solution.
    return A_1*Y_hat

def Verify(M, X, G, verbose = False): 
    k = len(G)
    n = len(X)
    v = n-k 
    q = G[0].base_ring().cardinality()
    FF = GF(q)
    RR = PolynomialRing(FF, 'x', n) 
    
    ver = True
    for e in range(k):
        if verbose :
            print(("G"+str(e)+"(X)= ", X*G[e]*X," M"+str(e)+"= ", M[e]))
        ver = ( X*G[e]*X==M[e])
        if not(ver):
            return False
    return True

def oracle(A, F) :
    """
    This function takes as input the secret key and returns a vector of the secret subspace O.
    This task is challenging without the secret key and is the basis of the security of UOV.
    """
    k = len(F)
    O = span(A.inverse().columns()[:k])
    #This is exactly the secret subspace.

    return O.random_element()