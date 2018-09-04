from Quasigroup import *
from SteinerTripleSystem import *
from itertools import *

'''
    Difficult instances for the STS isomorphism problem
    k: exponent
    returns: A STS of size n=2**k-1 and every cycle switch of length 4
'''
def all_pasches_sts(k):
    n = 2 ** k - 1
    B = set()
    for i in range(1, n // 2 + 1):
        for j in range(i + 1, n):
            B.add(tuple(sorted([i - 1, j - 1, i.__xor__(j) - 1])))
    return SteinerTripleSystem(n, B)

'''
'''
def quasigroup_from_sts(S):
    X = S.X[:]
    n = len(X)
    mult_table = {}
    for i in range(n):
        mult_table[(i, i)] = i
    for t in S.T:
        mult_table[(t[0], t[1])] = mult_table[(t[1], t[0])] = t[2]
        mult_table[(t[0], t[2])] = mult_table[(t[2], t[0])] = t[1]
        mult_table[(t[1], t[2])] = mult_table[(t[2], t[1])] = t[0]
    return Quasigroup(X, mult_table)

'''
'''
def find_quasigroup_isomorphism(Q1, Q2):
    # Check orders
    if Q1.order != Q2.order:
        return False
    # Find a set of generators for Q1
    G1 = find_quasigroup_generators(Q1)
    # For every set of m elements check if it is a well defined isomorphism
    p = permutations(Q2.X, len(G1))
    while True:
        G2 = p.next()
        iso = zip(G1, G2)
        # print(iso)
        d = verify_quasigroup_isomorphism(Q1, Q2, iso)
        if d:
            return list(map(lambda i: d[i], range(Q1.order)))
    return False

'''
'''
def verify_quasigroup_isomorphism(Q1, Q2, f):
    d = dict(f)
    R = [x[0] for x in f]
    for x in R:
        for y in R:
            z = Q1.operate(x, y)
            w = Q2.operate(d[x], d[y])
            if z in d.keys():
                if d[z] != w:
                    return False
            elif w in d.values():
                return False
            else:
                d[z] = w
                R.append(z)
    return d

'''
'''
def find_quasigroup_generators(Q):
    # Random pair to start
    R = set(Q.X[:])
    G = [np.random.choice(list(R))]
    R = R.difference(G)
    #    print(G)
    Q2 = Q.generate(G)
    while Q2.order != Q.order:
        #        print(G, Q2.X)
        R = R.difference(Q2.X)
        G.append(np.random.choice(list(R)))
        Q2 = Q.generate(G)
    return sorted(G)

'''
    Returns a permutation defining the isomorphism
'''
def miller_algorithm(S1, S2):
    Q1, Q2 = quasigroup_from_sts(S1), quasigroup_from_sts(S2)
    return find_quasigroup_isomorphism(Q1, Q2)
# import sts
# # All together
#
# '''
#     Tries one round of Alg.1 to complete a partial bijection.
#     input
#         S0, S1: Steiner triple systems of same length
#         partial_phi: partially defined bijection represented with a list. Not defined elements have value -1
#     returns
#         an extended bijection phi
#         FALSE if the bijection cannot be extended to an isomorphism between STSs
# '''
# def try_isomorphism_round(S0, S1, partial_phi):
#     n = len(S0.X)
#     phi = partial_phi[:]
#     for x in range(n):
#         for y in range(x+1, n):
#             # Verify if bijection is defined on pair (x,y), otherwise try with a new pair
#             u, v = phi[x], phi[y]
#             if u != -1 and v != -1:
#                 # Find the triplet containing (x, y)
#                 t = S0.get_triplet_with_pair(x, y)
#                 #print('xy found with {}, {} : {}'.format(x, y, z))
#                 # Get the third element
#                 z = list(filter(lambda x: x not in [x,y], t))[0]
#                 # if the third element is undefined, then define it with the triple (u, v, w) in S1
#                 if w == -1:
#                     t = S1.get_triplet_with_pair(u, v)
#                     #print('ab found with {}, {}: {}'.format(a, b, t))
#                     partial[z] = list(filter(lambda x: x not in [u,v], t))[0]
#                 else:
#                     # If phi was defined on w, verify that definition is an isomorphism
#                     if sorted([u, v, w]) not in S1.T:
#                         return False
#         if -1 not in partial:
#             return partial
#     return False
#
# '''
#     Tries to complete a bijection until:
#         - phi is completed
#         - phi extension remains the same
#         - phi cannot be extended to a STS isomorphism
#     input
#         S0, S1: Steiner triple systems of same length
#         partial_phi: partially defined bijection represented with a list. Not defined elements have value -1
#     returns
#         an extended bijection phi
#         FALSE if the bijection cannot be extended to an isomorphism between STSs
# '''
# def try_isomorphism(S0, S1, partial_phi):
#     phi = try_isomorphism_round(S0, S1, partial_phi)
#     while phi != partial_phi:
#         partial_phi = phi[:]
#         phi = try_isomorphism_round(S0, S1, phi)
#     return phi
#
# '''
#     Find all possible isomorphisms between S0 an S1 using pivots of the cycle graphs Swap0 and Swap1
#     Swap0 an Swap1 are triples of the form: (Cycles, a, b)
# '''
# def get_isomorphism_from_cycles(S0, S1, Swap0, Swap1):
#     # Verify if the cycle graphs are isomorphic (a simple count must work)
#     Gab, a, b = Swap0
#     Gcd, c, d = Swap1
#     lg = len(Gab)
#     n = len(S0.n)
#     isomorphisms = []
#     lens1, lens2 = list(map(len, Gab)), list(map(len, Gcd))
#     if sorted(lens1) != sorted(lens2):
#         return isomorphisms
#
#     idx1 = np.argsort(lens1)[::-1]
#     idx2 = np.argsort(lens2)[::-1]
#     # We start with longest cycle
#     for i in range(lg):
#         C0 = Gab[idx1[i]]
#         # Pivots of the first cycle
#         P0 = np.array(list(map(lambda t: t[1], C0)))
#         for j in range(lg):
#             C1 = Gcd[j]
#             cl = len(C1)
#             if cl != len(C0):
#                 break
#             # Try with every element in the cycle (forward and backward)
#             for k in range(cl):
#                 phi = [-1] * n
#                 # Pivots of the second cycle
#                 P1 = np.roll(list(map(lambda t: t[1], C1)), k)
#                 for l in range(len(P0)):
#                     phi[P0[l]] = P1[l]
#                 #print(phi)
#                 phi = try_isomorphism(S0, S1, phi)
#                 if phi:
#                     isomorphisms.append(tuple(phi))
#
#                 phi = [-1] * n
#                 # Pivots of the second cycle reversed
#                 P1 = P1[::-1]
#                 for l in range(len(P0)):
#                     phi[P0[l]] = P1[l]
#                 #print(phi)
#                 phi = try_isomorphism(S0, S1, phi)
#                 if phi:
#                     isomorphisms.append(tuple(phi))
#
#     return isomorphisms
#
# '''
#     Find all isomorphisms between Triples S0 and S1 by creating permutations using the cycle graphs
# '''
# def find_isomorphisms(S0, S1):
#     n = len(S0.X)
#     isomorphisms = set()
#     # For every graph
#     for a in range(n-1):
#         for b in range(a+1, n):
#             # Get the cycles of the corresponding pair
#             Gab = S0.cycle_graph(a, b)
#             # We are checking for every graph in the second system
#             for c in range(n-1):
#                 for d in range(c+1, n):
#                     Gcd = S1.cycle_graph(c, d)
#                     isomorphisms.update(set(get_isomorphism_from_cycles(S0, S1, (Gab, a, b), (Gcd, c, d))))
#     return sorted(list(isomorphisms))
#
# '''
#     Find any isomorphism. We use our method
# '''
# def find_isomorphisms(S0, S1):
#     n = len(S0.X)
#     isomorphisms = set()
#     '''
#         Look for a big cycle graph. Longest is of length n-3. Half that size might suffice
#         (which is the min required? Testing with (n-3))'''
#     lmax, l, r = 0, 0, (n-3) // 2
#     for a in range(n-1):
#         for b in range(a+1, n):
#             # Get the cycles of the corresponding pair
#             Gab = S0.cycle_graph(a, b)
#             l = max(map(len, Gab.pivots))
#             if l > lmax:
#                 lmax = l
#             # We are checking for every graph in the second system
#             for c in range(n-1):
#                 for d in range(c+1, n):
#                     Gcd = S1.cycle_graph(c, d)
#                     isomorphisms.update(set(get_isomorphism_from_cycles(S0, S1, (Gab, a, b), (Gcd, c, d))))
#     return sorted(list(isomorphisms))
#
# '''
#     Returns the first isomorphism found
# '''
# def find_an_isomorphism(S0, S1):
#     n = len(S0.X)
#     # For every graph
#     for a in range(n-1):
#         for b in range(a+1, n):
#             # Get the cycles of the corresponding pair
#             Gab = S0.cycle_graph(a, b)
#             # We are checking for every graph in the second system
#             for c in range(n-1):
#                 for d in range(c+1, n):
#                     Gcd = S1.cycle_graph(c, d)
#                     isomorphisms = get_isomorphism_from_cycles(S0, S1, (Gab, a, b), (Gcd, c, d))
#                     if len(isomorphisms) > 0:
#                         return isomorphisms[0]
#     return False