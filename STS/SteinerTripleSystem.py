import numpy as np, random

def is_sts(X, T):
    v = len(X)
    k = len(T)
    n = v*(v-1)//6
    if k != n:
        raise ValueError('The number of triples is incorrect: ' + str(k) + ' != ' + str(n))
    S = set()
    for t in T:
        S.add(frozenset([t[0], t[1]]))
        S.add(frozenset([t[0], t[2]]))
        S.add(frozenset([t[1], t[2]]))
    if len(S) != v*(v-1)//2:
        #print S
        return False
    return True

''' Cycle Graph 
'''
class CycleGraphComponent(object):
    ''' Initializes the cycle graph with the given pivots and element pair (a,b)
    '''
    def __init__(self, a, b, pivots, C):
        self.pair = (a, b)
        self.pivots = pivots[:]
        self.C = C
 
    def get_pivots(self):
        return self.pivots[:]
    
    def __len__(self):
        return len(self.pivots)

''' SteineTripleSystem class.
'''
class SteinerTripleSystem(object):
    ''' Initializes a STS of order n. If a set of triples T is given, a test is run to check if it defines a valid STS.
        If the triples are not given then they are generated
    '''
    def __init__(self, n, T = None):
        self.order = n
        self.X = range(n)
        if (T == None):
            self.__generate__()
        elif (is_sts(self.X, T)):
            self.T = T
        else:
            raise ValueError('Not a valid STS')

    ''' Generation of a triple system with a given algorithm
    '''
    def __generate__(self):
        n = self.order
        if n % 6 == 3:
            t = (n - 3) // 6
            Z = range(2 * t + 1)
            T = lambda x: x[0] + (2 * t + 1) * x[1]

            sts = [[(i, 0), (i, 1), (i, 2)] for i in Z] + [
                [(i, k), (j, k), (((t + 1) * (i + j)) % (2 * t + 1), (k + 1) % 3)] for k in range(3) for i in Z for j in
                Z if i != j]

        elif n % 6 == 1:
            t = (n - 1) / 6
            N = range(2 * t)
            T = lambda x, y: x + y * t * 2 if (x, y) != (-1, -1) else n - 1

            L1 = lambda i, j: (i + j) % (int((n - 1) / 3))
            L = lambda i, j: L1(i, j) / 2 if L1(i, j) % 2 == 0 else t + (L1(i, j) - 1) / 2

            sts = [[(i, 0), (i, 1), (i, 2)] for i in range(t)] + \
                  [[(-1, -1), (i, k), (i - t, (k + 1) % 3)] for i in range(t, 2 * t) for k in [0, 1, 2]] + \
                  [[(i, k), (j, k), (L(i, j), (k + 1) % 3)] for k in [0, 1, 2] for i in N for j in N if i < j]
        else:
            raise ValueError("Steiner triple systems only exist for n = 1 mod 6 or n = 3 mod 6")
        self.T = set(list(map(lambda x: tuple(sorted(set(map(T, x)))), sts)))

    ''' Checks stric equality (same elements)
    '''
    def __equals__(self, S):
        return self.order == S.n and self.T == S.T

    ''' Return a copy of the set of triples
    '''
    def get_triples(self):
        return [t[:] for t in self.T]

    ''' Finds the triple having the pair (a, b)
    '''
    def get_triplet_with_pair(self, a, b):
        return filter(lambda t: a in t and b in t, self.T)[0][:]
    
    '''  Creates the cycle switching graph with pair (a, b)    
    '''
    # Get cycle graph Gab
    def cycle_graph(self, a, b):
        Gab = []
        T1 = self.get_triples()
        
        # All triples having 'a' or 'b' but not both
        Tab = filter(lambda t: (a in t or b in t) and not (a in t and b in t), T1)
        # Continue until no more triples remain to be added
        while Tab != []:
            C = []
            s = a
            ab = a+b
            nxt_t = (t for t in Tab if (a in t)).next()
            p = (x for x in nxt_t if (x != a)).next()
            P = []
            # Verify if cycle is complete
            while p not in P:
                P.append(p)
                # Append triple and pivot
                C.append([nxt_t, p])
                # Toggle elements a and b
                s = ab-s
                nxt_t = (t for t in Tab if (p in t and s in t)).next()
                p = (x for x in nxt_t if (x != s and x != p)).next()
                Tab.remove(nxt_t) 
            Gab.append(CycleGraphComponent(a, b, P, C))
        return Gab
    
    ''' Get all cycle graphs
    '''
    def all_cycle_graphs(self):
        G = {}
        for i in range(self.order-1):
            for j in range(i+1, self.order):
                G[(i,j)] = self.cycle_graph(i,j)
        return G
    
    '''
        Applies permutation phi of elements to produce an isomorphic STS
        parameters
            phi: Permutation to be applied
        returns
            Isomorphic STS
    ''' 
    def permute(self, phi):
        T1 = self.get_triples()
        perm_sts = map(lambda t: map(lambda x: phi[x], t), T1)
        perm_sts = sorted(map(sorted, perm_sts))
        return SteinerTripleSystem(self.order, perm_sts)

    '''
        Applies a random permutation of elements to produce an isomorphic STS
        returns
            1. Random permutation phi
            2. Isomorphic STS
    '''
    def random_permute(self):
        phi = np.random.permutation(self.X)
        return phi, self.permute(phi)

    '''
        Applies the cycle switch operation to produce another (possibly non-isomorphic) STS.
        If first argument is provide, then the second is ignored 
        input
            C: Component to perform the cycle switch
            t: triple (a,b,p) for finding Cycle graph Gab an component with pivot p or
        returns
            2. New STS resulting from the Switch cycle operation
    '''
    def cycle_switch(self, Cabp=None, t = None):
        if Cabp is None:
            if t is None:
                raise ValueError('Triple or cycle component required')
            elif t not in self.T:
                raise ValueError('Triple not valid')
            else:
                a, b, p = t
                Gab = self.cycle_graph(a, b)
                for Cabp in Gab:
                    if not p in Cab.pivots:
                        continue
                    C = Cabp
        elif not type(Cabp) == CycleGraphComponent:
            raise ValueError('C must be a CycleGraphComponent')
        C = Cabp.C[:]
        #print('C: ' + str(C))
        a, b = Cabp.pair[:]
        t = self.get_triplet_with_pair(a, b)
        T1 = self.get_triples()
        T2 = [t1[:] for t1 in T1]
        #T2.remove(t)
        #print(T2)
        for t in [c[0] for c in C]:
            t = t[:]
            T2.remove(t)
            toadd = []
            if a in t:
                t.remove(a)
                toadd.append(b)
            elif b in t:
                t.remove(b)
                toadd.append(a)
            t.extend(toadd)
            t.sort()
            T2.append(t)
        #T2.extend(C)
        T2 = sorted([sorted(t) for t in T2])
        #print(T2)
        return SteinerTripleSystem(self.order, T2)
    
    '''
        Applies the cycle switch operation to a random pair
        inputs
            noisom: If true, a not isomorphic operation is required the resulting STS is non-isomorphic
                    to this. Otherwise the cycle switch could create an isomorphic STS (or not)
    '''
    def random_cycle_switch(self, notisom = False):
        if self.order < 10 and notisom:
            raise ValueError('A not-isomorphic cycle switch cannot be performed with this size')
        while True:
            # Random cycle graph
            a, b = random.sample(self.X, 2)
            Gab = self.cycle_graph(a, b)
            # If result is required to be non-isomorphic
            if notisom and len(Gab) == 1:
                continue
            # Select a random component
            i = random.choice(range(len(Gab)))
            S2 = self.cycle_switch(Gab[i])
            return S2

    def __eq__(self, other):
        return set(map(frozenset, self.T)) == set(map(frozenset, other.T))

'''
'''
def str_sts_decomp(TC, output = None):
    import string
    str0Z = '0123456789' + string.ascii_uppercase
    str_cycles = {}
    for k in sorted(TC.keys()):
        str_cycles[k] = [] 
        for c in (TC[k]):
            str_c = map(lambda t: '{}{}{}'.format(str0Z[t[0]],str0Z[t[1]],str0Z[t[2]]), map(sorted, c))
            str_cycle = str_c[0]
            for i in range(1, len(str_c)):
                str_cycle = str_cycle + '--{}--{}'.format((set(str_c[i-1]) & set(str_c[i])).pop(), str_c[i])
            str_cycles[k].append(str_cycle)
    if output== None:
        return str_cycles
    else:
        wf = open(output, 'w')
        for k in sorted(str_cycles.keys()):
            wf.write('{}: \t{}\n'.format(k, str_cycles[k]))
        wf.flush()
        wf.close()

