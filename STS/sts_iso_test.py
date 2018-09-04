from sts_isomorphism import *

S1 = SteinerTripleSystem(15)
# Permutation and isomorphic system
f1, S2 = S1.random_permute()
# Permutation found by isomorphism test (f1 and f2 are not equal in general)
f2 = miller_algorithm(S1, S2)
# Testing equality of triple systems after applying both isomorphisms
print(S1 == S2)
print(S1.permute(f1) == S2)
print(S1.permute(f2) == S2)

# Difficult instances for the isomorphism problem
S1 = all_pasches_sts(4)
f1, S2 = S1.random_permute()
f2 = miller_algorithm(S1, S2)
print(S1 == S2)
print(S1.permute(f2) == S2)
print(S1.permute(f2) == S2)