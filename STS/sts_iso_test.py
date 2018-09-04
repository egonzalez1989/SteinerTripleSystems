from sts_isomorphism import *

S1 = SteinerTripleSystem(15)
_, S2 = S1.random_permute()
print(find_stsq_isomorphism(quasigroup_from_sts(S1), quasigroup_from_sts(S2)))