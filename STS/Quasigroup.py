

class Quasigroup(object):
    def __init__(self, X_set, mult_dict):
        self.X = X_set
        self.order = len(X_set)
        if self.verify_quasi(mult_dict):
            self.op = mult_dict.copy()

    def verify_quasi(self, mult_table):
        # Trust by now. Only verify if it is latin square
        return True

    def operate(self, x, y):
        if not x in self.X:
            return ValueError('{} is not an element'.format(x))
        if not y in self.X:
            return ValueError('{} is not an element'.format(y))
        return self.op[(x, y)]

    def generate(self, Y):
        Y = Y[:]
        if not all(elem in self.X for elem in Y):
            raise ValueError('Not every member of Y in the set')
        n = len(Y)
        mtable = {}
        Ynew = Y[:]
        # Find new elements until closure
        while len(Ynew) > 0:
            Ynew2 = []
            for x in Y:
                for y in Ynew:
                    z = self.operate(x, y)
                    mtable[(x, y)] = z
                    mtable[(y, x)] = z
                    if z not in Y and z not in Ynew2:
                        Ynew2.append(z)
            Y.extend(Ynew2)
            Ynew = Ynew2
        return Quasigroup(Y, mtable)