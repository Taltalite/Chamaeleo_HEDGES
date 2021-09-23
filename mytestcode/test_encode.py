import numpy as np
import copy

class A:
    def __init__(self, seq):
        self.seq = seq

    def init_func(self, aa):
        self.seq = aa.seq + 1


class B:
    def __init__(self):
        self.arra = [A(seq) for seq in range(10)]

    def func(self):
        h1 = self.arra[1]
        h2 = self.arra[2]
        print(h1.seq)
        h1.init_func(h2)
        print(h1.seq)


b = B()
b.func()

bigval = 9999.99
ps = 10
ar = np.array([bigval for _ in range(ps)])
ar.resize((6,))
ar[2] = 0
print(1.e10)

hypo = {'111': 111, '222': 222, '333': 333}
hypostack = [copy.deepcopy(hypo) for i in range(10)]
# hypostack=[]
# for i in range(10):
#     hypostack.append(hypo)

hypostack[2].update({'444':4})

def fun(c_hypo):
    c_hypo['111'] = 222

fun(hypostack[0])
print(hypostack)

