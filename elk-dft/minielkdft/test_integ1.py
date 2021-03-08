# comparison with poly4i

import numpy as np

def f(x):
    return x**2 + 1.1

def integ_fx(x):
    return 1.0/3.0*x**3 + 1.1*x

xa = np.array([1.1, 2.0, 3.0, 4.0])
ya = f(xa)

print("xa = ", xa)
print("ya = ", ya)

x = 4.1
print(integ_fx(x) - integ_fx(xa[0]))