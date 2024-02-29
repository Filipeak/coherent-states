import math
import numpy as np


def derivative(f, a, n=1, h=1e-10, cache={}):
    if n == 0:
        return f(a)

    if (a, n) in cache:
        return cache[(a, n)]

    d1 = derivative(f, a + h, n - 1, h, cache)
    d2 = derivative(f, a - h, n - 1, h, cache)
    d = (d1 - d2) / (2 * h)

    cache[(a, n)] = d

    return d


def integrate_inf(f):
    a = -100000000
    b = +100000000
    n = 10000000

    sum = 0

    for k in range(1, n):
        sum += f(a + k * ((b - a) / n))

    return (b - a) / n * (f(a) / 2 + sum + f(b) / 2)


def phi(a, n, h=0.05, cache={}):
    if n == 0:
        return math.exp(-0.5 * (a**2))

    if (a, n) in cache:
        return cache[(a, n)]

    val = a * phi(a, n - 1, h) - (phi(a + h, n - 1, h) - phi(a - h, n - 1, h)) / (2 * h)

    cache[(a, n)] = val

    return val


N = 4
V = 1

print(phi(V, N) / math.sqrt(integrate_inf(lambda x: phi(x, N))))

import base
import sympy

x = sympy.Symbol("x")
print(base.calculate_phis(x, N)[N].subs(x, V).evalf())
