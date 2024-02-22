from sympy import *
import math
import sys
import os


sys.setrecursionlimit(10000)


def get_dir(path):
    filepath = os.path.join(os.path.dirname(__file__), path)

    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    return filepath


def fact_sqrt(n):
    """
    Faster method for calculation square root of factorial
    """

    result = 1

    for i in range(1, n + 1):
        result *= math.sqrt(i)

    return result


def normalize_wave_function(func, x):
    """
    Norming wave function
    """

    return func / sqrt(integrate(abs(func) ** 2, (x, -oo, oo)))


def calculate_average_of_wave_function(func, x):
    """
    Average of wave function
    """

    return re(integrate(x * abs(func) ** 2, (x, -oo, oo)).evalf())


def calculate_phis(x, count):
    phis = []

    # Calculating phi_0
    last_phi = normalize_wave_function(exp(-0.5 * (x**2)), x)

    print(f"phi_0 Calculated!")

    phis.append(last_phi)

    for i in range(count):
        # Calculating derivative using creation operator
        last_phi = normalize_wave_function(simplify(x * last_phi - diff(last_phi)), x)

        phis.append(last_phi)

        print(f"phi_{i + 1} Calculated!")

    print("All phis were calculated")

    return phis
