"""
Resources:
    - https://en.wikipedia.org/wiki/Coherent_states
    - https://en.wikipedia.org/wiki/Creation_and_annihilation_operators
    - https://pl.wikipedia.org/wiki/Stan_koherentny
"""

from sympy import *
from math import sqrt, pi
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import sys

sys.setrecursionlimit(10000)

MODE = "average"  # Modes: coherent, average, squezzed
SAVE_IMAGES = True  # Save images
PHI_COUNT = 20  # Count of phi functions (N)
ALPHA = 1  # Alpha
T_STEP = 0.25  # Time step

x = Symbol("x")  # Variable x
phis = []  # Phi list
avgs = []  # Averages of coherent states


def fact_sqrt(n):
    """
    Faster method for calculation square root of factorial
    """

    result = 1

    for i in range(1, n + 1):
        result *= sqrt(i)

    return result


def normalize_wave_function(phi):
    """
    Norming wave function
    """

    return phi / sqrt(integrate(abs(phi) ** 2, (x, -oo, oo)))


def calculate_average_of_wave_function(func):
    """
    Average of wave function
    """

    return re(integrate(x * abs(func) ** 2, (x, -oo, oo)).evalf())


def calculate_phis():
    global phis

    # Calculating phi_0
    last_phi = normalize_wave_function(exp(-0.5 * (x**2)))

    print(f"phi_0 Calculated!")

    phis.append(last_phi)

    for i in range(PHI_COUNT):
        # Calculating derivative using creation operator
        last_phi = normalize_wave_function(simplify(x * last_phi - diff(last_phi)))

        phis.append(last_phi)

        print(f"phi_{i + 1} Calculated!")


def calculate_coherent(t):
    """
    Normal method for calculation coherent state
    """

    sum = 0

    # Summing phis (and coefficents)
    for i in range(PHI_COUNT):
        # Add to sum
        sum += simplify(exp(-I * (i + 0.5) * t) * (ALPHA**i) / fact_sqrt(i) * phis[i])

    print(f"Calculating Coherent State (Normal), t={t}...")

    # Final result
    return exp(-0.5 * (abs(ALPHA) ** 2)) * sum


def plot_state(t):
    id = int(t * 100)

    print(f"Plotting state: id={id}...")

    if MODE == "coherent":
        state = calculate_coherent(t)
        graph_re = re(state)
        graph_im = im(state)
        graph_abs2 = abs(state) ** 2

        tmp_plt = plotting.plot(
            graph_re, graph_im, xlim=[-10, 10], ylim=[-1, 1], show=(not SAVE_IMAGES)
        )

        if SAVE_IMAGES:
            tmp_plt.save(f"./states_re_im/coherent_t_{id}.png")

        plt.close()

        tmp_plt = plotting.plot(
            graph_abs2, xlim=[-10, 10], ylim=[-1, 1], show=(not SAVE_IMAGES)
        )

        if SAVE_IMAGES:
            tmp_plt.save(f"./states_abs2/coherent_t_{id}.png")

        plt.close()
    elif MODE == "average":
        state = calculate_coherent(t)
        avg = calculate_average_of_wave_function(state)

        avgs.append(avg)

        c = Circle((avg, 0), 0.2)
        figure, axes = plt.subplots()
        axes.set_xlim([-3, 3])
        axes.set_ylim([-2, 2])
        axes.set_aspect(1)
        axes.add_artist(c)

        if SAVE_IMAGES:
            plt.savefig(f"./states_avg/coherent_t_{id}.png")
        else:
            plt.show()

        plt.close()
    elif MODE == "squezzed":
        pass

    print("Successfully plotted state!")


def main():
    calculate_phis()

    for i in range(0, int(4 * pi / T_STEP) + 1):
        plot_state(T_STEP * i)

    if MODE == "average":
        a = 0

        for x in avgs:
            a += x

        a /= len(avgs)

        print(f"Total average is: {a}")


main()
