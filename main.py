'''
Źródła:
    - https://en.wikipedia.org/wiki/Coherent_states
    - https://en.wikipedia.org/wiki/Creation_and_annihilation_operators
    - https://pl.wikipedia.org/wiki/Stan_koherentny
'''

from sympy import *
from math import sqrt, pi
import matplotlib.pyplot as plt
import sys

sys.setrecursionlimit(10000)

PHI_COUNT = 20  # Ilość funkcji phi (N)
ALPHA = 1  # Alfa
T_STEP = 0.05  # Krok czasowy

x = Symbol("x")  # Zmienna x
phis = []  # Lista phi


def fact_sqrt(n):
    '''
    Funkcja do szybszego obliczania pierwiastka z silni
    '''

    result = 1

    for i in range(1, n + 1):
        result *= sqrt(i)

    return result


def normalize_wave_function(phi):
    '''
    Funkcja do normowania funkcji falowej
    '''

    return phi / sqrt(integrate(abs(phi) ** 2, (x, -oo, oo)))


def calculate_phis():
    global phis

    # Obliczanie phi_0
    last_phi = normalize_wave_function(exp(-0.5 * x**2))

    print(f"phi_0 Calculated!")

    phis.append(last_phi)

    for i in range(PHI_COUNT):
        # Obliczanie pochodnej używając operatora kreacji
        last_phi = normalize_wave_function(
            simplify(x * last_phi - diff(last_phi)))

        phis.append(last_phi)

        print(f"phi_{i + 1} Calculated!")

    return phis


def calculate_avg_x(state):
    '''
    Obliczanie średniej wartości x
    '''

    return integrate(x * abs(state) ** 2, (x, -oo, oo))


def calculate_coherent(t):
    '''
    Normalna metoda obliczania stanu koherentnego
    '''

    sum = 0

    # Sumowanie phi (i współczynników)
    for i in range(PHI_COUNT):
        # Dodawanie do sumy
        sum += simplify(exp(-I * (i + 0.5) * t) *
                        (ALPHA ** i) / fact_sqrt(i) * phis[i])

    print(f"Calculating Coherent State (Normal), t={t}...")

    # Końcowy wynik
    return exp(-0.5 * (abs(ALPHA) ** 2)) * sum


def plot_state(t):
    state = calculate_coherent(t)
    graph_re = re(state)
    graph_im = im(state)
    graph_abs2 = abs(state) ** 2

    id = int(t * 100)
    
    print("Plotting state:", id, "...")

    plotting.plot(graph_re, graph_im, xlim=[-10, 10], ylim=[-1, 1], show=False).save(
        f"./states_re_im/coherent_t_{id}.png")
    plt.close()

    plotting.plot(graph_abs2, xlim=[-10, 10], ylim=[-1, 1], show=False).save(
        f"./states_abs2/coherent_t_{id}.png")
    plt.close()

    print("Successfully plotted state!")


def main():
    calculate_phis()

    for i in range(0, int(4 * pi / T_STEP) + 1):
        plot_state(T_STEP * i)


main()
