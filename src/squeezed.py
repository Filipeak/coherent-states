from sympy import *
from math import pi, cosh, sqrt
import matplotlib.pyplot as plt
import base


SAVE_IMAGES = True  # Save images
PHI_COUNT = 40  # Count of phi functions (N)
ALPHA = 1  # Alpha
BETA = 1  # Beta
T_STEP = 4 * pi / 20  # Time step


def calculate_squeezed(t, phis, x):
    """
    Normal method for calculation squeezed state
    """

    print(f"Calculating Squeezed State (Normal), t={t}...")

    sum = 0

    # Summing phis (and coefficents)
    for k in range(12):
        for l in range(12):
            for m in range(2 * k):
                fact_2k = base.fact_sqrt(2 * k)

                sum += (
                    ((BETA**k) * fact_2k * ((-ALPHA) ** m) * ((ALPHA) ** l))
                    / (
                        factorial(k)
                        * factorial(m)
                        * factorial(l)
                        * base.fact_sqrt(2 * k - m)
                    )
                    * (fact_2k * base.fact_sqrt(2 * k - m + l))
                    * exp(-I * (2 * k - m + l + 0.5) * t)
                    * phis[2 * k - m + l]
                )

    # Final result
    result = exp(-0.5 * (abs(ALPHA) ** 2)) * sqrt(1 / cosh(BETA)) * sum
    result = base.normalize_wave_function(result, x)

    print("Successfully calculated squeezed state!")

    return result


def plot_state(t, phis, x):
    id = int(t * 100)

    print(f"Plotting state: id={id}...")

    state = calculate_squeezed(t, phis, x)
    graph_re = re(state)
    graph_im = im(state)
    # graph_abs2 = abs(state) ** 2

    print("Computed graph functions")

    tmp_plt = plotting.plot(graph_re, graph_im, adaptive=False, nb_of_points=50, show=(not SAVE_IMAGES))

    if SAVE_IMAGES:
        tmp_plt.save(base.get_dir(f"../generated/squeezed_re_im/squezzed_t_{id}.png"))

    plt.close()

    # tmp_plt = plotting.plot(graph_abs2, adaptive=False, nb_of_points=50, show=(not SAVE_IMAGES))

    # if SAVE_IMAGES:
    #     tmp_plt.save(base.get_dir(f"../generated/squeezed_abs2/squezzed_t_{id}.png"))

    # plt.close()

    print("Successfully plotted state!")


def main():
    x = Symbol("x")
    phis = base.calculate_phis(x, PHI_COUNT)

    for i in range(0, int(4 * pi / T_STEP) + 1):
        plot_state(T_STEP * i, phis, x)


main()