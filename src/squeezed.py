from sympy import *
from math import pi, cosh, sqrt
import matplotlib.pyplot as plt
import base


SAVE_IMAGES = True  # Save images
PHI_COUNT = 20  # Count of phi functions (N)
ALPHA = 1  # Alpha
BETA = 1  # Beta
MAX_T = 4 * pi  # Max time
SAMPLES_COUNT = 25  # Samples count
GRAPH_RESOLUTION = 250  # Number of points


def calculate_squeezed(t, phis, x):
    """
    Method for calculation squeezed state
    """

    print(f"Calculating Squeezed State, t={t}...")

    sum = 0

    # for n in range(int(PHI_COUNT / 2)):
    #     sum += (
    #         exp(-I * (2 * n + 0.5) * t)
    #         * base.fact_sqrt(2 * n)
    #         / (2**n * factorial(n))
    #         * (-tanh(BETA) ** n)
    #         * phis[2 * n]
    #     )

    # result = sqrt(1 / cosh(BETA)) * sum

    for l in range(4):
        for k in range(4):
            for m in range(2 * k):
                sum += (
                    exp(-I * (2 * k - m + l + 0.5) * t)
                    * (-tanh(BETA) ** k)
                    * (ALPHA**l * ((-ALPHA) ** m))
                    * factorial(2 * k)
                    / (2**k * factorial(k) * factorial(l) * factorial(m))
                    / factorial(2 * k - m)
                    * base.fact_sqrt(2 * k - m + l)
                    * phis[2 * k - m + l]
                )

    result = exp(-0.5 * abs(ALPHA) ** 2) * sqrt(1 / cosh(BETA)) * sum

    result = base.normalize_wave_function(result, x)

    print("Successfully calculated squeezed state!")

    return result


def plot_state(t, phis, x):
    id = int(t * 100)

    print(f"Plotting state: id={id}...")

    state = calculate_squeezed(t, phis, x)
    graph_re = re(state)
    graph_im = im(state)
    graph_abs2 = graph_re**2 + graph_im**2

    print("Computed graph functions")

    tmp_plt = plotting.plot(
        graph_re,
        graph_im,
        xlim=[-10, 10],
        ylim=[-1.5, 1.5],
        adaptive=False,
        nb_of_points=GRAPH_RESOLUTION,
        show=(not SAVE_IMAGES),
    )

    if SAVE_IMAGES:
        tmp_plt.save(base.get_dir(f"../generated/squeezed_re_im/squezzed_t_{id}.png"))

    plt.close()

    tmp_plt = plotting.plot(
        graph_abs2,
        xlim=[-10, 10],
        ylim=[-1.5, 1.5],
        adaptive=False,
        nb_of_points=GRAPH_RESOLUTION,
        show=(not SAVE_IMAGES),
    )

    if SAVE_IMAGES:
        tmp_plt.save(base.get_dir(f"../generated/squeezed_abs2/squezzed_t_{id}.png"))

    plt.close()

    print("Successfully plotted state!")


def main():
    x = Symbol("x")
    phis = base.calculate_phis(x, PHI_COUNT)
    t_step = MAX_T / SAMPLES_COUNT

    for i in range(0, int(SAMPLES_COUNT) + 1):
        plot_state(t_step * i, phis, x)


main()
