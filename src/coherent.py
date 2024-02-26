from sympy import *
from math import pi
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import base


MODE = "coherent"  # Modes: coherent, average
SAVE_IMAGES = True  # Save images
PHI_COUNT = 20  # Count of phi functions (N)
ALPHA = 1  # Alpha
T_STEP = 0.25  # Time step


avgs = []  # Averages


def calculate_coherent(t, phis):
    """
    Normal method for calculation coherent state
    """

    sum = 0

    # Summing phis (and coefficents)
    for i in range(len(phis)):
        # Add to sum
        sum += simplify(
            exp(-I * (i + 0.5) * t) * (ALPHA**i) / base.fact_sqrt(i) * phis[i]
        )

    print(f"Calculating Coherent State (Normal), t={t}...")

    # Final result
    return exp(-0.5 * (abs(ALPHA) ** 2)) * sum


def plot_state(t, phis, x):
    id = int(t * 100)

    print(f"Plotting state: id={id}...")

    if MODE == "coherent":
        state = calculate_coherent(t, phis)
        graph_re = re(state)
        graph_im = im(state)
        graph_abs2 = graph_re**2 + graph_im**2

        tmp_plt = plotting.plot(
            graph_re, graph_im, xlim=[-10, 10], ylim=[-1, 1], show=(not SAVE_IMAGES)
        )

        if SAVE_IMAGES:
            tmp_plt.save(
                base.get_dir(f"../generated/coherent_re_im/coherent_t_{id}.png")
            )

        plt.close()

        tmp_plt = plotting.plot(
            graph_abs2, xlim=[-10, 10], ylim=[-1, 1], show=(not SAVE_IMAGES)
        )

        if SAVE_IMAGES:
            tmp_plt.save(
                base.get_dir(f"../generated/coherent_abs2/coherent_t_{id}.png")
            )

        plt.close()
    elif MODE == "average":
        state = calculate_coherent(t, phis)
        avg = base.calculate_average_of_wave_function(state, x)

        avgs.append(avg)

        c = Circle((avg, 0), 0.2)
        figure, axes = plt.subplots()
        axes.set_xlim([-3, 3])
        axes.set_ylim([-2, 2])
        axes.set_aspect(1)
        axes.add_artist(c)

        if SAVE_IMAGES:
            plt.savefig(base.get_dir(f"../generated/coherent_avg/coherent_t_{id}.png"))
        else:
            plt.show()

        plt.close()

    print("Successfully plotted state!")


def main():
    x = Symbol("x")
    phis = base.calculate_phis(x, PHI_COUNT)

    for i in range(0, int(4 * pi / T_STEP) + 1):
        plot_state(T_STEP * i, phis, x)

    if MODE == "average":
        a = 0

        for tmp in avgs:
            a += tmp

        a /= len(avgs)

        print(f"Total average is: {a}")


main()
