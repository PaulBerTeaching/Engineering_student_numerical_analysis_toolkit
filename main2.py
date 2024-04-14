# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import matplotlib.pyplot as plt
from numpy import log10
from numpy import polynomial as poly
import root_finding_methods_expert as tf
import root_finding_methods_expert2 as tf2

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    studied_polynom = poly.Polynomial((-10, 0, 4, 1))
    LOW_LIMIT_INTERVAL = 0
    HIGH_LIMIT_INTERVAL = 5
    ITER_MAX = 100
    EPSILON = 1e-9
    inter = studied_polynom.roots()

    res_bisect, it_bisect = tf.bisection_method(studied_polynom, LOW_LIMIT_INTERVAL, HIGH_LIMIT_INTERVAL,
                                                  ITER_MAX, EPSILON)

    res_trisect, it_trisect = tf.trisection_method(studied_polynom, LOW_LIMIT_INTERVAL,
                                                     HIGH_LIMIT_INTERVAL, ITER_MAX, EPSILON)

    res_fixed_point, it_fixed_point = tf.fixed_point_method(
        lambda x: pow(10 / (x + 4), 1 / 2), 0, ITER_MAX, EPSILON)

    res_newton, it_newton = tf.newton_method(studied_polynom, poly.Polynomial((0, 8, 3)), 5, ITER_MAX,
                                               EPSILON)

    res_secant, it_secant = tf.secant_method(studied_polynom, (0, 5), ITER_MAX, EPSILON)

    res_rf, it_rf = tf.regula_falsi_method(studied_polynom, (0, 5), ITER_MAX, EPSILON)

    res_stef, it_stef = tf.steffensen_method(lambda x: pow(10 / (x + 4), 1 / 2), 0, ITER_MAX, EPSILON)

    res_bisect2, it_bisect2 = tf2.BissectionSolverClass(f_x=studied_polynom, a_n=LOW_LIMIT_INTERVAL, b_n=HIGH_LIMIT_INTERVAL,
                                                  nb_iteration=ITER_MAX, tolerance=EPSILON).solve()

    res_trisect2, it_trisect2 = tf2.TrissectionSolverClass(f_x=studied_polynom, a_n=LOW_LIMIT_INTERVAL, b_n=HIGH_LIMIT_INTERVAL,
                                                  nb_iteration=ITER_MAX, tolerance=EPSILON).solve()

    # PROBLEME CONCEPTION ICI

    res_newton2, it_newton2 = tf2.NewtonSolverClass(f_x=studied_polynom, df_x=poly.Polynomial((0, 8, 3)), p_0=5, nb_iteration=ITER_MAX,
                                               tolerance=EPSILON).solve()

    res_secant2, it_secant2 = tf2.SecantSolverClass(f_x=studied_polynom, p_0=(0, 5), nb_iteration=ITER_MAX, tolerance=EPSILON).solve()

    res_rf2, it_rf2 = tf2.RegularFalsiSolverClass(f_x=studied_polynom, p_0=(0, 5), nb_iteration=ITER_MAX, tolerance=EPSILON).solve()

    res_stef2, it_stef2 = tf2.SteffensenSolverClass(f_x=lambda x: pow(10 / (x + 4), 1 / 2), p_0=0, nb_iteration=ITER_MAX, tolerance=EPSILON).solve()

    fig2, ax2 = plt.subplots()
    ax2.plot(list(range(1, it_bisect2 + 1)),
             log10([abs(inter[-1] - x) for x in res_bisect2]),
             label='Bisection method',
             linewidth=2.0)
    ax2.plot(list(range(1, it_bisect2 + 1)),
             [0] + [log10(1/pow(2, x)) for x in list(range(1, it_bisect2))],
             label='Bisection method decreased',
             linestyle='-.',
             color='lightblue',
             linewidth=2.0)
    ax2.plot(list(range(1, it_trisect2 + 1)),
             log10([abs(inter[-1] - x) for x in res_trisect2]),
             label='Trisection method',
             color='red',
             linewidth=2.0)
    ax2.plot(list(range(1, it_trisect2 + 1)),
             [0] + [log10(1 / pow(3, x)) for x in list(range(1, it_trisect2))],
             label='Trisection method decreased',
             linestyle='-.',
             color='salmon',
             linewidth=2.0)
    ax2.plot(list(range(1, it_newton2 + 1)),
             log10([abs(inter[-1] - x) for x in res_newton2]),
             label='Newton method',
             color='turquoise',
             linewidth=2.0)
    ax2.plot(list(range(1, it_secant2 + 1)),
             log10([abs(inter[-1] - x[1]) for x in res_secant2]),
             label='Secant method',
             color='purple',
             linewidth=2.0)
    ax2.plot(list(range(1, it_rf2 + 1)),
             log10([abs(inter[-1] - x[1]) for x in res_rf2]),
             label='Regula Falsi method',
             color='mediumorchid',
             linewidth=2.0)
    ax2.plot(list(range(1, it_stef2 + 1)),
             log10([abs(inter[-1] - x) for x in res_stef2]),
             label='Steffessen method',
             color='goldenrod',
             linewidth=2.0)

    ax2.legend()

    plt.show()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
