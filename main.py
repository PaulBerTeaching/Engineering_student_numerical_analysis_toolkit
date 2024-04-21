import matplotlib.pyplot as plt
from numpy import log10
from numpy import polynomial as poly
import root_finding_problem.root_finding_methods as tf2

if __name__ == '__main__':
    studied_polynom = poly.Polynomial((-10, 0, 4, 1))
    LOW_LIMIT_INTERVAL = 0
    HIGH_LIMIT_INTERVAL = 5
    ITER_MAX = 100
    EPSILON = 1e-9
    inter = studied_polynom.roots()

    res_bisect, it_bisect = tf2.BisectionSolverClass(f_x=studied_polynom, a_n=LOW_LIMIT_INTERVAL,
                                                     b_n=HIGH_LIMIT_INTERVAL,
                                                     nb_iteration=ITER_MAX, tolerance=EPSILON).solve()

    res_trisect, it_trisect = tf2.TrisectionSolverClass(f_x=studied_polynom, a_n=LOW_LIMIT_INTERVAL,
                                                        b_n=HIGH_LIMIT_INTERVAL,
                                                        nb_iteration=ITER_MAX, tolerance=EPSILON).solve()

    res_fixed_point, it_fixed_point = tf2.FixedPointSolverClass(f_x=lambda x: pow(10 / (x + 4), 1 / 2), p_0=0,
                                                                nb_iteration=ITER_MAX,
                                                                tolerance=EPSILON).solve()

    res_newton, it_newton = tf2.NewtonSolverClass(f_x=studied_polynom, df_x=poly.Polynomial((0, 8, 3)), p_0=5,
                                                  nb_iteration=ITER_MAX,
                                                  tolerance=EPSILON).solve()

    res_secant, it_secant = tf2.SecantSolverClass(f_x=studied_polynom, p_0=(0, 5), nb_iteration=ITER_MAX,
                                                  tolerance=EPSILON).solve()

    res_rf, it_rf = tf2.RegularFalsiSolverClass(f_x=studied_polynom, p_0=(0, 5), nb_iteration=ITER_MAX,
                                                tolerance=EPSILON).solve()

    res_stef, it_stef = tf2.SteffensenSolverClass(f_x=lambda x: pow(10 / (x + 4), 1 / 2), p_0=0,
                                                  nb_iteration=ITER_MAX, tolerance=EPSILON).solve()

    fig2, ax2 = plt.subplots()
    ax2.plot(list(range(1, it_bisect + 1)),
             log10([abs(inter[-1] - x) for x in res_bisect]),
             label='Bisection method',
             linewidth=2.0)
    ax2.plot(list(range(1, it_bisect + 1)),
             [0] + [log10(1 / pow(2, x)) for x in list(range(1, it_bisect))],
             label='Bisection method decreased',
             linestyle='-.',
             color='lightblue',
             linewidth=2.0)
    ax2.plot(list(range(1, it_trisect + 1)),
             log10([abs(inter[-1] - x) for x in res_trisect]),
             label='Trisection method',
             color='red',
             linewidth=2.0)
    ax2.plot(list(range(1, it_trisect + 1)),
             [0] + [log10(1 / pow(3, x)) for x in list(range(1, it_trisect))],
             label='Trisection method decreased',
             linestyle='-.',
             color='salmon',
             linewidth=2.0)
    ax2.plot(list(range(1, it_fixed_point + 1)),
             log10([abs(inter[-1] - x) for x in res_fixed_point]),
             label='Fixed point method',
             color='forestgreen',
             linewidth=2.0)
    ax2.plot(list(range(1, it_fixed_point + 1)),
             log10([(pow(0.197642353760524, x) / (1 - 0.197642353760524)) * 10 for x in
                    list(range(1, it_fixed_point + 1))]),
             label='Fixed point method decreased',
             linestyle='-.',
             color='lightgreen',
             linewidth=2.0)
    ax2.plot(list(range(1, it_newton + 1)),
             log10([abs(inter[-1] - x) for x in res_newton]),
             label='Newton method',
             color='turquoise',
             linewidth=2.0)
    ax2.plot(list(range(1, it_secant + 1)),
             log10([abs(inter[-1] - x[1]) for x in res_secant]),
             label='Secant method',
             color='purple',
             linewidth=2.0)
    ax2.plot(list(range(1, it_rf + 1)),
             log10([abs(inter[-1] - x[1]) for x in res_rf]),
             label='Regula Falsi method',
             color='mediumorchid',
             linewidth=2.0)
    ax2.plot(list(range(1, it_stef + 1)),
             log10([abs(inter[-1] - x) for x in res_stef]),
             label='Steffessen method',
             color='goldenrod',
             linewidth=2.0)

    ax2.legend()

    plt.show()
