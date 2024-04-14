import matplotlib.pyplot as plt
from numpy import log10
from numpy import polynomial as poly
import root_finding_methods as rfm
import root_finding_methods_expert as tf


if __name__ == '__main__':
    studied_polynom = poly.Polynomial((-10, 0, 4, 1))
    LOW_LIMIT_INTERVAL = 0
    HIGH_LIMIT_INTERVAL = 5
    ITER_MAX = 100
    EPSILON = 1e-9
    inter = studied_polynom.roots()

    res_bisect, it_bisect, error_bisect = rfm.bisection_method(studied_polynom, LOW_LIMIT_INTERVAL, HIGH_LIMIT_INTERVAL,
                                                               ITER_MAX, EPSILON, inter[-1])

    res_bisect2, it_bisect2 = tf.bisection_method(studied_polynom, LOW_LIMIT_INTERVAL, HIGH_LIMIT_INTERVAL,
                                                  ITER_MAX, EPSILON)

    res_trisect, it_trisect, error_trisect = rfm.trisection_method(studied_polynom, LOW_LIMIT_INTERVAL,
                                                                   HIGH_LIMIT_INTERVAL, ITER_MAX, EPSILON, inter[-1])

    res_trisect2, it_trisect2 = tf.trisection_method(studied_polynom, LOW_LIMIT_INTERVAL,
                                                     HIGH_LIMIT_INTERVAL, ITER_MAX, EPSILON)

    res_fixed_point, it_fixed_point, error_fixed_point = rfm.fixed_point_method(
        lambda x: 1 / 2 * pow(10 - pow(x, 3), 1 / 2), 0, ITER_MAX, EPSILON, inter[-1])

    #res_fixed_point2, it_fixed_point2 = tf.fixed_point_method(
    #    lambda x: 1 / 2 * pow(10 - pow(x, 3), 1 / 2), 0, ITER_MAX, EPSILON)

    res_fixed_point2, it_fixed_point2 = tf.fixed_point_method(
        lambda x: pow(10 / (x + 4), 1 / 2), 0, ITER_MAX, EPSILON)

    res_newton, it_newton, error_newton = rfm.newton_method(studied_polynom, poly.Polynomial((0, 8, 3)), 5, ITER_MAX,
                                                            EPSILON, inter[-1])

    res_newton2, it_newton2 = tf.newton_method(studied_polynom, poly.Polynomial((0, 8, 3)), 5, ITER_MAX,
                                               EPSILON)

    res_secant, it_secant, error_secant = rfm.secant_method(studied_polynom, 0, 5, ITER_MAX, EPSILON, inter[-1])

    res_secant2, it_secant2 = tf.secant_method(studied_polynom, (0, 5), ITER_MAX, EPSILON)

    res_rf, it_rf, error_rf = rfm.regula_falsi_method(studied_polynom, 0, 5, ITER_MAX, EPSILON, inter[-1])

    res_rf2, it_rf2 = tf.regula_falsi_method(studied_polynom, (0, 5), ITER_MAX, EPSILON)

    res_stef, it_stef, error_stef = rfm.steffensen_method(lambda x: pow(10 / (x + 4), 1 / 2), 0, ITER_MAX, EPSILON,
                                                          inter[-1])

    res_stef2, it_stef2 = tf.steffensen_method(lambda x: pow(10 / (x + 4), 1 / 2), 0, ITER_MAX, EPSILON)

    fig2, ax2 = plt.subplots()
    ax2.plot(list(range(1, it_bisect2 + 1)),
             log10([abs(inter[-1] - x) for x in res_bisect2]),
             label='Bisection method',
             linewidth=2.0)
    ax2.plot(list(range(1, it_bisect2 + 1)),
             log10([1/pow(2,x)for x in list(range(1, it_bisect2 + 1))]),
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
             log10([1 / pow(3, x) for x in list(range(1, it_trisect2 + 1))]),
             label='Trisection method decreased',
             linestyle='-.',
             color='salmon',
             linewidth=2.0)
    ax2.plot(list(range(1, it_fixed_point2 + 1)),
             log10([abs(inter[-1] - x) for x in res_fixed_point2]),
             label='Fixed point method',
             color='forestgreen',
             linewidth=2.0)
    ax2.plot(list(range(1, it_fixed_point2 + 1)),
             log10([(pow(0.197642353760524, x)/(1-0.197642353760524))*10 for x in list(range(1, it_fixed_point2 + 1))]),
             label='Fixed point method decreased',
             linestyle='-.',
             color='lightgreen',
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
