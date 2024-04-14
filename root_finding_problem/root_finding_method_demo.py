"""
This module shows some examples of root finding methods

Classes: IntervalError
Functions:
check_interval(f_x: Callable, a_0: float, b_0: float, tolerance: float)
bisection_method(f_x: Callable, a_n: float, b_n: float, nb_iteration: int, tolerance: float)
"""
from typing import Callable


class IntervalError(Exception):
    """This is a children class of exception focusing on interval errors for bisection and trisection methods

    :param a_n: Low limit of studied interval
    :type a_n: float
    :param b_n: High limit of studied interval
    :type b_n: float
    :param error_code: Type of error to display
    :type error_code: int
    """

    def __init__(self, a_n, b_n, error_code, *args):
        """Constructor method
        """
        super().__init__(args)
        self.a_n = a_n
        self.b_n = b_n

        # Error message thrown is saved in msg
        if error_code == 1:
            self.msg = 'No/multiple roots for stater interval [ %f, %f ]' % (a_n, b_n)
        elif error_code == 2:
            self.msg = 'Multiple roots in interval [ %f, %f ]' % (a_n, b_n)

    def __str__(self):
        """Print error message
        """
        return repr(self.msg)


def check_interval(f_x: Callable, a_0: float, b_0: float, tolerance: float):
    """This is the checking function for interval for bisection and trisection methods

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param a_0: Low limit of studied interval
    :type a_0: float
    :param b_0: High limit of studied interval
    :type b_0: float
    :param tolerance: Precision limit
    :type tolerance: float
    :raises IntervalError: In case of bad interval selected (multiple or no root.s)
    :return: low or high interval value if they are a sufficient approximation else return empty
    :rtype: list
    """

    if f_x(a_0) * f_x(b_0) > 0:
        raise IntervalError(a_0, b_0, 1)

    if abs(f_x(a_0)) <= tolerance:
        return [a_0]

    if abs(f_x(b_0)) <= tolerance:
        return [b_0]

    return []


def root_nearest_point(f_x: Callable, c_1: float, c_2: float):
    """This is a helper function to decide which point is the nearest to the root of a function f(x)

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param c_1: Low limit of studied interval
    :type c_1: float
    :param c_2: High limit of studied interval
    :type c_2: float
    :return: The best root between c_1 and c_2
    :rtype: float
    """
    return c_1 if min(abs(f_x(c_1)), abs(f_x(c_2))) == abs(f_x(c_1)) else c_2


def bisection_method_demo(f_x: Callable, a_n: float, b_n: float, nb_iteration: int, tolerance: float):
    """This is the bisection algorithm with all series in output in order to make a demonstration of the algorithm
    efficiency

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param a_n: Low limit of studied interval
    :type a_n: float
    :param b_n: High limit of studied interval
    :type b_n: float
    :param nb_iteration: Maximum iteration limit
    :type nb_iteration: int
    :param tolerance: Precision limit
    :type tolerance: float
    :raises IntervalError: In case of bad interval selected (multiple or no root.s)
    :return: Last value estimated, number of iteration, list with all values of a_n and b_n
    :rtype: list, int, list, list
    """

    c_n = check_interval(f_x, a_n, b_n, tolerance)
    a_n = [a_n]
    b_n = [b_n]

    if c_n:
        return c_n, 0, c_n, a_n, b_n

    for iteration_number in range(1, nb_iteration):
        c_n.append(a_n[-1] + (b_n[-1] - a_n[-1]) / 2)

        if abs(f_x(c_n[-1])) <= tolerance:
            break

        if f_x(a_n[-1]) * f_x(c_n[-1]) < 0:
            a_n.append(a_n[-1])
            b_n.append(c_n[-1])
        elif f_x(c_n[-1]) * f_x(b_n[-1]) < 0:
            a_n.append(c_n[-1])
            b_n.append(b_n[-1])
        else:
            raise IntervalError(a_n[-1], b_n[-1], 2)

    return c_n, iteration_number, a_n, b_n


def trisection_method_demo(f_x: Callable, a_n: float, b_n: float, nb_iteration: int, tolerance: float):
    """This is the trisection algorithm with all series in output in order to make a demonstration of the algorithm
    efficiency

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param a_n: Low limit of studied interval
    :type a_n: float
    :param b_n: High limit of studied interval
    :type b_n: float
    :param nb_iteration: Maximum iteration limit
    :type nb_iteration: int
    :param tolerance: Precision limit
    :type tolerance: float
    :raises IntervalError: In case of bad interval selected (multiple or no root.s)
    :return: Last value estimated, number of iteration, list with all values of c_1_n, c_2_n, a_n and b_n
    :rtype: list, int, list, list, list, list
    """

    best_c_n = check_interval(f_x, a_n, b_n, tolerance)

    a_n = [a_n]
    b_n = [b_n]
    c_1_n = []
    c_2_n = []

    if best_c_n:
        return best_c_n, 0, c_1_n, c_2_n, a_n, b_n

    for iteration_number in range(1, nb_iteration):
        c_1_n.append((2 * a_n[-1] + b_n[-1]) / 3)
        c_2_n.append((a_n[-1] + 2 * b_n[-1]) / 3)

        if abs(f_x(c_1_n[-1])) <= tolerance:
            best_c_n.append(c_1_n)
            break

        if abs(f_x(c_2_n[-1])) <= tolerance:
            best_c_n.append(c_2_n[-1])
            break

        best_c_n.append(root_nearest_point(f_x, c_1_n[-1], c_2_n[-1]))

        if f_x(a_n[-1]) * f_x(c_1_n[-1]) < 0:
            a_n.append(a_n[-1])
            b_n.append(c_1_n[-1])
        elif f_x(c_1_n[-1]) * f_x(c_2_n[-1]) < 0:
            a_n.append(c_1_n[-1])
            b_n.append(c_2_n[-1])
        elif f_x(c_2_n[-1]) * f_x(b_n[-1]) < 0:
            a_n.append(c_2_n)
        else:
            raise IntervalError(a_n, b_n, 2)

    return best_c_n, iteration_number, c_1_n, c_2_n, a_n, b_n
