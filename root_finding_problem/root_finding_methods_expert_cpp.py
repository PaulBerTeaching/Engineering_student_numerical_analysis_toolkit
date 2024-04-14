"""
This module shows some examples of root finding methods

Classes: IntervalError(Exception)
Functions:
check_interval(f_x: Callable, a_0: float, b_0: float, tolerance: float)
bisection_method(f_x: Callable, a_n: float, b_n: float, nb_iteration: int, tolerance: float)
root_nearest_point(f_x: Callable, c_1: float, c_2: float)
regula_falsi_function(f_x: Callable, x: tuple) -> tuple
"""
from typing import Callable, Union

# import the module
from ctypes import cdll
from numpy import polynomial as poly

# load the library
lib = cdll.LoadLibrary('./bissection.so')

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

    if abs(f_x(a_0)) < tolerance:
        return [a_0]

    if abs(f_x(b_0)) < tolerance:
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


def secant_function(f_x: Callable, x_ab: tuple) -> float:
    """This is a helper function to decide which point is the nearest to the root of a function f(x)

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param x_ab: Tuple values to update following the regula falsi condition
    :type x_ab: tuple
    :return: The new tuple with the regula falsi condition
    :rtype: tuple
    """
    if f_x(x_ab[1]) - f_x(x_ab[0]) == 0:
        raise ZeroDivisionError

    return x_ab[1] - f_x(x_ab[1]) * (x_ab[1] - x_ab[0]) / (f_x(x_ab[1]) - f_x(x_ab[0]))


def regula_falsi_function(f_x: Callable, x_ab: tuple) -> tuple:
    """This is a helper function to decide which point is the nearest to the root of a function f(x)

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param x_ab: Tuple values to update following the regula falsi condition
    :type x_ab: tuple
    :return: The new tuple with the regula falsi condition
    :rtype: tuple
    """

    c_n = secant_function(f_x, x_ab)

    return (x_ab[1], c_n) if f_x(c_n) * f_x(x_ab[1]) < 0 else (x_ab[0], c_n)


def steffensen_function(f_x: Callable, p_n: float) -> float:
    """This is a helper function to decide which point is the nearest to the root of a function f(x)

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param p_n: Value to compute with steffensen
    :type p_n: float
    :return: The new tuple with the regula falsi condition
    :rtype: tuple
    """

    if (f_x(f_x(p_n)) - 2 * f_x(p_n) + p_n) == 0:
        raise ZeroDivisionError

    return p_n - pow(f_x(p_n) - p_n, 2) / (f_x(f_x(p_n)) - 2 * f_x(p_n) + p_n)


def bisection_method(f_x: Callable, a_n: float, b_n: float, nb_iteration: int, tolerance: float):
    """This is the bisection algorithm

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
    :return: List of best estimated values, number of iteration
    :rtype: list, int
    """
    iteration_number = 0
    c_n = check_interval(f_x, a_n, b_n, tolerance)

    if c_n:
        return c_n, iteration_number

    for iteration_number in range(1, nb_iteration):
        c_n.append(a_n + (b_n - a_n) / 2)

        if abs(f_x(c_n[-1])) < tolerance:
            break

        if f_x(a_n) * f_x(c_n[-1]) < 0:
            b_n = c_n[-1]
        elif f_x(c_n[-1]) * f_x(b_n) < 0:
            a_n = c_n[-1]
        else:
            raise IntervalError(a_n, b_n, 2)

    return c_n, iteration_number

def eval_poly(poly_values, x):

    return poly.Polynomial(poly_values)

def bisection_method_cpp(f_x: Callable, a_n: float, b_n: float, nb_iteration: int, tolerance: float):
    """This is the bisection algorithm

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
    :return: List of best estimated values, number of iteration
    :rtype: list, int
    """
    iteration_number = 0
    c_n = check_interval(f_x, a_n, b_n, tolerance)

    if c_n:
        return c_n

    c_n = lib.bissection(f_x, a_n, b_n, nb_iteration, tolerance)

    return c_n