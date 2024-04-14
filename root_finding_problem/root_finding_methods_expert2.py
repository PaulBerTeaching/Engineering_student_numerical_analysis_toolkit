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
import numpy as np
import abc


class IterativeAlgorithm(abc.ABC):

    def __init__(self, nb_iteration: int = 100, tolerance: float = 1e-9):
        self.nb_iteration = nb_iteration
        self.tolerance = tolerance

    @abc.abstractmethod
    def solve(self):
        pass


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


class DividingSearchSolverClass(IterativeAlgorithm):
    """It's not quite a binary search since we divide the space n times, we can say a greedy algorithm but
    not quite sure

        :param f_x: Handle of studied function
        :type f_x: Callable
        :param nb_iteration: Maximum iteration limit
        :type nb_iteration: int
        :param tolerance: Precision limit
        :type tolerance: float
        :return: List of best estimate, number of iteration
        :rtype: list, int
        """
    number_of_division = None

    def __init__(self, f_x: Callable, a_n: float, b_n: float, **kwargs):
        if self.number_of_division < 2:
            raise Exception("Division number must be > 2")

        if f_x(a_n) * f_x(b_n) > 0:
            raise IntervalError(a_n, b_n, 1)

        self.f_x = f_x
        self.a_n = a_n
        self.b_n = b_n
        super().__init__(**kwargs)

    def check_interval(self):
        """This is the checking function for interval for bisection and trisection methods

        :raises IntervalError: In case of bad interval selected (multiple or no root.s)
        :return: low or high interval value if they are a sufficient approximation else return empty
        :rtype: list
        """

        if self.f_x(self.a_n) * self.f_x(self.b_n) > 0:
            raise IntervalError(self.a_n, self.b_n, 1)

        if abs(self.f_x(self.a_n)) < self.tolerance:
            return [self.a_n]
        elif abs(self.f_x(self.b_n)) < self.tolerance:
            return [self.b_n]

        return []

    def solve(self):
        iteration_number = 0
        best_c_n = self.check_interval()

        # Best solution already available ... LUCKY
        if best_c_n:
            return best_c_n, iteration_number

        # Sinon on est parti pour un tour
        for iteration_number in range(1, self.nb_iteration):
            c_n = [self.a_n] + \
                  [((self.number_of_division - n) * self.a_n + n * self.b_n) / self.number_of_division
                   for n in range(1, self.number_of_division)] + \
                  [self.b_n]
            f_c_n = [self.f_x(i) for i in c_n]

            best_c_n.append(c_n[np.argmin(f_c_n)])

            if abs(self.f_x(best_c_n[-1])) < self.tolerance:
                break

            evol_array = [(c_n[n], c_n[n+1]) for n in range(len(f_c_n) - 1) if f_c_n[n] * f_c_n[n + 1] < 0]

            if len(evol_array) > 1:
                raise IntervalError(self.a_n, self.b_n, 2)

            self.a_n = evol_array[0][0]
            self.b_n = evol_array[0][1]

        return best_c_n, iteration_number


class BissectionSolverClass(DividingSearchSolverClass):
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

    def __init__(self, **kwargs):
        self.number_of_division = 2
        super().__init__(**kwargs)


class TrissectionSolverClass(DividingSearchSolverClass):
    """This is the trisection algorithm

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
    :return: List of best estimate, number of iteration
    :rtype: list, int
    """

    def __init__(self, **kwargs):
        self.number_of_division = 3
        super().__init__(**kwargs)


class FixedPointSolverClass(IterativeAlgorithm):
    """This is the trisection algorithm

        :param p_0: Starting point of fixed point algorithm
        :type p_0: float
        :return: List of best estimate, number of iteration
        :rtype: list, int
        """
    f_x = None

    def __init__(self, p_0: Union[float, tuple], **kwargs):
        if self.f_x is None:
            raise Exception("Unknown function for fixed point problem")

        self.p_0 = p_0
        super().__init__(**kwargs)

    def solve(self):
        iteration_number = 0
        p_n = [self.p_0]

        for iteration_number in range(1, self.nb_iteration):
            try:
                p_n.append(self.f_x(p_n[-1]))

                if isinstance(p_n[-1], tuple):
                    if abs(p_n[-1][-1] - self.f_x(p_n[-1])[-1]) < self.tolerance:
                        break
                elif abs(p_n[-1] - self.f_x(p_n[-1])) < self.tolerance:
                    break

            except ZeroDivisionError:
                print('f(x) reach a 0 value - Maximum precision reached or error in selection')
                break

        return p_n, iteration_number + 1


class NewtonSolverClass(FixedPointSolverClass):
    """This is the trisection algorithm

        :param f_x: Handle of studied function
        :type f_x: Callable
        :param df_x: Handle of derivative of studied function
        :type f_x: Callable
        :param p_0: Starting point of fixed point algorithm
        :type p_0: float
        :param nb_iteration: Maximum iteration limit
        :type nb_iteration: int
        :param tolerance: Precision limit
        :type tolerance: float
        :raises ZeroDivisionError: In case of derivative being 0 at p_0
        :return: List of best estimate, number of iteration
        :rtype: list, int
    """

    def __init__(self, f_x: Callable, df_x: Callable, **kwargs):

        if df_x(kwargs['p_0']) == 0:
            raise ZeroDivisionError
        self.f_x = lambda x: x - f_x(x) / df_x(x)
        super().__init__(**kwargs)


class SecantSolverClass(FixedPointSolverClass):
    """This is the trisection algorithm

        :param f_x: Handle of studied function
        :type f_x: Callable
        :param df_x: Handle of derivative of studied function
        :type f_x: Callable
        :param p_0: Starting point of fixed point algorithm
        :type p_0: float
        :param nb_iteration: Maximum iteration limit
        :type nb_iteration: int
        :param tolerance: Precision limit
        :type tolerance: float
        :raises ZeroDivisionError: In case of derivative being 0 at p_0
        :return: List of best estimate, number of iteration
        :rtype: list, int
    """

    def __init__(self, f_x: Callable, **kwargs):
        if f_x(kwargs['p_0'][0]) - f_x(kwargs['p_0'][1]) == 0:
            raise ZeroDivisionError
        self.f_x = lambda x: (x[1], secant_function(f_x, x))
        super().__init__(**kwargs)


class RegularFalsiSolverClass(FixedPointSolverClass):
    """This is the trisection algorithm

        :param f_x: Handle of studied function
        :type f_x: Callable
        :param df_x: Handle of derivative of studied function
        :type f_x: Callable
        :param p_0: Starting point of fixed point algorithm
        :type p_0: float
        :param nb_iteration: Maximum iteration limit
        :type nb_iteration: int
        :param tolerance: Precision limit
        :type tolerance: float
        :raises ZeroDivisionError: In case of derivative being 0 at p_0
        :return: List of best estimate, number of iteration
        :rtype: list, int
    """

    def __init__(self, f_x: Callable, **kwargs):
        if f_x(kwargs['p_0'][0]) - f_x(kwargs['p_0'][1]) == 0:
            raise ZeroDivisionError
        self.f_x = lambda x: regula_falsi_function(f_x, x)
        super().__init__(**kwargs)


class SteffensenSolverClass(FixedPointSolverClass):
    """This is the trisection algorithm

        :param f_x: Handle of studied function
        :type f_x: Callable
        :param df_x: Handle of derivative of studied function
        :type f_x: Callable
        :param p_0: Starting point of fixed point algorithm
        :type p_0: float
        :param nb_iteration: Maximum iteration limit
        :type nb_iteration: int
        :param tolerance: Precision limit
        :type tolerance: float
        :raises ZeroDivisionError: In case of derivative being 0 at p_0
        :return: List of best estimate, number of iteration
        :rtype: list, int
    """

    @staticmethod
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

    def __init__(self, f_x: Callable, **kwargs):
        self.f_x = lambda x: self.steffensen_function(f_x, x)
        super().__init__(**kwargs)
