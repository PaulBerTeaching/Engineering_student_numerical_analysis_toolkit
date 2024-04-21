"""
This module shows some examples of root finding methods

Classes:
IntervalError(Exception)
IterativeAlgorithm(abc.ABC)
DividingSearchSolverClass(IterativeAlgorithm)
BisectionSolverClass(DividingSearchSolverClass)
TrisectionSolverClass(DividingSearchSolverClass)
FixedPointSolverClass(IterativeAlgorithm)
NewtonSolverClass(FixedPointSolverClass)
SecantSolverClass(FixedPointSolverClass)
RegularFalsiSolverClass(FixedPointSolverClass)
SteffensenSolverClass(FixedPointSolverClass)

Functions:
secant_function(f_x: Callable, x_ab: tuple) -> float
regula_falsi_function(f_x: Callable, x_ab: tuple) -> tuple
"""
from typing import Callable, Union
import numpy as np
import abc


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
        """Constructor method - Print an error message depending on the type of error
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


class IterativeAlgorithm(abc.ABC):
    """This is an abstractive class for iterative algorithms

        :param nb_iteration: Number of iterations (By default it is set to 100)
        :type nb_iteration: int, optional
        :param tolerance: Tolerance for root finding methods (By default it is set to 1e-9)
        :type tolerance: float, optional
        """

    def __init__(self, nb_iteration: int = 100, tolerance: float = 1e-9):
        self.nb_iteration = nb_iteration
        self.tolerance = tolerance

    @abc.abstractmethod
    def solve(self):
        """This is an abstract method for solving the root finding method, needs to be overwritten in children class
        """
        pass


def secant_function(f_x: Callable, x_ab: tuple) -> float:
    """This is a helper function to update the values of x_ab using the secant rule.

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param x_ab: Tuple values to update following the secant update rule
    :type x_ab: tuple
    :raises ZeroDivisionError: In case of update rule implies a division by zero
    :return: The new tuple with the secant update rule
    :rtype: tuple
    """
    if f_x(x_ab[1]) - f_x(x_ab[0]) == 0:
        raise ZeroDivisionError

    return x_ab[1] - f_x(x_ab[1]) * (x_ab[1] - x_ab[0]) / (f_x(x_ab[1]) - f_x(x_ab[0]))


def regula_falsi_function(f_x: Callable, x_ab: tuple) -> tuple:
    """This is a helper function to update the values of x_ab using the regula_falsi rule

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param x_ab: Tuple values to update following the regula_falsi update rule
    :type x_ab: tuple
    :return: The new tuple with the regula_falsi update rule
    :rtype: tuple
    """
    c_n = secant_function(f_x, x_ab)

    return (x_ab[1], c_n) if f_x(c_n) * f_x(x_ab[1]) < 0 else (x_ab[0], c_n)


class DividingSearchSolverClass(IterativeAlgorithm):
    """Class for algorithm of dividing search (greedy algorithm/divide and conquer type)

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param a_n: Low limit of studied interval
    :type a_n: float
    :param b_n: High limit of studied interval
    :type b_n: float
    :param number_of_division: Number of space division in order to find a solution (Default value = 2 - Bisection)
    :type number_of_division: int
    :param \**kwargs: Complementary keyword arguments for the optional argument of IterativeAlgorithm class
                   (see IterativeAlgorithm)
    :type kwargs: dict
    :raises Exception: If number_of_division is less than 2
    :raises IntervalError: In case of bad interval selected (multiple or no root.s)
    """

    def __init__(self, f_x: Callable, a_n: float, b_n: float, number_of_division: int = 2, **kwargs):
        if number_of_division < 2:
            raise Exception("Division number must be > 2")

        if f_x(a_n) * f_x(b_n) > 0:
            raise IntervalError(a_n, b_n, 1)

        self.number_of_division = number_of_division
        self.f_x = f_x
        self.a_n = a_n
        self.b_n = b_n
        super().__init__(**kwargs)

    def check_interval(self):
        """This is the checking function for interval for bisection and trisection methods

        :raises IntervalError: In case of bad interval selected (multiple or no root.s)
        :return: Low or high interval value if they are a sufficient approximation else return empty
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
        """ Method for solving root finding problem using Dividing Search

        :raises IntervalError: In case of bad interval selected (multiple or no root.s)
        :return: The solution in a first return and the number of iteration in a second
        :rtype: Union[list(float), int]
        """
        iteration_number = 0
        best_c_n = self.check_interval()

        # Best solution already available ... LUCKY
        if best_c_n:
            return best_c_n, iteration_number

        # Otherwise, we compute the solution
        for iteration_number in range(1, self.nb_iteration):
            c_n = [self.a_n] + \
                  [((self.number_of_division - n) * self.a_n + n * self.b_n) / self.number_of_division
                   for n in range(1, self.number_of_division)] + \
                  [self.b_n]
            f_c_n = [self.f_x(i) for i in c_n]

            best_c_n.append(c_n[np.argmin(f_c_n)])

            if abs(self.f_x(best_c_n[-1])) < self.tolerance:
                break

            evol_array = [(c_n[n], c_n[n + 1]) for n in range(len(f_c_n) - 1) if f_c_n[n] * f_c_n[n + 1] < 0]

            if len(evol_array) > 1:
                raise IntervalError(self.a_n, self.b_n, 2)

            self.a_n = evol_array[0][0]
            self.b_n = evol_array[0][1]

        return best_c_n, iteration_number


class BisectionSolverClass(DividingSearchSolverClass):
    """This is the bisection solver class - a type of dividing search algorithm with a space division of 2

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param a_n: Low limit of studied interval
    :type a_n: float
    :param b_n: High limit of studied interval
    :type b_n: float
    :param \**kwargs: Complementary keyword arguments for the optional argument of IterativeAlgorithm class
                   (see IterativeAlgorithm)
    :type \**kwargs: dict
    """

    def __init__(self, f_x: Callable, a_n: float, b_n: float, **kwargs):
        super().__init__(f_x, a_n, b_n, number_of_division=2, **kwargs)


class TrisectionSolverClass(DividingSearchSolverClass):
    """This is the trisection solver class - a type of dividing search algorithm with a space division of 3

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param a_n: Low limit of studied interval
    :type a_n: float
    :param b_n: High limit of studied interval
    :type b_n: float
    :param \**kwargs: Complementary keyword arguments for the optional argument of IterativeAlgorithm class
                   (see IterativeAlgorithm)
    :type \**kwargs: dict
    """

    def __init__(self, f_x: Callable, a_n: float, b_n: float, **kwargs):
        super().__init__(f_x, a_n, b_n, number_of_division=3, **kwargs)


class FixedPointSolverClass(IterativeAlgorithm):
    """This is the fixed point solver class - a type of iterative algorithm

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param p_0: Starting point of fixed point algorithm. Can be a tuple depending on the algorithm
    :type p_0: Union[float, tuple]
    :param \**kwargs: Complementary keyword arguments for the optional argument of IterativeAlgorithm class
                   (see IterativeAlgorithm)
    :type \**kwargs: dict
    """

    def __init__(self, f_x: Callable, p_0: [float, tuple], **kwargs):
        self.f_x = f_x
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
    """This is the Newton solver class - a type of fixed point algorithm

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param df_x: Handle of the derivative of the studied function
    :type df_x: Callable
    :param p_0: Starting point of the Newton algorithm
    :type p_0: float
    :param \**kwargs: Complementary keyword arguments for the optional argument of IterativeAlgorithm class
                   (see IterativeAlgorithm)
    :type \**kwargs: dict
    :raises ZeroDivisionError: If derivative at p_0 is nul
    """

    def __init__(self, f_x: Callable, df_x: Callable, p_0: float, **kwargs):
        if df_x(p_0) == 0:
            raise ZeroDivisionError
        f_x_super: Callable = lambda x: x - f_x(x) / df_x(x)
        super().__init__(f_x_super, p_0, **kwargs)


class SecantSolverClass(FixedPointSolverClass):
    """This is the Secant solver class - a type of fixed point algorithm

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param p_0: Starting point of the secant algorithm
    :type p_0: tuple
    :param \**kwargs: Complementary keyword arguments for the optional argument of IterativeAlgorithm class
                   (see IterativeAlgorithm)
    :type \**kwargs: dict
    :raises ZeroDivisionError: If estimated derivative at p_0 is nul
    """

    def __init__(self, f_x: Callable, p_0: tuple, **kwargs):
        if f_x(p_0[0]) - f_x(p_0[1]) == 0:
            raise ZeroDivisionError
        f_x_super: Callable = lambda x: (x[1], secant_function(f_x, x))
        super().__init__(f_x_super, p_0, **kwargs)


class RegularFalsiSolverClass(FixedPointSolverClass):
    """This is the Regula Falsi solver class - a type of fixed point algorithm

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param p_0: Starting point of the regula falsi algorithm
    :type p_0: tuple
    :param \**kwargs: Complementary keyword arguments for the optional argument of IterativeAlgorithm class
                   (see IterativeAlgorithm)
    :type \**kwargs: dict
    :raises ZeroDivisionError: If estimated derivative at p_0 is nul
    """

    def __init__(self, f_x: Callable, p_0: tuple, **kwargs):
        if f_x(p_0[0]) - f_x(p_0[1]) == 0:
            raise ZeroDivisionError
        f_x_super: Callable = lambda x: regula_falsi_function(f_x, x)
        super().__init__(f_x_super, p_0, **kwargs)


class SteffensenSolverClass(FixedPointSolverClass):
    """This is the Steffensen solver class - a type of fixed point algorithm

    :param f_x: Handle of studied function
    :type f_x: Callable
    :param p_0: Starting point of the steffensen algorithm
    :type p_0: float
    :param \**kwargs: Complementary keyword arguments for the optional argument of IterativeAlgorithm class
                   (see IterativeAlgorithm)
    :type \**kwargs: dict
    """

    def __init__(self, f_x: Callable, p_0: float, **kwargs):
        f_x_super: Callable = lambda x: self.steffensen_function(f_x, x)
        super().__init__(f_x_super, p_0, **kwargs)

    @staticmethod
    def steffensen_function(f_x: Callable, p_n: float) -> float:
        """This is a helper function to decide which point is the nearest to the root of a function f(x)

        :param f_x: Handle of studied function
        :type f_x: Callable
        :param p_n: Value to compute with steffensen
        :type p_n: float
        :raises ZeroDivisionError: If denominator of update rule is nul
        :return: The new tuple with the regula falsi condition
        :rtype: tuple
        """

        if (f_x(f_x(p_n)) - 2 * f_x(p_n) + p_n) == 0:
            raise ZeroDivisionError

        return p_n - pow(f_x(p_n) - p_n, 2) / (f_x(f_x(p_n)) - 2 * f_x(p_n) + p_n)
