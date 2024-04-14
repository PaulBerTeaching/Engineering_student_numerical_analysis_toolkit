"""
This module shows some examples of root finding methods
"""
from typing import Callable, Union


class IntervalError(Exception):
    """This is a children class of exception focusing on interval errors for bisection and trisection methods

    Attributes
    ----------
    msg : string
        Error message to display

    Methods
    ---------
    __str__()
        Display the message error
    """

    def __init__(self, a_n, b_n, error_code, *args):
        """Constructor method

        Parameters
        ----------
        a_n : float
            Low limit of studied interval.
        b_n : float
            High limit of studied interval.
        error_code : int
            Type of error to display
        """
        super().__init__(args)

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

    Parameters
    ----------
    f_x: Callable
        Handle of studied function
    a_0: float
        Low limit of studied interval
    b_0: float
        High limit of studied interval
    tolerance: float
        Precision limit

    Raises
    ------
    IntervalError
        In case of bad interval selected (multiple or no root.s)

    Returns
    ------
    list
        Low or high interval value if they are a sufficient approximation else return empty
    """

    if f_x(a_0) * f_x(b_0) > 0:
        raise IntervalError(a_0, b_0, 1)

    if abs(f_x(a_0)) < tolerance:
        return [a_0]

    return [b_0] if abs(f_x(b_0)) < tolerance else []


def root_nearest_point(f_x: Callable, c_1: float, c_2: float):
    """This is a helper function to decide which point is the nearest to the root of a function f(x) for the trisection
    method

    Parameters
    ----------
    f_x: Callable
        Handle of studied function
    c_1: float
        Low limit of studied interval
    c_2: float
        High limit of studied interval

    Returns
    ---------
    float
        The best root between c_1 and c_2
    """
    return c_1 if min(abs(f_x(c_1)), abs(f_x(c_2))) == abs(f_x(c_1)) else c_2


def secant_function(f_x: Callable, x_ab: tuple) -> float:
    """This is a helper function to update the values of x_ab using the secant rule

    Parameters
    ---------
    f_x: Callable
        Handle of studied function
    x_ab: tuple
        The tuple values to update following the secant update rule

    Returns
    -------
    tuple
        The new tuple with the secant update rule
    """
    if f_x(x_ab[1]) - f_x(x_ab[0]) == 0:
        raise ZeroDivisionError

    return x_ab[1] - f_x(x_ab[1]) * (x_ab[1] - x_ab[0]) / (f_x(x_ab[1]) - f_x(x_ab[0]))


def regula_falsi_function(f_x: Callable, x_ab: tuple) -> tuple:
    """This is a helper function to update the values of x_ab using the regula falsi rule

    Parameters
    ---------
    f_x: Callable
        Handle of studied function
    x_ab: tuple
        The tuple values to update following the regula falsi update rule

    Returns
    -------
    tuple
        The new tuple with the regula falsi update rule
    """

    c_n = secant_function(f_x, x_ab)

    return (x_ab[1], c_n) if f_x(c_n) * f_x(x_ab[1]) < 0 else (x_ab[0], c_n)


def steffensen_function(f_x: Callable, p_n: float) -> float:
    """This is a helper function to update the values of x_ab using the steffensen rule

    Parameters
    ---------
    f_x: Callable
        Handle of studied function
    p_n: float
        The float values to update following the steffensen update rule

    Returns
    -------
    float
        The new float with the steffensen update rule
    """

    if (f_x(f_x(p_n)) - 2 * f_x(p_n) + p_n) == 0:
        raise ZeroDivisionError

    return p_n - pow(f_x(p_n) - p_n, 2) / (f_x(f_x(p_n)) - 2 * f_x(p_n) + p_n)


def bisection_method(f_x: Callable, a_n: float, b_n: float, nb_iteration: int, tolerance: float):
    """This is the bisection algorithm

    Parameters
    ----------
    f_x: Callable
        Handle of studied function
    a_n: float
        Low limit of studied interval
    b_n: float
        High limit of studied interval
    nb_iteration: int
        Maximum iteration limit
    tolerance: float
        Precision limit

    Raises
    ----------
    IntervalError
        In case of bad interval selected (multiple or no root.s)

    Returns
    ----------
    list, int
        List of best estimated values, number of iteration
    """
    iteration_number = 0
    c_n = list(check_interval(f_x, a_n, b_n, tolerance))

    if c_n[0]:
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


def trisection_method(f_x: Callable, a_n: float, b_n: float, nb_iteration: int, tolerance: float):
    """This is the trisection algorithm

    Parameters
    ----------
    f_x: Callable
        Handle of studied function
    a_n: float
        Low limit of studied interval
    b_n: float
        High limit of studied interval
    nb_iteration: int
        Maximum iteration limit
    tolerance: float
        Precision limit

    Raises
    ----------
    IntervalError
        In case of bad interval selected (multiple or no root.s)

    Returns
    ----------
    list, int
        List of best estimated values, number of iteration
    """
    iteration_number = 0
    best_c_n = list(check_interval(f_x, a_n, b_n, tolerance))

    if best_c_n[0]:
        return best_c_n, iteration_number

    for iteration_number in range(1, nb_iteration):
        c_1_n: float = (2 * a_n + b_n) / 3
        c_2_n: float = (a_n + 2 * b_n) / 3

        if abs(f_x(c_1_n)) < tolerance:
            best_c_n.append(c_1_n)
            break

        if abs(f_x(c_2_n)) < tolerance:
            best_c_n.append(c_2_n)
            break

        best_c_n.append(root_nearest_point(f_x, c_1_n, c_2_n))

        if f_x(a_n) * f_x(c_1_n) < 0:
            b_n = c_1_n
        elif f_x(c_1_n) * f_x(c_2_n) < 0:
            a_n = c_1_n
            b_n = c_2_n
        elif f_x(c_2_n) * f_x(b_n) < 0:
            a_n = c_2_n
        else:
            raise IntervalError(a_n, b_n, 2)

    return best_c_n, iteration_number


def fixed_point_method(f_x: Callable, p_0: Union[float, tuple], nb_iteration: int, tolerance: float):
    """This is the fixed point algorithm

    Parameters
    ----------
    f_x: Callable
        Handle of studied function
    p_0: Union[float, tuple]
        Starting point of fixed point algorithm
    nb_iteration: int
        Maximum iteration limit
    tolerance: float
        Precision limit

    Returns
    --------
    list, int
        List of best estimate, number of iteration
    """
    iteration_number = 0
    p_n = [p_0]

    for iteration_number in range(1, nb_iteration):
        try:
            p_n.append(f_x(p_n[-1]))

            if isinstance(p_n[-1], tuple):
                if abs(p_n[-1][-1] - f_x(p_n[-1])[-1]) < tolerance:
                    break
            elif abs(p_n[-1] - f_x(p_n[-1])) < tolerance:
                break

        except ZeroDivisionError:
            print('f(x) reach a 0 value - Maximum precision reached or error in selection')
            break

    return p_n, iteration_number + 1


def newton_method(f_x: Callable, df_x: Callable, p_0: float, nb_iteration: int, tolerance: float):
    """This is the Newton algorithm

    Parameters
    ----------
    f_x: Callable
        Handle of studied function
    df_x: Callable
        Handle of derivative of studied function
    p_0: float
        Starting point of fixed point algorithm
    nb_iteration: int
        Maximum iteration limit
    tolerance: float
        Precision limit
    Raises
    -------
    ZeroDivisionError
        In case of derivative being 0 at p_0

    Returns
    ---------
    list, int
        List of best estimate, number of iteration
    """

    if df_x(p_0) == 0:
        raise ZeroDivisionError

    p_n, iteration_number = fixed_point_method(lambda x: x - f_x(x) / df_x(x), p_0, nb_iteration, tolerance)

    return p_n, iteration_number


def secant_method(f_x: Callable, p_0: tuple, nb_iteration: int, tolerance: float):
    """This is the Secant algorithm

    Parameters
    -----------
    f_x: Callable
        Handle of studied function
    p_0: tuple
        Starting tuple of fixed point algorithm - starting interval
    nb_iteration: int
        Maximum iteration limit
    tolerance: float
        Precision limit

    Raises
    -------
    ZeroDivisionError
        In case of the estimate of the derivative being 0 at p_0

    Returns
    ---------
    list, int
        List of best estimate, number of iteration
    """
    if f_x(p_0[0]) - f_x(p_0[1]) == 0:
        raise ZeroDivisionError

    p_n, iteration_number = fixed_point_method(lambda x: (x[1], secant_function(f_x, x)), p_0, nb_iteration, tolerance)

    return p_n, iteration_number


def regula_falsi_method(f_x: Callable, p_0: tuple, nb_iteration: int, tolerance: float):
    """This is the regula Falsi algorithm

    Parameters
    -----------
    f_x: Callable
        Handle of studied function
    p_0: tuple
        Starting tuple of fixed point algorithm - starting interval
    nb_iteration: int
        Maximum iteration limit
    tolerance: float
        Precision limit

    Raises
    -------
    ZeroDivisionError
        In case of the estimate of the derivative being 0 at p_0

    Returns
    ---------
    list, int
        List of best estimate, number of iteration
    """
    if f_x(p_0[0]) - f_x(p_0[1]) == 0:
        raise ZeroDivisionError

    p_n, iteration_number = fixed_point_method(lambda x: regula_falsi_function(f_x, x), p_0, nb_iteration, tolerance)

    return p_n, iteration_number


def steffensen_method(f_x: Callable, p_0: float, nb_iteration: int, tolerance: float):
    """This is the Steffensen algorithm

    Parameters
    -----------
    f_x: Callable
        Handle of studied function
    p_0: float
        Starting point of fixed point algorithm
    nb_iteration: int
        Maximum iteration limit
    tolerance: float
        Precision limit

    Returns
    -----------
    list, int
        List of best estimate, number of iteration
    """

    p_n, iteration_number = fixed_point_method(lambda x: steffensen_function(f_x, x), p_0, nb_iteration, tolerance)

    return p_n, iteration_number
