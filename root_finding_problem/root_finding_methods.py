from typing import Callable

class IntervalError(Exception):

    # Constructor or Initializer
    def __init__(self, a, b, error_code, *args):
        super().__init__(args)
        self.a = a
        self.b = b

        # Error message thrown is saved in msg
        if error_code == 1:
            self.msg = 'No/multiple roots for stater interval [ %f, %f ]' % (a, b)
        elif error_code == 2:
            self.msg = 'Multiple roots in interval [ %f, %f ]' % (a, b)

    # __str__ is to print() the value
    def __str__(self):
        return repr(self.msg)


def check_interval(f, a, b, tolerance):
    if f(a) * f(b) > 0:
        raise IntervalError(a, b, 1)

    if abs(f(a)) <= tolerance:
        return a
    elif abs(f(b)) <= tolerance:
        return b
    else:
        return None


def check_interval_value(f, a, b, tolerance):
    if f(a) * f(b) > 0:
        raise IntervalError(a, b, 1)

    if abs(f(a)) <= tolerance:
        return a
    elif abs(f(b)) <= tolerance:
        return b
    else:
        return []


def bisection_method(f: Callable, a: float, b: float, nb_iteration: int, tolerance: float, true_value: float):
    iteration_number = 1
    return_value = check_interval(f, a, b, tolerance)
    error_estimation = []

    while iteration_number < nb_iteration and return_value is None:
        c: float = a + (b - a) / 2

        error_estimation.append(min(abs(true_value - c),
                                    abs(true_value - a),
                                    abs(true_value - b)))

        if abs(f(c)) <= tolerance:
            return_value = c
            break

        if f(a) * f(c) < 0:
            b = c
        elif f(c) * f(b) < 0:
            a = c
        else:
            raise IntervalError(a, b, 2)

        iteration_number += 1

    if return_value is None:
        return_value = a if min(abs(f(a)), abs(f(b))) == abs(f(a)) else b

    return return_value, iteration_number, error_estimation


def trisection_method(f: Callable, a: float, b: float, nb_iteration: int, tolerance: float, true_value: float):
    iteration_number = 1
    return_value = check_interval(f, a, b, tolerance)
    error_estimation = []

    while iteration_number < nb_iteration and return_value is None:
        c1: float = (2 * a + b) / 3
        c2: float = (a + 2 * b) / 3

        error_estimation.append(min(abs(true_value - c1),
                                    abs(true_value - c2),
                                    abs(true_value - a),
                                    abs(true_value - b)))

        if abs(f(c1)) <= tolerance:
            return_value = c1
            break

        if abs(f(c2)) <= tolerance:
            return_value = c2
            break

        if f(a) * f(c1) < 0:
            b = c1
        elif f(c1) * f(c2) < 0:
            a = c1
            b = c2
        elif f(c2) * f(b) < 0:
            a = c2
        else:
            raise IntervalError(a, b, 2)

        iteration_number += 1

    if return_value is None:
        return_value = a if min(abs(f(a)), abs(f(b))) == abs(f(a)) else b

    return return_value, iteration_number, error_estimation


def fixed_point_method(f: Callable, p0: float, nb_iteration: int, tolerance: float, true_value: float):
    iteration_number = 1
    return_value = None
    error_estimation = []

    while iteration_number < nb_iteration and return_value is None:

        error_estimation.append(abs(true_value - f(p0)))

        if abs(p0 - f(p0)) <= tolerance:
            return_value = f(p0)
            break

        iteration_number += 1
        p0 = f(p0)

    if return_value is None:
        return_value = p0

    return return_value, iteration_number, error_estimation


def newton_method(f: Callable, df: Callable, p0: float, nb_iteration: int, tolerance: float, true_value: float):
    iteration_number = 1
    return_value = None
    error_estimation = []

    while iteration_number < nb_iteration and return_value is None:

        p = p0 - f(p0) / df(p0)

        error_estimation.append(abs(true_value - p))

        if abs(p0 - p) <= tolerance:
            return_value = p
            break

        iteration_number += 1
        p0 = p

    if return_value is None:
        return_value = p0

    return return_value, iteration_number, error_estimation


def secant_method(f: Callable, a: float, b: float, nb_iteration: int, tolerance: float, true_value: float):
    iteration_number = 1
    return_value = check_interval(f, a, b, tolerance)
    error_estimation = []
    q0 = f(a)
    q1 = f(b)

    while iteration_number < nb_iteration and return_value is None:

        p = b - q1 * (b - a) / (q1 - q0)
        error_estimation.append(abs(true_value - p))

        if abs(p - b)/abs(b) <= tolerance:
            return_value = p
            break

        a = b
        q0 = q1
        b = p
        q1 = f(p)

        iteration_number += 1

    if return_value is None:
        return_value = b - q1 * (b - a) / (q1 - q0)

    return return_value, iteration_number, error_estimation


def regula_falsi_method(f: Callable, a: float, b: float, nb_iteration: int, tolerance: float, true_value: float):
    iteration_number = 1
    return_value = check_interval(f, a, b, tolerance)
    error_estimation = []
    q0 = f(a)
    q1 = f(b)

    while iteration_number < nb_iteration and return_value is None:

        p: float = b - q1 * (b - a) / (q1 - q0)
        error_estimation.append(abs(true_value - p))

        if abs(p - b)/abs(b) <= tolerance:
            return_value = p
            break

        q = f(p)

        if q * q1 < 0:
            a = b
            q0 = q1

        b = p
        q1 = f(p)

        iteration_number += 1

    if return_value is None:
        return_value = b - q1 * (b - a) / (q1 - q0)

    return return_value, iteration_number, error_estimation


def steffensen_method(f: Callable, p0: float, nb_iteration: int, tolerance: float, true_value: float):
    iteration_number = 1
    return_value = None
    error_estimation = []

    while iteration_number < nb_iteration and return_value is None:
        p1 = f(p0)
        p2 = f(p1)

        try:
            p = p0 - pow(p1 - p0, 2) / (p2 - 2 * p1 + p0)
        except ZeroDivisionError:
            return_value = p0
            break

        error_estimation.append(abs(true_value - p))

        if abs(p0 - p)/abs(p) <= tolerance:
            return_value = p
            break

        iteration_number += 1
        p0 = p

    if return_value is None:
        return_value = p0

    return return_value, iteration_number, error_estimation
