

def dichotomie_2(f, a, b, nb_iteration):
    iteration_number = 1

    while iteration_number < nb_iteration:
        c = a + (b - a) / 2

        if f(a) * f(c) < 0:
            b = c
        elif f(c) * f(b) < 0:
            a = c

        iteration_number += 1

    return_value = c

    return return_value

def dichotomie(a, b, nb_iteration):
    iteration_number = 1

    while iteration_number < nb_iteration:

        c = a + (b - a) / 2

        fa = a ** 3 + 4 * a ** 2 - 10
        fb = b ** 3 + 4 * b ** 2 - 10
        fc = c ** 3 + 4 * c ** 2 - 10

        if fa * fc < 0:
            b = c
        elif fc* fb < 0:
            a = c

        iteration_number += 1

    return_value = c

    return return_value