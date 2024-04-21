"""
Created on Fri Apr 24 00:40:30 2020

@author: paulberraute
"""

import pytest
from .context import tf, poly

def test_char():
    res_bisect, it_bisect = tf.BisectionSolverClass(f_x=poly.Polynomial((-10, 0, 4, 1)), a_n=0,
                                                    b_n=5,
                                                    nb_iteration=100).solve()

    assert res_bisect == 1.3652300134140969


if __name__ == '__main__':
    try:
        test_char()
    except:
        print("Votre fonction de dichotomie ne retourne pas exactement la valeur attendu. Revérifiez s'il vous plaît.")